/**
 * @file main.cpp
 * @brief Command-line driver for wannier2sparse.
 *
 * Parses arguments, loads the Wannier90 model, builds requested operators, applies
 * optional Wigner-Seitz corrections and self-checks, and writes CSR output.
 */
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cctype>

#include "w2sp_arguments.hpp"
#include "tbmodel.hpp"
#include "hopping_list.hpp"
#include "operator_algebra.hpp"
#include "descriptor.hpp"
#include "gauge.hpp"
#include "local_operators.hpp"
#include "checks.hpp"
#include "bundle_writer.hpp"
#include "run_config.hpp"
#include "win_parser.hpp"
#include "qe_xml_parser.hpp"
#include "logger.hpp"
#include "run_stats.hpp"
#include "run_report.hpp"
#include <set>
#include <chrono>
#include <ctime>
#include <sys/stat.h>   // mkdir (POSIX)
#include <stdexcept>

using namespace std;

/// Tool version, echoed by --version and into the bundle manifest / .out receipt.
static const char* W2SP_VERSION = "1.0";

/**
 * @brief Local wall-clock timestamp "YYYY-MM-DD HH:MM:SS" for the run receipt.
 * @return formatted timestamp
 */
static string now_timestamp()
{
    const std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tmv;
    localtime_r(&t, &tmv);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tmv);
    return string(buf);
}

/**
 * @brief Reconstruct the invocation command line from argv for the receipt.
 * @param argc argument count
 * @param argv argument vector
 * @return space-joined command line (arguments containing spaces are quoted)
 */
static string join_command_line(int argc, char* argv[])
{
    string out;
    for (int i = 0; i < argc; ++i)
    {
        if (i) out += ' ';
        const string a = argv[i];
        if (a.find(' ') != string::npos) out += '"' + a + '"';
        else                             out += a;
    }
    return out;
}

/**
 * @brief Check whether a file exists and is readable.
 * @param path file path
 * @return true if the file can be opened
 */
static bool file_exists(const string& path)
{
    ifstream f(path.c_str());
    return f.good();
}

/**
 * @brief Create a directory if it does not already exist (no-op for "." / empty).
 * @param path directory path
 */
static void make_dir(const string& path)
{
    if (path.empty() || path == ".") return;
    ::mkdir(path.c_str(), 0755);   // ignore EEXIST; a real write failure surfaces later
}

/**
 * @brief Map a log-level name to a Logger::Level (defaults to INFO).
 * @param s "trace" | "debug" | "info" | "warn" | "error"
 * @return matching level
 */
static Logger::Level level_from_string(const string& s)
{
    if (s == "trace") return Logger::TRACE;
    if (s == "debug") return Logger::DEBUG;
    if (s == "warn")  return Logger::WARN;
    if (s == "error") return Logger::ERROR;
    return Logger::INFO;
}

/**
 * @brief Build a descriptor with the given fields.
 * @param observable observable name
 * @param component component code (e.g. "X", "XSZ")
 * @param units physical units
 * @param provenance how the operator was produced
 * @return populated OperatorDescriptor
 */
static OperatorDescriptor describe(const string& observable, const string& component,
                                   const string& units, const string& provenance)
{
    OperatorDescriptor d;
    d.observable = observable; d.component = component;
    d.units = units;           d.provenance = provenance;
    return d;
}

/**
 * @brief Descriptor for a generated operator code (VX, SZ, VXSZ, ...).
 * @param code operator code
 * @return populated OperatorDescriptor
 */
static OperatorDescriptor op_descriptor(const string& code)
{
    if (code.size() == 2 && code[0] == 'V')
        return describe("velocity", string(1, code[1]), "eV*Angstrom", "H_ij * (-i) dr");
    if (code.size() == 2 && code[0] == 'S')
        return describe("spin", string(1, code[1]), "hbar/2", "spin labels (_s+_/_s-_)");
    if (code.size() == 4 && code[0] == 'V' && code[2] == 'S')
        return describe("spin_current", code.substr(1), "eV*Angstrom*hbar/2", "label-based V*S");
    return describe("operator", code, "", "generated");
}

/**
 * @brief Build a requested operator as a primitive-cell hopping_list.
 * @param model tight-binding model
 * @param op operator code
 * @param out output hopping_list
 * @return true on success
 *
 * Dispatches the operator name to the matching tbmodel builder. Returns false for a
 * name the parser accepted but this build cannot produce.
 */
static bool build_operator(tbmodel& model, const string& op, hopping_list& out)
{
    if (op == "VX") { out = model.createHoppingCurrents_list(0); return true; }
    if (op == "VY") { out = model.createHoppingCurrents_list(1); return true; }
    if (op == "VZ") { out = model.createHoppingCurrents_list(2); return true; }

    if (op == "SX") { out = model.createHoppingSpinDensity_list('x'); return true; }
    if (op == "SY") { out = model.createHoppingSpinDensity_list('y'); return true; }
    if (op == "SZ") { out = model.createHoppingSpinDensity_list('z'); return true; }

    if (op.size() == 4 && op[0] == 'V' && op[2] == 'S')
    {
        const int dir = (op[1] == 'X') ? 0 : (op[1] == 'Y') ? 1 : (op[1] == 'Z') ? 2 : -1;
        const char sdir = static_cast<char>(std::tolower(static_cast<unsigned char>(op[3])));
        if (dir >= 0 && (sdir == 'x' || sdir == 'y' || sdir == 'z'))
        {
            out = model.createHoppingSpinCurrents_list(dir, sdir);
            return true;
        }
    }
    return false;
}

enum op_kind { K_GENERIC, K_HAM, K_SPIN, K_ORBITAL }; ///< operator category for checks

/**
 * @brief Evaluate requested self-checks on an assembled operator.
 * @param args parsed CLI arguments
 * @param opname operator name for the sidecar filename
 * @param hl assembled operator in hopping_list form
 * @param kind operator category
 *
 * Writes a <prefix>.<op>.check sidecar. Never touches the CSR.
 */
static void run_checks(const W2SP_arguments& args, const string& opname,
                       const hopping_list& hl, op_kind kind)
{
    if (args.check.empty()) return;
    const string& w = args.check;
    auto want = [&](const char* n){ return w == "all" || w == n; };
    std::vector<check_result> rs;
    if (want("hermiticity")) rs.push_back(checks::hermiticity(hl));
    if (want("sum_rules"))   rs.push_back(checks::trace_rule(hl, kind == K_SPIN || kind == K_ORBITAL));
    if (want("aliasing"))    rs.push_back(checks::aliasing(hl, args.cellDim));
    if (want("algebra") && kind == K_SPIN)    rs.push_back(checks::spin_algebra());
    if (want("algebra") && kind == K_ORBITAL) rs.push_back(checks::orbital_algebra());
    if (want("bounds") && kind == K_HAM)
    {
        double a = 0, b = 0; spectral_bounds(supercell_matrix(args.cellDim, hl), a, b);
        check_result r; r.name = "bounds"; r.pass = (a < b); r.residual = b - a;
        r.detail = "[a,b]=[" + std::to_string(a) + "," + std::to_string(b) + "]";
        rs.push_back(r);
    }
    if (!rs.empty())
        checks::write_report(rs, args.output_dir + "/" + args.label + "." + opname + ".check");
}

/**
 * @brief Collect the wannier sites (orbital positions + spin tag) for the manifest.
 * @param model loaded tight-binding model
 * @return wannier sites in orbital order
 */
static std::vector<WannierSite> collect_wannier_sites(tbmodel& model)
{
    std::vector<WannierSite> sites;
    const auto id2spin = model.map_id2spin();
    int idx = 0;
    for (const auto& orb : model.orbPos_list)
    {
        WannierSite s;
        s.index = idx;
        s.label = std::get<0>(orb);
        s.cart  = std::get<1>(orb);
        s.spin  = id2spin.count(idx) ? id2spin.at(idx) : 0;
        sites.push_back(s);
        ++idx;
    }
    return sites;
}

/**
 * @brief Bundle (provenance) output mode: emit primitive operators + manifest.
 * @param args parsed CLI arguments (mode == "bundle")
 * @return exit code
 *
 * Builds the same primitive-cell operators as the CSR path but, instead of
 * expanding them into a supercell, writes them as _hr.dat-format data files inside
 * a `<LABEL>.w2sp/` bundle alongside a JSON manifest. The supercell dimensions are
 * not used. lsquant rebuilds the Hamiltonian/supercell from the bundle.
 */
static int run_bundle_mode(const W2SP_arguments& args, Logger& log, RunStats& stats, RunReport& rep)
{
    const string in_prefix = args.input_prefix();
    const string f_uc  = in_prefix + ".uc";
    const string f_xyz = in_prefix + ".xyz";
    const string f_hr  = in_prefix + "_hr.dat";
    const string f_ws  = in_prefix + "_wsvec.dat";

    const bool need_positions = !args.operators.empty();

    bool missing = false;
    if (!file_exists(f_hr))
    {
        log.error("required input file '" + f_hr + "' not found");
        missing = true;
    }
    if (need_positions)
        for (const string& f : {f_uc, f_xyz})
            if (!file_exists(f))
            {
                log.error("operator generation requires '" + f + "'");
                missing = true;
            }
    for (const auto& nf : args.op_files)
        if (!file_exists(nf.second))
        {
            log.error("--op-file '" + nf.first + "': file '" + nf.second + "' not found");
            missing = true;
        }
    if (missing) return 1;

    log.info("bundle mode, label '" + args.label + "'");

    tbmodel model;
    {
        ScopedTimer t(log, stats, "load model");
        if (file_exists(f_xyz)) model.readOrbitalPositions(f_xyz);
        if (file_exists(f_uc))  model.readUnitCell(f_uc);
        model.readWannierModel(f_hr);
        log.at(Logger::INFO) << "loaded " << f_hr << " (num_wann=" << model.hl.WannierBasisSize() << ")";
    }

    const bool wsvec_applied = file_exists(f_ws);
    if (wsvec_applied)
    {
        log.info("applying Wigner-Seitz correction from " + f_ws);
        model.applyWsvec(f_ws);
    }

    SystemProvenance prov;
    prov.num_wann = model.hl.WannierBasisSize();
    prov.manual = args.manual;   // user-declared provenance from the input file
    if (file_exists(f_uc)) { prov.has_lattice = true; prov.lattice = model.lat_vecs; }
    if (file_exists(f_xyz)) prov.wannier_sites = collect_wannier_sites(model);

    // Wannier90 provenance from the .win file (config-specified path, else
    // <prefix>.win). Absent or unparsable => the manifest block stays null.
    const string f_win_prov = !args.win_path.empty() ? args.win_path : (in_prefix + ".win");
    if (file_exists(f_win_prov))
    {
        ScopedTimer t(log, stats, "parse .win");
        try { prov.wann = parse_win_provenance(f_win_prov); log.info("parsed Wannier provenance from " + f_win_prov); }
        catch (const std::exception& e)
        {
            log.warn(string("could not parse .win provenance (") + f_win_prov + "): " + e.what());
        }
    }
    else if (!args.win_path.empty())
        log.warn(".win provenance file '" + args.win_path + "' not found");

    // DFT (Quantum ESPRESSO) provenance from data-file-schema.xml. Only parsed
    // when a path was given (no default location); absent/unparsable => null.
    if (!args.qe_xml_path.empty())
    {
        if (file_exists(args.qe_xml_path))
        {
            ScopedTimer t(log, stats, "parse QE XML");
            try { prov.dft = parse_qe_xml(args.qe_xml_path); log.info("parsed DFT provenance from " + args.qe_xml_path); }
            catch (const std::exception& e)
            {
                log.warn(string("could not parse QE XML provenance (") + args.qe_xml_path + "): " + e.what());
            }
        }
        else
            log.warn("QE XML provenance file '" + args.qe_xml_path + "' not found");
    }

    std::vector<BundleOperator> ops;
    {
    ScopedTimer t_build(log, stats, "build operators");

    // Hamiltonian (always).
    {
        BundleOperator op;
        op.name = "HAM";
        op.desc = describe("hamiltonian", "", "eV", "wannier90 _hr.dat");
        if (args.emit_descriptor)
        {
            double emin = 0.0, emax = 0.0;
            if (read_eig_bounds(in_prefix + ".eig", emin, emax))
            {
                op.desc.has_bounds = true; op.desc.a = emin; op.desc.b = emax;
                op.desc.provenance += "; bounds from .eig";
            }
            else
                log.warn("bundle mode skips Lanczos bounds (they need k-sampling); provide "
                         + in_prefix + ".eig for HAM bounds");
        }
        op.hl = model.hl;
        ops.push_back(op);
    }

    // Requested generated operators (velocity / spin / spin-current).
    for (const auto& opname : args.operators)
    {
        hopping_list h;
        if (!build_operator(model, opname, h))
        {
            log.error("operator '" + opname + "' is not supported by this build");
            return 1;
        }
        BundleOperator op;
        op.name = opname;
        op.desc = op_descriptor(opname);
        op.hl   = h;
        ops.push_back(op);
    }

    // External operators ingested verbatim from _hr.dat-format files.
    for (const auto& nf : args.op_files)
    {
        BundleOperator op;
        op.name = nf.first;
        op.desc = describe("external", nf.first, "", "ingested from " + nf.second);
        op.hl   = model.readOperatorModel(nf.second);
        ops.push_back(op);
    }

    try {
        // Exact spin operators via the gauge transform (same primitive builder as the CSR path).
        if (args.exact_spin)
        {
            const string f_spn = in_prefix + ".spn";
            const string f_umat = in_prefix + "_u.mat";
            if (!file_exists(f_spn) || !file_exists(f_umat))
            {
                log.error("--exact-spin requires " + f_spn + " and " + f_umat);
                return 1;
            }
            gauge_data g = read_gauge(in_prefix);
            hopping_list rawH = create_hopping_list(read_wannier_file(f_hr));
            if (g.num_wann != rawH.WannierBasisSize())
            {
                log.at(Logger::ERROR) << "num_wann mismatch between gauge (" << g.num_wann
                                      << ") and Hamiltonian (" << rawH.WannierBasisSize() << ")";
                return 1;
            }
            std::set<hopping_list::cellID_t> Rs;
            for (const auto& h : rawH.hoppings) Rs.insert(get<0>(h));
            const std::vector<hopping_list::cellID_t> Rset(Rs.begin(), Rs.end());

            const char* names[3] = {"SXexact", "SYexact", "SZexact"};
            for (int alpha = 0; alpha < 3; ++alpha)
            {
                hopping_list s = exact_spin_operator(g, alpha, Rset);
                if (wsvec_applied) s = apply_wsvec(s, read_wsvec(f_ws));
                BundleOperator op;
                op.name = names[alpha];
                op.desc = describe("spin", string(1, "XYZ"[alpha]), "hbar/2",
                                   "exact gauge transform V^dag S_B V from .spn");
                op.hl   = s;
                ops.push_back(op);
            }
        }

        // Orbital angular momentum via the projector route.
        if (args.orbital_l)
        {
            const string f_amn = in_prefix + ".amn";
            const string f_win = in_prefix + ".win";
            const string f_umat = in_prefix + "_u.mat";
            if (!file_exists(f_amn) || !file_exists(f_umat) || !file_exists(f_win))
            {
                log.error("--orbital-L requires " + f_amn + ", " + f_umat + " and " + f_win);
                return 1;
            }
            const std::vector<int> shells = parse_projection_shells(f_win);
            gauge_data g = read_gauge(in_prefix);
            amn_data   A = read_amn(f_amn);
            hopping_list rawH = create_hopping_list(read_wannier_file(f_hr));
            std::set<hopping_list::cellID_t> Rs;
            for (const auto& h : rawH.hoppings) Rs.insert(get<0>(h));
            const std::vector<hopping_list::cellID_t> Rset(Rs.begin(), Rs.end());

            const char* names[3] = {"LX", "LY", "LZ"};
            for (int alpha = 0; alpha < 3; ++alpha)
            {
                hopping_list L = orbital_L_operator(g, A, shells, alpha, Rset);
                if (wsvec_applied) L = apply_wsvec(L, read_wsvec(f_ws));
                BundleOperator op;
                op.name = names[alpha];
                op.desc = describe("orbital_L", string(1, "XYZ"[alpha]), "hbar",
                                   "projector route C=A^dag V, C^dag L_local C");
                op.hl   = L;
                ops.push_back(op);
            }
        }
    } catch (const std::exception& e) {
        log.error(e.what());
        return 1;
    }
    } // build operators

    BundleSpec spec;
    spec.label = args.label;
    spec.truncation_threshold = args.has_truncation
        ? args.truncation_threshold : std::numeric_limits<double>::epsilon();
    spec.ndegen_applied = true;
    spec.wsvec_applied = wsvec_applied;
    if (wsvec_applied) spec.wsvec_src = f_ws;   // copied verbatim into the bundle

    string bundle_dir;
    {
        ScopedTimer t(log, stats, "write bundle");
        bundle_dir = write_bundle(spec, prov, ops, args.output_dir);
    }
    log.at(Logger::INFO) << "wrote bundle with " << ops.size() << " operator(s) -> " << bundle_dir;
    log.info("manifest: " + bundle_dir + "/manifest.json");

    // Output ledger for the .out receipt: each primitive operator's data file plus
    // the input files that produced it (the provenance link).
    for (const auto& op : ops)
    {
        std::vector<string> s;
        if      (op.name == "HAM") s.push_back(f_hr);
        else if (op.name == "SXexact" || op.name == "SYexact" || op.name == "SZexact")
            { s.push_back(in_prefix + ".spn"); s.push_back(in_prefix + "_u.mat"); s.push_back(f_hr); }
        else if (op.name == "LX" || op.name == "LY" || op.name == "LZ")
            { s.push_back(in_prefix + ".amn"); s.push_back(in_prefix + "_u.mat"); s.push_back(in_prefix + ".win"); s.push_back(f_hr); }
        else
        {
            bool from_file = false;
            for (const auto& nf : args.op_files) if (nf.first == op.name) { s.push_back(nf.second); from_file = true; break; }
            if (!from_file) { s.push_back(f_hr); s.push_back(f_uc); s.push_back(f_xyz); }   // generated velocity/spin/current
        }
        if (wsvec_applied) s.push_back(f_ws);
        rep.add_output(op.name, bundle_dir + "/operators/" + op.name + ".hr.dat",
                       op.desc.observable, op.desc.component, op.desc.units, op.desc.provenance, s);
    }
    return 0;
}

/**
 * @brief Sparse (supercell CSR) output mode: expand and export operators as CSR.
 * @param args parsed + config-applied arguments
 * @param log  run logger
 * @param stats per-stage timing collector
 * @return exit code (0 on success)
 *
 * This is the original positional-CLI export path, unchanged in what it writes
 * (the CSR output is byte-stable); only the progress/error messages now route
 * through the logger and the heavy stages are wrapped in timers.
 */
static int run_sparse_mode(const W2SP_arguments& args, Logger& log, RunStats& stats, RunReport& rep)
{
    const string in_prefix = args.input_prefix();   // <project_dir>/<seed>, defaults to LABEL
    const string f_uc  = in_prefix + ".uc";
    const string f_xyz = in_prefix + ".xyz";
    const string f_hr  = in_prefix + "_hr.dat";

    // Every generated operator (velocity / spin / spin-current) needs the
    // orbital positions (.xyz) and lattice vectors (.uc); a plain Hamiltonian or
    // an externally ingested --op-file operator needs only its _hr.dat. So .uc
    // and .xyz are required only when at least one such operator is requested.
    const bool need_positions = !args.operators.empty() || !args.spin_currents.empty();

    bool missing = false;
    if (!file_exists(f_hr))
    {
        log.error("required input file '" + f_hr + "' not found");
        missing = true;
    }
    if (need_positions)
        for (const string& f : {f_uc, f_xyz})
            if (!file_exists(f))
            {
                log.error("operator generation requires '" + f + "'");
                missing = true;
            }
    for (const auto& nf : args.op_files)
        if (!file_exists(nf.second))
        {
            log.error("--op-file '" + nf.first + "': file '" + nf.second + "' not found");
            missing = true;
        }
    if (missing) return 1;

    log.info("sparse mode, label '" + args.label + "'");

    tbmodel model;
    {
        ScopedTimer t(log, stats, "load model");
        if (file_exists(f_xyz)) model.readOrbitalPositions(f_xyz);
        if (file_exists(f_uc))  model.readUnitCell(f_uc);
        model.readWannierModel(f_hr);
        log.at(Logger::INFO) << "loaded " << f_hr << " (num_wann=" << model.hl.WannierBasisSize() << ")";
    }

    // Optional Wigner-Seitz minimum-image correction (Plan 5): applied to every
    // operator since it acts on the shared Hamiltonian hopping list.
    const string f_ws = in_prefix + "_wsvec.dat";
    const bool ws = file_exists(f_ws);
    if (ws)
    {
        log.info("applying Wigner-Seitz correction from " + f_ws);
        model.applyWsvec(f_ws);
    }

    // Source lists for the .out ledger: a generated operator is built from the
    // Hamiltonian (+ .uc/.xyz for the geometric factor); a plain HAM only needs
    // its _hr.dat. The wsvec correction, when present, is folded into every one.
    auto with_ws = [&](std::vector<string> s) { if (ws) s.push_back(f_ws); return s; };
    const std::vector<string> hr_src  = with_ws({f_hr});
    const std::vector<string> gen_src = with_ws({f_hr, f_uc, f_xyz});

    const string prefix = args.output_dir + "/" + args.label;
    make_dir(args.output_dir);   // ensure the CSR destination exists

    log.at(Logger::INFO) << "creating the supercell (" << args.cellDim[0] << ","
                         << args.cellDim[1] << "," << args.cellDim[2] << ")";

    try {

    // Plan 10D: refuse to emit corrupted operators. Every generated operator
    // shares H's R-range (velocity/spin/current/exact-spin/orbital-L are all
    // built on the H Wigner-Seitz set), so guarding the (post-wsvec) H covers
    // them; external op-files are guarded individually below.
    guard_minimum_image(model.hl, args.cellDim);

    {
        ScopedTimer t(log, stats, "write HAM");
        log.info("writing Hamiltonian -> " + prefix + ".HAM.CSR");
        save_supercell_as_csr(args.cellDim, model.hl, prefix + ".HAM.CSR");
        run_checks(args, "HAM", model.hl, K_HAM);
        rep.add_output("HAM", prefix + ".HAM.CSR", "hamiltonian", "", "eV", "wannier90 _hr.dat", hr_src);

        if (args.emit_descriptor)
        {
            OperatorDescriptor d = describe("hamiltonian", "", "eV", "wannier90 _hr.dat");
            double emin = 0.0, emax = 0.0;
            if (read_eig_bounds(in_prefix + ".eig", emin, emax))
                d.provenance += "; bounds from .eig";
            else
            {
                spectral_bounds(supercell_matrix(args.cellDim, model.hl), emin, emax);
                d.provenance += "; bounds from Lanczos";
            }
            d.has_bounds = true; d.a = emin; d.b = emax;
            write_descriptor(d, prefix + ".HAM.desc");
            log.at(Logger::INFO) << "  spectral bounds [a,b] = [" << emin << ", " << emax << "]";
        }
    }

    if (!args.operators.empty())
    {
        ScopedTimer t(log, stats, "write operators");
        for (const auto& op : args.operators)
        {
            hopping_list h;
            if (!build_operator(model, op, h))
            {
                log.error("operator '" + op + "' is not supported by this build");
                return 1;
            }
            log.info("writing operator " + op + " -> " + prefix + "." + op + ".CSR");
            save_supercell_as_csr(args.cellDim, h, prefix + "." + op + ".CSR");
            if (args.emit_descriptor)
                write_descriptor(op_descriptor(op), prefix + "." + op + ".desc");
            run_checks(args, op, h, (op.size()==2 && op[0]=='S') ? K_SPIN : K_GENERIC);
            const OperatorDescriptor od = op_descriptor(op);
            rep.add_output(op, prefix + "." + op + ".CSR", od.observable, od.component, od.units, od.provenance, gen_src);
        }
    }

    // External operators ingested verbatim from _hr.dat-format files and
    // expanded through the same engine (Plan 2).
    if (!args.op_files.empty())
    {
        ScopedTimer t(log, stats, "write op-files");
        for (const auto& nf : args.op_files)
        {
            log.info("ingesting operator " + nf.first + " from " + nf.second
                     + " -> " + prefix + "." + nf.first + ".CSR");
            hopping_list h = model.readOperatorModel(nf.second);
            guard_minimum_image(h, args.cellDim);           // external op may have its own range
            save_supercell_as_csr(args.cellDim, h, prefix + "." + nf.first + ".CSR");
            if (args.emit_descriptor)
                write_descriptor(describe("external", nf.first, "", "ingested from " + nf.second),
                                 prefix + "." + nf.first + ".desc");
            rep.add_output(nf.first, prefix + "." + nf.first + ".CSR", "external", nf.first, "",
                           "ingested from " + nf.second, std::vector<string>{nf.second});
        }
    }

    // Derived spin currents J = 1/2{V,S}, formed by sparse matrix product on the
    // expanded operators (Plan 3).
    if (!args.spin_currents.empty())
    {
        ScopedTimer t(log, stats, "write spin-currents");
        for (const auto& vs : args.spin_currents)
        {
            const int vdir = (vs.first == 'X') ? 0 : (vs.first == 'Y') ? 1 : 2;
            const char sdir = static_cast<char>(std::tolower(static_cast<unsigned char>(vs.second)));
            const string name = string("J") + vs.first + "S" + vs.second;

            log.at(Logger::INFO) << "writing spin current " << name << " = 1/2{V" << vs.first
                                 << ",S" << vs.second << "} -> " << prefix << "." << name << ".CSR";

            const SparseMatrix_t V = supercell_matrix(args.cellDim, model.createHoppingCurrents_list(vdir));
            const SparseMatrix_t S = supercell_matrix(args.cellDim, model.createHoppingSpinDensity_list(sdir));
            write_csr(anticommutator(V, S), prefix + "." + name + ".CSR");
            if (args.emit_descriptor)
                write_descriptor(describe("spin_current", string(1, vs.first) + "S" + vs.second,
                                          "eV*Angstrom*hbar/2", "anticommutator 1/2{V,S}"),
                                 prefix + "." + name + ".desc");
            rep.add_output(name, prefix + "." + name + ".CSR", "spin_current",
                           string(1, vs.first) + "S" + vs.second, "eV*Angstrom*hbar/2",
                           "anticommutator 1/2{V,S}", gen_src);
        }
    }

    // Exact spin operators via the gauge transform (Plan 7). Built on the raw
    // Wigner-Seitz R-set of H, then given the SAME wsvec correction as H (the
    // correction is operator-agnostic) so S and H stay on one grid/gauge.
    if (args.exact_spin)
    {
        ScopedTimer t(log, stats, "exact-spin");
        const string f_spn = in_prefix + ".spn";
        const string f_umat = in_prefix + "_u.mat";
        if (!file_exists(f_spn) || !file_exists(f_umat))
        {
            log.error("--exact-spin requires " + f_spn + " and " + f_umat);
            return 1;
        }
        gauge_data g = read_gauge(in_prefix);

        hopping_list rawH = create_hopping_list(read_wannier_file(f_hr));   // raw WS R-set
        if (g.num_wann != rawH.WannierBasisSize())
        {
            log.at(Logger::ERROR) << "num_wann mismatch between gauge (" << g.num_wann
                                  << ") and Hamiltonian (" << rawH.WannierBasisSize() << ")";
            return 1;
        }
        std::set<hopping_list::cellID_t> Rs;
        for (const auto& h : rawH.hoppings) Rs.insert(get<0>(h));
        const std::vector<hopping_list::cellID_t> Rset(Rs.begin(), Rs.end());

        const char* names[3] = {"SXexact", "SYexact", "SZexact"};
        for (int alpha = 0; alpha < 3; ++alpha)
        {
            hopping_list s = exact_spin_operator(g, alpha, Rset);
            if (file_exists(f_ws)) s = apply_wsvec(s, read_wsvec(f_ws));  // same correction as H
            log.info(string("writing exact spin ") + names[alpha] + " -> " + prefix + "." + names[alpha] + ".CSR");
            save_supercell_as_csr(args.cellDim, s, prefix + "." + names[alpha] + ".CSR");
            run_checks(args, names[alpha], s, K_SPIN);
            if (args.emit_descriptor)
                write_descriptor(describe("spin", string(1, "XYZ"[alpha]), "hbar/2",
                                          "exact gauge transform V^dag S_B V from .spn"),
                                 prefix + "." + names[alpha] + ".desc");
            rep.add_output(names[alpha], prefix + "." + names[alpha] + ".CSR", "spin",
                           string(1, "XYZ"[alpha]), "hbar/2", "exact gauge transform V^dag S_B V from .spn",
                           with_ws({f_spn, f_umat, f_hr}));
        }
    }

    // Orbital angular momentum via the projector route (Plan 8). Pure p/d shells
    // only; hybrids (e.g. Fe's sp3d2) raise a clear error.
    if (args.orbital_l)
    {
        ScopedTimer t(log, stats, "orbital-L");
        const string f_amn = in_prefix + ".amn";
        const string f_win = in_prefix + ".win";
        const string f_umat = in_prefix + "_u.mat";
        if (!file_exists(f_amn) || !file_exists(f_umat) || !file_exists(f_win))
        {
            log.error("--orbital-L requires " + f_amn + ", " + f_umat + " and " + f_win);
            return 1;
        }
        try
        {
            const std::vector<int> shells = parse_projection_shells(f_win);
            gauge_data g = read_gauge(in_prefix);
            amn_data   A = read_amn(f_amn);

            hopping_list rawH = create_hopping_list(read_wannier_file(f_hr));
            std::set<hopping_list::cellID_t> Rs;
            for (const auto& h : rawH.hoppings) Rs.insert(get<0>(h));
            const std::vector<hopping_list::cellID_t> Rset(Rs.begin(), Rs.end());

            const char* names[3] = {"LX", "LY", "LZ"};
            for (int alpha = 0; alpha < 3; ++alpha)
            {
                hopping_list L = orbital_L_operator(g, A, shells, alpha, Rset);
                if (file_exists(f_ws)) L = apply_wsvec(L, read_wsvec(f_ws));
                log.info(string("writing orbital L ") + names[alpha] + " -> " + prefix + "." + names[alpha] + ".CSR");
                save_supercell_as_csr(args.cellDim, L, prefix + "." + names[alpha] + ".CSR");
                run_checks(args, names[alpha], L, K_ORBITAL);
                if (args.emit_descriptor)
                    write_descriptor(describe("orbital_L", string(1, "XYZ"[alpha]), "hbar",
                                              "projector route C=A^dag V, C^dag L_local C"),
                                     prefix + "." + names[alpha] + ".CSR.desc");
                rep.add_output(names[alpha], prefix + "." + names[alpha] + ".CSR", "orbital_L",
                               string(1, "XYZ"[alpha]), "hbar", "projector route C=A^dag V, C^dag L_local C",
                               with_ws({f_amn, f_umat, f_win, f_hr}));
            }
        }
        catch (const std::exception& e)
        {
            log.error(e.what());
            return 1;
        }
    }

    log.info("supercells created successfully");

    } catch (const std::exception& e) {
        log.error(e.what());
        return 1;
    }
    return 0;
}

/**
 * @brief Dispatch to the requested output mode.
 * @param args parsed + config-applied arguments
 * @param log  run logger
 * @param stats per-stage timing collector
 * @return exit code (0 on success)
 */
static int run_app(const W2SP_arguments& args, Logger& log, RunStats& stats, RunReport& rep)
{
    if (args.mode == "bundle")
        return run_bundle_mode(args, log, stats, rep);
    return run_sparse_mode(args, log, stats, rep);
}

/**
 * @brief Main entry point: parse CLI / input file, set up logging + timing, run,
 *        and always emit the end-of-run summary (timings, peak memory, warn/err).
 * @param argc argument count
 * @param argv argument vector
 * @return exit code (non-zero on any error)
 */
int main(int argc, char* argv[])
{
    W2SP_arguments args;
    switch (args.parse(argc, argv))
    {
        case W2SP_arguments::EXIT_OK:    return 0;  // help / version / list-operators
        case W2SP_arguments::EXIT_ERROR: return 1;  // bad arguments
        case W2SP_arguments::PROCEED:    break;
    }

    // --create-template [NAME]: write a fully-annotated .w2s template and exit. The
    // scaffold is itself a valid input file (the reader tolerates its // comments).
    if (args.create_template)
    {
        const string stem = args.create_inp.empty() ? string("template") : args.create_inp;
        const string path = stem + ".w2s";
        std::ofstream f(path.c_str());
        if (!f.good()) { cerr << args.program_name << ": error: cannot write template '" << path << "'\n"; return 1; }
        write_run_config_template(f);
        cout << "Wrote template -> " << path << "\n";
        cout << "Edit it, then run:  " << args.program_name << " " << path << "\n";
        return 0;
    }

    // --create "LABEL N1 N2 N3 [OP...] [options]" [-inp NAME]: parse the quoted
    // command line and serialize the equivalent .w2s, instead of running. The target
    // stem is -inp NAME, else the LABEL from the command line, else "template".
    if (args.create_from_string)
    {
        std::vector<string> toks;
        { std::istringstream iss(args.create_string); string t; while (iss >> t) toks.push_back(t); }
        std::vector<char*> cargv;
        cargv.push_back(const_cast<char*>(args.program_name.c_str()));
        for (auto& t : toks) cargv.push_back(const_cast<char*>(t.c_str()));

        W2SP_arguments sub;
        if (sub.parse((int)cargv.size(), cargv.data()) != W2SP_arguments::PROCEED)
        {
            cerr << args.program_name << ": error: could not parse the --create command line\n";
            return 1;
        }
        const string stem = !args.create_inp.empty() ? args.create_inp
                          : (!sub.label.empty() ? sub.label : string("template"));
        const string path = stem + ".w2s";
        std::ofstream f(path.c_str());
        if (!f.good()) { cerr << args.program_name << ": error: cannot write input file '" << path << "'\n"; return 1; }
        write_run_config(sub, f);
        cout << "Wrote input file -> " << path << "\n";
        cout << "Edit it, then run:  " << args.program_name << " " << path << "\n";
        return 0;
    }

    // Input file (--input / --config): a JSON file supplies label / mode /
    // supercell / operators / output / provenance instead of the positional CLI.
    // It overrides only the keys it sets; the positional CLI keeps working.
    if (!args.config_path.empty())
    {
        RunConfig cfg;
        try { cfg = read_run_config(args.config_path); }
        catch (const std::exception& e)
        {
            cerr << args.program_name << ": error: " << e.what() << "\n";
            return 1;
        }
        cfg.apply_to(args);
        if (!cfg.has_mode) args.mode = "bundle";        // --run/--input/--config defaults to bundle
        if (args.label.empty())
        {
            cerr << args.program_name << ": error: input file must set a non-empty \"label\"\n";
            return 1;
        }
    }

    // --write: scaffold the equivalent JSON input file from the CLI and exit,
    // instead of running. The positional CLI is thus a generator for the input
    // file (edit it, then run it with --run), not an execution path of its own.
    if (args.write_config)
    {
        make_dir(args.output_dir);
        const string path = (args.output_dir.empty() ? string(".") : args.output_dir)
                            + "/" + args.label + ".w2s";
        std::ofstream f(path.c_str());
        if (!f.good())
        {
            cerr << args.program_name << ": error: cannot write input file '" << path << "'\n";
            return 1;
        }
        write_run_config(args, f);
        cout << "Wrote input file -> " << path << "\n";
        cout << "Edit it, then run:  " << args.program_name << " " << path << "\n";
        return 0;
    }

    // Logger: console at the requested verbosity, plus a complete log file (unless
    // disabled) located beside the outputs for traceability.
    Logger log(level_from_string(args.log_level));
    if (!args.no_log_file)
    {
        string logpath = args.log_file;
        if (logpath.empty())
        {
            make_dir(args.output_dir);
            logpath = (args.output_dir.empty() ? string(".") : args.output_dir)
                      + "/" + args.label + ".run.log";
        }
        if (log.open_file(logpath)) log.info("logging to " + logpath);
        else                        log.warn("could not open log file '" + logpath + "'; console only");
    }

    // The .out run receipt: invocation + per-output provenance ledger, filled by the
    // run modes and serialized below. Scalar facts are captured before the run.
    RunReport rep;
    rep.label        = args.label;
    rep.mode         = args.mode;
    rep.input_file   = args.config_path;                 // "" for the positional CLI
    rep.output_dir   = args.output_dir;
    rep.command_line = join_command_line(argc, argv);
    rep.started_at   = now_timestamp();
    rep.version      = W2SP_VERSION;
    rep.manual       = args.manual;

    RunStats stats;
    const auto t_start = std::chrono::steady_clock::now();
    int rc = 0;
    try
    {
        rc = run_app(args, log, stats, rep);
    }
    catch (const std::exception& e)
    {
        log.error(string("fatal: ") + e.what());
        rc = 1;
    }
    const double total_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    report_run(log, stats, total_s);

    if (log.errors() > 0 && rc == 0) rc = 1;   // surface any logged error as a failure
    rep.success = (rc == 0);

    // Write the .out receipt next to the outputs. It is a run deliverable (the
    // machine-readable provenance record), distinct from the human .run.log, so it
    // is written even when --no-log-file silenced the log file.
    {
        make_dir(args.output_dir);
        const string outpath = (args.output_dir.empty() ? string(".") : args.output_dir)
                             + "/" + args.label + ".out";
        std::ofstream of(outpath.c_str());
        if (of.good())
        {
            write_run_report(of, rep, stats, log.warnings(), log.errors(), total_s);
            log.info("run receipt: " + outpath);
        }
        else
            log.warn("could not write run receipt '" + outpath + "'");
    }

    return rc;
}
