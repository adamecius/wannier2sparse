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
#include <cctype>
#include <functional>
#include <memory>
#include <array>

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
#include "input_file.hpp"
#include "win_parser.hpp"
#include "qe_xml_parser.hpp"
#include <set>
#include <stdexcept>

using namespace std;

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
/**
 * @brief Build the velocity builder for the active velocity_mode (the velocity ladder).
 * @param model loaded tight-binding model
 * @param args  parsed arguments (velocity_mode, optional r_dat path)
 * @param in_prefix resolved input stem (for the default <seed>_r.dat)
 * @return functor dir -> velocity hopping_list
 *
 * bare/berry_connection need only H + .xyz; covariant additionally reads the
 * position matrix `_r.dat` once and forms v = -i(R.lat)H - i[H,A]. Throws if the
 * covariant mode is requested without a readable `_r.dat`.
 */
static std::function<hopping_list(int)>
make_velocity_builder(tbmodel& model, const W2SP_arguments& args, const string& in_prefix)
{
    if (args.velocity_mode == "covariant")
    {
        const string rpath = args.r_dat_path.empty() ? (in_prefix + "_r.dat") : args.r_dat_path;
        if (!file_exists(rpath))
            throw std::runtime_error("velocity_mode=covariant needs the position matrix '" + rpath +
                "' (Wannier90 'write_rmn' _r.dat). Pass --r-dat PATH, set r_dat in the input file, "
                "or use velocity_mode=berry_connection.");
        auto A = std::make_shared<std::array<hopping_list,3> >(model.readPositionMatrix(rpath));
        cout << "Velocity mode: covariant (Berry connection from " << rpath << ")\n";
        return [&model, A](int dir){ return model.createCovariantVelocity_list(dir, (*A)[dir]); };
    }
    const bool bare = (args.velocity_mode == "bare");
    cout << "Velocity mode: " << args.velocity_mode << "\n";
    return [&model, bare](int dir){ return model.createHoppingCurrents_list(dir, bare); };
}

static bool build_operator(tbmodel& model, const string& op, hopping_list& out,
                           const std::function<hopping_list(int)>& make_velocity)
{
    if (op == "VX") { out = make_velocity(0); return true; }
    if (op == "VY") { out = make_velocity(1); return true; }
    if (op == "VZ") { out = make_velocity(2); return true; }

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
static int run_bundle_mode(const W2SP_arguments& args)
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
        cerr << args.program_name << ": error: required input file '" << f_hr << "' not found\n";
        missing = true;
    }
    if (need_positions)
        for (const string& f : {f_uc, f_xyz})
            if (!file_exists(f))
            {
                cerr << args.program_name << ": error: operator generation requires '" << f << "'\n";
                missing = true;
            }
    for (const auto& nf : args.op_files)
        if (!file_exists(nf.second))
        {
            cerr << args.program_name << ": error: --op-file '" << nf.first
                 << "': file '" << nf.second << "' not found\n";
            missing = true;
        }
    if (missing) return 1;

    cout << "Using " << args.label << " as the system's identification label (bundle mode)\n";

    tbmodel model;
    if (file_exists(f_xyz)) model.readOrbitalPositions(f_xyz);
    if (file_exists(f_uc))  model.readUnitCell(f_uc);
    model.readWannierModel(f_hr);

    const bool wsvec_applied = file_exists(f_ws);
    if (wsvec_applied)
    {
        cout << "Applying Wigner-Seitz correction from " << f_ws << "\n";
        model.applyWsvec(f_ws);
    }

    SystemProvenance prov;
    prov.num_wann = model.hl.WannierBasisSize();
    if (file_exists(f_uc)) { prov.has_lattice = true; prov.lattice = model.lat_vecs; }
    if (file_exists(f_xyz)) prov.wannier_sites = collect_wannier_sites(model);

    // Wannier90 provenance from the .win file (config-specified path, else
    // <prefix>.win). Absent or unparsable => the manifest block stays null.
    const string f_win_prov = !args.win_path.empty() ? args.win_path : (in_prefix + ".win");
    if (file_exists(f_win_prov))
    {
        try { prov.wann = parse_win_provenance(f_win_prov); }
        catch (const std::exception& e)
        {
            cerr << args.program_name << ": warning: could not parse .win provenance ("
                 << f_win_prov << "): " << e.what() << "\n";
        }
    }

    // DFT (Quantum ESPRESSO) provenance from data-file-schema.xml. Only parsed
    // when a path was given (no default location); absent/unparsable => null.
    if (!args.qe_xml_path.empty())
    {
        if (file_exists(args.qe_xml_path))
        {
            try { prov.dft = parse_qe_xml(args.qe_xml_path); }
            catch (const std::exception& e)
            {
                cerr << args.program_name << ": warning: could not parse QE XML provenance ("
                     << args.qe_xml_path << "): " << e.what() << "\n";
            }
        }
        else
            cerr << args.program_name << ": warning: QE XML provenance file '"
                 << args.qe_xml_path << "' not found\n";
    }

    std::vector<BundleOperator> ops;

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
                cerr << args.program_name << ": note: bundle mode skips Lanczos bounds "
                        "(they need k-sampling); provide " << in_prefix << ".eig for HAM bounds\n";
        }
        op.hl = model.hl;
        ops.push_back(op);
    }

    // Requested generated operators (velocity / spin / spin-current).
    std::function<hopping_list(int)> make_velocity;
    try { make_velocity = make_velocity_builder(model, args, in_prefix); }
    catch (const std::exception& e) { cerr << args.program_name << ": error: " << e.what() << "\n"; return 1; }
    for (const auto& opname : args.operators)
    {
        hopping_list h;
        if (!build_operator(model, opname, h, make_velocity))
        {
            cerr << args.program_name << ": error: operator '" << opname << "' is not supported by this build\n";
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
                cerr << args.program_name << ": error: --exact-spin requires " << f_spn
                     << " and " << f_umat << "\n";
                return 1;
            }
            gauge_data g = read_gauge(in_prefix);
            hopping_list rawH = create_hopping_list(read_wannier_file(f_hr));
            if (g.num_wann != rawH.WannierBasisSize())
            {
                cerr << args.program_name << ": error: num_wann mismatch between gauge ("
                     << g.num_wann << ") and Hamiltonian (" << rawH.WannierBasisSize() << ")\n";
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
                cerr << args.program_name << ": error: --orbital-L requires " << f_amn
                     << ", " << f_umat << " and " << f_win << "\n";
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
        cerr << args.program_name << ": error: " << e.what() << "\n";
        return 1;
    }

    BundleSpec spec;
    spec.label = args.label;
    spec.truncation_threshold = args.has_truncation
        ? args.truncation_threshold : std::numeric_limits<double>::epsilon();
    spec.ndegen_applied = true;
    spec.wsvec_applied = wsvec_applied;
    if (wsvec_applied) spec.wsvec_src = f_ws;   // copied verbatim into the bundle

    const string bundle_dir = write_bundle(spec, prov, ops, args.output_dir);
    cout << "Wrote bundle with " << ops.size() << " operator(s) -> " << bundle_dir << "\n";
    cout << "  manifest: " << bundle_dir << "/manifest.json\n";
    return 0;
}

/**
 * @brief Main entry point.
 * @param argc argument count
 * @param argv argument vector
 * @return exit code
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

    // ---- input-file CLI: --create (template) / --write (set option) / --run -----
    const std::string default_inp = "w2sp.inp";
    if (args.do_create)
    {
        const std::string p = args.inp_path.empty() ? default_inp : args.inp_path;
        try { w2sp::input_file_create(p); }
        catch (const std::exception& e) { cerr << args.program_name << ": error: " << e.what() << "\n"; return 1; }
        cout << "Wrote template input file -> " << p << "\n  edit it, or: " << args.program_name
             << " --write KEY [VALUE...] -inp " << p << "\n  then run:    " << args.program_name << " --run " << p << "\n";
        return 0;
    }
    if (args.do_write)
    {
        const std::string p = args.inp_path.empty() ? default_inp : args.inp_path;
        if (args.write_tokens.empty()) { cerr << args.program_name << ": error: --write needs an option, e.g. --write covariant-velocity -inp " << p << "\n"; return 1; }
        const std::string t0 = args.write_tokens[0];
        std::string rest;
        for (std::size_t i = 1; i < args.write_tokens.size(); ++i) rest += (i > 1 ? " " : "") + args.write_tokens[i];
        std::string norm = t0; for (char& c : norm) if (c == '-') c = '_';   // flag names use '-' on CLI
        std::string key, value; bool append = false;
        if (norm == "covariant_velocity" || norm == "exact_spin" || norm == "bounds" ||
            norm == "orbital_L" || norm == "orbital_l")
            { key = norm; value = rest.empty() ? "true" : rest; }
        else if (W2SP_arguments::is_valid_operator(t0))      { key = "operators"; value = t0; append = true; }
        else if (norm == "spin_current")                     { key = "spin_current"; value = rest; }
        else if (norm == "op_file")                          { key = "op_file"; value = rest; }
        else if (t0.find('=') != std::string::npos)          { key = t0.substr(0, t0.find('=')); value = t0.substr(t0.find('=') + 1); if (!rest.empty()) value += " " + rest; }
        else                                                 { key = norm; value = rest; }
        try { w2sp::input_file_set(p, key, value, append); }
        catch (const std::exception& e) { cerr << args.program_name << ": error: " << e.what() << "\n"; return 1; }
        cout << "Set '" << key << " = " << value << "' in " << p << "\n";
        return 0;
    }
    if (!args.run_path.empty())
    {
        try { w2sp::input_file_apply(args.run_path, args); }
        catch (const std::exception& e) { cerr << args.program_name << ": error: " << e.what() << "\n"; return 1; }
        if (args.label.empty()) { cerr << args.program_name << ": error: input file must set a non-empty 'label'\n"; return 1; }
    }
    // Config-driven bundle runs: a run.json supplies label/operators/output/
    // provenance instead of the positional CLI. It overrides only the keys it
    // sets, then drives the same run_bundle_mode path below.
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
        if (!cfg.has_mode) args.mode = "bundle";        // --config drives bundle mode by default
        if (args.label.empty())
        {
            cerr << args.program_name << ": error: run.json must set a non-empty \"label\"\n";
            return 1;
        }
    }

    if (args.mode == "bundle")
    {
        const int rc = run_bundle_mode(args);
        if (rc == 0 && !args.run_path.empty())
            w2sp::input_file_write_output(args.output_dir + "/" + args.label + ".w2sp.out", args);
        return rc;
    }

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
        cerr << args.program_name << ": error: required input file '" << f_hr << "' not found\n";
        missing = true;
    }
    if (need_positions)
        for (const string& f : {f_uc, f_xyz})
            if (!file_exists(f))
            {
                cerr << args.program_name << ": error: operator generation requires '" << f << "'\n";
                missing = true;
            }
    for (const auto& nf : args.op_files)
        if (!file_exists(nf.second))
        {
            cerr << args.program_name << ": error: --op-file '" << nf.first
                 << "': file '" << nf.second << "' not found\n";
            missing = true;
        }
    if (missing) return 1;

    cout << "Using " << args.label << " as the system's identification label\n";

    tbmodel model;
    if (file_exists(f_xyz)) model.readOrbitalPositions(f_xyz);
    if (file_exists(f_uc))  model.readUnitCell(f_uc);
    model.readWannierModel(f_hr);

    // Optional Wigner-Seitz minimum-image correction (Plan 5): applied to every
    // operator since it acts on the shared Hamiltonian hopping list.
    const string f_ws = in_prefix + "_wsvec.dat";
    if (file_exists(f_ws))
    {
        cout << "Applying Wigner-Seitz correction from " << f_ws << "\n";
        model.applyWsvec(f_ws);
    }

    const string prefix = args.output_dir + "/" + args.label;

    cout << "Creating the supercell (" << args.cellDim[0] << ","
         << args.cellDim[1] << "," << args.cellDim[2] << ")\n";

    try {

    // Plan 10D: refuse to emit corrupted operators. Every generated operator
    // shares H's R-range (velocity/spin/current/exact-spin/orbital-L are all
    // built on the H Wigner-Seitz set), so guarding the (post-wsvec) H covers
    // them; external op-files are guarded individually below.
    guard_minimum_image(model.hl, args.cellDim);

    cout << "Writing Hamiltonian -> " << prefix << ".HAM.CSR\n";
    save_supercell_as_csr(args.cellDim, model.hl, prefix + ".HAM.CSR");
    run_checks(args, "HAM", model.hl, K_HAM);

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
        cout << "  spectral bounds [a,b] = [" << emin << ", " << emax << "]\n";
    }

    const std::function<hopping_list(int)> make_velocity = make_velocity_builder(model, args, in_prefix);

    for (const auto& op : args.operators)
    {
        hopping_list h;
        if (!build_operator(model, op, h, make_velocity))
        {
            cerr << args.program_name << ": error: operator '" << op << "' is not supported by this build\n";
            return 1;
        }
        cout << "Writing operator " << op << " -> " << prefix << "." << op << ".CSR\n";
        save_supercell_as_csr(args.cellDim, h, prefix + "." + op + ".CSR");
        if (args.emit_descriptor)
            write_descriptor(op_descriptor(op), prefix + "." + op + ".desc");
        run_checks(args, op, h, (op.size()==2 && op[0]=='S') ? K_SPIN : K_GENERIC);
    }

    // External operators ingested verbatim from _hr.dat-format files and
    // expanded through the same engine (Plan 2).
    for (const auto& nf : args.op_files)
    {
        cout << "Ingesting operator " << nf.first << " from " << nf.second
             << " -> " << prefix << "." << nf.first << ".CSR\n";
        hopping_list h = model.readOperatorModel(nf.second);
        guard_minimum_image(h, args.cellDim);           // external op may have its own range
        save_supercell_as_csr(args.cellDim, h, prefix + "." + nf.first + ".CSR");
        if (args.emit_descriptor)
            write_descriptor(describe("external", nf.first, "", "ingested from " + nf.second),
                             prefix + "." + nf.first + ".desc");
    }

    // Derived spin currents J = 1/2{V,S}, formed by sparse matrix product on the
    // expanded operators (Plan 3).
    for (const auto& vs : args.spin_currents)
    {
        const int vdir = (vs.first == 'X') ? 0 : (vs.first == 'Y') ? 1 : 2;
        const char sdir = static_cast<char>(std::tolower(static_cast<unsigned char>(vs.second)));
        const string name = string("J") + vs.first + "S" + vs.second;

        cout << "Writing spin current " << name << " = 1/2{V" << vs.first
             << ",S" << vs.second << "} -> " << prefix << "." << name << ".CSR\n";

        const SparseMatrix_t V = supercell_matrix(args.cellDim, make_velocity(vdir));
        const SparseMatrix_t S = supercell_matrix(args.cellDim, model.createHoppingSpinDensity_list(sdir));
        write_csr(anticommutator(V, S), prefix + "." + name + ".CSR");
        if (args.emit_descriptor)
            write_descriptor(describe("spin_current", string(1, vs.first) + "S" + vs.second,
                                      "eV*Angstrom*hbar/2", "anticommutator 1/2{V,S} (" + args.velocity_mode + " velocity)"),
                             prefix + "." + name + ".desc");
    }

    // Exact spin operators via the gauge transform (Plan 7). Built on the raw
    // Wigner-Seitz R-set of H, then given the SAME wsvec correction as H (the
    // correction is operator-agnostic) so S and H stay on one grid/gauge.
    if (args.exact_spin)
    {
        const string f_spn = in_prefix + ".spn";
        const string f_umat = in_prefix + "_u.mat";
        if (!file_exists(f_spn) || !file_exists(f_umat))
        {
            cerr << args.program_name << ": error: --exact-spin requires " << f_spn
                 << " and " << f_umat << "\n";
            return 1;
        }
        gauge_data g = read_gauge(in_prefix);

        hopping_list rawH = create_hopping_list(read_wannier_file(f_hr));   // raw WS R-set
        if (g.num_wann != rawH.WannierBasisSize())
        {
            cerr << args.program_name << ": error: num_wann mismatch between gauge ("
                 << g.num_wann << ") and Hamiltonian (" << rawH.WannierBasisSize() << ")\n";
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
            cout << "Writing exact spin " << names[alpha] << " -> " << prefix << "." << names[alpha] << ".CSR\n";
            save_supercell_as_csr(args.cellDim, s, prefix + "." + names[alpha] + ".CSR");
            run_checks(args, names[alpha], s, K_SPIN);
            if (args.emit_descriptor)
                write_descriptor(describe("spin", string(1, "XYZ"[alpha]), "hbar/2",
                                          "exact gauge transform V^dag S_B V from .spn"),
                                 prefix + "." + names[alpha] + ".desc");
        }
    }

    // Orbital angular momentum via the projector route (Plan 8). Pure p/d shells
    // only; hybrids (e.g. Fe's sp3d2) raise a clear error.
    if (args.orbital_l)
    {
        const string f_amn = in_prefix + ".amn";
        const string f_win = in_prefix + ".win";
        const string f_umat = in_prefix + "_u.mat";
        if (!file_exists(f_amn) || !file_exists(f_umat) || !file_exists(f_win))
        {
            cerr << args.program_name << ": error: --orbital-L requires " << f_amn
                 << ", " << f_umat << " and " << f_win << "\n";
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
                cout << "Writing orbital L " << names[alpha] << " -> " << prefix << "." << names[alpha] << ".CSR\n";
                save_supercell_as_csr(args.cellDim, L, prefix + "." + names[alpha] + ".CSR");
                run_checks(args, names[alpha], L, K_ORBITAL);
                if (args.emit_descriptor)
                    write_descriptor(describe("orbital_L", string(1, "XYZ"[alpha]), "hbar",
                                              "projector route C=A^dag V, C^dag L_local C"),
                                     prefix + "." + names[alpha] + ".CSR.desc");
            }
        }
        catch (const std::exception& e)
        {
            cerr << args.program_name << ": error: " << e.what() << "\n";
            return 1;
        }
    }

    cout << "Supercells created successfully\n";

    } catch (const std::exception& e) {
        cerr << args.program_name << ": error: " << e.what() << "\n";
        return 1;
    }

    if (!args.run_path.empty())   // --run: write the lsquant-style output summary
        w2sp::input_file_write_output(args.output_dir + "/" + args.label + ".w2sp.out", args);
    return 0;
}
