#include <string>
#include <fstream>
#include <iostream>
#include <cctype>

#include "w2sp_arguments.hpp"
#include "tbmodel.hpp"
#include "hopping_list.hpp"
#include "operator_algebra.hpp"
#include "descriptor.hpp"
#include "gauge.hpp"
#include <set>

using namespace std;

static bool file_exists(const string& path)
{
    ifstream f(path.c_str());
    return f.good();
}

static OperatorDescriptor describe(const string& observable, const string& component,
                                   const string& units, const string& provenance)
{
    OperatorDescriptor d;
    d.observable = observable; d.component = component;
    d.units = units;           d.provenance = provenance;
    return d;
}

// Metadata for a generated operator code (VX, SZ, VXSZ, ...).
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

// Build the requested operator as a (primitive-cell) hopping_list, dispatching
// the operator name to the matching tbmodel builder.  Returns false for a name
// the parser accepted but this build cannot produce (should not happen, since
// the parser validates against the same operator set).
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

int main(int argc, char* argv[])
{
    W2SP_arguments args;
    switch (args.parse(argc, argv))
    {
        case W2SP_arguments::EXIT_OK:    return 0;  // help / version / list-operators
        case W2SP_arguments::EXIT_ERROR: return 1;  // bad arguments
        case W2SP_arguments::PROCEED:    break;
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

    cout << "Writing Hamiltonian -> " << prefix << ".HAM.CSR\n";
    save_supercell_as_csr(args.cellDim, model.hl, prefix + ".HAM.CSR");

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

    for (const auto& op : args.operators)
    {
        hopping_list h;
        if (!build_operator(model, op, h))
        {
            cerr << args.program_name << ": error: operator '" << op << "' is not supported by this build\n";
            return 1;
        }
        cout << "Writing operator " << op << " -> " << prefix << "." << op << ".CSR\n";
        save_supercell_as_csr(args.cellDim, h, prefix + "." + op + ".CSR");
        if (args.emit_descriptor)
            write_descriptor(op_descriptor(op), prefix + "." + op + ".desc");
    }

    // External operators ingested verbatim from _hr.dat-format files and
    // expanded through the same engine (Plan 2).
    for (const auto& nf : args.op_files)
    {
        cout << "Ingesting operator " << nf.first << " from " << nf.second
             << " -> " << prefix << "." << nf.first << ".CSR\n";
        hopping_list h = model.readOperatorModel(nf.second);
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

        const SparseMatrix_t V = supercell_matrix(args.cellDim, model.createHoppingCurrents_list(vdir));
        const SparseMatrix_t S = supercell_matrix(args.cellDim, model.createHoppingSpinDensity_list(sdir));
        write_csr(anticommutator(V, S), prefix + "." + name + ".CSR");
        if (args.emit_descriptor)
            write_descriptor(describe("spin_current", string(1, vs.first) + "S" + vs.second,
                                      "eV*Angstrom*hbar/2", "anticommutator 1/2{V,S}"),
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
            if (args.emit_descriptor)
                write_descriptor(describe("spin", string(1, "XYZ"[alpha]), "hbar/2",
                                          "exact gauge transform V^dag S_B V from .spn"),
                                 prefix + "." + names[alpha] + ".desc");
        }
    }

    cout << "Supercells created successfully\n";
    return 0;
}
