#include <string>
#include <fstream>
#include <iostream>
#include <cctype>

#include "w2sp_arguments.hpp"
#include "tbmodel.hpp"
#include "hopping_list.hpp"
#include "operator_algebra.hpp"

using namespace std;

static bool file_exists(const string& path)
{
    ifstream f(path.c_str());
    return f.good();
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

    const string f_uc  = args.label + ".uc";
    const string f_xyz = args.label + ".xyz";
    const string f_hr  = args.label + "_hr.dat";

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

    const string prefix = args.output_dir + "/" + args.label;

    cout << "Creating the supercell (" << args.cellDim[0] << ","
         << args.cellDim[1] << "," << args.cellDim[2] << ")\n";

    cout << "Writing Hamiltonian -> " << prefix << ".HAM.CSR\n";
    save_supercell_as_csr(args.cellDim, model.hl, prefix + ".HAM.CSR");

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
    }

    // External operators ingested verbatim from _hr.dat-format files and
    // expanded through the same engine (Plan 2).
    for (const auto& nf : args.op_files)
    {
        cout << "Ingesting operator " << nf.first << " from " << nf.second
             << " -> " << prefix << "." << nf.first << ".CSR\n";
        hopping_list h = model.readOperatorModel(nf.second);
        save_supercell_as_csr(args.cellDim, h, prefix + "." + nf.first + ".CSR");
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
    }

    cout << "Supercells created successfully\n";
    return 0;
}
