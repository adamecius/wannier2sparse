#include <string>
#include <fstream>
#include <iostream>
#include <cctype>

#include "w2sp_arguments.hpp"
#include "tbmodel.hpp"
#include "hopping_list.hpp"

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

    bool missing = false;
    for (const string& f : {f_uc, f_xyz, f_hr})
        if (!file_exists(f))
        {
            cerr << args.program_name << ": error: required input file '" << f << "' not found\n";
            missing = true;
        }
    if (missing) return 1;

    cout << "Using " << args.label << " as the system's identification label\n";

    tbmodel model;
    model.readOrbitalPositions(f_xyz);
    model.readUnitCell(f_uc);
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

    cout << "Supercells created successfully\n";
    return 0;
}
