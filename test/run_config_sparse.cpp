// One input file drives the sparse (supercell CSR) mode too, not just the bundle.
// This test verifies the generalized run.json schema parses the sparse-mode keys
// (supercell, spin_currents, op_files, checks, log settings) into a RunConfig with
// the right presence flags and values, and that apply_to() projects them onto a
// W2SP_arguments so the existing sparse path is driven exactly as the positional
// CLI would be.
#include <cassert>
#include <fstream>
#include <string>
#include <iostream>
#include "run_config.hpp"
#include "w2sp_arguments.hpp"

using namespace std;

int main()
{
    ofstream("sparse.json") <<
        "{\n"
        "  \"label\": \"graphene\", \"seed\": \"gr\", \"mode\": \"sparse\",\n"
        "  \"output_dir\": \"out\",\n"
        "  \"supercell\": [80, 80, 1],\n"
        "  \"operators\": [\"VX\", \"SZ\"],\n"
        "  \"spin_currents\": [[\"X\", \"Z\"], [\"Y\", \"X\"]],\n"
        "  \"op_files\": [{\"name\": \"OX\", \"path\": \"ox_hr.dat\"}],\n"
        "  \"checks\": \"hermiticity\", \"emit_bounds\": true,\n"
        "  \"log_level\": \"debug\", \"log_file\": \"out/run.log\"\n"
        "}\n";

    RunConfig cfg = read_run_config("sparse.json");
    assert(cfg.has_mode && cfg.mode == "sparse");
    assert(cfg.has_supercell && cfg.supercell[0] == 80 && cfg.supercell[1] == 80 && cfg.supercell[2] == 1);
    assert(cfg.has_operators && cfg.operators.size() == 2);
    assert(cfg.has_spin_currents && cfg.spin_currents.size() == 2);
    assert(cfg.spin_currents[0].first == 'X' && cfg.spin_currents[0].second == 'Z');
    assert(cfg.spin_currents[1].first == 'Y' && cfg.spin_currents[1].second == 'X');
    assert(cfg.has_op_files && cfg.op_files.size() == 1);
    assert(cfg.op_files[0].first == "OX" && cfg.op_files[0].second == "ox_hr.dat");
    assert(cfg.has_checks && cfg.checks == "hermiticity");
    assert(cfg.has_emit_bounds && cfg.emit_bounds);
    assert(cfg.has_log_level && cfg.log_level == "debug");
    assert(cfg.has_log_file && cfg.log_file == "out/run.log");

    // apply_to projects every set key onto the arguments.
    W2SP_arguments a;
    cfg.apply_to(a);
    assert(a.mode == "sparse");
    assert(a.cellDim[0] == 80 && a.cellDim[1] == 80 && a.cellDim[2] == 1);
    assert(a.operators.size() == 2);
    assert(a.spin_currents.size() == 2 && a.spin_currents[0].first == 'X' && a.spin_currents[0].second == 'Z');
    assert(a.op_files.size() == 1 && a.op_files[0].first == "OX");
    assert(a.check == "hermiticity");
    assert(a.emit_descriptor == true);
    assert(a.log_level == "debug");
    assert(a.log_file == "out/run.log");
    assert(a.input_prefix() == "gr");   // seed overrides label for inputs

    // A bad supercell (dimension < 1) is rejected.
    ofstream("badsc.json") << "{ \"label\":\"x\", \"supercell\":[0,1,1] }\n";
    bool threw = false;
    try { read_run_config("badsc.json"); } catch (const std::exception&) { threw = true; }
    assert(threw && "supercell dimension < 1 must throw");

    // write_run_config (--write) round-trips: CLI args -> JSON -> read -> args.
    {
        W2SP_arguments src;
        src.label = "gr"; src.seed = "seedname"; src.mode = "sparse"; src.output_dir = "o";
        src.cellDim = {{4, 5, 6}};
        src.operators.push_back("VX"); src.operators.push_back("SZ");
        src.spin_currents.push_back(make_pair('X', 'Z'));
        src.op_files.push_back(make_pair(string("OX"), string("ox_hr.dat")));
        src.emit_descriptor = true; src.check = "all"; src.exact_spin = true;
        src.qe_xml_path = "scf.xml"; src.win_path = "gr.win";

        { ofstream f("written.json"); write_run_config(src, f); }

        RunConfig back = read_run_config("written.json");
        W2SP_arguments dst;
        back.apply_to(dst);
        assert(dst.label == "gr" && dst.seed == "seedname" && dst.mode == "sparse");
        assert(dst.output_dir == "o");
        assert(dst.cellDim[0] == 4 && dst.cellDim[1] == 5 && dst.cellDim[2] == 6);
        assert(dst.operators.size() == 2 && dst.operators[0] == "VX" && dst.operators[1] == "SZ");
        assert(dst.spin_currents.size() == 1 && dst.spin_currents[0].first == 'X' && dst.spin_currents[0].second == 'Z');
        assert(dst.op_files.size() == 1 && dst.op_files[0].first == "OX" && dst.op_files[0].second == "ox_hr.dat");
        assert(dst.emit_descriptor && dst.check == "all" && dst.exact_spin);
        assert(dst.qe_xml_path == "scf.xml" && dst.win_path == "gr.win");
    }

    cout << "RUN CONFIG SPARSE TEST PASSED" << endl;
    return 0;
}
