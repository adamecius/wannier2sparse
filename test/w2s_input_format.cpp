// .w2s input-format test. The `.w2s` file is the single, validated entry point:
// JSON that tolerates // and /* */ comments, carries an optional user-declared
// `provenance.manual` block, and round-trips through write_run_config /
// write_run_config_template. This test exercises those three contracts directly on
// the reader/writer (no model load), so a regression in the input surface is caught
// without needing fixtures.
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include "run_config.hpp"
#include "w2sp_arguments.hpp"

using namespace std;

static string slurp(const string& p) { ifstream f(p.c_str()); stringstream ss; ss << f.rdbuf(); return ss.str(); }
static bool has(const string& hay, const string& needle) { return hay.find(needle) != string::npos; }

int main()
{
    // 1) Comments (// line and /* block */) are tolerated, and the manual
    //    provenance block parses into the right fields.
    ofstream("commented.w2s") <<
        "{\n"
        "  // a line comment\n"
        "  \"label\": \"gr\", /* inline block */ \"mode\": \"sparse\",\n"
        "  \"supercell\": [2, 2, 1],          // dims\n"
        "  \"operators\": [\"VX\", \"SZ\"],\n"
        "  \"provenance\": {\n"
        "    \"manual\": {\n"
        "      \"code\": \"QE 7.2\",\n"
        "      \"xc_functional\": \"PBE\",\n"
        "      \"basis\": \"pz\",\n"
        "      \"ecutwfc_Ry\": 90,\n"
        "      \"kpoint_grid\": [12, 12, 4],\n"
        "      \"pseudopotentials\": [{ \"species\": \"Fe\", \"file\": \"Fe.UPF\", \"z_valence\": 16 }],\n"
        "      \"orbitals\": [{ \"index\": 0, \"label\": \"Fe-d\", \"element\": \"Fe\", \"xyz\": [0.0, 0.5, 1.0] }]\n"
        "    }\n"
        "  }\n"
        "}\n";

    RunConfig cfg = read_run_config("commented.w2s");
    assert(cfg.has_label && cfg.label == "gr");
    assert(cfg.has_supercell && cfg.supercell[0] == 2);
    assert(cfg.has_operators && cfg.operators.size() == 2);
    assert(cfg.has_manual && cfg.manual.present);
    assert(cfg.manual.code == "QE 7.2");
    assert(cfg.manual.xc_functional == "PBE");
    assert(cfg.manual.basis == "pz");
    assert(cfg.manual.has_ecutwfc && cfg.manual.ecutwfc_Ry == 90.0);
    assert(cfg.manual.has_kpoint_grid && cfg.manual.kpoint_grid[0] == 12 && cfg.manual.kpoint_grid[2] == 4);
    assert(cfg.manual.pseudopotentials.size() == 1);
    assert(cfg.manual.pseudopotentials[0].species == "Fe" && cfg.manual.pseudopotentials[0].file == "Fe.UPF");
    assert(cfg.manual.pseudopotentials[0].has_z_valence && cfg.manual.pseudopotentials[0].z_valence == 16.0);
    assert(cfg.manual.orbitals.size() == 1);
    assert(cfg.manual.orbitals[0].index == 0 && cfg.manual.orbitals[0].label == "Fe-d");
    assert(cfg.manual.orbitals[0].xyz[1] == 0.5 && cfg.manual.orbitals[0].xyz[2] == 1.0);

    // apply_to carries the manual block onto the resolved arguments.
    W2SP_arguments a; cfg.apply_to(a);
    assert(a.manual.present && a.manual.code == "QE 7.2" && a.manual.orbitals.size() == 1);

    // 1b) velocity_mode round-trips and an invalid value is rejected.
    {
        ofstream("vmode.w2s") << "{ \"label\":\"x\", \"velocity_mode\":\"covariant\", \"r_dat\":\"x_r.dat\" }\n";
        RunConfig vc = read_run_config("vmode.w2s");
        assert(vc.has_velocity_mode && vc.velocity_mode == "covariant");
        assert(vc.has_r_dat && vc.r_dat == "x_r.dat");
        W2SP_arguments va; vc.apply_to(va);
        assert(va.velocity_mode == "covariant" && va.covariant_velocity && va.r_dat_path == "x_r.dat");

        ofstream("vbad.w2s") << "{ \"label\":\"x\", \"velocity_mode\":\"nope\" }\n";
        bool threw = false;
        try { read_run_config("vbad.w2s"); } catch (const std::exception&) { threw = true; }
        assert(threw && "invalid velocity_mode must throw");
    }

    // 2) An unterminated block comment is a clear error, not silent acceptance.
    {
        ofstream("badcomment.w2s") << "{ \"label\": \"x\" /* never closed }\n";
        bool threw = false;
        try { read_run_config("badcomment.w2s"); } catch (const std::exception&) { threw = true; }
        assert(threw && "unterminated /* */ must throw");
    }

    // 3) The --create-template scaffold is itself a valid .w2s: it parses back,
    //    with its manual block present (the reader skips its documentation comments).
    {
        ofstream tf("template.w2s");
        write_run_config_template(tf);
        tf.close();
        RunConfig t = read_run_config("template.w2s");
        assert(t.has_label && t.has_mode && t.has_operators);
        assert(t.has_manual && t.manual.present);
        assert(t.manual.has_kpoint_grid && t.manual.pseudopotentials.size() == 1);
    }

    // 4) write_run_config (--create / --write) round-trips the manual block.
    {
        W2SP_arguments src;
        src.label = "gr"; src.mode = "sparse"; src.output_dir = ".";
        src.cellDim = {{3, 3, 1}};
        src.operators.push_back("VX");
        src.manual.present = true;
        src.manual.code = "manual code";
        src.manual.has_kpoint_grid = true; src.manual.kpoint_grid = {{6, 6, 6}};
        ManualOrbital o; o.index = 2; o.label = "pz"; o.element = "C"; o.xyz = {{1.0, 2.0, 3.0}};
        src.manual.orbitals.push_back(o);

        { ofstream f("written.w2s"); write_run_config(src, f); }
        const string body = slurp("written.w2s");
        assert(has(body, "\"manual\"") && has(body, "manual code") && has(body, "\"kpoint_grid\""));

        RunConfig back = read_run_config("written.w2s");
        assert(back.has_manual && back.manual.present);
        assert(back.manual.code == "manual code");
        assert(back.manual.has_kpoint_grid && back.manual.kpoint_grid[0] == 6);
        assert(back.manual.orbitals.size() == 1 && back.manual.orbitals[0].index == 2);
        assert(back.manual.orbitals[0].xyz[2] == 3.0);
    }

    cout << "W2S INPUT FORMAT TEST PASSED" << endl;
    return 0;
}
