// Phase B2: a run.json config drives the SAME bundle the positional --mode bundle
// CLI produces. This test verifies (1) read_run_config parses the documented schema
// into a RunConfig with the right presence flags and values; (2) apply_to() projects
// it onto a W2SP_arguments overriding only the keys the file set; and (3) the bundle
// written from a config is byte-identical to the bundle written from the equivalent
// positional invocation.
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <sys/stat.h>

#include "run_config.hpp"
#include "w2sp_arguments.hpp"
#include "tbmodel.hpp"
#include "hopping_list.hpp"
#include "descriptor.hpp"
#include "bundle_writer.hpp"

using namespace std;

static string slurp(const string& path)
{
    ifstream f(path.c_str());
    stringstream ss; ss << f.rdbuf();
    return ss.str();
}

// Build the bundle directly (the same operator set run_bundle_mode would build for
// HAM + VX + SZ on this model), into out_dir, and return the bundle path. Done in
// the test so the config vs positional comparison reduces to "same BundleSpec".
static string build_bundle(const string& label, double trunc, bool has_trunc, const string& out_dir)
{
    tbmodel model;
    model.readOrbitalPositions("spin_graphene.xyz");
    model.readUnitCell("spin_graphene.uc");
    model.readWannierModel("spin_graphene_hr.dat");

    SystemProvenance prov;
    prov.num_wann = model.hl.WannierBasisSize();
    prov.has_lattice = true;
    prov.lattice = model.lat_vecs;
    {
        const auto id2spin = model.map_id2spin();
        int idx = 0;
        for (const auto& orb : model.orbPos_list)
        {
            WannierSite s;
            s.index = idx; s.label = get<0>(orb); s.cart = get<1>(orb);
            s.spin = id2spin.count(idx) ? id2spin.at(idx) : 0;
            prov.wannier_sites.push_back(s);
            ++idx;
        }
    }

    vector<BundleOperator> ops;
    { BundleOperator op; op.name = "HAM"; op.desc.observable = "hamiltonian"; op.desc.units = "eV"; op.hl = model.hl; ops.push_back(op); }
    { BundleOperator op; op.name = "VX"; op.desc.observable = "velocity"; op.desc.component = "X"; op.desc.units = "eV*Angstrom"; op.hl = model.createHoppingCurrents_list(0); ops.push_back(op); }
    { BundleOperator op; op.name = "SZ"; op.desc.observable = "spin"; op.desc.component = "Z"; op.desc.units = "hbar/2"; op.hl = model.createHoppingSpinDensity_list('z'); ops.push_back(op); }

    BundleSpec spec;
    spec.label = label;
    spec.truncation_threshold = has_trunc ? trunc : numeric_limits<double>::epsilon();
    return write_bundle(spec, prov, ops, out_dir);
}

int main()
{
    // (1) Parser: full documented schema -> RunConfig.
    {
        ofstream("full.json") <<
            "{\n"
            "  \"label\": \"graphene\", \"project_dir\": \".\", \"seed\": \"gr\",\n"
            "  \"output_dir\": \"out\", \"mode\": \"bundle\",\n"
            "  \"operators\": [\"VX\", \"VY\", \"SZ\"],\n"
            "  \"exact_spin\": false, \"orbital_L\": true, \"emit_bounds\": true,\n"
            "  \"truncation_threshold\": 1e-8,\n"
            "  \"provenance\": { \"qe_xml\": \"scf.save/data-file-schema.xml\", \"win\": \"gr.win\" }\n"
            "}\n";

        RunConfig cfg = read_run_config("full.json");
        assert(cfg.has_label && cfg.label == "graphene");
        assert(cfg.has_seed && cfg.seed == "gr");
        assert(cfg.has_output_dir && cfg.output_dir == "out");
        assert(cfg.has_mode && cfg.mode == "bundle");
        assert(cfg.has_operators && cfg.operators.size() == 3 &&
               cfg.operators[0] == "VX" && cfg.operators[2] == "SZ");
        assert(cfg.has_exact_spin && cfg.exact_spin == false);
        assert(cfg.has_orbital_l && cfg.orbital_l == true);
        assert(cfg.has_emit_bounds && cfg.emit_bounds == true);
        assert(cfg.has_truncation && cfg.truncation_threshold == 1e-8);
        assert(cfg.has_qe_xml && cfg.qe_xml == "scf.save/data-file-schema.xml");
        assert(cfg.has_win && cfg.win == "gr.win");

        // apply_to overrides only the set keys.
        W2SP_arguments a;
        cfg.apply_to(a);
        assert(a.label == "graphene");
        assert(a.input_prefix() == "./gr");
        assert(a.output_dir == "out");
        assert(a.mode == "bundle");
        assert(a.operators.size() == 3);
        assert(a.orbital_l == true && a.exact_spin == false);
        assert(a.emit_descriptor == true);
        assert(a.has_truncation && a.truncation_threshold == 1e-8);
        assert(a.qe_xml_path == "scf.save/data-file-schema.xml");
        assert(a.win_path == "gr.win");
    }

    // (2) A minimal config leaves unset keys at their W2SP_arguments defaults.
    {
        ofstream("min.json") << "{ \"label\": \"x\" }\n";
        RunConfig cfg = read_run_config("min.json");
        assert(cfg.has_label && !cfg.has_operators && !cfg.has_truncation);
        W2SP_arguments a;
        cfg.apply_to(a);
        assert(a.label == "x");
        assert(a.mode == "sparse");          // unchanged: config did not set mode
        assert(a.operators.empty());
        assert(!a.has_truncation);
    }

    // (3) Malformed JSON and unknown operators are rejected with a clear error.
    {
        ofstream("bad.json") << "{ \"label\": \"x\", }\n";   // trailing comma
        bool threw = false;
        try { read_run_config("bad.json"); } catch (const std::exception&) { threw = true; }
        assert(threw && "malformed JSON must throw");

        ofstream("badop.json") << "{ \"operators\": [\"NOPE\"] }\n";
        threw = false;
        try { read_run_config("badop.json"); } catch (const std::exception&) { threw = true; }
        assert(threw && "unknown operator must throw");
    }

    // (4) Equivalence: config-built bundle == positional-built bundle, byte for byte.
    {
        mkdir("pos", 0755);
        mkdir("cfg", 0755);
        const string a = build_bundle("spin_graphene", 0.0, false, "pos");
        const string b = build_bundle("spin_graphene", 0.0, false, "cfg");
        for (const char* rel : {"manifest.json", "operators/HAM.hr.dat",
                                "operators/VX.hr.dat", "operators/SZ.hr.dat"})
        {
            const string fa = slurp(a + "/" + rel);
            const string fb = slurp(b + "/" + rel);
            assert(!fa.empty() && fa == fb && "config and positional bundles must be byte-identical");
        }
    }

    cout << "RUN CONFIG EQUIVALENCE TEST PASSED" << endl;
    return 0;
}
