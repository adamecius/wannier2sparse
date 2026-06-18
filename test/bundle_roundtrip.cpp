// Phase B1: the bundle path ships the PRIMITIVE operator O_ij(R) as _hr.dat-shaped
// text. This test verifies the two guarantees that make the bundle usable by
// lsquant: (1) write_hopping_list_as_hr round-trips through the existing
// create_hopping_list(read_wannier_file(...)) parser, byte-of-meaning identical to
// the in-memory operator; and (2) write_bundle lays out a well-formed bundle whose
// manifest indexes every emitted operator data file.
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include "tbmodel.hpp"
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
#include "bundle_writer.hpp"

using namespace std;

// Write hl to a _hr.dat file, read it back, and assert equality with the original.
static void assert_roundtrip(hopping_list hl, const string& tmp)
{
    {
        ofstream f(tmp.c_str());
        write_hopping_list_as_hr(hl, f, "roundtrip");
    }
    hopping_list back = create_hopping_list(read_wannier_file(tmp));
    assert(back == hl && "primitive operator must round-trip through the _hr.dat writer/parser");
}

static string slurp(const string& path)
{
    ifstream f(path.c_str());
    stringstream ss; ss << f.rdbuf();
    return ss.str();
}

int main()
{
    tbmodel model;
    model.readOrbitalPositions("spin_graphene.xyz");
    model.readUnitCell("spin_graphene.uc");
    model.readWannierModel("spin_graphene_hr.dat");

    // Round-trip the Hamiltonian and two derived primitive operators.
    assert_roundtrip(model.hl, "rt_HAM.hr.dat");
    assert_roundtrip(model.createHoppingCurrents_list(0), "rt_VX.hr.dat");
    assert_roundtrip(model.createHoppingSpinDensity_list('z'), "rt_SZ.hr.dat");

    // Build a bundle and check its layout + manifest indexing.
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
    { BundleOperator op; op.name = "VX";  op.desc.observable = "velocity"; op.desc.component = "X"; op.desc.units = "eV*Angstrom"; op.hl = model.createHoppingCurrents_list(0); ops.push_back(op); }
    { BundleOperator op; op.name = "SZ";  op.desc.observable = "spin"; op.desc.component = "Z"; op.desc.units = "hbar/2"; op.hl = model.createHoppingSpinDensity_list('z'); ops.push_back(op); }

    BundleSpec spec;
    spec.label = "spin_graphene";
    spec.truncation_threshold = 1e-16;

    const string dir = write_bundle(spec, prov, ops, ".");
    assert(dir == "./spin_graphene.w2sp");

    // Manifest and every operator data file must exist.
    const string manifest = slurp(dir + "/manifest.json");
    assert(!manifest.empty() && "manifest.json must be written");
    for (const char* name : {"HAM", "VX", "SZ"})
    {
        ifstream df((dir + "/operators/" + name + ".hr.dat").c_str());
        assert(df.good() && "each operator must have a data file");
        // The manifest operator index must reference the operator and its file.
        assert(manifest.find(string("\"name\": \"") + name + "\"") != string::npos);
        assert(manifest.find(string("operators/") + name + ".hr.dat") != string::npos);
    }
    // Structure block must carry the lattice and the wannier sites.
    assert(manifest.find("lattice_vectors") != string::npos);
    assert(manifest.find("wannier_sites") != string::npos);
    assert(manifest.find("\"num_wann\": 4") != string::npos);
    // Unparsed provenance degrades to null, and the bundle is flagged incomplete.
    assert(manifest.find("\"dft_provenance\": null") != string::npos);
    assert(manifest.find("\"wannier_provenance\": null") != string::npos);
    assert(manifest.find("\"provenance_complete\": false") != string::npos);

    // The data file a consumer re-ingests must reproduce the in-memory operator.
    hopping_list ham_back = create_hopping_list(read_wannier_file(dir + "/operators/HAM.hr.dat"));
    assert(ham_back == model.hl);

    cout << "BUNDLE ROUNDTRIP TEST PASSED" << endl;
    return 0;
}
