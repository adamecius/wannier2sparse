// Phase B4: parse a QE data-file-schema.xml into DftProvenance. Two synthetic
// snippets exercise the QE 7.x layout (ecutwfc + spin flags in <output>) and the
// QE 6.x layout (ecutwfc in <input>/basis, spin in <input>/spin). The test asserts
// the lattice (Bohr->Angstrom), atoms (species + fractional coords), one symmetry
// op (Fortran column-major rotation + fractional translation, with ops beyond nsym
// skipped), the k-mesh, SOC/noncolin, ecutwfc (Hartree->Ry) and pseudo names, plus
// graceful degradation when nodes are absent.
#include <cassert>
#include <fstream>
#include <cmath>
#include <string>
#include <iostream>
#include "qe_xml_parser.hpp"

using namespace std;

static const double BOHR = 0.529177210903;
static bool close(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; }

int main()
{
    // ---- QE 7.x flavor: structure + symmetry + DFT in <output> --------------
    {
        ofstream("qe7.xml") <<
"<?xml version=\"1.0\"?>\n"
"<qes:espresso xmlns:qes=\"http://www.quantum-espresso.org/ns/qes/qes-1.0\">\n"
"  <general_info><creator NAME=\"PWSCF\" VERSION=\"7.2\"/></general_info>\n"
"  <input>\n"
"    <k_points_IBZ>\n"
"      <monkhorst_pack nk1=\"4\" nk2=\"4\" nk3=\"4\" k1=\"0\" k2=\"0\" k3=\"1\">MP</monkhorst_pack>\n"
"    </k_points_IBZ>\n"
"  </input>\n"
"  <output>\n"
"    <atomic_species ntyp=\"1\">\n"
"      <species name=\"Si\"><mass>28.0</mass><pseudo_file>Si.pbe-n.UPF</pseudo_file></species>\n"
"    </atomic_species>\n"
"    <atomic_structure nat=\"2\" alat=\"5.0\">\n"
"      <atomic_positions>\n"
"        <atom name=\"Si\" index=\"1\">0.0 0.0 0.0</atom>\n"
"        <atom name=\"Si\" index=\"2\">2.5 2.5 2.5</atom>\n"
"      </atomic_positions>\n"
"      <cell>\n"
"        <a1>5.0 0.0 0.0</a1>\n"
"        <a2>0.0 5.0 0.0</a2>\n"
"        <a3>0.0 0.0 5.0</a3>\n"
"      </cell>\n"
"    </atomic_structure>\n"
"    <symmetries>\n"
"      <nsym>1</nsym>\n"
"      <nrot>2</nrot>\n"
"      <symmetry>\n"
"        <info name=\"r4\">crystal_symmetry</info>\n"
"        <rotation rank=\"2\" dims=\"3 3\" order=\"F\">0 1 0 -1 0 0 0 0 1</rotation>\n"
"        <fractional_translation>0.0 0.0 0.5</fractional_translation>\n"
"      </symmetry>\n"
"      <symmetry>\n"
"        <info name=\"lattice_only\">not_a_crystal_symmetry</info>\n"
"        <rotation rank=\"2\" dims=\"3 3\" order=\"F\">1 0 0 0 1 0 0 0 1</rotation>\n"
"        <fractional_translation>0.0 0.0 0.0</fractional_translation>\n"
"      </symmetry>\n"
"    </symmetries>\n"
"    <basis_set><ecutwfc>30.0</ecutwfc></basis_set>\n"
"    <dft><functional>PBE</functional></dft>\n"
"    <band_structure><spinorbit>true</spinorbit><noncolin>true</noncolin></band_structure>\n"
"  </output>\n"
"</qes:espresso>\n";

        DftProvenance d = parse_qe_xml("qe7.xml");
        assert(d.present);
        assert(d.code == "Quantum ESPRESSO");
        assert(d.version == "7.2");

        // Lattice: 5 Bohr -> Angstrom along the diagonal.
        assert(d.has_structure);
        assert(close(d.lattice[0][0], 5.0 * BOHR) && close(d.lattice[1][1], 5.0 * BOHR));
        assert(close(d.lattice[2][2], 5.0 * BOHR));

        // Atoms: species + Cartesian (Bohr->Ang) + fractional (cubic => cart/5).
        assert(d.atoms.size() == 2);
        assert(d.atoms[0].species == "Si" && d.atoms[1].species == "Si");
        assert(close(d.atoms[0].frac[0], 0.0) && close(d.atoms[1].frac[0], 0.5) &&
               close(d.atoms[1].frac[2], 0.5));
        assert(close(d.atoms[1].cart[0], 2.5 * BOHR));

        // Symmetry: only the first nsym=1 op; Fortran column-major rotation.
        assert(d.symmetry.size() == 1);
        assert(close(d.symmetry[0].rotation[0][1], -1.0));
        assert(close(d.symmetry[0].rotation[1][0],  1.0));
        assert(close(d.symmetry[0].rotation[2][2],  1.0));
        assert(close(d.symmetry[0].rotation[0][0],  0.0));
        assert(close(d.symmetry[0].translation[2], 0.5));
        assert(d.symmetry[0].time_reversal == false);

        // DFT conditions.
        assert(d.xc_functional == "PBE");
        assert(d.spin_orbit == true && d.noncolin == true);
        assert(d.has_ecutwfc && close(d.ecutwfc_Ry, 60.0));   // 30 Hartree -> 60 Ry
        assert(d.has_k_mesh && d.k_mesh[0] == 4 && d.k_mesh[2] == 4);
        assert(d.k_mesh_shift[2] == 1);
        assert(d.pseudopotentials.size() == 1);
        assert(d.pseudopotentials[0].species == "Si" &&
               d.pseudopotentials[0].file == "Si.pbe-n.UPF");
    }

    // ---- QE 6.x flavor: ecutwfc in <input>/basis, spin in <input>/spin ------
    {
        ofstream("qe6.xml") <<
"<?xml version=\"1.0\"?>\n"
"<qes:espresso xmlns:qes=\"http://www.quantum-espresso.org/ns/qes/qes-1.0\">\n"
"  <general_info><creator NAME=\"PWSCF\" VERSION=\"6.7MaX\"/></general_info>\n"
"  <input>\n"
"    <basis><ecutwfc>25.0</ecutwfc></basis>\n"
"    <spin><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit></spin>\n"
"    <k_points_IBZ>\n"
"      <monkhorst_pack nk1=\"8\" nk2=\"8\" nk3=\"1\" k1=\"0\" k2=\"0\" k3=\"0\">MP</monkhorst_pack>\n"
"    </k_points_IBZ>\n"
"    <atomic_species ntyp=\"1\">\n"
"      <species name=\"C\"><pseudo_file>C.UPF</pseudo_file></species>\n"
"    </atomic_species>\n"
"  </input>\n"
"  <output>\n"
"    <atomic_structure nat=\"1\">\n"
"      <atomic_positions><atom name=\"C\">0.0 0.0 0.0</atom></atomic_positions>\n"
"      <cell><a1>4.0 0.0 0.0</a1><a2>0.0 4.0 0.0</a2><a3>0.0 0.0 4.0</a3></cell>\n"
"    </atomic_structure>\n"
"    <dft><functional>LDA</functional></dft>\n"
"  </output>\n"
"</qes:espresso>\n";

        DftProvenance d = parse_qe_xml("qe6.xml");
        assert(d.present && d.version == "6.7MaX");
        assert(d.has_structure && close(d.lattice[0][0], 4.0 * BOHR));
        assert(d.atoms.size() == 1 && d.atoms[0].species == "C");
        assert(d.xc_functional == "LDA");
        assert(d.spin_orbit == false && d.noncolin == false);
        assert(d.has_ecutwfc && close(d.ecutwfc_Ry, 50.0));   // 25 Hartree -> 50 Ry (input fallback)
        assert(d.has_k_mesh && d.k_mesh[0] == 8 && d.k_mesh[2] == 1);
        assert(d.pseudopotentials.size() == 1 && d.pseudopotentials[0].file == "C.UPF");
        // No <symmetries> node -> empty symmetry list, but structure still present.
        assert(d.symmetry.empty());
    }

    // ---- Graceful degradation: malformed XML and missing file throw ---------
    {
        ofstream("broken.xml") << "<qes:espresso><output><cell>\n";   // unterminated
        bool threw = false;
        try { parse_qe_xml("broken.xml"); } catch (const std::exception&) { threw = true; }
        assert(threw && "malformed XML must throw");

        threw = false;
        try { parse_qe_xml("nope.xml"); } catch (const std::exception&) { threw = true; }
        assert(threw && "missing file must throw");
    }

    cout << "QE PROVENANCE TEST PASSED" << endl;
    return 0;
}
