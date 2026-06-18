/**
 * @file qe_xml_parser.cpp
 * @brief Implementation of the QE data-file-schema.xml provenance parser.
 *
 * Navigation is by unqualified child names off the document root (only QE's root
 * element, `qes:espresso`, carries the namespace prefix; its children do not), so
 * the same code reads the QE 6.x and 7.x layouts. Each field is looked up
 * independently and missing nodes are tolerated.
 */
#include "qe_xml_parser.hpp"
#include "third_party/pugixml/pugixml.hpp"

#include <stdexcept>
#include <sstream>
#include <string>
#include <cmath>
#include <cctype>

namespace {

// CODATA Bohr radius in Angstrom; QE XML stores structure in Bohr (atomic units).
const double BOHR_TO_ANG  = 0.529177210903;
// QE XML stores energies in Hartree; the manifest reports ecutwfc in Rydberg.
const double HARTREE_TO_RY = 2.0;

std::string trim(const std::string& s)
{
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

bool parse_bool(const std::string& s)
{
    const std::string t = trim(s);
    return (t == "true" || t == "TRUE" || t == "T" || t == "true." ||
            t == ".true." || t == "1");
}

std::array<double,3> parse_vec3(const std::string& text)
{
    std::array<double,3> v{{0,0,0}};
    std::stringstream ss(text);
    ss >> v[0] >> v[1] >> v[2];
    return v;
}

// Fractional coordinates of a Cartesian point given lattice rows (3x3 inverse).
std::array<double,3> cart_to_frac(const std::array<double,3>& cart,
                                  const std::array<std::array<double,3>,3>& a, bool& ok)
{
    // Lattice vector i has components a[i][k]; solve (sum_i f_i a[i]) = cart, i.e.
    // M f = cart with M[k][i] = a[i][k].
    double M[3][3];
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            M[k][i] = a[i][k];
    const double det =
        M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
      - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
      + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    std::array<double,3> f{{0,0,0}};
    if (std::fabs(det) < 1e-300) { ok = false; return f; }
    ok = true;
    double inv[3][3];
    inv[0][0] =  (M[1][1]*M[2][2]-M[1][2]*M[2][1])/det;
    inv[0][1] = -(M[0][1]*M[2][2]-M[0][2]*M[2][1])/det;
    inv[0][2] =  (M[0][1]*M[1][2]-M[0][2]*M[1][1])/det;
    inv[1][0] = -(M[1][0]*M[2][2]-M[1][2]*M[2][0])/det;
    inv[1][1] =  (M[0][0]*M[2][2]-M[0][2]*M[2][0])/det;
    inv[1][2] = -(M[0][0]*M[1][2]-M[0][2]*M[1][0])/det;
    inv[2][0] =  (M[1][0]*M[2][1]-M[1][1]*M[2][0])/det;
    inv[2][1] = -(M[0][0]*M[2][1]-M[0][1]*M[2][0])/det;
    inv[2][2] =  (M[0][0]*M[1][1]-M[0][1]*M[1][0])/det;
    for (int i = 0; i < 3; ++i)
        f[i] = inv[i][0]*cart[0] + inv[i][1]*cart[1] + inv[i][2]*cart[2];
    return f;
}

// First non-empty child node of the same name among the given parents.
pugi::xml_node first_child(const pugi::xml_node& a, const pugi::xml_node& b, const char* name)
{
    pugi::xml_node n = a.child(name);
    return n ? n : b.child(name);
}

} // namespace

DftProvenance parse_qe_xml(const std::string& xml_path)
{
    pugi::xml_document doc;
    const pugi::xml_parse_result res = doc.load_file(xml_path.c_str());
    if (!res)
        throw std::runtime_error("parse_qe_xml: cannot read '" + xml_path + "': " +
                                 res.description());

    const pugi::xml_node root = doc.document_element();      // qes:espresso
    if (!root)
        throw std::runtime_error("parse_qe_xml: '" + xml_path + "' has no root element");

    DftProvenance d;
    d.present = true;
    d.source_file = xml_path;
    d.code = "Quantum ESPRESSO";

    // Code version (general_info/creator VERSION).
    pugi::xml_node creator = root.child("general_info").child("creator");
    if (creator) d.version = creator.attribute("VERSION").value();

    const pugi::xml_node out = root.child("output");
    const pugi::xml_node in  = root.child("input");

    // ---- structure: lattice (Bohr->Angstrom) + atoms ------------------------
    pugi::xml_node astruct = first_child(out, in, "atomic_structure");
    if (astruct)
    {
        pugi::xml_node cell = astruct.child("cell");
        if (cell)
        {
            const char* names[3] = {"a1", "a2", "a3"};
            bool ok = true;
            for (int i = 0; i < 3 && ok; ++i)
            {
                pugi::xml_node a = cell.child(names[i]);
                if (!a) { ok = false; break; }
                const std::array<double,3> v = parse_vec3(a.child_value());
                for (int k = 0; k < 3; ++k) d.lattice[i][k] = v[k] * BOHR_TO_ANG;
            }
            if (ok) d.has_structure = true;
        }

        pugi::xml_node apos = astruct.child("atomic_positions");
        if (apos)
            for (pugi::xml_node atom = apos.child("atom"); atom; atom = atom.next_sibling("atom"))
            {
                AtomSite s;
                s.species = atom.attribute("name").value();
                const std::array<double,3> c = parse_vec3(atom.child_value());
                for (int k = 0; k < 3; ++k) s.cart[k] = c[k] * BOHR_TO_ANG;
                if (d.has_structure)
                {
                    bool ok = false;
                    s.frac = cart_to_frac(s.cart, d.lattice, ok);
                }
                d.atoms.push_back(s);
            }
    }

    // ---- symmetry: the first nsym crystal operations ------------------------
    pugi::xml_node syms = out.child("symmetries");
    if (syms && d.has_structure)
    {
        const int nsym = syms.child("nsym").text().as_int(-1);
        int count = 0;
        for (pugi::xml_node sym = syms.child("symmetry"); sym; sym = sym.next_sibling("symmetry"))
        {
            if (nsym >= 0 && count >= nsym) break;   // skip lattice-only rotations
            pugi::xml_node rot = sym.child("rotation");
            if (!rot) continue;

            std::vector<double> r;
            { std::stringstream ss(rot.child_value()); double x; while (ss >> x) r.push_back(x); }
            if (r.size() < 9) continue;

            // QE writes matrices in Fortran (column-major) order unless told "C".
            const std::string order = rot.attribute("order").value();
            const bool col_major = (order.find('C') == std::string::npos &&
                                    order.find('c') == std::string::npos);

            SymmetryOp op;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    op.rotation[i][j] = col_major ? r[i + 3*j] : r[3*i + j];

            pugi::xml_node ft = sym.child("fractional_translation");
            if (ft) op.translation = parse_vec3(ft.child_value());

            pugi::xml_node info = sym.child("info");
            if (info)
            {
                const std::string tr = info.attribute("time_reversal").value();
                if (!tr.empty()) op.time_reversal = parse_bool(tr);
            }

            d.symmetry.push_back(op);
            ++count;
        }
    }

    // ---- DFT conditions -----------------------------------------------------
    pugi::xml_node functional = out.child("dft").child("functional");
    if (functional) d.xc_functional = trim(functional.child_value());

    // spin-orbit / noncollinear: output band_structure, else input spin block.
    pugi::xml_node bs   = out.child("band_structure");
    pugi::xml_node spin = in.child("spin");
    {
        pugi::xml_node so = bs.child("spinorbit");
        if (!so) so = spin.child("spinorbit");
        if (so) d.spin_orbit = parse_bool(so.child_value());
        pugi::xml_node nc = bs.child("noncolin");
        if (!nc) nc = spin.child("noncolin");
        if (nc) d.noncolin = parse_bool(nc.child_value());
    }

    // ecutwfc (Hartree in the XML -> Rydberg in the manifest).
    pugi::xml_node ecut = out.child("basis_set").child("ecutwfc");
    if (!ecut) ecut = in.child("basis").child("ecutwfc");
    if (ecut)
    {
        d.ecutwfc_Ry = ecut.text().as_double() * HARTREE_TO_RY;
        d.has_ecutwfc = true;
    }

    // k-mesh (Monkhorst-Pack grid + shift).
    pugi::xml_node mp = in.child("k_points_IBZ").child("monkhorst_pack");
    if (!mp) mp = out.child("k_points_IBZ").child("monkhorst_pack");
    if (mp)
    {
        d.k_mesh       = {{ mp.attribute("nk1").as_int(), mp.attribute("nk2").as_int(), mp.attribute("nk3").as_int() }};
        d.k_mesh_shift = {{ mp.attribute("k1").as_int(),  mp.attribute("k2").as_int(),  mp.attribute("k3").as_int()  }};
        d.has_k_mesh = true;
    }

    // pseudopotentials (one per species).
    pugi::xml_node asp = first_child(out, in, "atomic_species");
    if (asp)
        for (pugi::xml_node sp = asp.child("species"); sp; sp = sp.next_sibling("species"))
        {
            PseudoInfo p;
            p.species = sp.attribute("name").value();
            p.file    = trim(sp.child("pseudo_file").child_value());
            d.pseudopotentials.push_back(p);
        }

    return d;
}
