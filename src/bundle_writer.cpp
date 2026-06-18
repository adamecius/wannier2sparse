#include "bundle_writer.hpp"
#include "json_writer.hpp"

#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h>   // mkdir (POSIX; the tool targets Linux/macOS)

namespace {

const char* W2SP_BUNDLE_SCHEMA_VERSION = "1.0";
const char* W2SP_BUNDLE_GENERATOR_VERSION = "1.0.0";

/// Create a directory if it does not already exist (no-op if present).
void make_dir(const std::string& path)
{
    if (path.empty() || path == ".") return;
    ::mkdir(path.c_str(), 0755);   // ignore EEXIST and other races; the writes below will surface real errors
}

/// Reciprocal vectors (rows) for lattice vectors given as rows; 0 if degenerate.
std::array<std::array<double,3>,3> reciprocal(const std::array<std::array<double,3>,3>& a)
{
    auto cross = [](const std::array<double,3>& u, const std::array<double,3>& v) {
        return std::array<double,3>{{ u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] }};
    };
    auto dot = [](const std::array<double,3>& u, const std::array<double,3>& v) {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    };
    const auto a1xa2 = cross(a[1], a[2]);
    const double V = dot(a[0], a1xa2);
    std::array<std::array<double,3>,3> b = {{{{0,0,0}},{{0,0,0}},{{0,0,0}}}};
    if (std::fabs(V) < 1e-300) return b;
    const double two_pi = 6.283185307179586476925286766559;
    const std::array<std::array<double,3>,3> num = {{ a1xa2, cross(a[2], a[0]), cross(a[0], a[1]) }};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            b[i][k] = two_pi * num[i][k] / V;
    return b;
}

/// Fractional coordinates of a Cartesian point given lattice rows (3x3 inverse).
std::array<double,3> to_fractional(const std::array<double,3>& cart,
                                   const std::array<std::array<double,3>,3>& a, bool& ok)
{
    // Columns of A^T are the lattice vectors; solve A^T f = cart for f.
    // Build matrix M with M[k][i] = a[i][k] (lattice vector i, component k).
    double M[3][3];
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            M[k][i] = a[i][k];
    const double det =
        M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
      - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
      + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    std::array<double,3> f = {{0,0,0}};
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

void emit_structure(JsonWriter& w, const SystemProvenance& prov)
{
    w.key("structure");
    w.begin_object();
        if (prov.has_lattice)
        {
            w.key("lattice_vectors");   w.matrix3(prov.lattice);
            w.key("lattice_units");     w.str("Angstrom");
            w.key("reciprocal_vectors");w.matrix3(reciprocal(prov.lattice));
        }
        else
        {
            w.member_null("lattice_vectors");
        }

        // DFT atoms (from QE). Null until the QE parser populates them.
        if (prov.dft.present && prov.dft.has_structure)
        {
            w.key("atoms");
            w.begin_array();
            for (const auto& at : prov.dft.atoms)
            {
                w.begin_object();
                    w.member("species", at.species);
                    w.key("frac"); w.array_d(at.frac);
                    w.key("cart"); w.array_d(at.cart);
                w.end_object();
            }
            w.end_array();
        }
        else
        {
            w.member_null("atoms");
        }

        // Wannier sites (from .xyz), always present when positions were loaded.
        w.key("wannier_sites");
        w.begin_array();
        for (const auto& s : prov.wannier_sites)
        {
            w.begin_object();
                w.member("index", (long long)s.index);
                w.member("label", s.label);
                w.key("cart"); w.array_d(s.cart);
                if (prov.has_lattice)
                {
                    bool ok = false;
                    auto f = to_fractional(s.cart, prov.lattice, ok);
                    if (ok) { w.key("frac"); w.array_d(f); } else { w.member_null("frac"); }
                }
                else { w.member_null("frac"); }
                w.member("spin", (long long)s.spin);
            w.end_object();
        }
        w.end_array();

        w.member("num_wann", (long long)prov.num_wann);
    w.end_object();
}

void emit_symmetry(JsonWriter& w, const SystemProvenance& prov)
{
    w.key("symmetry");
    if (!(prov.dft.present && prov.dft.has_structure && !prov.dft.symmetry.empty()))
    {
        w.null();
        return;
    }
    w.begin_object();
        w.member("n_sym", (long long)prov.dft.symmetry.size());
        w.member("translation_convention", std::string("fractional"));
        w.member("rotation_convention", std::string("crystal"));
        w.key("operators");
        w.begin_array();
        for (const auto& s : prov.dft.symmetry)
        {
            w.begin_object();
                w.key("rotation"); w.matrix3(s.rotation);
                w.key("translation"); w.array_d(s.translation);
                w.member("time_reversal", s.time_reversal);
            w.end_object();
        }
        w.end_array();
    w.end_object();
}

void emit_dft(JsonWriter& w, const DftProvenance& d)
{
    w.key("dft_provenance");
    if (!d.present) { w.null(); return; }
    w.begin_object();
        w.member("code", d.code);
        w.member("version", d.version);
        w.member("source_file", d.source_file);
        if (!d.xc_functional.empty()) w.member("xc_functional", d.xc_functional);
        else                          w.member_null("xc_functional");
        w.member("spin_orbit", d.spin_orbit);
        w.member("noncolin", d.noncolin);
        if (d.has_ecutwfc) w.member("ecutwfc_Ry", d.ecutwfc_Ry); else w.member_null("ecutwfc_Ry");
        if (d.has_k_mesh)
        {
            w.key("k_mesh");       w.begin_array(); for (int x : d.k_mesh)       w.inum(x); w.end_array();
            w.key("k_mesh_shift"); w.begin_array(); for (int x : d.k_mesh_shift) w.inum(x); w.end_array();
        }
        else { w.member_null("k_mesh"); w.member_null("k_mesh_shift"); }
        w.key("pseudopotentials");
        w.begin_array();
        for (const auto& p : d.pseudopotentials)
        {
            w.begin_object();
                w.member("species", p.species);
                w.member("file", p.file);
                if (p.has_z_valence) w.member("z_valence", p.z_valence);
            w.end_object();
        }
        w.end_array();
    w.end_object();
}

void emit_wannier(JsonWriter& w, const WannierProvenance& v)
{
    w.key("wannier_provenance");
    if (!v.present) { w.null(); return; }
    w.begin_object();
        w.member("source_file", v.source_file);
        if (v.has_num_wann)  w.member("num_wann", (long long)v.num_wann);   else w.member_null("num_wann");
        if (v.has_num_bands) w.member("num_bands", (long long)v.num_bands); else w.member_null("num_bands");
        w.key("exclude_bands"); w.begin_array(); for (int b : v.exclude_bands) w.inum(b); w.end_array();
        if (v.has_mp_grid) { w.key("mp_grid"); w.begin_array(); for (int x : v.mp_grid) w.inum(x); w.end_array(); }
        else               w.member_null("mp_grid");
        w.key("projections");
        w.begin_array();
        for (const auto& p : v.projections)
        {
            w.begin_object();
                w.member("site", p.site);
                w.member("orbital", p.orbital);
                w.member("raw", p.raw);
            w.end_object();
        }
        w.end_array();
        w.key("disentanglement");
        w.begin_object();
            w.member("enabled", v.dis_enabled);
            if (v.has_dis_win)  { w.member("dis_win_min", v.dis_win_min);   w.member("dis_win_max", v.dis_win_max); }
            if (v.has_dis_froz) { w.member("dis_froz_min", v.dis_froz_min); w.member("dis_froz_max", v.dis_froz_max); }
        w.end_object();
        w.member("use_ws_distance", v.use_ws_distance);
    w.end_object();
}

} // namespace

std::string write_bundle(const BundleSpec& spec, const SystemProvenance& prov,
                         const std::vector<BundleOperator>& ops, const std::string& out_dir)
{
    const std::string bundle_dir = (out_dir.empty() ? std::string(".") : out_dir) + "/" + spec.label + ".w2sp";
    const std::string ops_dir    = bundle_dir + "/operators";
    make_dir(out_dir);
    make_dir(bundle_dir);
    make_dir(ops_dir);

    // Operator data files (primitive O_ij(R) in _hr.dat shape) + remember n_rpts.
    std::vector<size_t> n_rpts(ops.size(), 0);
    for (size_t k = 0; k < ops.size(); ++k)
    {
        const std::string rel = "operators/" + ops[k].name + ".hr.dat";
        std::ofstream f((bundle_dir + "/" + rel).c_str());
        std::ostringstream hdr;
        hdr << "wannier2sparse bundle operator " << ops[k].name
            << " (" << ops[k].desc.observable << ")";
        // Count distinct R from the serialized output by capturing it once.
        write_hopping_list_as_hr(ops[k].hl, f, hdr.str());
        f.close();
    }

    // Copy _wsvec.dat verbatim if it was supplied.
    if (!spec.wsvec_src.empty())
    {
        std::ifstream src(spec.wsvec_src.c_str(), std::ios::binary);
        if (src.good())
        {
            std::ofstream dst((bundle_dir + "/wsvec.dat").c_str(), std::ios::binary);
            dst << src.rdbuf();
        }
    }

    // manifest.json
    std::ofstream mf((bundle_dir + "/manifest.json").c_str());
    JsonWriter w(mf);
    w.begin_object();
        w.member("schema_version", std::string(W2SP_BUNDLE_SCHEMA_VERSION));
        w.key("generator");
        w.begin_object();
            w.member("tool", std::string("wannier2sparse"));
            w.member("version", std::string(W2SP_BUNDLE_GENERATOR_VERSION));
        w.end_object();
        w.member("label", spec.label);
        w.member("provenance_complete", prov.dft.present && prov.wann.present);

        w.key("units");
        w.begin_object();
            w.member("length", std::string("Angstrom"));
            w.member("energy", std::string("eV"));
            w.member("spin", std::string("hbar/2"));
            w.member("orbital_L", std::string("hbar"));
        w.end_object();

        emit_structure(w, prov);
        emit_symmetry(w, prov);
        emit_dft(w, prov.dft);
        emit_wannier(w, prov.wann);

        w.key("normalization");
        w.begin_object();
            w.member("ndegen_applied", spec.ndegen_applied);
            w.member("wsvec_applied", spec.wsvec_applied);
            w.member("truncation_threshold", spec.truncation_threshold);
        w.end_object();

        w.key("operators");
        w.begin_array();
        for (const auto& op : ops)
        {
            w.begin_object();
                w.member("name", op.name);
                w.member("observable", op.desc.observable);
                w.member("component", op.desc.component);
                w.member("units", op.desc.units);
                w.member("provenance", op.desc.provenance);
                w.member("data_file", std::string("operators/") + op.name + ".hr.dat");
                w.member("data_format", std::string("wannier90_hr"));
                w.member("num_wann", (long long)op.hl.num_wann);
                if (op.desc.has_bounds)
                {
                    w.member("has_bounds", true);
                    w.member("spectral_min", op.desc.a);
                    w.member("spectral_max", op.desc.b);
                }
                else { w.member("has_bounds", false); }
            w.end_object();
        }
        w.end_array();
    w.end_object();
    mf << "\n";
    mf.close();

    return bundle_dir;
}
