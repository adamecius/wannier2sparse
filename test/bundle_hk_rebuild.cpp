// Phase B5: end-to-end "lsquant could rebuild" acceptance test.
//
// The bundle ships the PRIMITIVE operator O_ij(R) plus a manifest; the whole point
// is that a consumer can rebuild the Bloch Hamiltonian H(k) = sum_R e^{ik.R} H(R)
// from it, with no help from this tool's supercell engine. This test proves that:
//
//   1. Build a physical nearest-neighbour graphene model (Hermitian -> real bands).
//   2. Write the bundle and re-ingest operators/HAM.hr.dat exactly as a consumer
//      would (create_hopping_list(read_wannier_file(...))).
//   3. Rebuild H(k) at the commensurate k-grid purely from the re-ingested H(R)
//      (integer phases 2*pi*m.R/N), diagonalise each small dense block.
//   4. Assert the multiset of H(k) eigenvalues equals the spectrum of the SAME
//      model expanded by the (sacred) supercell engine -- i.e. the bundle alone is
//      sufficient to reconstruct what the engine produces.
//   5. Cross-check the physics: Dirac touching (E=0) at K and |E|=3t at Gamma.
//   6. Determinism: regenerating the bundle is byte-identical.
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sys/stat.h>
#include <Eigen/Dense>

#include "hopping_list.hpp"
#include "wannier_parser.hpp"
#include "bundle_writer.hpp"
#include "sparse_matrix.hpp"

using namespace std;
typedef complex<double> cd;

// Nearest-neighbour graphene: 2 sublattices, hopping t on the three A-B bonds
// {R=(0,0), (-1,0), (0,-1)} and their Hermitian conjugates. Real bands +/-|f(k)|.
static hopping_list graphene_nn(double t)
{
    hopping_list hl;
    hl.SetWannierBasisSize(2);
    auto add = [&](array<int,3> R, int i, int j, double re)
    { hl.hoppings.push_back(make_tuple(R, cd(re, 0.0), array<int,2>{{i,j}})); };
    add({{0,0,0}},  0, 1, t);
    add({{-1,0,0}}, 0, 1, t);
    add({{0,-1,0}}, 0, 1, t);
    add({{0,0,0}},  1, 0, t);     // Hermitian conjugates (real t)
    add({{1,0,0}},  1, 0, t);
    add({{0,1,0}},  1, 0, t);
    return hl;
}

static vector<double> sorted_eigs(const Eigen::MatrixXcd& H)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H);
    vector<double> e(es.eigenvalues().data(), es.eigenvalues().data() + es.eigenvalues().size());
    sort(e.begin(), e.end());
    return e;
}

// Rebuild H(k) at a commensurate k = (m1/N1, m2/N2, m3/N3) from a primitive H(R),
// using only integer phases 2*pi*(m.R)/N -- no lattice vectors needed.
static Eigen::MatrixXcd Hk(const hopping_list& hl, int m1, int m2, int m3,
                           int N1, int N2, int N3)
{
    const int nw = hl.WannierBasisSize();
    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(nw, nw);
    const double tau = 2.0 * acos(-1.0);
    for (const auto& h : hl.hoppings)
    {
        const auto& R = get<0>(h);
        const auto& e = get<2>(h);
        const double theta = tau * ((double)m1*R[0]/N1 + (double)m2*R[1]/N2 + (double)m3*R[2]/N3);
        H(e[0], e[1]) += get<1>(h) * polar(1.0, theta);
    }
    return H;
}

static string slurp(const string& p) { ifstream f(p.c_str()); stringstream ss; ss << f.rdbuf(); return ss.str(); }

static string write_gr_bundle(const hopping_list& hl, const string& out_dir)
{
    SystemProvenance prov;
    prov.num_wann = hl.WannierBasisSize();
    prov.has_lattice = true;
    prov.lattice = {{ {{0.5, 0.8660254, 0.0}}, {{0.5, -0.8660254, 0.0}}, {{0.0, 0.0, 1.0}} }};
    { WannierSite s; s.index = 0; s.label = "A"; s.cart = {{0,0,0}}; prov.wannier_sites.push_back(s); }
    { WannierSite s; s.index = 1; s.label = "B"; s.cart = {{0.5,0.288675,0}}; prov.wannier_sites.push_back(s); }

    vector<BundleOperator> ops;
    BundleOperator op; op.name = "HAM"; op.desc.observable = "hamiltonian"; op.desc.units = "eV"; op.hl = hl;
    ops.push_back(op);

    BundleSpec spec; spec.label = "gr";
    return write_bundle(spec, prov, ops, out_dir);
}

int main()
{
    const double t = 1.0;
    hopping_list hl = graphene_nn(t);

    // --- write bundle twice and re-ingest the primitive H(R) -----------------
    mkdir("A", 0755); mkdir("B", 0755);
    const string dirA = write_gr_bundle(hl, "A");
    const string dirB = write_gr_bundle(hl, "B");

    // Determinism: both runs byte-identical (manifest + operator data file).
    assert(slurp(dirA + "/manifest.json")          == slurp(dirB + "/manifest.json"));
    assert(slurp(dirA + "/operators/HAM.hr.dat")   == slurp(dirB + "/operators/HAM.hr.dat"));

    hopping_list back = create_hopping_list(read_wannier_file(dirA + "/operators/HAM.hr.dat"));
    assert(back == hl && "re-ingested HAM must equal the original primitive operator");

    // --- (1) bundle-rebuilt H(k) spectrum == engine supercell spectrum -------
    const int N = 3;   // commensurate grid; also satisfies N >= 2*range+1 = 3
    vector<double> from_bundle;
    for (int m1 = 0; m1 < N; ++m1)
        for (int m2 = 0; m2 < N; ++m2)
        {
            vector<double> e = sorted_eigs(Hk(back, m1, m2, 0, N, N, 1));
            from_bundle.insert(from_bundle.end(), e.begin(), e.end());
        }
    sort(from_bundle.begin(), from_bundle.end());

    const Eigen::MatrixXcd super(supercell_matrix({{N, N, 1}}, hl));   // engine expansion (dense)
    vector<double> from_engine = sorted_eigs(super);

    assert(from_bundle.size() == from_engine.size());
    for (size_t i = 0; i < from_bundle.size(); ++i)
        assert(fabs(from_bundle[i] - from_engine[i]) < 1e-9 &&
               "bundle-rebuilt H(k) spectrum must match the supercell engine");

    // --- (2) physics cross-check: Gamma and Dirac ---------------------------
    {
        vector<double> g = sorted_eigs(Hk(back, 0, 0, 0, N, N, 1));   // Gamma: +/-3t
        assert(fabs(g.front() + 3.0*t) < 1e-9 && fabs(g.back() - 3.0*t) < 1e-9);

        // K point lives at (m1,m2)=(1,2) on the 3x3 grid: 1+e^{-i2pi/3}+e^{-i4pi/3}=0.
        vector<double> k = sorted_eigs(Hk(back, 1, 2, 0, N, N, 1));
        assert(fabs(k.front()) < 1e-9 && fabs(k.back()) < 1e-9 && "Dirac touching at K");
    }

    cout << "BUNDLE H(k) REBUILD TEST PASSED" << endl;
    return 0;
}
