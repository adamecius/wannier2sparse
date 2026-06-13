#include <cassert>
#include <cstdlib>
#include <fstream>
#include <map>
#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include "gauge.hpp"
#include "local_operators.hpp"
#include "wannier_parser.hpp"
#include "hopping_list.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;

// Plan 11, Test 2 -- band-texture cross-check (docs/conventions.md sec 7).
// Compares the band-resolved expectation <O_alpha>_{nk} = <psi_nk| O_alpha(k) |psi_nk>
// reconstructed by w2s (build O_alpha(k), diagonalize H(k), project) against a
// committed WannierBerri golden. This is GAUGE-INVARIANT -> more robust than
// Test 1 and the independent-physics closure. Degenerate multiplets make a
// per-band <O> ill-defined (Fe+SOC), so the comparator sums <O> over each
// degenerate block (subspace trace), using the golden's block tags as the shared
// partition. Reads only the .ref (no WannierBerri). Skip-clean.

#ifndef W2S_GOLDEN_DIR
#define W2S_GOLDEN_DIR "."
#endif
#ifndef W2S_WS_CONVENTION
#define W2S_WS_CONVENTION "II"
#endif
static const int SKIP = 77;

static Eigen::MatrixXcd at(const hopping_list& hl, const cellID_t& R, int nw)
{
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nw, nw);
    for (const auto& h : hl.hoppings)
        if (get<0>(h) == R) { const auto e = get<2>(h); M(e[0],e[1]) += get<1>(h); }
    return M;
}

static string basename_of(const string& p)
{ const size_t s = p.find_last_of('/'); return s == string::npos ? p : p.substr(s+1); }

static string golden_dir()
{ const char* e = getenv("W2SP_GOLDEN_DIR"); return e ? e : W2S_GOLDEN_DIR; }

static string expected_ws()
{ const char* e = getenv("W2SP_WS_CONVENTION"); return e ? e : W2S_WS_CONVENTION; }

static double texture_tol()
{ const char* e = getenv("W2SP_TEXTURE_TOL"); return e ? atof(e) : 1e-5; }

static bool run_one(const string& prefix, const string& op_tag, bool spin, bool& ok)
{
    const string seed   = basename_of(prefix);
    const string golden = golden_dir() + "/" + seed + "_" + op_tag + "_texture.ref";
    const string probe  = prefix + (spin ? ".spn" : ".amn");
    { ifstream g(golden.c_str()), x(probe.c_str());
      if (!g.good() || !x.good()) return false; }

    ifstream f(golden.c_str());
    string h_seed, h_op, h_ws, h_ver; int nk = 0, nw = 0; double deg_tol = 1e-5;
    { string line; getline(f, line); istringstream hs(line);
      hs >> h_seed >> h_op >> nk >> nw >> h_ws >> h_ver >> deg_tol; }
    if (h_op != op_tag) { cerr << "texture golden " << golden << ": operator '" << h_op
                               << "' != expected '" << op_tag << "'\n"; ok = false; return true; }
    if (h_ws != expected_ws()) { cerr << "texture golden " << golden << ": WS_convention '"
                                      << h_ws << "' != expected '" << expected_ws()
                                      << "' (docs/conventions.md sec 6/7)\n"; ok = false; return true; }

    // ref_val[ik][ibnd][alpha], block[ik][ibnd] (the shared degenerate partition).
    vector<vector<array<double,3> > > refv(nk, vector<array<double,3> >(nw));
    vector<vector<int> > block(nk, vector<int>(nw, 0));
    { int ik, ibnd, a, blk; double v;
      while (f >> ik >> ibnd >> a >> blk >> v) { refv[ik][ibnd][a] = v; block[ik][ibnd] = blk; } }

    gauge_data g = read_gauge(prefix);
    assert(g.num_wann == nw); assert(g.num_kpts == nk);

    auto wd = read_wannier_file(prefix + "_hr.dat");
    const vector<int>&    ndeg  = std::get<1>(wd);
    const vector<string>& lines = std::get<2>(wd);
    const long blk = (long)nw * nw;
    vector<cellID_t> Rset; vector<double> wgt;
    for (size_t b = 0; b < ndeg.size(); ++b)
    { stringstream ss(lines[b*blk]); cellID_t R; ss>>R[0]>>R[1]>>R[2];
      Rset.push_back(R); wgt.push_back(1.0/ndeg[b]); }

    // H(R): create_hopping_list already folds 1/ndegen (Plan 1), so H(k) carries no
    // extra weight. The operator engines do NOT fold ndegen -> apply wgt below.
    const hopping_list H = create_hopping_list(wd);
    map<cellID_t, Eigen::MatrixXcd> Hm; for (const auto& R : Rset) Hm[R] = at(H, R, nw);

    vector<int> shells; amn_data A;
    if (!spin) { shells = parse_projection_shells(prefix + ".win"); A = read_amn(prefix + ".amn"); }
    array<map<cellID_t, Eigen::MatrixXcd>,3> Om;
    for (int a = 0; a < 3; ++a)
    { const hopping_list O = spin ? exact_spin_operator(g, a, Rset)
                                  : orbital_L_operator(g, A, shells, a, Rset);
      for (const auto& R : Rset) Om[a][R] = at(O, R, nw); }

    const double twopi = 2.0*M_PI; double maxerr = 0.0;
    for (int k = 0; k < nk; ++k)
    {
        Eigen::MatrixXcd Hk = Eigen::MatrixXcd::Zero(nw, nw);
        array<Eigen::MatrixXcd,3> Ok = { Eigen::MatrixXcd::Zero(nw,nw),
                                         Eigen::MatrixXcd::Zero(nw,nw),
                                         Eigen::MatrixXcd::Zero(nw,nw) };
        for (size_t b = 0; b < Rset.size(); ++b)
        { const auto& R = Rset[b];
          const double kr = g.kpt[k][0]*R[0]+g.kpt[k][1]*R[1]+g.kpt[k][2]*R[2];
          const complex<double> ph = std::exp(complex<double>(0.0,+twopi*kr));
          Hk += ph * Hm[R];
          for (int a = 0; a < 3; ++a) Ok[a] += wgt[b] * ph * Om[a][R]; }
        Hk = 0.5*(Hk + Hk.adjoint().eval());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
        const Eigen::MatrixXcd& U = es.eigenvectors();   // columns, ascending energy

        // Subspace trace per (golden) degenerate block: sum <O_alpha> over its bands.
        map<int, array<double,3> > w2s_sum, ref_sum;
        for (int ibnd = 0; ibnd < nw; ++ibnd)
        {
            const int bl = block[k][ibnd];
            const Eigen::VectorXcd psi = U.col(ibnd);
            for (int a = 0; a < 3; ++a)
            { const double ev = (psi.adjoint()*Ok[a]*psi).value().real();
              w2s_sum[bl][a] += ev; ref_sum[bl][a] += refv[k][ibnd][a]; }
        }
        for (const auto& kv : w2s_sum) for (int a = 0; a < 3; ++a)
            maxerr = std::max(maxerr, std::abs(kv.second[a] - ref_sum[kv.first][a]));
    }
    cout << "wberri_texture_crosscheck[" << seed << "_" << op_tag
         << "]: max|<O>_block(w2s) - <O>_block(ref)| = " << maxerr
         << " (deg_tol " << deg_tol << ")\n";
    if (maxerr >= texture_tol()) { cerr << "  -> exceeds " << texture_tol() << " tolerance\n"; ok = false; }
    return true;
}

int main()
{
    const char* sp = getenv("W2SP_SPIN_FIXTURE");
    const char* or_ = getenv("W2SP_ORBITAL_FIXTURE");
    const string spin_fx = sp  ? sp  : "/tmp/fix/Fe/Fe";
    const string orb_fx  = or_ ? or_ : "/tmp/fix/copper/copper";

    bool ok = true, ran = false;
    ran |= run_one(spin_fx, "S", true,  ok);
    ran |= run_one(orb_fx,  "L", false, ok);

    if (!ran) { cout << "wberri_texture_crosscheck: SKIP (no golden+fixture pair found)\n"; return SKIP; }
    assert(ok);
    cout << "wberri_texture_crosscheck: PASS\n";
    return 0;
}
