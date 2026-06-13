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

// Plan 11, Test 1 -- matrix-element cross-check (docs/conventions.md sec 7).
// Compares the operator in the WANNIER GAUGE at every mp_grid k, element by
// element, against a committed WannierBerri golden:
//   O_W(k) = sum_R (1/ndegen_R) e^{+i 2pi k.R} O_W(R)   vs   ref O_W(k)
// This is gauge-dependent -> the STRICTEST check (FT-sign, .spn packing,
// real-harmonic order, disentanglement). Needs the same gauge files + WS
// convention on both sides; asserts the golden header's WS_convention matches.
// Reads only the .ref file (no WannierBerri). Skip-clean: missing golden OR
// fixture -> SKIP_RETURN_CODE. Goldens live in test/golden/ (W2S_GOLDEN_DIR).

#ifndef W2S_GOLDEN_DIR
#define W2S_GOLDEN_DIR "."
#endif
#ifndef W2S_WS_CONVENTION
#define W2S_WS_CONVENTION "II"          // use_ws_distance applied; override via env
#endif
static const int SKIP = 77;             // matches SKIP_RETURN_CODE in CMakeLists

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

// One (fixture, operator) candidate. Returns true if it actually ran a comparison;
// sets ok=false on a real mismatch. Leaves ran=false (returns false) when the
// golden or fixture is absent so the caller can decide to skip.
static bool run_one(const string& prefix, const string& op_tag, bool spin, bool& ok)
{
    const string seed   = basename_of(prefix);
    const string golden = golden_dir() + "/" + seed + "_" + op_tag + "_matrix.ref";
    const string probe  = prefix + (spin ? ".spn" : ".amn");
    { ifstream g(golden.c_str()), x(probe.c_str());
      if (!g.good() || !x.good()) return false; }

    ifstream f(golden.c_str());
    string h_seed, h_op, h_ws, h_ver; int nk = 0, nw = 0;
    { string line; getline(f, line); istringstream hs(line);
      hs >> h_seed >> h_op >> nk >> nw >> h_ws >> h_ver; }
    if (h_op != op_tag) { cerr << "matrix golden " << golden << ": operator '" << h_op
                               << "' != expected '" << op_tag << "'\n"; ok = false; return true; }
    if (h_ws != expected_ws()) { cerr << "matrix golden " << golden << ": WS_convention '"
                                      << h_ws << "' != expected '" << expected_ws()
                                      << "' (docs/conventions.md sec 6/7)\n"; ok = false; return true; }

    vector<array<Eigen::MatrixXcd,3> > refk(nk);
    for (int k = 0; k < nk; ++k) for (int a = 0; a < 3; ++a)
        refk[k][a] = Eigen::MatrixXcd::Zero(nw, nw);
    { int ik, m, n, a; double re, im;
      while (f >> ik >> m >> n >> a >> re >> im) refk[ik][a](m,n) = complex<double>(re, im); }

    gauge_data g = read_gauge(prefix);
    assert(g.num_wann == nw); assert(g.num_kpts == nk);

    // Complete WS R-set + ndegen weights from the _hr.dat (one R per nw*nw block).
    auto wd = read_wannier_file(prefix + "_hr.dat");
    const vector<int>&    ndeg  = std::get<1>(wd);
    const vector<string>& lines = std::get<2>(wd);
    const long blk = (long)nw * nw;
    vector<cellID_t> Rset; vector<double> wgt;
    for (size_t b = 0; b < ndeg.size(); ++b)
    { stringstream ss(lines[b*blk]); cellID_t R; ss>>R[0]>>R[1]>>R[2];
      Rset.push_back(R); wgt.push_back(1.0/ndeg[b]); }

    // Per-alpha operator R-matrices via the same engines the production path uses.
    vector<int> shells; amn_data A;
    if (!spin) { shells = parse_projection_shells(prefix + ".win"); A = read_amn(prefix + ".amn"); }

    const double twopi = 2.0*M_PI; double maxerr = 0.0;
    for (int alpha = 0; alpha < 3; ++alpha)
    {
        const hopping_list O = spin ? exact_spin_operator(g, alpha, Rset)
                                    : orbital_L_operator(g, A, shells, alpha, Rset);
        map<cellID_t, Eigen::MatrixXcd> Om; for (const auto& R : Rset) Om[R] = at(O, R, nw);
        for (int k = 0; k < nk; ++k)
        {
            Eigen::MatrixXcd recon = Eigen::MatrixXcd::Zero(nw, nw);
            for (size_t b = 0; b < Rset.size(); ++b)
            { const auto& R = Rset[b];
              const double kr = g.kpt[k][0]*R[0]+g.kpt[k][1]*R[1]+g.kpt[k][2]*R[2];
              recon += wgt[b]*std::exp(complex<double>(0.0,+twopi*kr))*Om[R]; }
            maxerr = std::max(maxerr, (recon - refk[k][alpha]).cwiseAbs().maxCoeff());
        }
    }
    cout << "wberri_matrix_crosscheck[" << seed << "_" << op_tag
         << "]: max|O_W(k)_w2s - O_W(k)_ref| = " << maxerr << "\n";
    if (maxerr >= 1e-6) { cerr << "  -> exceeds 1e-6 tolerance\n"; ok = false; }
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

    if (!ran) { cout << "wberri_matrix_crosscheck: SKIP (no golden+fixture pair found)\n"; return SKIP; }
    assert(ok);
    cout << "wberri_matrix_crosscheck: PASS\n";
    return 0;
}
