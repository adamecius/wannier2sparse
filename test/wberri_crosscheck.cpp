#include <cassert>
#include <cstdlib>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include "gauge.hpp"
#include "wannier_parser.hpp"
#include "hopping_list.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;

// Plan 10C (independent path): compare the wannier2sparse spin operator, in the
// Wannier gauge on the mp_grid, against a WannierBerri-generated reference
// (test/fixtures/wberri_reference.py -> ref_Ok.dat). OFF by default
// (-DW2SP_WBERRI_CHECK=ON, label `wberri`); skips if inputs are missing.
//
// WS-convention trap: the reference MUST use the same use_ws_distance convention
// as <seed>_hr.dat, else O_W(k) differs by trivial phases (docs/conventions.md sec 6).
static Eigen::MatrixXcd at(const hopping_list& hl, const cellID_t& R, int nw)
{
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nw, nw);
    for (const auto& h : hl.hoppings)
        if (get<0>(h) == R) { const auto e = get<2>(h); M(e[0],e[1]) += get<1>(h); }
    return M;
}

int main()
{
    const char* fx = std::getenv("W2SP_SPIN_FIXTURE");
    const string p = fx ? fx : "/tmp/fix/Fe/Fe";
    const char* rf = std::getenv("W2SP_WBERRI_REF");
    const string ref = rf ? rf : "ref_Ok.dat";

    { ifstream a((p + ".spn").c_str()), b(ref.c_str());
      if (!a.good() || !b.good())
      { std::cout << "wberri_crosscheck: SKIP (need fixture " << p << " and reference " << ref << ")\n"; return 0; } }

    // Reference O_W(k): num_wann num_kpts ncomp, then Re Im (i fast, j, comp, k).
    ifstream f(ref.c_str());
    int nw, nk, nc; f >> nw >> nk >> nc;
    std::vector<std::array<Eigen::MatrixXcd,3> > refk(nk);
    for (int k = 0; k < nk; ++k) for (int c = 0; c < nc; ++c)
    {
        refk[k][c] = Eigen::MatrixXcd(nw, nw);
        for (int j = 0; j < nw; ++j) for (int i = 0; i < nw; ++i)
        { double re, im; f >> re >> im; refk[k][c](i,j) = complex<double>(re, im); }
    }

    gauge_data g = read_gauge(p);
    auto wd = read_wannier_file(p + "_hr.dat");
    const std::vector<int>& ndeg = std::get<1>(wd);
    const std::vector<string>& lines = std::get<2>(wd);
    const long blk = (long)nw * nw;
    std::vector<cellID_t> Rset; std::vector<double> wgt;
    for (size_t b = 0; b < ndeg.size(); ++b)
    { stringstream ss(lines[b*blk]); cellID_t R; ss>>R[0]>>R[1]>>R[2]; Rset.push_back(R); wgt.push_back(1.0/ndeg[b]); }

    const double twopi = 2.0*M_PI; double maxerr = 0.0;
    for (int alpha = 0; alpha < 3; ++alpha)
    {
        hopping_list S = exact_spin_operator(g, alpha, Rset);
        std::map<cellID_t, Eigen::MatrixXcd> Smat; for (const auto& R : Rset) Smat[R] = at(S, R, nw);
        for (int k = 0; k < nk; ++k)
        {
            Eigen::MatrixXcd recon = Eigen::MatrixXcd::Zero(nw, nw);
            for (size_t b = 0; b < Rset.size(); ++b)
            { const auto& R = Rset[b]; const double kr = g.kpt[k][0]*R[0]+g.kpt[k][1]*R[1]+g.kpt[k][2]*R[2];
              recon += wgt[b]*std::exp(complex<double>(0.0,+twopi*kr))*Smat[R]; }
            maxerr = std::max(maxerr, (recon - refk[k][alpha]).cwiseAbs().maxCoeff());
        }
    }
    std::cout << "wberri_crosscheck: max|S_W(k)_w2s - S_W(k)_ref| = " << maxerr << "\n";
    assert(maxerr < 1e-6);
    std::cout << "wberri_crosscheck: PASS\n";
    return 0;
}
