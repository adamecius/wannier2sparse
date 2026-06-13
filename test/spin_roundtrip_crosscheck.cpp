#include <cassert>
#include <cstdlib>
#include <fstream>
#include <map>
#include <array>
#include <vector>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include "gauge.hpp"
#include "wannier_parser.hpp"
#include "hopping_list.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;

// Plan 10C Level 1 (implementation cross-check, no external package): the exact
// spin operator written by wannier2sparse must, when inverse-transformed with
// the W90 convention S_W(k)=sum_R (1/ndegen) e^{+i2pi k.R} S_W(R), reproduce the
// direct V(k)^dag S_B(k) V(k) at every mesh k. Same definition by two paths ->
// decisive for FT-sign / gauge convention bugs. Skips if the fixture is absent.
static Eigen::MatrixXcd at(const hopping_list& hl, const cellID_t& R, int nw)
{
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nw, nw);
    for (const auto& h : hl.hoppings)
        if (get<0>(h) == R) { const auto e = get<2>(h); M(e[0],e[1]) += get<1>(h); }
    return M;
}

int main()
{
    const char* env = std::getenv("W2SP_SPIN_FIXTURE");
    const string p = env ? env : "/tmp/fix/Fe/Fe";
    { ifstream f((p + ".spn").c_str()); if (!f.good())
      { std::cout << "spin_roundtrip_crosscheck: SKIP (no fixture at " << p << ")\n"; return 0; } }

    gauge_data g = read_gauge(p);

    // Complete WS R-set + ndegen from the _hr.dat (R per num_wann^2 block).
    auto wd = read_wannier_file(p + "_hr.dat");
    const int nw = std::get<0>(wd);
    const std::vector<int>&    ndeg  = std::get<1>(wd);
    const std::vector<string>& lines = std::get<2>(wd);
    const long blk = (long)nw * nw;
    std::vector<cellID_t> Rset; std::vector<double> wgt;
    for (size_t b = 0; b < ndeg.size(); ++b)
    {
        std::stringstream ss(lines[b*blk]); cellID_t R; ss >> R[0] >> R[1] >> R[2];
        Rset.push_back(R); wgt.push_back(1.0 / ndeg[b]);
    }

    const int alpha = 2;                                  // S_z
    hopping_list Sr = exact_spin_operator(g, alpha, Rset);
    std::map<cellID_t, Eigen::MatrixXcd> Smat;
    for (const auto& R : Rset) Smat[R] = at(Sr, R, nw);

    const double twopi = 2.0 * M_PI;
    double maxerr = 0.0;
    for (int k = 0; k < g.num_kpts; ++k)
    {
        Eigen::MatrixXcd recon = Eigen::MatrixXcd::Zero(nw, nw);
        for (size_t b = 0; b < Rset.size(); ++b)
        {
            const auto& R = Rset[b];
            const double kr = g.kpt[k][0]*R[0] + g.kpt[k][1]*R[1] + g.kpt[k][2]*R[2];
            recon += wgt[b] * std::exp(complex<double>(0.0, +twopi*kr)) * Smat[R];
        }
        const Eigen::MatrixXcd V = g.V(k);
        const Eigen::MatrixXcd direct = V.adjoint() * g.Sb[k][alpha] * V;
        maxerr = std::max(maxerr, (recon - direct).cwiseAbs().maxCoeff());
    }
    std::cout << "spin_roundtrip_crosscheck: max|S_W(k)_fromR - V^dag S_B V| = " << maxerr << "\n";
    assert(maxerr < 1e-6);
    std::cout << "spin_roundtrip_crosscheck: PASS\n";
    return 0;
}
