#include <cassert>
#include <fstream>
#include <vector>
#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include "gauge.hpp"
#include "hopping_list.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;

// Extract the num_wann x num_wann matrix at a given R from a hopping_list.
static Eigen::MatrixXcd at(const hopping_list& hl, const cellID_t& R, int nw)
{
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nw, nw);
    for (const auto& h : hl.hoppings)
        if (get<0>(h) == R)
        {
            const auto e = get<2>(h);
            M(e[0], e[1]) += get<1>(h);
        }
    return M;
}

int main()
{
    const complex<double> I(0, 1);

    // ---- engine: S_W(R) must be Hermitian-conjugate across R <-> -R --------
    // 2-orbital, non-entangled (V = U = I), 4-point 1D k-mesh, a k-dependent
    // Hermitian S_B (non-diagonal) so the check is non-trivial.
    {
        gauge_data g;
        g.num_bands = 2; g.num_wann = 2; g.num_kpts = 4; g.disentangled = false;
        g.kpt.resize(4); g.U.resize(4); g.Sb.resize(4);
        for (int k = 0; k < 4; ++k)
        {
            g.kpt[k] = {0.25 * k, 0.0, 0.0};
            g.U[k] = Eigen::MatrixXcd::Identity(2, 2);
            const complex<double> a = 0.5 * std::exp(I * (2.0 * M_PI * k / 4.0));
            for (int s = 0; s < 3; ++s) g.Sb[k][s] = Eigen::MatrixXcd::Zero(2, 2);
            // use the x-component as a k-dependent Hermitian matrix
            g.Sb[k][0](0,0) = 1.0; g.Sb[k][0](1,1) = -1.0;
            g.Sb[k][0](0,1) = a;   g.Sb[k][0](1,0) = conj(a);
        }
        std::vector<cellID_t> Rset = { {0,0,0},{1,0,0},{-1,0,0},{2,0,0},{-2,0,0} };
        hopping_list S = exact_spin_operator(g, 0, Rset);

        // S(R) == S(-R)^dagger
        assert((at(S, cellID_t({1,0,0}),  2) - at(S, cellID_t({-1,0,0}), 2).adjoint()).cwiseAbs().maxCoeff() < 1e-9);
        assert((at(S, cellID_t({2,0,0}),  2) - at(S, cellID_t({-2,0,0}), 2).adjoint()).cwiseAbs().maxCoeff() < 1e-9);
        // S(R=0) Hermitian
        Eigen::MatrixXcd S0 = at(S, cellID_t({0,0,0}), 2);
        assert((S0 - S0.adjoint()).cwiseAbs().maxCoeff() < 1e-9);
        // non-trivial (engine actually produced something)
        assert(S.hoppings.size() > 0);
    }

    // ---- readers: parse a tiny _u.mat + .spn (no _u_dis.mat -> V = U) ------
    {
        { ofstream f("syn_u.mat");
          f << "header\n"
               "    1    2    2\n"
               "\n"
               "  0.0000000000  +0.0000000000  +0.0000000000\n"
               // u(i,j), i fast, j slow: identity
               "  1.0000000000  +0.0000000000\n"   // (0,0)
               "  0.0000000000  +0.0000000000\n"   // (1,0)
               "  0.0000000000  +0.0000000000\n"   // (0,1)
               "  1.0000000000  +0.0000000000\n"; } // (1,1)
        { ofstream f("syn.spn");
          f << "header\n"
               "    2    1\n"
               // counter (m=0,n=0): sx sy sz
               "  0.0E+00  0.0E+00\n  0.0E+00  0.0E+00\n  1.0E+00  0.0E+00\n"
               // counter (m=1,n=0): sx sy sz
               "  0.0E+00  0.0E+00\n  0.0E+00  0.0E+00\n  0.0E+00  0.0E+00\n"
               // counter (m=1,n=1): sx sy sz
               "  0.0E+00  0.0E+00\n  0.0E+00  0.0E+00\n -1.0E+00  0.0E+00\n"; }

        gauge_data g = read_gauge("syn");
        assert(g.num_wann == 2 && g.num_bands == 2 && g.num_kpts == 1 && !g.disentangled);
        assert((g.U[0] - Eigen::MatrixXcd::Identity(2,2)).cwiseAbs().maxCoeff() < 1e-12);
        // Sb_z = diag(1,-1), reconstructed Hermitian
        assert(abs(g.Sb[0][2](0,0) - complex<double>( 1,0)) < 1e-12);
        assert(abs(g.Sb[0][2](1,1) - complex<double>(-1,0)) < 1e-12);
        assert((g.Sb[0][2] - g.Sb[0][2].adjoint()).cwiseAbs().maxCoeff() < 1e-12);
    }

    std::cout << "GAUGE SPIN TEST PASSED" << std::endl;
    return 0;
}
