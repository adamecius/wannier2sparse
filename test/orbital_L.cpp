#include <cassert>
#include <fstream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <iostream>
#include <Eigen/Dense>
#include "local_operators.hpp"
#include "gauge.hpp"
#include "hopping_list.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;
static const complex<double> I(0,1);

static Eigen::MatrixXcd at(const hopping_list& hl, const cellID_t& R, int nw)
{
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nw, nw);
    for (const auto& h : hl.hoppings)
        if (get<0>(h) == R) { const auto e = get<2>(h); M(e[0], e[1]) += get<1>(h); }
    return M;
}

static void check_shell(int l, const vector<double>& expect_lz)
{
    auto L = shell_L(l);                       // {Lx, Ly, Lz}
    const int n = 2*l + 1;
    for (int a = 0; a < 3; ++a)                // Hermitian
        assert((L[a] - L[a].adjoint()).cwiseAbs().maxCoeff() < 1e-12);
    // [Lx, Ly] = i Lz  (catches ordering/phase errors)
    assert(((L[0]*L[1] - L[1]*L[0]) - I*L[2]).cwiseAbs().maxCoeff() < 1e-12);
    for (int a = 0; a < 3; ++a)                // traceless
        assert(abs(L[a].trace()) < 1e-12);
    // eig(Lz) = expect_lz
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(L[2]);
    for (int i = 0; i < n; ++i)
        assert(abs(es.eigenvalues()(i) - expect_lz[i]) < 1e-9);
}

static vector<int> shells_of(const string& proj)
{
    ofstream("t.win") << "begin projections\n" << proj << "\nend projections\n";
    return parse_projection_shells("t.win");
}

int main()
{
    // 1) L algebra for p and d
    check_shell(1, {-1, 0, 1});
    check_shell(2, {-2, -1, 0, 1, 2});

    // 2) projection parsing: pure shells accepted, everything else errors by name
    assert((shells_of("X: s")   == vector<int>{0}));    // s -> trivial L=0 (allowed)
    assert((shells_of("X: p")   == vector<int>{1}));
    assert((shells_of("X: d")   == vector<int>{2}));
    assert((shells_of("X: d;s;s") == vector<int>{2,0,0}));   // copper-like
    assert((shells_of("X: p;d") == vector<int>{1,2}));
    for (const char* bad : {"X: sp3d2", "X: sp3", "X: dxy", "X: f"})
    {
        bool threw = false;
        try { shells_of(bad); } catch (const std::exception&) { threw = true; }
        assert(threw && "hybrid/incomplete/unsupported projection must throw");
    }

    // 3) projector route: V=I, A=I, one p shell -> L(R=0) == L_local; Hermitian across +-R
    {
        gauge_data g;
        g.num_bands = 3; g.num_wann = 3; g.num_kpts = 2; g.disentangled = false;
        g.kpt = {{0.0,0.0,0.0}, {0.5,0.0,0.0}};
        g.U = {Eigen::MatrixXcd::Identity(3,3), Eigen::MatrixXcd::Identity(3,3)};
        amn_data A; A.num_bands = 3; A.num_proj = 3; A.num_kpts = 2;
        A.A = {Eigen::MatrixXcd::Identity(3,3), Eigen::MatrixXcd::Identity(3,3)};
        vector<int> shells = {1};
        vector<cellID_t> Rset = {{0,0,0},{1,0,0},{-1,0,0}};

        auto Lp = shell_L(1);
        for (int a = 0; a < 3; ++a)
        {
            hopping_list L = orbital_L_operator(g, A, shells, a, Rset);
            assert((at(L, cellID_t({0,0,0}), 3) - Lp[a]).cwiseAbs().maxCoeff() < 1e-9);
            assert((at(L, cellID_t({1,0,0}), 3) - at(L, cellID_t({-1,0,0}), 3).adjoint()).cwiseAbs().maxCoeff() < 1e-9);
        }
        // shell/num_proj mismatch must throw
        bool threw = false;
        try { orbital_L_operator(g, A, vector<int>{2}, 0, Rset); } catch (const std::exception&) { threw = true; }
        assert(threw);
    }

    std::cout << "ORBITAL L TEST PASSED" << std::endl;
    return 0;
}
