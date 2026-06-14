// Algebra check on the ASSEMBLED spin-density operators (Phase 2c).
//
// checks::spin_algebra() only validates the static 2x2 Pauli matrices. This test
// validates the operators that createHoppingSpinDensity_list actually builds (support
// from the spin doubling) on a collinear doubled model: each S_alpha must be
// Hermitian, traceless, involutory (S_alpha^2 = I on the fully-paired basis), and
// satisfy [S_x, S_y] = 2 i S_z. A model whose onsite block has same-spin orbital bonds
// would have failed this before the support-from-doubling fix.
#include <cassert>
#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include "hopping_list.hpp"
#include "tbmodel.hpp"
using namespace std;

int main()
{
    tbmodel model;
    model.readOrbitalPositions("spin_algebra_chain.xyz");
    model.readUnitCell("spin_algebra_chain.uc");
    model.readWannierModel("spin_algebra_chain_hr.dat");

    const hopping_list::cellID_t dim({4, 1, 1});
    const Eigen::MatrixXcd X(supercell_matrix(dim, model.createHoppingSpinDensity_list('x')));
    const Eigen::MatrixXcd Y(supercell_matrix(dim, model.createHoppingSpinDensity_list('y')));
    const Eigen::MatrixXcd Z(supercell_matrix(dim, model.createHoppingSpinDensity_list('z')));

    const int n = (int)X.rows();
    const Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(n, n);
    const complex<double> i_(0, 1);
    const double tol = 1e-12;
    auto mx = [](const Eigen::MatrixXcd& M){ return M.cwiseAbs().maxCoeff(); };

    const double hx = mx(X - X.adjoint()), hy = mx(Y - Y.adjoint()), hz = mx(Z - Z.adjoint());
    cout << "Hermiticity residual: Sx " << hx << " Sy " << hy << " Sz " << hz << endl;
    cout << "trace: Sx " << abs(X.trace()) << " Sy " << abs(Y.trace()) << " Sz " << abs(Z.trace()) << endl;
    cout << "involutory residual: Sx^2 " << mx(X*X - Id) << " Sy^2 " << mx(Y*Y - Id) << " Sz^2 " << mx(Z*Z - Id) << endl;
    cout << "[Sx,Sy]-2i Sz residual: " << mx(X*Y - Y*X - 2.0*i_*Z) << endl;

    assert(hx < tol && hy < tol && hz < tol);                                  // Hermitian
    assert(abs(X.trace()) < tol && abs(Y.trace()) < tol && abs(Z.trace()) < tol); // traceless
    assert(mx(X*X - Id) < tol && mx(Y*Y - Id) < tol && mx(Z*Z - Id) < tol);    // S_a^2 = I
    assert(mx(X*Y - Y*X - 2.0*i_*Z) < tol);                                    // [Sx,Sy] = 2i Sz

    cout << "PASS: assembled spin operators satisfy the su(2) algebra" << endl;
    return 0;
}
