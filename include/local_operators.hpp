/**
 * @file local_operators.hpp
 * @brief Local (on-site) orbital angular momentum generators and projection parsing.
 *
 * The orbital angular momentum operators are built symbolically in the complex
 * \f$|l,m\rangle\f$ basis and rotated to the Wannier90 real-harmonic ordering. Pure
 * s, p, and d shells are supported; hybrids and incomplete shells are rejected.
 */
#ifndef LOCAL_OPERATORS_HPP
#define LOCAL_OPERATORS_HPP

#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <Eigen/Dense>

/**
 * @brief Real-basis orbital angular momentum \f$L_x,L_y,L_z\f$ (units \f$\hbar\f$)
 *        for a complete shell of angular momentum l.
 *
 * l=0 (s) gives the trivial L=0; l=1 (p) and l=2 (d) are the non-trivial cases.
 * The matrices are built symbolically from \f$L_\pm\f$ in the \f$|l,m\rangle\f$
 * basis and rotated to the W90 real-harmonic order via C. Callers/tests assert
 * \f$L=L^\dagger\f$, \f$[L_x,L_y]=iL_z\f$, \f$\mathrm{Tr}\,L=0\f$, and integer
 * eigenvalues of \f$L_z\f$.
 *
 * W90 real-harmonic order (parameters.F90 / projections.tex):
 *   - p: pz, px, py
 *   - d: dz2, dxz, dyz, dx2-y2, dxy
 *
 * @param l angular momentum (0, 1, or 2)
 * @return array of three \f$2l+1\f$ x \f$2l+1\f$ complex matrices
 */
inline std::array<Eigen::MatrixXcd, 3> shell_L(int l)
{
    typedef std::complex<double> cd;
    const cd I(0.0, 1.0);
    const int n = 2 * l + 1;

    Eigen::MatrixXcd Lz = Eigen::MatrixXcd::Zero(n, n);
    Eigen::MatrixXcd Lp = Eigen::MatrixXcd::Zero(n, n);          // raising
    for (int idx = 0; idx < n; ++idx) { const int m = idx - l; Lz(idx, idx) = cd(m, 0); }
    for (int idx = 0; idx < n; ++idx)
    {
        const int m = idx - l;
        if (m < l) Lp(idx + 1, idx) = cd(std::sqrt(double(l)*(l+1) - double(m)*(m+1)), 0);
    }
    const Eigen::MatrixXcd Lm = Lp.adjoint();
    const Eigen::MatrixXcd Lx = 0.5 * (Lp + Lm);
    const Eigen::MatrixXcd Ly = (Lp - Lm) / cd(0.0, 2.0);

    // C: columns = W90 real harmonics expressed in |l,m> (Condon-Shortley).
    Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(n, n);
    const double s = 1.0 / std::sqrt(2.0);
    auto put = [&](int col, int m, cd coef) { C(m + l, col) += coef; };
    if (l == 0) {
        return { Eigen::MatrixXcd::Zero(1,1), Eigen::MatrixXcd::Zero(1,1), Eigen::MatrixXcd::Zero(1,1) };
    } else if (l == 1) {
        put(0, 0, cd(1,0));                       // pz = |0>
        put(1, -1, cd(s,0)); put(1, 1, cd(-s,0)); // px = (|-1> - |+1>)/√2
        put(2, -1, I*s);     put(2, 1, I*s);      // py = i(|-1> + |+1>)/√2
    } else if (l == 2) {
        put(0, 0, cd(1,0));                        // dz2
        put(1, -1, cd(s,0)); put(1, 1, cd(-s,0)); // dxz
        put(2, -1, I*s);     put(2, 1, I*s);      // dyz
        put(3, -2, cd(s,0)); put(3, 2, cd(s,0));  // dx2-y2 = (|-2> + |+2>)/√2
        put(4, -2, I*s);     put(4, 2, -I*s);     // dxy   = i(|-2> - |+2>)/√2
    } else {
        // only p,d in this cut; n!=3,5 never reached via parse_projection_shells
        return {Eigen::MatrixXcd::Zero(n,n), Eigen::MatrixXcd::Zero(n,n), Eigen::MatrixXcd::Zero(n,n)};
    }

    return { C.adjoint()*Lx*C, C.adjoint()*Ly*C, C.adjoint()*Lz*C };
}

/**
 * @brief Parse a `.win` "begin projections ... end projections" block.
 *
 * Returns the list of shell angular momenta (0 for s, 1 for p, 2 for d), in
 * projector-column order. Throws std::runtime_error on anything that is not a
 * complete pure s/p/d shell (hybrids, incomplete shells, f, ...).
 *
 * @param winfile path to `.win`
 * @return vector of angular momenta
 */
std::vector<int> parse_projection_shells(const std::string& winfile);

#endif
