#ifndef LOCAL_OPERATORS_HPP
#define LOCAL_OPERATORS_HPP

#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <Eigen/Dense>

/**
 * Real-basis orbital angular momentum Lx, Ly, Lz (units ħ) for a complete shell
 * of angular momentum l. l=0 (s) gives the trivial L=0; l=1 (p) and l=2 (d) are
 * the non-trivial cases (see docs/conventions.md sec 5). Built symbolically from
 * L± in the |l,m> basis,
 * rotated to the W90 real-harmonic order via C (never transcribed); callers/tests
 * assert L=L†, [Lx,Ly]=iLz, Tr=0, eig(Lz)={-1,0,1}/{-2..2}.
 *
 * W90 real-harmonic order (parameters.F90 / projections.tex):
 *   p: pz, px, py            d: dz2, dxz, dyz, dx2-y2, dxy
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

// Parse a .win "begin projections ... end projections" block into the list of
// shell angular momenta (0 for s, 1 for p, 2 for d), in projector-column order.
// Throws a std::runtime_error naming the offending token on anything that is not
// a complete pure s/p/d shell (hybrids, incomplete shells, f, ...).
std::vector<int> parse_projection_shells(const std::string& winfile);

#endif
