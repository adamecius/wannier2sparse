/**
 * @file checks.hpp
 * @brief Self-verification invariants for exported operators.
 *
 * Each check returns a pass/fail result plus a numeric residual. They operate on an
 * assembled operator's hopping_list O(R) (gauge- and supercell-independent) or on
 * the local generators. The check accompanies the operator as a sidecar; it never
 * alters the CSR.
 */
#ifndef CHECKS_HPP
#define CHECKS_HPP

#include <string>
#include <vector>
#include <map>
#include <array>
#include <complex>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include "hopping_list.hpp"
#include "local_operators.hpp"

/**
 * @brief Result of a single self-check.
 */
struct check_result
{
    std::string name;     ///< check name
    bool        pass;     ///< true if the residual is within tolerance
    double      residual; ///< numerical residual
    std::string detail;   ///< human-readable description
};

/**
 * @brief Self-verification invariants.
 */
namespace checks {

/**
 * @brief Reduce a hopping_list to a map keyed by (R, i, j).
 * @param hl hopping list
 * @return map from \f$(R_x,R_y,R_z,i,j)\f$ to summed complex value
 */
inline std::map<std::array<int,5>, std::complex<double> > reduce(const hopping_list& hl)
{
    std::map<std::array<int,5>, std::complex<double> > m;
    for (const auto& h : hl.hoppings)
    {
        const auto R = std::get<0>(h); const auto e = std::get<2>(h);
        m[{R[0],R[1],R[2],e[0],e[1]}] += std::get<1>(h);
    }
    return m;
}

/**
 * @brief Hermiticity check: \f$O_{ij}(R) = O_{ji}(-R)^*\f$.
 * @param hl hopping list
 * @param tol tolerance (default 1e-8)
 * @return check result
 */
inline check_result hermiticity(const hopping_list& hl, double tol = 1e-8)
{
    const auto m = reduce(hl);
    double res = 0.0;
    for (const auto& kv : m)
    {
        const auto& k = kv.first;
        const std::array<int,5> mirror = {-k[0],-k[1],-k[2], k[4], k[3]};
        auto it = m.find(mirror);
        const std::complex<double> other = (it == m.end()) ? std::complex<double>(0,0) : it->second;
        res = std::max(res, std::abs(kv.second - std::conj(other)));
    }
    return { "hermiticity", res <= tol, res, "max|O_ij(R)-conj(O_ji(-R))|" };
}

/**
 * @brief Trace of \f$O(R=0)\f$.
 *
 * For spin/orbital operators the trace should vanish (traceless). When
 * expect_traceless is false the check is informational (always passes).
 *
 * @param hl hopping list
 * @param expect_traceless true for spin/orbital operators
 * @param tol tolerance (default 1e-8)
 * @return check result
 */
inline check_result trace_rule(const hopping_list& hl, bool expect_traceless, double tol = 1e-8)
{
    const auto m = reduce(hl);
    std::complex<double> tr(0,0);
    for (const auto& kv : m)
        if (kv.first[0]==0 && kv.first[1]==0 && kv.first[2]==0 && kv.first[3]==kv.first[4])
            tr += kv.second;
    const double res = std::abs(tr);
    return { "sum_rules", expect_traceless ? (res <= tol) : true, res,
             "Tr O(R=0)" + std::string(expect_traceless ? " (expect 0)" : "") };
}

/**
 * @brief Minimum-image aliasing check.
 *
 * Detects distinct R vectors that collapse onto the same supercell bond under
 * periodic wrapping. Empty result means the supercell is safe.
 *
 * @param hl hopping list
 * @param cellDim supercell dimensions
 * @return check result
 */
inline check_result aliasing(const hopping_list& hl, const hopping_list::cellID_t& cellDim)
{
    const auto c = ::aliasing_collisions(hl, cellDim);
    return { "aliasing", c.empty(), double(c.size()),
             c.empty() ? "no minimum-image collisions" : (std::to_string(c.size()) + " colliding pairs: " + c.front()) };
}

/**
 * @brief Orbital angular-momentum algebra check.
 *
 * Verifies \f$[L_x,L_y]=iL_z\f$ and integer eigenvalues \f$\{-l..l\}\f$ for p and d
 * shells.
 *
 * @return check result
 */
inline check_result orbital_algebra()
{
    const std::complex<double> I(0,1);
    double res = 0.0;
    for (int l : {1, 2})
    {
        auto L = shell_L(l);
        res = std::max(res, ((L[0]*L[1]-L[1]*L[0]) - I*L[2]).cwiseAbs().maxCoeff());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(L[2]);
        for (int i = 0; i < 2*l+1; ++i) res = std::max(res, std::abs(es.eigenvalues()(i) - double(i - l)));
    }
    return { "algebra", res <= 1e-9, res, "[Lx,Ly]=iLz and eig(Lz)={-l..l} (p,d)" };
}

/**
 * @brief Spin-\f$\frac12\f$ algebra check.
 *
 * Verifies \f$[S_x,S_y]=iS_z\f$ and eigenvalues \f$\pm\frac12\f$.
 *
 * @return check result
 */
inline check_result spin_algebra()
{
    typedef std::complex<double> cd; const cd I(0,1);
    Eigen::MatrixXcd sx(2,2), sy(2,2), sz(2,2);
    sx << 0,1, 1,0;  sy << 0,-I, I,0;  sz << 1,0, 0,-1;
    const Eigen::MatrixXcd Sx=0.5*sx, Sy=0.5*sy, Sz=0.5*sz;   // su(2): eig=±1/2
    double res = ((Sx*Sy-Sy*Sx) - I*Sz).cwiseAbs().maxCoeff();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Sz);
    res = std::max(res, std::abs(es.eigenvalues()(0)+0.5));
    res = std::max(res, std::abs(es.eigenvalues()(1)-0.5));
    return { "algebra", res <= 1e-9, res, "[Sx,Sy]=iSz and eig(Sz)={-1/2,1/2}" };
}

/**
 * @brief Write a collection of check results to a sidecar file.
 * @param rs vector of results
 * @param filename output path
 */
inline void write_report(const std::vector<check_result>& rs, const std::string& filename)
{
    std::ofstream f(filename.c_str());
    f.precision(6);
    for (const auto& r : rs)
        f << r.name << ": " << (r.pass ? "PASS" : "FAIL")
          << "  residual=" << std::scientific << r.residual << "  [" << r.detail << "]\n";
    f.close();
}

} // namespace checks
#endif
