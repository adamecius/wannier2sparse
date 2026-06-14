/**
 * @file gauge.hpp
 * @brief Gauge data and exact operators from Wannier90 / pw2wannier90 output.
 *
 * All conventions are documented and source-cited in docs/conventions.md; this
 * code follows that doc verbatim.
 */
#ifndef GAUGE_HPP
#define GAUGE_HPP

#include <array>
#include <string>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "hopping_list.hpp"

/**
 * @brief Gauge data read from a Wannier90 / pw2wannier90 run, used to build exact
 *        operators.
 *
 * The conventions are:
 * \f[
 *   V(k)    = U_{dis}(k) \, U(k)        \quad (\text{num_bands} \times \text{num_wann})
 * \f]
 * \f[
 *   O_W(k)  = V(k)^\dagger \, O_B(k) \, V(k)
 * \f]
 * \f[
 *   O_W(R)  = \frac{1}{N_k} \sum_k e^{-i 2\pi k\cdot R} O_W(k)
 * \f]
 * matching the Wannier90 Fourier transform convention.
 */
struct gauge_data
{
    int num_bands = 0, num_wann = 0, num_kpts = 0; ///< matrix/k-space dimensions
    bool disentangled = false;                       ///< true if Udis is present

    std::vector< std::array<double,3> >  kpt;   ///< fractional k-point coordinates, size num_kpts
    std::vector< Eigen::MatrixXcd >      U;      ///< per k: num_wann x num_wann
    std::vector< Eigen::MatrixXcd >      Udis;   ///< per k: num_bands x num_wann (empty if not disentangled)
    std::vector< std::array<Eigen::MatrixXcd,3> > Sb;  ///< per k: 3 Pauli components, num_bands x num_bands

    /**
     * @brief Composite gauge matrix V(k).
     * @param k k-point index
     * @return Udis[k]*U[k] if disentangled, else U[k]
     */
    Eigen::MatrixXcd V(int k) const { return disentangled ? (Udis[k] * U[k]) : U[k]; }
};

/**
 * @brief Read `_u.mat`, `_u_dis.mat` (optional), and `.spn` files.
 * @param prefix file prefix, e.g. "/path/Fe"
 * @return populated gauge_data
 */
gauge_data read_gauge(const std::string& prefix);

/**
 * @brief Exact spin operator \f$S_\alpha(R)\f$ (alpha = 0,1,2 for x,y,z).
 *
 * Values are in units of \f$\hbar/2\f$ (bare Pauli). The operator is evaluated on
 * the supplied R vectors, typically the Hamiltonian Wigner-Seitz set.
 *
 * @param g gauge data
 * @param alpha spin component index (0=x, 1=y, 2=z)
 * @param Rset set of R vectors
 * @return hopping_list representation of the spin operator
 */
hopping_list exact_spin_operator(const gauge_data& g, int alpha,
                                 const std::vector<hopping_list::cellID_t>& Rset);

/**
 * @brief Projection overlaps \f$A_{m,\alpha}(k)\f$ from a `.amn` file.
 *
 * Per k the matrix is num_bands x num_proj.
 */
struct amn_data
{
    int num_bands = 0, num_proj = 0, num_kpts = 0; ///< dimensions
    std::vector< Eigen::MatrixXcd > A;             ///< per k: num_bands x num_proj
};

/**
 * @brief Read a Wannier90 `.amn` file.
 * @param amn_file path to `.amn`
 * @return populated amn_data
 */
amn_data read_amn(const std::string& amn_file);

/**
 * @brief Orbital angular momentum \f$L_\alpha(R)\f$ via the projector route.
 *
 * Computes
 * \f[
 *   C(k) = A(k)^\dagger V(k), \quad
 *   L_W(k) = C(k)^\dagger L_{local}^\alpha C(k),
 * \f]
 * and Fourier transforms to real space. shell_ls gives the per-shell angular
 * momentum in projector-column order; sum(2l+1) must equal num_proj. Units are
 * \f$\hbar\f$.
 *
 * @param g gauge data
 * @param a projection overlaps
 * @param shell_ls per-shell angular momenta
 * @param alpha component index (0=x, 1=y, 2=z)
 * @param Rset set of R vectors
 * @return hopping_list representation of L_alpha
 */
hopping_list orbital_L_operator(const gauge_data& g, const amn_data& a,
                                const std::vector<int>& shell_ls, int alpha,
                                const std::vector<hopping_list::cellID_t>& Rset);

#endif
