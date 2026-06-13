#ifndef GAUGE_HPP
#define GAUGE_HPP

#include <array>
#include <string>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "hopping_list.hpp"

/**
 * Gauge data read from a Wannier90 / pw2wannier90 run, used to build exact
 * operators (Plan 7). All conventions are documented and source-cited in
 * docs/conventions.md; this code follows that doc verbatim.
 *
 *   V(k)    = U_dis(k) * U(k)        (num_bands x num_wann), or U(k) if not disentangled
 *   O_W(k)  = V(k)^dagger O_B(k) V(k)   (num_wann x num_wann)
 *   O_W(R)  = (1/N_k) sum_k exp(-i 2pi k.R) O_W(k)     [same FT as W90 H(R)]
 */
struct gauge_data
{
    int num_bands = 0, num_wann = 0, num_kpts = 0;
    bool disentangled = false;

    std::vector< std::array<double,3> >  kpt;   // fractional kpt_latt, size num_kpts
    std::vector< Eigen::MatrixXcd >      U;      // per k: num_wann x num_wann
    std::vector< Eigen::MatrixXcd >      Udis;   // per k: num_bands x num_wann (empty if not disentangled)
    std::vector< std::array<Eigen::MatrixXcd,3> > Sb;  // per k: 3 Pauli components, num_bands x num_bands

    // V(k) per the disentanglement convention (docs/conventions.md sec 2).
    Eigen::MatrixXcd V(int k) const { return disentangled ? (Udis[k] * U[k]) : U[k]; }
};

// Read <prefix>_u.mat, <prefix>_u_dis.mat (optional), <prefix>.spn (formatted).
// Asserts per-k dimensional consistency. <prefix> is e.g. "/path/Fe".
gauge_data read_gauge(const std::string& prefix);

// Exact spin operator S_alpha(R) (alpha = 0,1,2 for x,y,z) as a hopping_list,
// evaluated on the given set of R vectors (typically the H(R) WS set). Values
// are in units of hbar/2 (bare Pauli; see docs/conventions.md sec 3).
hopping_list exact_spin_operator(const gauge_data& g, int alpha,
                                 const std::vector<hopping_list::cellID_t>& Rset);

#endif
