#ifndef DESCRIPTOR
#define DESCRIPTOR

#include <string>
#include <fstream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>
#include <limits>
#include <Eigen/Dense>
#include "sparse_matrix.hpp"

/**
 * Physical metadata "sidecar" for an exported operator.
 *
 * Invariant: this descriptor lives BESIDE the numerical CSR and never enters the
 * Chebyshev/KPM recursion. It carries provenance and, for the Hamiltonian, the
 * spectral bounds (a,b) used to rescale H into [-1,1] downstream.
 */
struct OperatorDescriptor
{
    std::string observable;   // "hamiltonian", "velocity", "spin", "spin_current", ...
    std::string component;    // "", "X", "Z", "XSZ", ...
    std::string units;        // "eV", "eV*Angstrom", "hbar/2", "dimensionless", ...
    std::string provenance;   // how it was produced
    bool        has_bounds;   // (a,b) meaningful (typically only for H)
    double      a, b;         // spectral_min, spectral_max

    OperatorDescriptor() : has_bounds(false), a(0.0), b(0.0) {}
};

inline void write_descriptor(const OperatorDescriptor& d, const std::string& filename)
{
    std::ofstream f(filename.c_str());
    f.precision(std::numeric_limits<double>::digits10 + 2);
    f << "observable: "  << d.observable  << "\n";
    f << "component: "   << d.component   << "\n";
    f << "units: "       << d.units       << "\n";
    f << "provenance: "  << d.provenance  << "\n";
    if (d.has_bounds)
    {
        f << "spectral_min: " << d.a << "\n";
        f << "spectral_max: " << d.b << "\n";
    }
    f.close();
}

/**
 * Estimate the extremal eigenvalues (a = min, b = max) of a Hermitian sparse
 * operator with `iters` Lanczos steps. The extreme Ritz values of the small
 * tridiagonal matrix converge to the spectral edges, which is exactly what KPM
 * needs to rescale H. The start vector is fixed-seed random so the result is
 * deterministic and not accidentally orthogonal to an edge eigenvector.
 */
inline void spectral_bounds(const SparseMatrix_t& H, double& a, double& b, int iters = 40)
{
    const int n = static_cast<int>(H.rows());
    if (n == 0) { a = b = 0.0; return; }
    iters = std::min(iters, n);

    std::mt19937 rng(12345u);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    Eigen::VectorXcd v(n);
    for (int i = 0; i < n; ++i) v(i) = std::complex<double>(u(rng), u(rng));
    v.normalize();

    Eigen::VectorXcd vprev = Eigen::VectorXcd::Zero(n);
    std::vector<double> alpha, beta;
    double beta_prev = 0.0;
    for (int j = 0; j < iters; ++j)
    {
        Eigen::VectorXcd w = H * v - beta_prev * vprev;
        const double aj = std::real(v.dot(w));    // v^H w
        w -= aj * v;
        alpha.push_back(aj);
        const double bj = w.norm();
        if (bj < 1e-12) break;                    // invariant subspace exhausted
        vprev = v;
        v = w / bj;
        beta_prev = bj;
        beta.push_back(bj);
    }

    const int m = static_cast<int>(alpha.size());
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m, m);
    for (int i = 0; i < m; ++i)     T(i, i) = alpha[i];
    for (int i = 0; i + 1 < m; ++i) { T(i, i + 1) = beta[i]; T(i + 1, i) = beta[i]; }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
    a = es.eigenvalues().minCoeff();
    b = es.eigenvalues().maxCoeff();
}

#endif
