/**
 * @file operator_algebra.hpp
 * @brief Sparse-matrix operator algebra helpers.
 */
#ifndef OPERATOR_ALGEBRA
#define OPERATOR_ALGEBRA

#include "sparse_matrix.hpp"

/**
 * @brief Symmetrized product (anticommutator) of two sparse operators.
 *
 * Computes
 * \f[
 *     \{A,B\}/2 = \frac12 (AB + BA)
 * \f]
 * and returns the resulting sparse matrix. Used to form the conventional spin
 * current \f$J = \frac12 (VS + SV)\f$ from already-expanded velocity and spin
 * operators. If A and B are Hermitian the result is Hermitian. Fill-in is bounded
 * by the combined range of A and B.
 *
 * @param A first operator
 * @param B second operator
 * @return \f$\frac12(AB + BA)\f$
 */
inline SparseMatrix_t anticommutator(const SparseMatrix_t& A, const SparseMatrix_t& B)
{
    SparseMatrix_t J = (A * B + B * A);
    J *= 0.5;
    J.makeCompressed();
    return J;
}

#endif
