#ifndef OPERATOR_ALGEBRA
#define OPERATOR_ALGEBRA

#include "sparse_matrix.hpp"

/**
 * Symmetrized product (anticommutator) of two sparse operators:
 *
 *     {A,B}/2 = 1/2 (A*B + B*A)
 *
 * Used to form the conventional spin current J = 1/2 (V*S + S*V) from the
 * already-expanded velocity and spin operators. If A and B are Hermitian the
 * result is Hermitian. Fill-in is bounded by the combined range of A and B.
 */
inline SparseMatrix_t anticommutator(const SparseMatrix_t& A, const SparseMatrix_t& B)
{
    SparseMatrix_t J = (A * B + B * A);
    J *= 0.5;
    J.makeCompressed();
    return J;
}

#endif
