/**
 * @file sparse_matrix.hpp
 * @brief Sparse matrix aliases used throughout wannier2sparse.
 *
 * All CSR output and internal sparse algebra is built on Eigen. The matrix is
 * stored row-major so the CSR text layout (values / column indices / row pointers)
 * matches the file format directly. The FULL Hermitian matrix is written (both
 * triangles); the consumer folds it as H = 0.5*(M + M^dagger), an identity for a
 * Hermitian operator.
 */
#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <Eigen/Sparse>
#include <complex>

typedef Eigen::SparseMatrix< std::complex<double>, Eigen::RowMajor > SparseMatrix_t; ///< Row-major complex sparse matrix.
typedef Eigen::Triplet< std::complex<double> > Triplet_t;                           ///< (row, col, value) triplet for assembly.


#endif
