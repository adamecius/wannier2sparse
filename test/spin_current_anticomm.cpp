#include <cassert>
#include <fstream>
#include <complex>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
#include "tbmodel.hpp"
#include "operator_algebra.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;

static double max_abs(const SparseMatrix_t& a, const SparseMatrix_t& b)
{
    Eigen::MatrixXcd d = Eigen::MatrixXcd(a) - Eigen::MatrixXcd(b);
    return d.cwiseAbs().maxCoeff();
}

static SparseMatrix_t mat2(complex<double> a, complex<double> b,
                           complex<double> c, complex<double> d)
{
    std::vector<Triplet_t> t;
    t.push_back(Triplet_t(0,0,a)); t.push_back(Triplet_t(0,1,b));
    t.push_back(Triplet_t(1,0,c)); t.push_back(Triplet_t(1,1,d));
    SparseMatrix_t m(2,2); m.setFromTriplets(t.begin(),t.end()); m.makeCompressed();
    return m;
}

int main()
{
    const complex<double> I(0,1), Z(0,0), U(1,0);

    // --- anticommutator algebra: known Pauli identities -------------------
    const SparseMatrix_t sx = mat2(Z, U, U, Z);
    const SparseMatrix_t sy = mat2(Z,-I, I, Z);
    const SparseMatrix_t sz = mat2(U, Z, Z,-U);
    const SparseMatrix_t id = mat2(U, Z, Z, U);
    const SparseMatrix_t zero2(2,2);

    assert(max_abs(anticommutator(sx, sz), zero2) < 1e-12);   // {sx,sz} = 0
    assert(max_abs(anticommutator(sx, sy), zero2) < 1e-12);   // {sx,sy} = 0
    assert(max_abs(anticommutator(sx, sx), id)    < 1e-12);   // 1/2{sx,sx} = I
    assert(max_abs(anticommutator(sx, sz),
                   anticommutator(sz, sx))         < 1e-12);   // symmetric in A,B

    // Hermiticity of 1/2{A,B} for Hermitian A,B.
    const SparseMatrix_t A = mat2(U, I, -I, complex<double>(2,0)); // [[1,i],[-i,2]]
    const SparseMatrix_t Jab = anticommutator(A, sx);
    assert(max_abs(Jab, SparseMatrix_t(Jab.adjoint())) < 1e-12);

    // --- supercell_matrix: assemble a known 1D nearest-neighbor ring ------
    {
        ofstream f("ring_hr.dat");
        f << "ring\n1\n3\n1 1 1\n"
             "-1 0 0 1 1 1.0 0.0\n"
             " 0 0 0 1 1 0.0 0.0\n"
             " 1 0 0 1 1 1.0 0.0\n";
    }
    {
        hopping_list hl = create_hopping_list(read_wannier_file("ring_hr.dat"));
        const SparseMatrix_t M = supercell_matrix(cellID_t({3,1,1}), hl);
        Eigen::MatrixXcd D = Eigen::MatrixXcd(M);
        Eigen::MatrixXcd expect(3,3);
        expect << 0,1,1, 1,0,1, 1,1,0;          // 3-site ring adjacency
        assert((D - expect).cwiseAbs().maxCoeff() < 1e-12);
    }

    // --- full pipeline: J = 1/2{V,S} on a model is Hermitian and non-trivial
    {
        { ofstream f("col.xyz");
          f << "4\nA_s+_ 0.0 0.0 0.0\nB_s+_ 0.5 0.0 0.0\nA_s-_ 0.0 0.0 0.0\nB_s-_ 0.5 0.0 0.0\n"; }
        { ofstream f("col.uc"); f << "1 0 0\n0 1 0\n0 0 1\n"; }
        { const int R[3] = {-1,0,1}; const double t = 1.0;
          ofstream f("col_hr.dat"); f << "collinear spin chain\n4\n3\n1 1 1\n";
          for (int r=0;r<3;++r) for (int i=1;i<=4;++i) for (int j=1;j<=4;++j) {
              double v=0.0;
              if (R[r]== 0 && ((i==1&&j==2)||(i==2&&j==1)||(i==3&&j==4)||(i==4&&j==3))) v=t;
              if (R[r]== 1 && ((i==2&&j==1)||(i==4&&j==3))) v=t;
              if (R[r]==-1 && ((i==1&&j==2)||(i==3&&j==4))) v=t;
              f << R[r] << " 0 0 " << i << " " << j << " " << v << " 0.0\n";
          } }

        tbmodel model;
        model.readOrbitalPositions("col.xyz");
        model.readUnitCell("col.uc");
        model.readWannierModel("col_hr.dat");

        const cellID_t dim({3,1,1});
        const SparseMatrix_t V = supercell_matrix(dim, model.createHoppingCurrents_list(0));
        const SparseMatrix_t S = supercell_matrix(dim, model.createHoppingSpinDensity_list('z'));
        const SparseMatrix_t J = anticommutator(V, S);

        assert(J.nonZeros() > 0);
        assert(max_abs(J, SparseMatrix_t(J.adjoint())) < 1e-9);   // Hermitian
    }

    std::cout << "SPIN CURRENT ANTICOMMUTATOR TEST PASSED" << std::endl;
    return 0;
}
