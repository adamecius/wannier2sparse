#include <cassert>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
#include "descriptor.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;

int main()
{
    // Analytic check: a 1D nearest-neighbor ring H(k) = 2t cos(k). For an even
    // number of sites the allowed k include 0 and pi, so the finite spectrum hits
    // the band edges exactly: [a,b] = [-2t, 2t]. Use t = 1, N = 8.
    {
        ofstream("ring_hr.dat")
            << "ring\n1\n3\n1 1 1\n"
               "-1 0 0 1 1 1.0 0.0\n"
               " 0 0 0 1 1 0.0 0.0\n"
               " 1 0 0 1 1 1.0 0.0\n";
        hopping_list hl = create_hopping_list(read_wannier_file("ring_hr.dat"));
        SparseMatrix_t H = supercell_matrix(cellID_t({8,1,1}), hl);

        double a = 0.0, b = 0.0;
        spectral_bounds(H, a, b);
        std::cout << "ring bounds [a,b] = [" << a << ", " << b << "]\n";
        assert(std::fabs(a - (-2.0)) < 1e-9);
        assert(std::fabs(b -  ( 2.0)) < 1e-9);
    }

    // read_eig_bounds: standard W90 .eig format "band kpt energy".
    {
        ofstream("model.eig")
            << "1 1 -3.5\n"
               "2 1  1.2\n"
               "1 2 -0.7\n"
               "2 2  4.0\n";
        double emin = 0.0, emax = 0.0;
        assert(read_eig_bounds("model.eig", emin, emax));
        assert(std::fabs(emin - (-3.5)) < 1e-12);
        assert(std::fabs(emax -  ( 4.0)) < 1e-12);
        // Missing file -> false, outputs untouched.
        double x = 7.0, y = 9.0;
        assert(!read_eig_bounds("does_not_exist.eig", x, y));
        assert(x == 7.0 && y == 9.0);
    }

    // Descriptor sidecar round-trip: written fields are present and parseable.
    {
        OperatorDescriptor d;
        d.observable = "hamiltonian"; d.units = "eV"; d.provenance = "test";
        d.has_bounds = true; d.a = -2.0; d.b = 2.0;
        write_descriptor(d, "H.desc");

        ifstream f("H.desc");
        string key; bool saw_min = false, saw_max = false; double v = 0.0;
        while (f >> key)
        {
            if (key == "spectral_min:") { f >> v; saw_min = (std::fabs(v + 2.0) < 1e-12); }
            else if (key == "spectral_max:") { f >> v; saw_max = (std::fabs(v - 2.0) < 1e-12); }
            else { string rest; getline(f, rest); }
        }
        assert(saw_min && saw_max);
    }

    std::cout << "SPECTRAL BOUNDS / DESCRIPTOR TEST PASSED" << std::endl;
    return 0;
}
