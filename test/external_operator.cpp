#include <cassert>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <iostream>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
using namespace std;

typedef hopping_list::value_t  value_t;
typedef hopping_list::cellID_t cellID_t;
typedef hopping_list::edge_t   edge_t;
typedef hopping_list::hopping_t hopping_t;

static void write_hr(const string& fname, int num_wann,
                     const vector<int>& ndegen, const vector<string>& hoppings)
{
    ofstream f(fname.c_str());
    f << "hand-built hr.dat\n" << num_wann << "\n" << ndegen.size() << "\n";
    for (size_t i = 0; i < ndegen.size(); ++i)
    {
        f << ndegen[i];
        if ((i + 1) % 15 == 0) f << "\n"; else f << " ";
    }
    if (ndegen.size() % 15 != 0) f << "\n";
    for (const auto& h : hoppings) f << h << "\n";
    f.close();
}

int main()
{
    // A hand-built 2-band model over two R-points (Plan 2: ingest an arbitrary
    // operator in _hr.dat format through the same parse chain as H). Zero terms
    // in the second block must be dropped, leaving 5 nonzero hoppings.
    write_hr("ext_hr.dat", 2, {1, 1},
             { "0 0 0 1 1  0.5  0.0",
               "0 0 0 1 2  1.0 -0.2",
               "0 0 0 2 1  1.0  0.2",
               "0 0 0 2 2 -0.5  0.0",
               "1 0 0 1 1  0.3  0.0",
               "1 0 0 1 2  0.0  0.0",
               "1 0 0 2 1  0.0  0.0",
               "1 0 0 2 2  0.0  0.0" });

    hopping_list got = create_hopping_list(read_wannier_file("ext_hr.dat"));

    hopping_list want;
    want.num_wann = 2;
    want.hoppings.push_back(hopping_t(cellID_t({0,0,0}), value_t( 0.5,  0.0), edge_t({0,0})));
    want.hoppings.push_back(hopping_t(cellID_t({0,0,0}), value_t( 1.0, -0.2), edge_t({0,1})));
    want.hoppings.push_back(hopping_t(cellID_t({0,0,0}), value_t( 1.0,  0.2), edge_t({1,0})));
    want.hoppings.push_back(hopping_t(cellID_t({0,0,0}), value_t(-0.5,  0.0), edge_t({1,1})));
    want.hoppings.push_back(hopping_t(cellID_t({1,0,0}), value_t( 0.3,  0.0), edge_t({0,0})));

    assert(got.hoppings.size() == 5);
    assert(got == want);

    // The ingested operator must expand and export through the existing engine
    // exactly like a Hamiltonian (provenance-agnostic). Expand to a 2x1x1
    // supercell and confirm a well-formed CSR (dim = num_wann * 2).
    save_supercell_as_csr(cellID_t({2,1,1}), got, "ext.CSR");
    {
        ifstream f("ext.CSR");
        size_t dim = 0, nnz = 0; f >> dim >> nnz;
        assert(dim == 4);
        assert(nnz > 0);
    }

    std::cout << "EXTERNAL OPERATOR INGESTION TEST PASSED" << std::endl;
    return 0;
}
