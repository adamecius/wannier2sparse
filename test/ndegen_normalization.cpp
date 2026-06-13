#include <cassert>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <iostream>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
using namespace std;

// Write a minimal Wannier90 _hr.dat: comment, num_wann, nrpts, the ndegen block
// (15 integers per line), then the hopping lines verbatim.
static void write_hr(const string& fname, int num_wann,
                     const vector<int>& ndegen, const vector<string>& hoppings)
{
    ofstream f(fname.c_str());
    f << "synthetic hr.dat\n";
    f << num_wann << "\n";
    f << ndegen.size() << "\n";
    for (size_t i = 0; i < ndegen.size(); ++i)
    {
        f << ndegen[i];
        if ((i + 1) % 15 == 0) f << "\n"; else f << " ";
    }
    if (ndegen.size() % 15 != 0) f << "\n";   // close a partial line (no spurious blank line)
    for (const auto& h : hoppings) f << h << "\n";
    f.close();
}

int main()
{
    // (2) A block with ndegen = 2 must halve the hopping value.
    {
        write_hr("t2_hr.dat", 1, {2}, { "0 0 0 1 1 4.0 0.0" });
        hopping_list hl = create_hopping_list(read_wannier_file("t2_hr.dat"));
        assert(hl.hoppings.size() == 1);
        assert(abs(get<1>(hl.hoppings[0]) - hopping_list::value_t(2.0, 0.0)) < 1e-12);
    }

    // (3) Block index is taken from RAW line order: a filtered ~0 hopping in an
    //     earlier block must not shift later blocks. num_wann=1, ndegen=[1,4];
    //     block0 is zero (filtered), block1 value 8 -> /4 = 2 (not /1 = 8).
    {
        write_hr("t3_hr.dat", 1, {1, 4},
                 { "0 0 0 1 1 0.0 0.0", "1 0 0 1 1 8.0 0.0" });
        hopping_list hl = create_hopping_list(read_wannier_file("t3_hr.dat"));
        assert(hl.hoppings.size() == 1);
        assert(abs(get<1>(hl.hoppings[0]) - hopping_list::value_t(2.0, 0.0)) < 1e-12);
    }

    // (4) nrpts a multiple of 15 must NOT consume the first hopping line as a
    //     degeneracy line (the old nrpts/15 + 1 line count did).
    {
        vector<int> deg(15, 1);
        vector<string> hops;
        for (int r = 0; r < 15; ++r) hops.push_back(to_string(r) + " 0 0 1 1 1.0 0.0");
        write_hr("t4_hr.dat", 1, deg, hops);
        hopping_list hl = create_hopping_list(read_wannier_file("t4_hr.dat"));
        assert(hl.hoppings.size() == 15);   // all 15 preserved; the bug would yield 14
    }

    // (1) No-op regression: every ndegen == 1 leaves values unchanged.
    {
        write_hr("t1_hr.dat", 2, {1},
                 { "0 0 0 1 1 1.5 0.0",
                   "0 0 0 1 2 0.0 0.0",
                   "0 0 0 2 1 0.0 0.0",
                   "0 0 0 2 2 -2.5 0.0" });
        hopping_list hl = create_hopping_list(read_wannier_file("t1_hr.dat"));
        assert(hl.hoppings.size() == 2);
        bool found_p = false, found_n = false;
        for (const auto& h : hl.hoppings)
        {
            if (abs(get<1>(h) - hopping_list::value_t( 1.5, 0.0)) < 1e-12) found_p = true;
            if (abs(get<1>(h) - hopping_list::value_t(-2.5, 0.0)) < 1e-12) found_n = true;
        }
        assert(found_p && found_n);
    }

    std::cout << "NDEGEN NORMALIZATION TEST PASSED" << std::endl;
    return 0;
}
