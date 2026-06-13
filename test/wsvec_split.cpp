#include <cassert>
#include <fstream>
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

int main()
{
    // A 2-orbital model with two hoppings.
    hopping_list hl;
    hl.num_wann = 2;
    hl.hoppings.push_back(hopping_t(cellID_t({0,0,0}), value_t(4.0, 0.0), edge_t({0,1})));
    hl.hoppings.push_back(hopping_t(cellID_t({1,0,0}), value_t(2.0, 0.0), edge_t({0,0})));

    // wsvec: (0,0,0;0,1) stays put (nT=1, T=0); (1,0,0;0,0) splits into two images.
    vector<wsvec_entry> ws(2);
    ws[0].R = {0,0,0}; ws[0].i = 0; ws[0].j = 1; ws[0].T = { {0,0,0} };
    ws[1].R = {1,0,0}; ws[1].i = 0; ws[1].j = 0; ws[1].T = { {0,0,0}, {0,-1,0} };

    hopping_list got = apply_wsvec(hl, ws);

    hopping_list want;
    want.num_wann = 2;
    want.hoppings.push_back(hopping_t(cellID_t({0, 0,0}), value_t(4.0, 0.0), edge_t({0,1})));
    want.hoppings.push_back(hopping_t(cellID_t({1, 0,0}), value_t(1.0, 0.0), edge_t({0,0})));  // 2.0/2
    want.hoppings.push_back(hopping_t(cellID_t({1,-1,0}), value_t(1.0, 0.0), edge_t({0,0})));  // shifted image
    assert(got == want);

    // Empty wsvec (no use_ws_distance) is a pass-through.
    vector<wsvec_entry> none;
    hopping_list same = apply_wsvec(hl, none);
    assert(same == hl);

    // read_wsvec parses the W90 format (1-based i,j -> 0-based) and the T lists.
    {
        ofstream("w_wsvec.dat")
            << "## written with use_ws_distance=.true.\n"
               "    0    0    0    1    2\n"
               "    1\n"
               "    0    0    0\n"
               "    1    0    0    1    1\n"
               "    2\n"
               "    0    0    0\n"
               "    0   -1    0\n";
        vector<wsvec_entry> r = read_wsvec("w_wsvec.dat");
        assert(r.size() == 2);
        assert(r[0].i == 0 && r[0].j == 1 && r[0].T.size() == 1);
        assert(r[1].i == 0 && r[1].j == 0 && r[1].T.size() == 2);
        assert((r[1].R == cellID_t({1,0,0})));
        assert((r[1].T[1] == array<int,3>({0,-1,0})));
    }

    std::cout << "WSVEC SPLIT TEST PASSED" << std::endl;
    return 0;
}
