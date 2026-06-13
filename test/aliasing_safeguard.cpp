#include <cassert>
#include <stdexcept>
#include <iostream>
#include "hopping_list.hpp"
using namespace std;

typedef hopping_list::value_t  value_t;
typedef hopping_list::cellID_t cellID_t;
typedef hopping_list::edge_t   edge_t;
typedef hopping_list::hopping_t hopping_t;

static bool guards(const hopping_list& hl, const cellID_t& N)
{
    try { guard_minimum_image(hl, N); return false; }    // no throw -> safe
    catch (const std::exception&) { return true; }       // threw -> guarded
}

int main()
{
    // 1D chain with 4th-nearest neighbours: R = -4..+4 on the single orbital.
    hopping_list chain; chain.num_wann = 1;
    for (int R = -4; R <= 4; ++R)
        chain.hoppings.push_back(hopping_t(cellID_t({R,0,0}), value_t(1,0), edge_t({0,0})));

    // range 4 -> need N >= 9; N<=8 aliases (e.g. +4 and -4 mod 8).
    for (int N = 1; N <= 8; ++N)
        assert(guards(chain, cellID_t({N,1,1})) && "4-NN chain must be guarded for N<=8");
    for (int N = 9; N <= 12; ++N)
        assert(!guards(chain, cellID_t({N,1,1})) && "4-NN chain must pass for N>=9");

    // range-1 model (graphene/Haldane-like): R = -1,0,1 -> need N >= 3.
    hopping_list r1; r1.num_wann = 1;
    for (int R = -1; R <= 1; ++R)
        r1.hoppings.push_back(hopping_t(cellID_t({R,0,0}), value_t(1,0), edge_t({0,0})));
    assert(guards(r1, cellID_t({2,1,1})));      // N=2 aliases +1,-1
    assert(!guards(r1, cellID_t({3,1,1})));     // N=3 ok
    assert(guards(r1, cellID_t({1,1,1})));      // N0=1 with range 1 -> alias

    std::cout << "ALIASING SAFEGUARD TEST PASSED" << std::endl;
    return 0;
}
