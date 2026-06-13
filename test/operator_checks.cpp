#include <cassert>
#include <iostream>
#include "hopping_list.hpp"
#include "checks.hpp"
using namespace std;

typedef hopping_list::value_t  value_t;
typedef hopping_list::cellID_t cellID_t;
typedef hopping_list::edge_t   edge_t;
typedef hopping_list::hopping_t hopping_t;

int main()
{
    // Hermitian operator: O_ij(R) = conj(O_ji(-R)).
    hopping_list H; H.num_wann = 2;
    H.hoppings.push_back(hopping_t(cellID_t({ 1,0,0}), value_t(1.0, 0.5), edge_t({0,1})));
    H.hoppings.push_back(hopping_t(cellID_t({-1,0,0}), value_t(1.0,-0.5), edge_t({1,0})));  // mirror
    H.hoppings.push_back(hopping_t(cellID_t({ 0,0,0}), value_t(2.0, 0.0), edge_t({0,0})));
    assert(checks::hermiticity(H).pass);

    // Corrupted: an entry with no Hermitian mirror -> FAIL with residual > tol.
    hopping_list B = H;
    B.hoppings.push_back(hopping_t(cellID_t({2,0,0}), value_t(0.0, 1.0), edge_t({0,1})));
    { auto r = checks::hermiticity(B); assert(!r.pass && r.residual > 1e-8); }

    // sum_rules: H has Tr O(0)=2 -> traceless check FAILs; an actually traceless op PASSes.
    assert(!checks::trace_rule(H, /*expect_traceless=*/true).pass);
    assert( checks::trace_rule(H, /*expect_traceless=*/false).pass);   // informational
    hopping_list T; T.num_wann = 2;
    T.hoppings.push_back(hopping_t(cellID_t({0,0,0}), value_t( 1,0), edge_t({0,0})));
    T.hoppings.push_back(hopping_t(cellID_t({0,0,0}), value_t(-1,0), edge_t({1,1})));
    assert(checks::trace_rule(T, true).pass);

    // aliasing: +4 and -4 collide mod 8 (both -> 4), are distinct mod 9 (4,5).
    hopping_list C; C.num_wann = 1;
    C.hoppings.push_back(hopping_t(cellID_t({ 4,0,0}), value_t(1,0), edge_t({0,0})));
    C.hoppings.push_back(hopping_t(cellID_t({-4,0,0}), value_t(1,0), edge_t({0,0})));
    assert(!checks::aliasing(C, cellID_t({8,1,1})).pass);
    assert( checks::aliasing(C, cellID_t({9,1,1})).pass);

    // generator algebra
    assert(checks::spin_algebra().pass);
    assert(checks::orbital_algebra().pass);

    std::cout << "OPERATOR CHECKS TEST PASSED" << std::endl;
    return 0;
}
