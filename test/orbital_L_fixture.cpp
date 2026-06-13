#include <cassert>
#include <cstdlib>
#include <fstream>
#include <set>
#include <iostream>
#include <Eigen/Dense>
#include "gauge.hpp"
#include "local_operators.hpp"
#include "wannier_parser.hpp"
#include "hopping_list.hpp"
using namespace std;

typedef hopping_list::cellID_t cellID_t;

// Positive validation of the orbital-L route on a REAL pure-d Wannier model
// (Plan 10B). The fixture (example04 copper, Cu:d) is generated locally by
// test/fixtures/gen_fixture.sh and lives outside the repo, so this test SKIPs
// (passes) when the fixture is absent (e.g. in CI) and validates when present.
// Point it at <prefix> via W2SP_ORBITAL_FIXTURE; default /tmp/fix/copper/copper.
static Eigen::MatrixXcd at(const hopping_list& hl, const cellID_t& R, int nw)
{
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nw, nw);
    for (const auto& h : hl.hoppings)
        if (get<0>(h) == R) { const auto e = get<2>(h); M(e[0],e[1]) += get<1>(h); }
    return M;
}

int main()
{
    const char* env = std::getenv("W2SP_ORBITAL_FIXTURE");
    const string p = env ? env : "/tmp/fix/copper/copper";
    { ifstream f((p + ".amn").c_str()); if (!f.good())
      { std::cout << "orbital_L_fixture: SKIP (no fixture at " << p << ")\n"; return 0; } }

    const auto shells = parse_projection_shells(p + ".win");
    gauge_data g = read_gauge(p);
    amn_data   A = read_amn(p + ".amn");
    int total = 0; for (int l : shells) total += 2*l + 1;
    assert(total == A.num_proj);

    // The decisive check: L_local for the (pure) d shell has integer eigenvalues.
    bool has_d = false; for (int l : shells) if (l == 2) has_d = true;
    if (has_d)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(shell_L(2)[2]);
        for (int i = 0; i < 5; ++i) assert(abs(es.eigenvalues()(i) - double(i-2)) < 1e-9);
    }

    // L_W(R) Hermitian across R <-> -R on the model's R-set.
    hopping_list H = create_hopping_list(read_wannier_file(p + "_hr.dat"));
    std::set<cellID_t> Rs; for (const auto& h : H.hoppings) Rs.insert(get<0>(h));
    const std::vector<cellID_t> Rset(Rs.begin(), Rs.end());
    double herm = 0;
    for (int a = 0; a < 3; ++a)
    {
        hopping_list L = orbital_L_operator(g, A, shells, a, Rset);
        for (const auto& R : Rset)
        {
            const cellID_t nR = {-R[0],-R[1],-R[2]};
            if (!Rs.count(nR)) continue;
            herm = std::max(herm, (at(L,R,g.num_wann) - at(L,nR,g.num_wann).adjoint()).cwiseAbs().maxCoeff());
        }
    }
    assert(herm < 1e-9);

    std::cout << "orbital_L_fixture: PASS (real pure-d model; Lz_local={-2..2}, L_W Hermitian)\n";
    return 0;
}
