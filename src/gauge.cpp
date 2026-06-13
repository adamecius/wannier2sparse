#include "gauge.hpp"
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
using namespace std;

// All readers use operator>> for the numeric payload: it skips the blank lines
// and the per-k k-point lines transparently, and parses both the sign-glued
// f15.10 (u.mat) and the es26.16 (.spn) layouts. Header text lines are consumed
// with getline. Formats verified in docs/conventions.md (sec 2 and 3).

static void read_umat(const string& file, gauge_data& g, bool is_dis)
{
    ifstream f(file.c_str());
    assert(f.is_open());
    string header; getline(f, header);
    int nk, a, b; f >> nk >> a >> b;     // u.mat: nk nw nw ; u_dis.mat: nk nw nb

    if (!is_dis) { g.num_kpts = nk; g.num_wann = a; g.kpt.resize(nk); g.U.resize(nk); }
    else         { g.num_bands = b; g.disentangled = true; g.Udis.resize(nk); }

    const int rows = is_dis ? g.num_bands : g.num_wann;   // data block is rows x num_wann
    const int cols = g.num_wann;
    for (int k = 0; k < nk; ++k)
    {
        double kx, ky, kz; f >> kx >> ky >> kz;           // per-k k-point line
        if (!is_dis) g.kpt[k] = {kx, ky, kz};
        Eigen::MatrixXcd M(rows, cols);
        for (int j = 0; j < cols; ++j)                    // column-major: i fastest
            for (int i = 0; i < rows; ++i)
            {
                double re, im; f >> re >> im;
                M(i, j) = complex<double>(re, im);
            }
        if (is_dis) g.Udis[k] = M; else g.U[k] = M;
    }
}

static void read_spn(const string& file, gauge_data& g)
{
    ifstream f(file.c_str());
    assert(f.is_open());
    string header; getline(f, header);
    int nb, nk; f >> nb >> nk;
    if (g.num_bands == 0) g.num_bands = nb;
    assert(nb == g.num_bands);
    assert(g.num_kpts == 0 || nk == g.num_kpts);

    g.Sb.resize(nk);
    for (int k = 0; k < nk; ++k)
    {
        for (int s = 0; s < 3; ++s)
            g.Sb[k][s] = Eigen::MatrixXcd::Zero(nb, nb);
        // packing: m outer (col), n<=m inner (row); per pair sx, sy, sz.
        for (int m = 0; m < nb; ++m)
            for (int n = 0; n <= m; ++n)
                for (int s = 0; s < 3; ++s)
                {
                    double re, im; f >> re >> im;
                    const complex<double> v(re, im);
                    g.Sb[k][s](n, m) = v;                 // upper triangle (row n <= col m)
                    g.Sb[k][s](m, n) = conj(v);           // Hermiticity
                }
    }
}

gauge_data read_gauge(const string& prefix)
{
    gauge_data g;
    read_umat(prefix + "_u.mat", g, /*is_dis=*/false);

    ifstream probe((prefix + "_u_dis.mat").c_str());
    if (probe.good()) { probe.close(); read_umat(prefix + "_u_dis.mat", g, /*is_dis=*/true); }
    else              { g.disentangled = false; g.num_bands = g.num_wann; }  // V = U

    read_spn(prefix + ".spn", g);

    // Per-k dimensional consistency (docs/conventions.md sec 2).
    for (int k = 0; k < g.num_kpts; ++k)
    {
        assert(g.U[k].rows() == g.num_wann && g.U[k].cols() == g.num_wann);
        if (g.disentangled)
            assert(g.Udis[k].rows() == g.num_bands && g.Udis[k].cols() == g.num_wann);
        for (int s = 0; s < 3; ++s)
            assert(g.Sb[k][s].rows() == g.num_bands && g.Sb[k][s].cols() == g.num_bands);
    }
    return g;
}

hopping_list exact_spin_operator(const gauge_data& g, int alpha,
                                 const std::vector<hopping_list::cellID_t>& Rset)
{
    const int nw = g.num_wann;
    const double twopi = 2.0 * M_PI;

    // S_W^alpha(k) = V(k)^dagger S_B^alpha(k) V(k)   (num_wann x num_wann)
    std::vector<Eigen::MatrixXcd> Sw(g.num_kpts);
    for (int k = 0; k < g.num_kpts; ++k)
    {
        const Eigen::MatrixXcd V = g.V(k);
        Sw[k] = V.adjoint() * g.Sb[k][alpha] * V;
    }

    hopping_list out;
    out.num_wann = nw;
    for (const auto& R : Rset)
    {
        Eigen::MatrixXcd acc = Eigen::MatrixXcd::Zero(nw, nw);
        for (int k = 0; k < g.num_kpts; ++k)
        {
            const double kr = g.kpt[k][0]*R[0] + g.kpt[k][1]*R[1] + g.kpt[k][2]*R[2];
            const complex<double> phase = std::exp(complex<double>(0.0, -twopi * kr));
            acc += phase * Sw[k];
        }
        acc /= static_cast<double>(g.num_kpts);                 // 1/N_k, no ndegen (sec 1)

        for (int i = 0; i < nw; ++i)
            for (int j = 0; j < nw; ++j)
            {
                const complex<double> v = acc(i, j);
                if (std::abs(v) > 1e-10)
                    out.hoppings.push_back(hopping_list::hopping_t(R, v, hopping_list::edge_t({i, j})));
            }
    }
    return out;
}
