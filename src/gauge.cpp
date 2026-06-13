#include "gauge.hpp"
#include "local_operators.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <algorithm>
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

    ifstream spnprobe((prefix + ".spn").c_str());          // optional: only spin needs it
    if (spnprobe.good()) { spnprobe.close(); read_spn(prefix + ".spn", g); }

    // Per-k dimensional consistency (docs/conventions.md sec 2).
    for (int k = 0; k < g.num_kpts; ++k)
    {
        assert(g.U[k].rows() == g.num_wann && g.U[k].cols() == g.num_wann);
        if (g.disentangled)
            assert(g.Udis[k].rows() == g.num_bands && g.Udis[k].cols() == g.num_wann);
        if (!g.Sb.empty())
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


amn_data read_amn(const std::string& amn_file)
{
    amn_data a;
    ifstream f(amn_file.c_str());
    assert(f.is_open());
    string header; getline(f, header);
    f >> a.num_bands >> a.num_kpts >> a.num_proj;     // header: nbnd nkpt nwann
    a.A.assign(a.num_kpts, Eigen::MatrixXcd::Zero(a.num_bands, a.num_proj));
    int m, alpha, k; double re, im;
    while (f >> m >> alpha >> k >> re >> im)           // lines: m alpha k Re Im (self-indexed)
        a.A[k-1](m-1, alpha-1) = complex<double>(re, im);
    return a;
}

static string lower_trim(const string& s)
{
    string t;
    for (char c : s) if (!isspace((unsigned char)c)) t += (char)tolower((unsigned char)c);
    return t;
}

std::vector<int> parse_projection_shells(const std::string& winfile)
{
    ifstream f(winfile.c_str());
    if (!f.is_open()) throw std::runtime_error("parse_projection_shells: cannot open " + winfile);

    std::vector<int> shells;
    string line; bool in_block = false;
    while (getline(f, line))
    {
        const string low = lower_trim(line);
        if (low.find("beginprojections") != string::npos) { in_block = true; continue; }
        if (low.find("endprojections")   != string::npos) break;
        if (!in_block) continue;

        const size_t colon = line.find(':');
        if (colon == string::npos) continue;            // unit/"random"/blank line

        // orbital list = between the first ':' and the next ':' (axis specs, z=/x=)
        const size_t colon2 = line.find(':', colon + 1);
        string orbitals = line.substr(colon + 1,
                                      (colon2 == string::npos) ? string::npos : colon2 - colon - 1);

        std::stringstream ss(orbitals);
        string tok;
        while (getline(ss, tok, ';'))
        {
            const string t = lower_trim(tok);
            if (t.empty()) continue;
            if      (t == "s") shells.push_back(0);   // complete l=0 shell -> L=0 (trivial)
            else if (t == "p") shells.push_back(1);
            else if (t == "d") shells.push_back(2);
            else throw std::runtime_error(
                "orbital L: unsupported projection '" + t +
                "' -- on-site L needs a complete pure s, p or d shell "
                "(hybrids and partial shells are out of scope; see docs/conventions.md sec 5)");
        }
    }
    return shells;
}

hopping_list orbital_L_operator(const gauge_data& g, const amn_data& a,
                                const std::vector<int>& shell_ls, int alpha,
                                const std::vector<hopping_list::cellID_t>& Rset)
{
    assert(g.num_bands == a.num_bands && g.num_kpts == a.num_kpts);

    // Validate shells vs num_proj BEFORE assembling (avoid out-of-range blocks).
    int total = 0;
    for (int l : shell_ls) total += 2*l + 1;
    if (total != a.num_proj)
        throw std::runtime_error("orbital L: projector count from .win shells (" +
                                 std::to_string(total) + ") != .amn num_proj (" +
                                 std::to_string(a.num_proj) + ")");

    // Assemble L_local^alpha (num_proj x num_proj), block-diagonal over shells.
    Eigen::MatrixXcd Lloc = Eigen::MatrixXcd::Zero(a.num_proj, a.num_proj);
    int off = 0;
    for (int l : shell_ls)
    {
        const int d = 2*l + 1;
        Lloc.block(off, off, d, d) = shell_L(l)[alpha];
        off += d;
    }

    const int nw = g.num_wann;
    const double twopi = 2.0 * M_PI;
    std::vector<Eigen::MatrixXcd> Lk(g.num_kpts);
    for (int k = 0; k < g.num_kpts; ++k)
    {
        const Eigen::MatrixXcd C = a.A[k].adjoint() * g.V(k);   // num_proj x num_wann
        Lk[k] = C.adjoint() * Lloc * C;                         // num_wann x num_wann
    }

    hopping_list out;
    out.num_wann = nw;
    for (const auto& R : Rset)
    {
        Eigen::MatrixXcd acc = Eigen::MatrixXcd::Zero(nw, nw);
        for (int k = 0; k < g.num_kpts; ++k)
        {
            const double kr = g.kpt[k][0]*R[0] + g.kpt[k][1]*R[1] + g.kpt[k][2]*R[2];
            acc += std::exp(complex<double>(0.0, -twopi*kr)) * Lk[k];
        }
        acc /= static_cast<double>(g.num_kpts);
        for (int i = 0; i < nw; ++i)
            for (int j = 0; j < nw; ++j)
                if (std::abs(acc(i,j)) > 1e-10)
                    out.hoppings.push_back(hopping_list::hopping_t(R, acc(i,j), hopping_list::edge_t({i,j})));
    }
    return out;
}
