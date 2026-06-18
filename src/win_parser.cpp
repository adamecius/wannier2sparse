/**
 * @file win_parser.cpp
 * @brief Implementation of the Wannier90 `.win` provenance parser.
 */
#include "win_parser.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace {

/// Trim leading/trailing ASCII whitespace.
std::string trim(const std::string& s)
{
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

/// Lower-case copy.
std::string lower(const std::string& s)
{
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(),
                   [](unsigned char c){ return (char)std::tolower(c); });
    return t;
}

/// Collapse to lower-case with all whitespace removed (for keyword matching).
std::string squash(const std::string& s)
{
    std::string t;
    for (char c : s) if (!std::isspace((unsigned char)c)) t += (char)std::tolower((unsigned char)c);
    return t;
}

/// Strip a trailing `!`/`#` comment.
std::string strip_comment(const std::string& s)
{
    const size_t p = s.find_first_of("!#");
    return (p == std::string::npos) ? s : s.substr(0, p);
}

/// Parse a Fortran/Wannier90 logical value (.true./.false./T/F/true/false).
bool parse_logical(const std::string& v)
{
    const std::string t = squash(v);
    if (t == "t" || t == ".true." || t == "true" || t == "1") return true;
    return false;   // anything else (incl. .false./f/false/0) is false
}

/// Parse Wannier90 `exclude_bands` (e.g. "1-4,6,9-10") into a sorted band list.
std::vector<int> parse_band_list(const std::string& v)
{
    std::vector<int> out;
    std::stringstream ss(v);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        tok = trim(tok);
        if (tok.empty()) continue;
        const size_t dash = tok.find('-');
        try {
            if (dash == std::string::npos)
                out.push_back(std::stoi(tok));
            else
            {
                const int lo = std::stoi(tok.substr(0, dash));
                const int hi = std::stoi(tok.substr(dash + 1));
                for (int b = lo; b <= hi; ++b) out.push_back(b);
            }
        } catch (const std::exception&) { /* ignore malformed token */ }
    }
    return out;
}

/// Split a `.win` statement into (key, value) on the first '=' or ':'.
bool split_kv(const std::string& line, std::string& key, std::string& val)
{
    const size_t eq = line.find('=');
    const size_t co = line.find(':');
    size_t sep = std::string::npos;
    if (eq != std::string::npos && co != std::string::npos) sep = std::min(eq, co);
    else if (eq != std::string::npos) sep = eq;
    else if (co != std::string::npos) sep = co;
    if (sep == std::string::npos) return false;
    key = lower(trim(line.substr(0, sep)));
    val = trim(line.substr(sep + 1));
    return !key.empty();
}

} // namespace

WannierProvenance parse_win_provenance(const std::string& winfile)
{
    std::ifstream f(winfile.c_str());
    if (!f.is_open()) throw std::runtime_error("parse_win_provenance: cannot open " + winfile);

    WannierProvenance w;
    w.present = true;
    w.source_file = winfile;

    std::string raw;
    std::string current_block;   // non-empty while inside a begin/end block
    while (std::getline(f, raw))
    {
        const std::string line = trim(strip_comment(raw));
        if (line.empty()) continue;

        const std::string sq = squash(line);

        // Block boundaries. We only capture the projections block; other blocks
        // (unit_cell_cart, atoms_frac, kpoints, ...) are skipped.
        if (sq.rfind("begin", 0) == 0)
        {
            current_block = squash(line.substr(line.find_first_not_of(" \t") + 5));
            continue;
        }
        if (sq.rfind("end", 0) == 0)
        {
            current_block.clear();
            continue;
        }

        if (current_block == "projections")
        {
            // Skip option lines like "random" or a leading units specifier.
            const size_t colon = line.find(':');
            if (colon == std::string::npos) continue;
            Projection p;
            p.raw     = line;
            p.site    = trim(line.substr(0, colon));
            // orbital = text up to the next ':' (which would start axis specs)
            const size_t colon2 = line.find(':', colon + 1);
            p.orbital = trim(line.substr(colon + 1,
                             (colon2 == std::string::npos) ? std::string::npos : colon2 - colon - 1));
            w.projections.push_back(p);
            continue;
        }
        if (!current_block.empty()) continue;   // inside some other block: skip

        std::string key, val;
        if (!split_kv(line, key, val)) continue;

        if (key == "num_wann")  { try { w.num_wann  = std::stoi(val); w.has_num_wann  = true; } catch (...) {} }
        else if (key == "num_bands") { try { w.num_bands = std::stoi(val); w.has_num_bands = true; } catch (...) {} }
        else if (key == "exclude_bands") { w.exclude_bands = parse_band_list(val); }
        else if (key == "mp_grid")
        {
            std::stringstream ss(val);
            std::array<int,3> g{{0,0,0}};
            if (ss >> g[0] >> g[1] >> g[2]) { w.mp_grid = g; w.has_mp_grid = true; }
        }
        else if (key == "dis_win_min")  { try { w.dis_win_min  = std::stod(val); w.has_dis_win  = true; } catch (...) {} }
        else if (key == "dis_win_max")  { try { w.dis_win_max  = std::stod(val); w.has_dis_win  = true; } catch (...) {} }
        else if (key == "dis_froz_min") { try { w.dis_froz_min = std::stod(val); w.has_dis_froz = true; } catch (...) {} }
        else if (key == "dis_froz_max") { try { w.dis_froz_max = std::stod(val); w.has_dis_froz = true; } catch (...) {} }
        else if (key == "use_ws_distance") { w.use_ws_distance = parse_logical(val); }
    }

    // Disentanglement is in play when there are more bands than wannier functions,
    // or when an explicit window was given.
    w.dis_enabled = w.has_dis_win || w.has_dis_froz ||
                    (w.has_num_bands && w.has_num_wann && w.num_bands > w.num_wann);

    return w;
}
