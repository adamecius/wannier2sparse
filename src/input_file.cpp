/**
 * @file input_file.cpp
 * @brief The --create / --write / --run input-file workflow (see include/input_file.hpp).
 *
 * A wannier2sparse run is described by a plain `key = value` text file whose keys
 * are exactly the old CLI options. This keeps the established formats intact (the
 * file is just their serialization) while giving a self-documenting, reproducible,
 * incrementally-editable run description, and an lsquant-style output summary.
 */
#include "input_file.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace w2sp {

// --- small string helpers --------------------------------------------------
static std::string strip_comment(const std::string& s)
{
    const std::size_t h = s.find('#');
    return h == std::string::npos ? s : s.substr(0, h);
}
static std::string trim(const std::string& s)
{
    const std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}
// Split "key = value" (first '='). Returns false if the line has no key.
static bool split_kv(const std::string& raw, std::string& key, std::string& val)
{
    const std::string line = strip_comment(raw);
    const std::size_t eq = line.find('=');
    if (eq == std::string::npos) { key = trim(line); val = ""; return !key.empty(); }
    key = trim(line.substr(0, eq)); val = trim(line.substr(eq + 1));
    return !key.empty();
}
static std::vector<std::string> tokens(const std::string& s)
{
    std::vector<std::string> out; std::istringstream is(s); std::string t;
    while (is >> t) out.push_back(t);
    return out;
}
static bool as_bool(const std::string& v)
{
    std::string s = v; std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s == "true" || s == "1" || s == "yes" || s == "on";
}
// Is `key` a repeatable key (may appear on several lines)?
static bool repeatable(const std::string& k) { return k == "spin_current" || k == "op_file"; }

// --- template --------------------------------------------------------------
void input_file_create(const std::string& path)
{
    std::ofstream f(path.c_str());
    if (!f) throw std::runtime_error("input_file_create: cannot write '" + path + "'");
    f <<
"# wannier2sparse input file\n"
"#   run :  wannier2sparse --run " << path << "\n"
"#   edit:  wannier2sparse --write KEY [VALUE...] -inp " << path << "\n"
"# Keys are the CLI options as `key = value`; '#' starts a comment.\n"
"\n"
"label        = MODEL          # system label; reads <project_dir>/<seed|label>_hr.dat\n"
"project_dir  = .              # directory holding the input files\n"
"seed         =                # input seedname (default: label)\n"
"output_dir   = .              # where the .CSR / bundle is written\n"
"mode         = sparse         # sparse (supercell CSR) | bundle (primitive O(R) + manifest)\n"
"supercell    = 1 1 1          # N1 N2 N3 (sparse mode only)\n"
"operators    =                # space list: VX VY VZ SX SY SZ VXSZ ...\n"
"spin_current =                # repeatable; \"Vdir Sdir\", e.g.  X Z   (1/2{V,S}, CSR anticommutator)\n"
"exact_spin   = false          # gauge-route spin from .spn + _u.mat (SOC/noncollinear)\n"
"orbital_L    = false          # orbital angular momentum from .amn + _u.mat\n"
"velocity_mode = berry_connection # velocity ladder: bare | berry_connection (default) | covariant\n"
"r_dat        =                # position matrix _r.dat for velocity_mode=covariant (default <seed>_r.dat)\n"
"bounds       = false          # write .desc spectral-bound sidecars\n"
"check        =                # self-checks: all|hermiticity|sum_rules|algebra|aliasing|bounds\n"
"op_file      =                # repeatable; \"NAME PATH\" external _hr.dat operator (e.g. JXSZ model_JXSZ_hr.dat)\n";
}

// --- set / append one key --------------------------------------------------
void input_file_set(const std::string& path, const std::string& key,
                    const std::string& value, bool append)
{
    std::vector<std::string> lines;
    {
        std::ifstream in(path.c_str());
        if (!in) { input_file_create(path); std::ifstream in2(path.c_str()); }
        std::ifstream in3(path.c_str()); std::string l;
        while (std::getline(in3, l)) lines.push_back(l);
    }
    bool done = false;
    if (!repeatable(key))
    {
        for (auto& l : lines)
        {
            std::string k, v;
            if (split_kv(l, k, v) && k == key)
            {
                const std::string nv = (append && !v.empty()) ? (v + " " + value) : value;
                l = key + " = " + nv;
                done = true; break;
            }
        }
    }
    if (!done) lines.push_back(key + " = " + value);     // repeatable -> always a new line

    std::ofstream out(path.c_str());
    if (!out) throw std::runtime_error("input_file_set: cannot write '" + path + "'");
    for (const auto& l : lines) out << l << "\n";
}

// --- parse -> args ---------------------------------------------------------
void input_file_apply(const std::string& path, W2SP_arguments& args)
{
    std::ifstream in(path.c_str());
    if (!in) throw std::runtime_error("input_file_apply: cannot open '" + path + "'");
    std::string raw;
    while (std::getline(in, raw))
    {
        std::string key, val;
        if (!split_kv(raw, key, val)) continue;          // blank / comment-only
        if (val.empty() && key.find('=') == std::string::npos &&
            strip_comment(raw).find('=') == std::string::npos)
            continue;                                    // bare key with no '=' and no value -> skip
        if      (key == "label")        args.label = val;
        else if (key == "project_dir")  args.project_dir = val;
        else if (key == "seed")         args.seed = val;
        else if (key == "output_dir")   { if (!val.empty()) args.output_dir = val; }
        else if (key == "mode")         { if (!val.empty()) args.mode = val; }
        else if (key == "operators")    { for (auto& o : tokens(val)) args.operators.push_back(o); }
        else if (key == "supercell")
        {
            auto t = tokens(val);
            if (!t.empty())
            {
                if (t.size() != 3) throw std::runtime_error("input file: supercell needs 3 integers, got '" + val + "'");
                for (int i = 0; i < 3; ++i) args.cellDim[i] = std::stoi(t[i]);
            }
        }
        else if (key == "spin_current")
        {
            auto t = tokens(val);
            if (!t.empty())
            {
                if (t.size() != 2) throw std::runtime_error("input file: spin_current needs 'Vdir Sdir', got '" + val + "'");
                args.spin_currents.push_back(std::make_pair(t[0][0], t[1][0]));
            }
        }
        else if (key == "op_file")
        {
            auto t = tokens(val);
            if (!t.empty())
            {
                if (t.size() != 2) throw std::runtime_error("input file: op_file needs 'NAME PATH', got '" + val + "'");
                args.op_files.push_back(std::make_pair(t[0], t[1]));
            }
        }
        else if (key == "exact_spin")          args.exact_spin = as_bool(val);
        else if (key == "orbital_L" || key == "orbital_l") args.orbital_l = as_bool(val);
        else if (key == "covariant_velocity")  { if (as_bool(val)) { args.velocity_mode = "covariant"; args.covariant_velocity = true; } }
        else if (key == "velocity_mode")
        {
            if (!val.empty())
            {
                if (val != "bare" && val != "berry_connection" && val != "covariant")
                    throw std::runtime_error("input file: velocity_mode must be bare|berry_connection|covariant, got '" + val + "'");
                args.velocity_mode = val;
                args.covariant_velocity = (val == "covariant");
            }
        }
        else if (key == "r_dat")               { if (!val.empty()) args.r_dat_path = val; }
        else if (key == "bounds")              args.emit_descriptor = as_bool(val);
        else if (key == "check")               { if (!val.empty()) args.check = val; }
        else throw std::runtime_error("input file: unknown key '" + key + "' in " + path);
    }
}

// --- output summary --------------------------------------------------------
static long file_size(const std::string& p)
{
    std::ifstream f(p.c_str(), std::ios::binary | std::ios::ate);
    return f ? static_cast<long>(f.tellg()) : -1;
}
void input_file_write_output(const std::string& outpath, const W2SP_arguments& args)
{
    std::ofstream f(outpath.c_str());
    if (!f) return;
    f << "# wannier2sparse run summary\n";
    f << "label        = " << args.label << "\n";
    f << "mode         = " << args.mode << "\n";
    f << "supercell    = " << args.cellDim[0] << " " << args.cellDim[1] << " " << args.cellDim[2] << "\n";
    f << "output_dir   = " << args.output_dir << "\n";
    f << "exact_spin   = " << (args.exact_spin ? "true" : "false") << "\n";
    f << "orbital_L    = " << (args.orbital_l ? "true" : "false") << "\n";
    f << "velocity_mode = " << args.velocity_mode << "\n";
    f << "bounds       = " << (args.emit_descriptor ? "true" : "false") << "\n";
    f << "operators    =";
    for (const auto& o : args.operators) f << " " << o;
    f << "\n";
    for (const auto& sc : args.spin_currents) f << "spin_current = " << sc.first << " " << sc.second << "\n";
    for (const auto& of : args.op_files)      f << "op_file      = " << of.first << " " << of.second << "\n";
    // Produced operator files (name + bytes), the lsquant-style results part.
    const std::string prefix = args.output_dir + "/" + args.label + ".";
    const std::string suffix = (args.mode == "bundle") ? "" : ".CSR";
    f << "# produced operator files:\n";
    std::vector<std::string> names = {"HAM"};
    for (const auto& o : args.operators) names.push_back(o);
    for (const auto& of : args.op_files) names.push_back(of.first);
    for (const auto& sc : args.spin_currents) { std::string n = "J"; n += sc.first; n += 'S'; n += sc.second; names.push_back(n); }
    for (const auto& n : names)
    {
        const std::string p = prefix + n + suffix;
        const long sz = file_size(p);
        if (sz >= 0) f << "produced     = " << n << " " << p << " " << sz << " bytes\n";
    }
}

} // namespace w2sp
