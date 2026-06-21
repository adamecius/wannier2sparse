/**
 * @file run_config.cpp
 * @brief Minimal, dependency-free JSON reader for the bundle `run.json` config.
 *
 * A small recursive-descent JSON parser (objects, arrays, strings, numbers,
 * booleans, null) is used instead of a heavy JSON dependency: run.json is flat
 * and small, and keeping the reader hand-rolled means the config format stays
 * under our control. The parser is strict enough to reject malformed input with
 * a clear message but otherwise tolerant of whitespace and key order.
 */
#include "run_config.hpp"
#include "json_writer.hpp"
#include "provenance_writer.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <memory>
#include <cctype>
#include <cstdlib>

namespace {

/// A parsed JSON value (null / bool / number / string / array / object).
struct JsonValue
{
    enum Type { NUL, BOOL, NUM, STR, ARR, OBJ } type;
    bool                              b;
    double                            num;
    std::string                       str;
    std::vector<JsonValue>            arr;
    std::map<std::string, JsonValue>  obj;
    JsonValue() : type(NUL), b(false), num(0.0) {}
};

/// Recursive-descent JSON parser over an in-memory string.
class JsonParser
{
public:
    explicit JsonParser(const std::string& text) : s(text), pos(0) {}

    JsonValue parse()
    {
        JsonValue v = parse_value();
        skip_ws();
        if (pos != s.size()) fail("trailing characters after JSON value");
        return v;
    }

private:
    const std::string& s;
    size_t             pos;

    [[noreturn]] void fail(const std::string& msg) const
    {
        throw std::runtime_error("run.json: " + msg + " (at offset " + std::to_string(pos) + ")");
    }

    void skip_ws()
    {
        // Whitespace plus `//` line and `/* ... */` block comments. Strict JSON has
        // no comments, but the `.w2s` input is meant to be hand-edited and the
        // --create-template scaffold is self-documenting, so the reader tolerates
        // them. A bare `_hr.dat`-only run.json (no comments) is unaffected.
        while (pos < s.size())
        {
            const char c = s[pos];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') { ++pos; continue; }
            if (c == '/' && pos + 1 < s.size())
            {
                if (s[pos + 1] == '/')                       // line comment to EOL
                {
                    pos += 2;
                    while (pos < s.size() && s[pos] != '\n') ++pos;
                    continue;
                }
                if (s[pos + 1] == '*')                       // block comment to */
                {
                    pos += 2;
                    while (pos + 1 < s.size() && !(s[pos] == '*' && s[pos + 1] == '/')) ++pos;
                    if (pos + 1 >= s.size()) fail("unterminated /* */ comment");
                    pos += 2;
                    continue;
                }
            }
            break;
        }
    }

    char peek()
    {
        skip_ws();
        if (pos >= s.size()) fail("unexpected end of input");
        return s[pos];
    }

    JsonValue parse_value()
    {
        const char c = peek();
        switch (c)
        {
            case '{': return parse_object();
            case '[': return parse_array();
            case '"': { JsonValue v; v.type = JsonValue::STR; v.str = parse_string(); return v; }
            case 't': case 'f': return parse_bool();
            case 'n': return parse_null();
            default:  return parse_number();
        }
    }

    JsonValue parse_object()
    {
        JsonValue v; v.type = JsonValue::OBJ;
        ++pos;                                   // consume '{'
        if (peek() == '}') { ++pos; return v; }
        for (;;)
        {
            if (peek() != '"') fail("expected string key in object");
            const std::string k = parse_string();
            if (peek() != ':') fail("expected ':' after object key");
            ++pos;                               // consume ':'
            v.obj[k] = parse_value();
            const char nc = peek();
            if (nc == ',') { ++pos; continue; }
            if (nc == '}') { ++pos; break; }
            fail("expected ',' or '}' in object");
        }
        return v;
    }

    JsonValue parse_array()
    {
        JsonValue v; v.type = JsonValue::ARR;
        ++pos;                                   // consume '['
        if (peek() == ']') { ++pos; return v; }
        for (;;)
        {
            v.arr.push_back(parse_value());
            const char nc = peek();
            if (nc == ',') { ++pos; continue; }
            if (nc == ']') { ++pos; break; }
            fail("expected ',' or ']' in array");
        }
        return v;
    }

    std::string parse_string()
    {
        ++pos;                                   // consume opening '"'
        std::string out;
        while (pos < s.size())
        {
            const char c = s[pos++];
            if (c == '"') return out;
            if (c == '\\')
            {
                if (pos >= s.size()) fail("unterminated escape in string");
                const char e = s[pos++];
                switch (e)
                {
                    case '"':  out += '"';  break;
                    case '\\': out += '\\'; break;
                    case '/':  out += '/';  break;
                    case 'n':  out += '\n'; break;
                    case 't':  out += '\t'; break;
                    case 'r':  out += '\r'; break;
                    case 'b':  out += '\b'; break;
                    case 'f':  out += '\f'; break;
                    case 'u':  fail("\\u escapes are not supported in run.json");
                    default:   fail("invalid escape in string");
                }
            }
            else out += c;
        }
        fail("unterminated string");
    }

    JsonValue parse_bool()
    {
        JsonValue v; v.type = JsonValue::BOOL;
        if (s.compare(pos, 4, "true") == 0)       { v.b = true;  pos += 4; }
        else if (s.compare(pos, 5, "false") == 0) { v.b = false; pos += 5; }
        else fail("invalid literal");
        return v;
    }

    JsonValue parse_null()
    {
        if (s.compare(pos, 4, "null") != 0) fail("invalid literal");
        pos += 4;
        return JsonValue();                      // NUL
    }

    JsonValue parse_number()
    {
        const size_t start = pos;
        if (pos < s.size() && (s[pos] == '-' || s[pos] == '+')) ++pos;
        bool any = false;
        while (pos < s.size())
        {
            const char c = s[pos];
            if ((c >= '0' && c <= '9') || c == '.' || c == 'e' || c == 'E' ||
                c == '+' || c == '-') { ++pos; any = true; }
            else break;
        }
        if (!any) fail("invalid value");
        JsonValue v; v.type = JsonValue::NUM;
        v.num = std::strtod(s.c_str() + start, nullptr);
        return v;
    }
};

const JsonValue* find(const JsonValue& obj, const std::string& key)
{
    if (obj.type != JsonValue::OBJ) return nullptr;
    auto it = obj.obj.find(key);
    return (it == obj.obj.end()) ? nullptr : &it->second;
}

std::string expect_string(const JsonValue& v, const std::string& key)
{
    if (v.type != JsonValue::STR) throw std::runtime_error("run.json: key '" + key + "' must be a string");
    return v.str;
}

bool expect_bool(const JsonValue& v, const std::string& key)
{
    if (v.type != JsonValue::BOOL) throw std::runtime_error("run.json: key '" + key + "' must be a boolean");
    return v.b;
}

double expect_number(const JsonValue& v, const std::string& key)
{
    if (v.type != JsonValue::NUM) throw std::runtime_error("run.json: key '" + key + "' must be a number");
    return v.num;
}

/// Read a JSON array of exactly three integers into an array (e.g. a k-point grid).
std::array<int,3> expect_int3(const JsonValue& v, const std::string& key)
{
    if (v.type != JsonValue::ARR || v.arr.size() != 3)
        throw std::runtime_error("run.json: '" + key + "' must be an array of three integers");
    std::array<int,3> out = {{0,0,0}};
    for (int k = 0; k < 3; ++k) out[k] = (int)expect_number(v.arr[k], key + "[]");
    return out;
}

/// Parse the optional `provenance.manual` block into a ManualProvenance.
ManualProvenance parse_manual_provenance(const JsonValue& m)
{
    if (m.type != JsonValue::OBJ)
        throw std::runtime_error("run.json: 'provenance.manual' must be an object");
    ManualProvenance mp;
    mp.present = true;

    if (const JsonValue* v = find(m, "code"))  mp.code  = expect_string(*v, "provenance.manual.code");
    if (const JsonValue* v = find(m, "basis")) mp.basis = expect_string(*v, "provenance.manual.basis");
    if (const JsonValue* v = find(m, "xc_functional")) mp.xc_functional = expect_string(*v, "provenance.manual.xc_functional");
    if (const JsonValue* v = find(m, "notes")) mp.notes = expect_string(*v, "provenance.manual.notes");
    if (const JsonValue* v = find(m, "ecutwfc_Ry")) { mp.has_ecutwfc = true; mp.ecutwfc_Ry = expect_number(*v, "provenance.manual.ecutwfc_Ry"); }
    if (const JsonValue* v = find(m, "kpoint_grid")) { mp.has_kpoint_grid = true; mp.kpoint_grid = expect_int3(*v, "provenance.manual.kpoint_grid"); }

    if (const JsonValue* v = find(m, "pseudopotentials"))
    {
        if (v->type != JsonValue::ARR)
            throw std::runtime_error("run.json: 'provenance.manual.pseudopotentials' must be an array of {species,file} objects");
        for (const JsonValue& e : v->arr)
        {
            const JsonValue* sp = find(e, "species");
            const JsonValue* fl = find(e, "file");
            if (!sp || !fl)
                throw std::runtime_error("run.json: each 'pseudopotentials' entry needs \"species\" and \"file\"");
            PseudoInfo p;
            p.species = expect_string(*sp, "pseudopotentials[].species");
            p.file    = expect_string(*fl, "pseudopotentials[].file");
            if (const JsonValue* z = find(e, "z_valence")) { p.has_z_valence = true; p.z_valence = expect_number(*z, "pseudopotentials[].z_valence"); }
            mp.pseudopotentials.push_back(p);
        }
    }

    if (const JsonValue* v = find(m, "orbitals"))
    {
        if (v->type != JsonValue::ARR)
            throw std::runtime_error("run.json: 'provenance.manual.orbitals' must be an array of orbital objects");
        for (const JsonValue& e : v->arr)
        {
            ManualOrbital o;
            if (const JsonValue* iv = find(e, "index"))   o.index   = (int)expect_number(*iv, "orbitals[].index");
            if (const JsonValue* lv = find(e, "label"))   o.label   = expect_string(*lv, "orbitals[].label");
            if (const JsonValue* ev = find(e, "element")) o.element = expect_string(*ev, "orbitals[].element");
            if (const JsonValue* xv = find(e, "xyz"))
            {
                if (xv->type != JsonValue::ARR || xv->arr.size() != 3)
                    throw std::runtime_error("run.json: 'orbitals[].xyz' must be an array of three numbers");
                for (int k = 0; k < 3; ++k) o.xyz[k] = expect_number(xv->arr[k], "orbitals[].xyz[]");
            }
            mp.orbitals.push_back(o);
        }
    }
    return mp;
}

} // namespace

RunConfig read_run_config(const std::string& path)
{
    std::ifstream f(path.c_str());
    if (!f.good()) throw std::runtime_error("run.json: cannot open config file '" + path + "'");
    std::stringstream ss; ss << f.rdbuf();
    const std::string text = ss.str();

    JsonParser parser(text);
    const JsonValue root = parser.parse();
    if (root.type != JsonValue::OBJ)
        throw std::runtime_error("run.json: top-level value must be an object");

    RunConfig cfg;

    if (const JsonValue* v = find(root, "label"))        { cfg.has_label = true;       cfg.label = expect_string(*v, "label"); }
    if (const JsonValue* v = find(root, "project_dir"))  { cfg.has_project_dir = true; cfg.project_dir = expect_string(*v, "project_dir"); }
    if (const JsonValue* v = find(root, "seed"))         { cfg.has_seed = true;        cfg.seed = expect_string(*v, "seed"); }
    if (const JsonValue* v = find(root, "output_dir"))   { cfg.has_output_dir = true;  cfg.output_dir = expect_string(*v, "output_dir"); }
    if (const JsonValue* v = find(root, "mode"))         { cfg.has_mode = true;        cfg.mode = expect_string(*v, "mode"); }

    if (const JsonValue* v = find(root, "supercell"))
    {
        if (v->type != JsonValue::ARR || v->arr.size() != 3)
            throw std::runtime_error("run.json: 'supercell' must be an array of three integers");
        for (int k = 0; k < 3; ++k)
        {
            const double d = expect_number(v->arr[k], "supercell[]");
            const int n = (int)d;
            if (n < 1) throw std::runtime_error("run.json: 'supercell' dimensions must be >= 1");
            cfg.supercell[k] = n;
        }
        cfg.has_supercell = true;
    }

    if (const JsonValue* v = find(root, "operators"))
    {
        if (v->type != JsonValue::ARR) throw std::runtime_error("run.json: 'operators' must be an array of strings");
        cfg.has_operators = true;
        for (const JsonValue& e : v->arr)
        {
            const std::string op = expect_string(e, "operators[]");
            if (!W2SP_arguments::is_valid_operator(op))
                throw std::runtime_error("run.json: unknown operator '" + op + "'");
            cfg.operators.push_back(op);
        }
    }

    if (const JsonValue* v = find(root, "spin_currents"))
    {
        if (v->type != JsonValue::ARR)
            throw std::runtime_error("run.json: 'spin_currents' must be an array of [V,S] pairs");
        cfg.has_spin_currents = true;
        for (const JsonValue& e : v->arr)
        {
            if (e.type != JsonValue::ARR || e.arr.size() != 2)
                throw std::runtime_error("run.json: each 'spin_currents' entry must be a [V,S] pair, e.g. [\"X\",\"Z\"]");
            const char vd = W2SP_arguments::to_axis(expect_string(e.arr[0], "spin_currents[].V"));
            const char sd = W2SP_arguments::to_axis(expect_string(e.arr[1], "spin_currents[].S"));
            if (vd == '?' || sd == '?')
                throw std::runtime_error("run.json: 'spin_currents' axes must be X, Y or Z");
            cfg.spin_currents.push_back(std::make_pair(vd, sd));
        }
    }

    if (const JsonValue* v = find(root, "op_files"))
    {
        if (v->type != JsonValue::ARR)
            throw std::runtime_error("run.json: 'op_files' must be an array of {name,path} objects");
        cfg.has_op_files = true;
        for (const JsonValue& e : v->arr)
        {
            const JsonValue* nm = find(e, "name");
            const JsonValue* pt = find(e, "path");
            if (!nm || !pt)
                throw std::runtime_error("run.json: each 'op_files' entry needs \"name\" and \"path\"");
            cfg.op_files.push_back(std::make_pair(expect_string(*nm, "op_files[].name"),
                                                  expect_string(*pt, "op_files[].path")));
        }
    }

    if (const JsonValue* v = find(root, "exact_spin")) { cfg.has_exact_spin = true; cfg.exact_spin = expect_bool(*v, "exact_spin"); }
    // Accept both "orbital_L" (matches the CLI flag spelling) and "orbital_l".
    if (const JsonValue* v = find(root, "orbital_L")) { cfg.has_orbital_l = true; cfg.orbital_l = expect_bool(*v, "orbital_L"); }
    if (const JsonValue* v = find(root, "orbital_l")) { cfg.has_orbital_l = true; cfg.orbital_l = expect_bool(*v, "orbital_l"); }
    if (const JsonValue* v = find(root, "emit_bounds")) { cfg.has_emit_bounds = true; cfg.emit_bounds = expect_bool(*v, "emit_bounds"); }
    // "checks"/"check" selects the self-checks (same selector as --check).
    if (const JsonValue* v = find(root, "checks")) { cfg.has_checks = true; cfg.checks = expect_string(*v, "checks"); }
    else if (const JsonValue* v2 = find(root, "check")) { cfg.has_checks = true; cfg.checks = expect_string(*v2, "check"); }
    if (const JsonValue* v = find(root, "truncation_threshold")) { cfg.has_truncation = true; cfg.truncation_threshold = expect_number(*v, "truncation_threshold"); }
    if (const JsonValue* v = find(root, "log_level")) { cfg.has_log_level = true; cfg.log_level = expect_string(*v, "log_level"); }
    if (const JsonValue* v = find(root, "log_file"))  { cfg.has_log_file = true;  cfg.log_file = expect_string(*v, "log_file"); }

    if (const JsonValue* prov = find(root, "provenance"))
    {
        if (prov->type != JsonValue::OBJ) throw std::runtime_error("run.json: 'provenance' must be an object");
        if (const JsonValue* v = find(*prov, "qe_xml")) { cfg.has_qe_xml = true; cfg.qe_xml = expect_string(*v, "provenance.qe_xml"); }
        if (const JsonValue* v = find(*prov, "win"))    { cfg.has_win = true;    cfg.win = expect_string(*v, "provenance.win"); }
        if (const JsonValue* v = find(*prov, "manual")) { cfg.has_manual = true; cfg.manual = parse_manual_provenance(*v); }
    }

    return cfg;
}

void write_run_config(const W2SP_arguments& a, std::ostream& os)
{
    JsonWriter w(os);
    w.begin_object();
        w.member("label", a.label);
        w.member("mode", a.mode);
        if (!a.project_dir.empty()) w.member("project_dir", a.project_dir);
        if (!a.seed.empty())        w.member("seed", a.seed);
        w.member("output_dir", a.output_dir);

        // Supercell is meaningful in sparse mode; harmless (ignored) in bundle mode.
        w.key("supercell");
        w.begin_array();
            w.inum((long long)a.cellDim[0]);
            w.inum((long long)a.cellDim[1]);
            w.inum((long long)a.cellDim[2]);
        w.end_array();

        if (!a.operators.empty())
        {
            w.key("operators");
            w.begin_array();
            for (size_t i = 0; i < a.operators.size(); ++i) w.str(a.operators[i]);
            w.end_array();
        }

        if (!a.spin_currents.empty())
        {
            w.key("spin_currents");
            w.begin_array();
            for (size_t i = 0; i < a.spin_currents.size(); ++i)
            {
                w.begin_array();
                    w.str(std::string(1, a.spin_currents[i].first));
                    w.str(std::string(1, a.spin_currents[i].second));
                w.end_array();
            }
            w.end_array();
        }

        if (!a.op_files.empty())
        {
            w.key("op_files");
            w.begin_array();
            for (size_t i = 0; i < a.op_files.size(); ++i)
            {
                w.begin_object();
                    w.member("name", a.op_files[i].first);
                    w.member("path", a.op_files[i].second);
                w.end_object();
            }
            w.end_array();
        }

        if (a.exact_spin)      w.member("exact_spin", true);
        if (a.orbital_l)       w.member("orbital_L", true);
        if (a.emit_descriptor) w.member("emit_bounds", true);
        if (!a.check.empty())  w.member("checks", a.check);
        if (a.has_truncation)  w.member("truncation_threshold", a.truncation_threshold);

        if (!a.qe_xml_path.empty() || !a.win_path.empty() || a.manual.present)
        {
            w.key("provenance");
            w.begin_object();
                if (!a.qe_xml_path.empty()) w.member("qe_xml", a.qe_xml_path);
                if (!a.win_path.empty())    w.member("win", a.win_path);
                if (a.manual.present)     { w.key("manual"); emit_manual_provenance(w, a.manual); }
            w.end_object();
        }

        if (a.log_level != "info") w.member("log_level", a.log_level);
        if (!a.log_file.empty())   w.member("log_file", a.log_file);
    w.end_object();
    os << "\n";
}

void write_run_config_template(std::ostream& os)
{
    // Hand-written (not via JsonWriter) so it can carry `//` documentation. The
    // reader skips comments, so this scaffold parses back unchanged. Values are
    // valid defaults: edit them, delete the keys you do not need, then run with
    //   wannier2sparse <thisfile>.w2s
    os <<
"{\n"
"  // wannier2sparse input file (.w2s). JSON with // and /* */ comments allowed.\n"
"  // Run it with:  wannier2sparse <thisfile>.w2s   (or  --run <thisfile>)\n"
"\n"
"  \"label\": \"LABEL\",          // system label; default seed for LABEL_hr.dat, LABEL.uc, LABEL.xyz\n"
"  \"mode\": \"sparse\",          // \"sparse\" -> expand to supercell CSR; \"bundle\" -> primitive ops + manifest\n"
"  \"project_dir\": \".\",        // directory holding the input files\n"
"  \"seed\": \"LABEL\",           // seedname of the input files (defaults to label)\n"
"  \"output_dir\": \"out\",       // where the .CSR / bundle / .out are written\n"
"\n"
"  \"supercell\": [1, 1, 1],    // N1 N2 N3 (sparse mode only; ignored in bundle mode)\n"
"\n"
"  // Operators to build. Valid: VX VY VZ  SX SY SZ  and spin currents VXSX..VZSZ.\n"
"  // The Hamiltonian (HAM) is always written; leave this empty for HAM only.\n"
"  \"operators\": [\"VX\", \"VY\", \"SZ\"],\n"
"\n"
"  // Derived spin currents J = 1/2{V,S}, each [velocity axis, spin axis].\n"
"  \"spin_currents\": [[\"X\", \"Z\"]],\n"
"\n"
"  // External operators ingested from _hr.dat-format files (written as <LABEL>.NAME.CSR).\n"
"  \"op_files\": [{ \"name\": \"OX\", \"path\": \"ox_hr.dat\" }],\n"
"\n"
"  \"exact_spin\": false,        // build SXexact/SYexact/SZexact from <seed>.spn + <seed>_u.mat\n"
"  \"orbital_L\": false,         // build LX/LY/LZ from <seed>.amn + <seed>_u.mat + <seed>.win\n"
"  \"emit_bounds\": false,       // also write .desc sidecars with spectral bounds for HAM\n"
"  \"checks\": \"all\",            // self-checks: all|hermiticity|sum_rules|algebra|aliasing|bounds\n"
"  \"truncation_threshold\": 1e-8, // bundle-mode amplitude truncation (echoed into the manifest)\n"
"\n"
"  \"log_level\": \"info\",        // trace|debug|info|warn|error\n"
"  \"log_file\": \"\",             // explicit log path (\"\" => <output_dir>/<LABEL>.run.log)\n"
"\n"
"  \"provenance\": {\n"
"    \"qe_xml\": \"\",             // Quantum ESPRESSO data-file-schema.xml (auto-parsed DFT provenance)\n"
"    \"win\": \"\",                // Wannier90 .win (auto-parsed Wannier provenance)\n"
"\n"
"    // Manual provenance: declare here what no side-file provides. Use this when\n"
"    // the model arrives as a bare _hr.dat. Documentation only; never affects numerics.\n"
"    \"manual\": {\n"
"      \"code\": \"Quantum ESPRESSO 7.2 + Wannier90 3.1.0\",\n"
"      \"xc_functional\": \"PBE\",\n"
"      \"basis\": \"plane waves\",\n"
"      \"ecutwfc_Ry\": 90,        // plane-wave cutoff (Ry)\n"
"      \"kpoint_grid\": [12, 12, 12], // SCF k-point mesh\n"
"      \"pseudopotentials\": [\n"
"        { \"species\": \"Fe\", \"file\": \"Fe.pbe-spn-rrkjus_psl.1.0.0.UPF\", \"z_valence\": 16 }\n"
"      ],\n"
"      \"orbitals\": [            // per-index basis identity + Cartesian position (Angstrom)\n"
"        { \"index\": 0, \"label\": \"Fe-d_xy\", \"element\": \"Fe\", \"xyz\": [0.0, 0.0, 0.0] }\n"
"      ],\n"
"      \"notes\": \"\"\n"
"    }\n"
"  }\n"
"}\n";
}
