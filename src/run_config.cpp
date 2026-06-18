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
        while (pos < s.size())
        {
            const char c = s[pos];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') ++pos;
            else break;
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

    if (const JsonValue* v = find(root, "exact_spin")) { cfg.has_exact_spin = true; cfg.exact_spin = expect_bool(*v, "exact_spin"); }
    // Accept both "orbital_L" (matches the CLI flag spelling) and "orbital_l".
    if (const JsonValue* v = find(root, "orbital_L")) { cfg.has_orbital_l = true; cfg.orbital_l = expect_bool(*v, "orbital_L"); }
    if (const JsonValue* v = find(root, "orbital_l")) { cfg.has_orbital_l = true; cfg.orbital_l = expect_bool(*v, "orbital_l"); }
    if (const JsonValue* v = find(root, "emit_bounds")) { cfg.has_emit_bounds = true; cfg.emit_bounds = expect_bool(*v, "emit_bounds"); }
    if (const JsonValue* v = find(root, "truncation_threshold")) { cfg.has_truncation = true; cfg.truncation_threshold = expect_number(*v, "truncation_threshold"); }

    if (const JsonValue* prov = find(root, "provenance"))
    {
        if (prov->type != JsonValue::OBJ) throw std::runtime_error("run.json: 'provenance' must be an object");
        if (const JsonValue* v = find(*prov, "qe_xml")) { cfg.has_qe_xml = true; cfg.qe_xml = expect_string(*v, "provenance.qe_xml"); }
        if (const JsonValue* v = find(*prov, "win"))    { cfg.has_win = true;    cfg.win = expect_string(*v, "provenance.win"); }
    }

    return cfg;
}
