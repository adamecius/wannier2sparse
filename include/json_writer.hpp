/**
 * @file json_writer.hpp
 * @brief Tiny dependency-free, deterministic JSON emitter.
 *
 * A minimal streaming writer used to produce the bundle `manifest.json`. It is
 * deliberately hand-rolled (rather than pulling a heavy JSON dependency) so the
 * output is byte-deterministic: the caller controls key order, and floating-point
 * values are printed at a fixed precision. The writer tracks indentation and
 * comma placement; the caller is responsible for balancing begin/end calls.
 *
 * Usage:
 *   JsonWriter w(os);
 *   w.begin_object();
 *     w.key("label"); w.str("graphene");
 *     w.key("dims");  w.begin_array(); w.inum(1); w.inum(2); w.end_array();
 *   w.end_object();
 */
#ifndef JSON_WRITER
#define JSON_WRITER

#include <ostream>
#include <string>
#include <vector>
#include <array>
#include <limits>

/**
 * @brief Streaming JSON writer with deterministic layout.
 */
class JsonWriter
{
public:
    /**
     * @brief Construct around an output stream.
     * @param o destination stream (must outlive the writer)
     */
    explicit JsonWriter(std::ostream& o) : os(o), pending(false)
    {
        os.precision(std::numeric_limits<double>::digits10 + 2);
    }

    void begin_object() { begin_value(); os << "{"; push('{'); }   ///< open an object `{`
    void end_object()   { pop('}'); }                              ///< close an object `}`
    void begin_array()  { begin_value(); os << "["; push('['); }   ///< open an array `[`
    void end_array()    { pop(']'); }                              ///< close an array `]`

    /**
     * @brief Emit an object member key. The next value call fills it.
     * @param k member name
     */
    void key(const std::string& k)
    {
        Frame& f = stack.back();
        if (f.count++ > 0) os << ",";
        os << "\n";
        indent(stack.size());
        os << "\"" << escape(k) << "\": ";
        pending = true;
    }

    void str(const std::string& v) { begin_value(); os << "\"" << escape(v) << "\""; } ///< string value
    void num(double v)             { begin_value(); os << v; }                         ///< number value
    void inum(long long v)         { begin_value(); os << v; }                         ///< integer value
    void boolean(bool v)           { begin_value(); os << (v ? "true" : "false"); }    ///< boolean value
    void null()                    { begin_value(); os << "null"; }                    ///< null value

    /// Convenience: write an object member in one call.
    void member(const std::string& k, const std::string& v) { key(k); str(v); }
    void member(const std::string& k, long long v)          { key(k); inum(v); }
    void member(const std::string& k, double v)             { key(k); num(v); }
    void member(const std::string& k, bool v)               { key(k); boolean(v); }
    void member_null(const std::string& k)                  { key(k); null(); }

    /// Convenience: write an array of doubles as a JSON array value.
    void array_d(const std::array<double, 3>& v)
    {
        begin_array(); for (double x : v) num(x); end_array();
    }
    /// Convenience: write a 3x3 matrix as a JSON array of arrays.
    void matrix3(const std::array<std::array<double, 3>, 3>& m)
    {
        begin_array(); for (const auto& row : m) array_d(row); end_array();
    }

private:
    struct Frame { char type; int count; };
    std::ostream&      os;
    std::vector<Frame> stack;
    bool               pending;   ///< true after key(): the next value fills the member

    void begin_value()
    {
        if (pending) { pending = false; return; }   // object member value: key() did the layout
        if (!stack.empty())                         // array element: emit our own separator
        {
            Frame& f = stack.back();
            if (f.count++ > 0) os << ",";
            os << "\n";
            indent(stack.size());
        }
    }

    void push(char t) { stack.push_back(Frame{t, 0}); }

    void pop(char close)
    {
        const Frame f = stack.back();
        stack.pop_back();
        if (f.count > 0) { os << "\n"; indent(stack.size()); }   // non-empty: close on its own line
        os << close;
    }

    void indent(size_t depth) { for (size_t i = 0; i < depth; ++i) os << "  "; }

    static std::string escape(const std::string& s)
    {
        std::string out;
        out.reserve(s.size());
        for (char c : s)
        {
            switch (c)
            {
                case '"':  out += "\\\""; break;
                case '\\': out += "\\\\"; break;
                case '\n': out += "\\n";  break;
                case '\t': out += "\\t";  break;
                case '\r': out += "\\r";  break;
                default:   out += c;      break;
            }
        }
        return out;
    }
};

#endif
