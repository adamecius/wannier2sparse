/**
 * @file logger.hpp
 * @brief Minimal, dependency-free leveled logger for tracking a run.
 *
 * The tool needs a single, traceable record of what a simulation did: which files
 * it loaded, which operators it built, and, above all, every WARNING and ERROR.
 * Rather than pull a heavy logging dependency (spdlog and friends bundle fmt and a
 * few hundred KB of headers), this is a small hand-rolled logger in the spirit of
 * the repo's json_writer: leveled messages with timestamps, a console sink
 * (stdout for INFO/below, stderr for WARN/ERROR) and an optional file sink, and a
 * running count of warnings and errors so the end-of-run summary can report them.
 *
 * Design notes:
 *  - Levels TRACE < DEBUG < INFO < WARN < ERROR < SILENT. The console threshold is
 *    configurable (--quiet / --verbose); the file sink records everything from its
 *    own (lower) threshold so the on-disk trace is complete even when the console
 *    is quiet.
 *  - WARN and ERROR are ALWAYS counted, regardless of the console threshold, so a
 *    run that prints nothing still reports "warnings: N  errors: M" at the end.
 *  - Two call styles: log.info("msg") for fixed strings, and a streaming proxy
 *    log.at(Logger::INFO) << "x=" << x for built messages. The proxy accumulates
 *    into a std::string (not a moved std::ostringstream) so it builds on the
 *    GCC 4.8 baseline this project targets.
 *  - Single-threaded tool: no locking. POSIX isatty/localtime_r are used (the tool
 *    targets Linux/macOS, matching the rest of the codebase).
 */
#ifndef W2SP_LOGGER
#define W2SP_LOGGER

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <unistd.h>   // isatty, fileno (POSIX)

class Logger
{
public:
    enum Level { TRACE = 0, DEBUG = 1, INFO = 2, WARN = 3, ERROR = 4, SILENT = 5 };

    /**
     * @brief Construct with a console verbosity threshold.
     * @param console_level messages at or above this level reach the console
     */
    explicit Logger(Level console_level = INFO)
        : console_level_(console_level), file_level_(TRACE),
          n_warn_(0), n_err_(0), color_(::isatty(::fileno(stderr)) != 0) {}

    ~Logger() { if (file_.is_open()) file_.close(); }

    void set_console_level(Level l) { console_level_ = l; }

    /**
     * @brief Open (truncating) a log file sink.
     * @param path destination path
     * @return true if the file was opened
     *
     * The file sink records from file_level_ (TRACE by default) so the on-disk log
     * is the complete, traceable record even when the console is quiet.
     */
    bool open_file(const std::string& path)
    {
        file_.open(path.c_str(), std::ios::out | std::ios::trunc);
        if (file_.good()) { file_path_ = path; return true; }
        return false;
    }

    const std::string& file_path() const { return file_path_; }
    int warnings() const { return n_warn_; }
    int errors()   const { return n_err_; }

    /**
     * @brief Emit one message at a level. Counts WARN/ERROR unconditionally.
     */
    void log(Level lvl, const std::string& msg)
    {
        if (lvl == WARN)  ++n_warn_;
        if (lvl == ERROR) ++n_err_;
        const std::string line = format(lvl, msg);
        if (lvl >= console_level_)
        {
            std::ostream& os = (lvl >= WARN) ? std::cerr : std::cout;
            if (color_) os << color_for(lvl) << line << "\033[0m\n";
            else        os << line << "\n";
        }
        if (file_.is_open() && lvl >= file_level_) { file_ << line << "\n"; file_.flush(); }
    }

    void trace(const std::string& m) { log(TRACE, m); }
    void debug(const std::string& m) { log(DEBUG, m); }
    void info (const std::string& m) { log(INFO,  m); }
    void warn (const std::string& m) { log(WARN,  m); }
    void error(const std::string& m) { log(ERROR, m); }

    /**
     * @brief Streaming proxy: log.at(INFO) << "x=" << x; flushes on destruction.
     *
     * Each operator<< stringifies its argument independently and appends to a
     * std::string buffer (chosen over a moved std::ostringstream for GCC 4.8
     * portability). Consequence: iomanip manipulators (std::setw, std::setprecision,
     * std::fixed) do NOT persist across <<. For formatted numbers, build the line in
     * a local std::ostringstream and pass the string to info()/warn()/etc.
     */
    class LineProxy
    {
    public:
        LineProxy(Logger* lg, Level lvl) : lg_(lg), lvl_(lvl) {}
        LineProxy(LineProxy&& o) : lg_(o.lg_), lvl_(o.lvl_), buf_(std::move(o.buf_)) { o.lg_ = 0; }
        ~LineProxy() { if (lg_) lg_->log(lvl_, buf_); }
        template <class T> LineProxy& operator<<(const T& v)
        {
            std::ostringstream o; o << v; buf_ += o.str(); return *this;
        }
    private:
        LineProxy(const LineProxy&);            // non-copyable
        LineProxy& operator=(const LineProxy&);
        Logger*     lg_;
        Level       lvl_;
        std::string buf_;
    };

    LineProxy at(Level l) { return LineProxy(this, l); }

private:
    Level         console_level_;
    Level         file_level_;
    int           n_warn_, n_err_;
    bool          color_;
    std::ofstream file_;
    std::string   file_path_;

    static const char* level_name(Level l)
    {
        switch (l)
        {
            case TRACE: return "TRACE";
            case DEBUG: return "DEBUG";
            case INFO:  return "INFO ";
            case WARN:  return "WARN ";
            case ERROR: return "ERROR";
            default:    return "?????";
        }
    }

    static const char* color_for(Level l)
    {
        switch (l)
        {
            case WARN:  return "\033[33m";   // yellow
            case ERROR: return "\033[31m";   // red
            case DEBUG: return "\033[2m";    // dim
            case TRACE: return "\033[2m";    // dim
            default:    return "\033[0m";    // default
        }
    }

    static std::string timestamp()
    {
        using namespace std::chrono;
        const auto now = system_clock::now();
        const std::time_t t = system_clock::to_time_t(now);
        const long ms = (long)(duration_cast<milliseconds>(now.time_since_epoch()).count() % 1000);
        std::tm tmv;
        localtime_r(&t, &tmv);
        char buf[32];
        std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tmv);
        std::ostringstream o;
        o << buf << "." << std::setfill('0') << std::setw(3) << ms;
        return o.str();
    }

    static std::string format(Level lvl, const std::string& msg)
    {
        std::ostringstream o;
        o << "[" << timestamp() << "] [" << level_name(lvl) << "] " << msg;
        return o.str();
    }
};

#endif
