/**
 * @file run_stats.hpp
 * @brief Wall-clock stage timing + peak-memory accounting + end-of-run summary.
 *
 * A run should not finish silently: the user wants to see where the time went and
 * how much memory it cost. RunStats collects named (stage -> seconds) entries; the
 * ScopedTimer measures a stage by RAII (it records on destruction, so an early
 * return or an exception still books the elapsed time) and logs start/stop at
 * DEBUG. peak_rss_mb() reads the high-water resident set from getrusage. report()
 * prints the summary at INFO so it always reaches the console and the log file.
 *
 * Units/portability: getrusage's ru_maxrss is kilobytes on Linux and bytes on
 * macOS; peak_rss_mb() normalizes both to MiB. Timing uses steady_clock (monotonic,
 * immune to wall-clock adjustments). Single-threaded, no locking.
 */
#ifndef W2SP_RUN_STATS
#define W2SP_RUN_STATS

#include <string>
#include <vector>
#include <utility>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <sys/resource.h>   // getrusage (POSIX)

#include "logger.hpp"

/**
 * @brief Ordered collection of per-stage wall-clock times for the run summary.
 */
struct RunStats
{
    std::vector<std::pair<std::string, double> > stages;   ///< (name, seconds), in order

    void record(const std::string& name, double seconds) { stages.push_back(std::make_pair(name, seconds)); }

    double total() const
    {
        double s = 0.0;
        for (size_t i = 0; i < stages.size(); ++i) s += stages[i].second;
        return s;
    }

    /// Peak resident set size in MiB (high-water mark of the whole process).
    static double peak_rss_mb()
    {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
        return (double)ru.ru_maxrss / (1024.0 * 1024.0);   // bytes on macOS
#else
        return (double)ru.ru_maxrss / 1024.0;              // kilobytes on Linux
#endif
    }
};

/**
 * @brief RAII stage timer: records elapsed time into RunStats on destruction.
 *
 * Measures even when the scope is left early or unwound by an exception, so the
 * summary's per-stage times stay honest. Logs "start: <name>" / "done: <name>
 * (<s> s)" at DEBUG.
 */
class ScopedTimer
{
public:
    ScopedTimer(Logger& log, RunStats& stats, const std::string& name)
        : log_(&log), stats_(&stats), name_(name), stopped_(false),
          t0_(std::chrono::steady_clock::now())
    {
        log_->at(Logger::DEBUG) << "start: " << name_;
    }

    /// Stop and record (idempotent). Returns elapsed seconds.
    double stop()
    {
        if (stopped_) return 0.0;
        stopped_ = true;
        const double s = elapsed();
        stats_->record(name_, s);
        // Build the whole line in one stream: iomanip state does not persist
        // across the Logger LineProxy's operator<<, so format numbers here.
        std::ostringstream o;
        o << "done:  " << name_ << " (" << std::fixed << std::setprecision(3) << s << " s)";
        log_->debug(o.str());
        return s;
    }

    ~ScopedTimer() { stop(); }

    double elapsed() const
    {
        const auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(now - t0_).count();
    }

private:
    ScopedTimer(const ScopedTimer&);
    ScopedTimer& operator=(const ScopedTimer&);
    Logger*   log_;
    RunStats* stats_;
    std::string name_;
    bool      stopped_;
    std::chrono::steady_clock::time_point t0_;
};

/**
 * @brief Emit the end-of-run summary (timings, peak memory, warning/error tally).
 * @param log      logger (summary is at INFO so it always surfaces)
 * @param stats    collected per-stage times
 * @param total_s  total wall time of the run (overall timer, not the stage sum)
 */
inline void report_run(Logger& log, const RunStats& stats, double total_s)
{
    // Each line is formatted in a single ostringstream: iomanip state (setw,
    // setprecision) does not carry across the Logger LineProxy's operator<<.
    log.info("==== run summary ====");
    {
        std::ostringstream o;
        o << "total wall time: " << std::fixed << std::setprecision(3) << total_s << " s";
        log.info(o.str());
    }
    for (size_t i = 0; i < stats.stages.size(); ++i)
    {
        std::ostringstream o;
        o << "  " << std::left << std::setw(22) << stats.stages[i].first << std::right
          << std::fixed << std::setprecision(3) << stats.stages[i].second << " s";
        log.info(o.str());
    }
    {
        std::ostringstream o;
        o << "peak memory: " << std::fixed << std::setprecision(1) << RunStats::peak_rss_mb() << " MiB";
        log.info(o.str());
    }
    {
        std::ostringstream o;
        o << "warnings: " << log.warnings() << "   errors: " << log.errors();
        log.info(o.str());
    }
}

#endif
