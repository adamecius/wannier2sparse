// Logger + run-stats unit test. Verifies the leveled logger's contract that the
// driver relies on: warnings and errors are ALWAYS counted (independent of the
// console threshold), the file sink records the full trace with level tags, the
// streaming proxy concatenates correctly, and the timer/peak-memory/report plumbing
// runs. The file sink is the assertion surface (console output is not captured).
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include "logger.hpp"
#include "run_stats.hpp"

using namespace std;

static string slurp(const string& p) { ifstream f(p.c_str()); stringstream ss; ss << f.rdbuf(); return ss.str(); }
static bool has(const string& hay, const string& needle) { return hay.find(needle) != string::npos; }

int main()
{
    const string logpath = "test.log";

    {
        // Console set to WARN (info suppressed on console) but the file sink still
        // records everything, and counts are independent of the console threshold.
        Logger log(Logger::WARN);
        assert(log.open_file(logpath));
        log.info("loaded model");
        log.warn("missing provenance");
        log.error("bad supercell");
        log.at(Logger::INFO) << "x=" << 42 << " y=" << 3;

        assert(log.warnings() == 1 && "exactly one warning counted");
        assert(log.errors() == 1 && "exactly one error counted");
    }

    const string body = slurp(logpath);
    assert(has(body, "loaded model"));        // INFO reached the file even though console is WARN
    assert(has(body, "missing provenance"));
    assert(has(body, "bad supercell"));
    assert(has(body, "[INFO ]") && has(body, "[WARN ]") && has(body, "[ERROR]"));
    assert(has(body, "x=42 y=3") && "LineProxy concatenates arguments");

    // Timer records a named stage; report runs end to end; peak memory is positive.
    {
        Logger log(Logger::ERROR);            // keep console quiet during the test
        RunStats stats;
        {
            ScopedTimer t(log, stats, "phase one");
            volatile double acc = 0.0;
            for (int i = 0; i < 100000; ++i) acc += i;   // a little measurable work
            (void)acc;
        }
        assert(stats.stages.size() == 1 && stats.stages[0].first == "phase one");
        assert(stats.stages[0].second >= 0.0);
        assert(RunStats::peak_rss_mb() > 0.0);
        report_run(log, stats, stats.total());           // must not throw
    }

    cout << "LOGGER AND STATS TEST PASSED" << endl;
    return 0;
}
