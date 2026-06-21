/**
 * @file run_report.hpp
 * @brief The `.out` run receipt: a machine-readable record of one simulation.
 *
 * A run should leave behind more than scattered CSR files: it should leave a single
 * document that says *what was produced, how, and at what cost*. `<LABEL>.out` is
 * that document — a JSON receipt written next to the outputs at the end of every
 * run. It complements (does not replace) the human console / `.run.log` summary
 * produced by run_stats.hpp.
 *
 * The receipt records four things:
 *   - **invocation**: the resolved label/mode, the input file that drove the run
 *     (when run from a `.w2s`), the output directory, the command line, the tool
 *     version, and a start timestamp;
 *   - **steps**: the per-stage wall-clock times collected by RunStats, in order.
 *     Every stage that began (and unwound, even on error) is listed, so a failed
 *     run still shows how far it got;
 *   - **resources**: total wall time, peak resident memory (MiB), and the
 *     warning / error tally from the logger;
 *   - **outputs**: the ledger — one entry per operator written, carrying its label,
 *     on-disk path, observable / component / units, the build provenance string,
 *     and the list of *input files that produced it*. This is the provenance link
 *     the maintainer asked for: given an operator, the receipt names exactly which
 *     `_hr.dat` / `.spn` / `.amn` / op-file went into it.
 *
 * Header-only and inline, matching logger.hpp / run_stats.hpp, so it adds no
 * translation unit. The JSON layout is deterministic (fixed key order) so two runs
 * with the same inputs differ only in the timing / memory / timestamp fields.
 */
#ifndef W2SP_RUN_REPORT
#define W2SP_RUN_REPORT

#include <string>
#include <vector>
#include <ostream>

#include "json_writer.hpp"
#include "run_stats.hpp"
#include "system_provenance.hpp"
#include "provenance_writer.hpp"

/**
 * @brief One written operator, with the input files that produced it.
 */
struct OutputRecord
{
    std::string              name;        ///< operator name, e.g. "HAM", "VX", "SXexact"
    std::string              path;        ///< on-disk path (CSR file, or the bundle dir)
    std::string              observable;  ///< "hamiltonian", "velocity", "spin", ...
    std::string              component;   ///< "X", "XSZ", ...
    std::string              units;       ///< "eV", "eV*Angstrom", "hbar/2", ...
    std::string              provenance;  ///< how it was built (free text)
    std::vector<std::string> sources;     ///< input files consumed to build it
};

/**
 * @brief Accumulated run-level facts for the `.out` receipt.
 *
 * Populated incrementally by main as it resolves arguments and writes operators,
 * then serialized once by write_run_report().
 */
struct RunReport
{
    std::string label;
    std::string mode;
    std::string input_file;     ///< the `.w2s` that drove the run, or "" for the positional CLI
    std::string output_dir;
    std::string command_line;
    std::string started_at;     ///< local timestamp, "YYYY-MM-DD HH:MM:SS"
    std::string version;
    bool        success;
    ManualProvenance manual;    ///< user-declared provenance, echoed for traceability
    std::vector<OutputRecord> outputs;

    RunReport() : success(false) {}

    /// Append one output operator with its consumed input files.
    void add_output(const std::string& name, const std::string& path,
                    const std::string& observable, const std::string& component,
                    const std::string& units, const std::string& provenance,
                    const std::vector<std::string>& sources)
    {
        OutputRecord r;
        r.name = name; r.path = path; r.observable = observable;
        r.component = component; r.units = units; r.provenance = provenance;
        r.sources = sources;
        outputs.push_back(r);
    }
};

/**
 * @brief Serialize the run receipt as JSON to a stream.
 * @param os       destination stream
 * @param rep      run-level facts + output ledger
 * @param stats    per-stage timings (the "steps" array)
 * @param warnings logger warning count
 * @param errors   logger error count
 * @param total_s  total wall time of the run
 */
inline void write_run_report(std::ostream& os, const RunReport& rep, const RunStats& stats,
                             int warnings, int errors, double total_s)
{
    JsonWriter w(os);
    w.begin_object();
        w.member("schema_version", std::string("1.0"));
        w.key("generator");
        w.begin_object();
            w.member("tool", std::string("wannier2sparse"));
            w.member("version", rep.version);
        w.end_object();

        w.key("invocation");
        w.begin_object();
            w.member("label", rep.label);
            w.member("mode", rep.mode);
            if (!rep.input_file.empty()) w.member("input_file", rep.input_file);
            else                         w.member_null("input_file");
            w.member("output_dir", rep.output_dir);
            w.member("command_line", rep.command_line);
            w.member("started_at", rep.started_at);
        w.end_object();

        w.member("status", std::string(rep.success ? "success" : "failed"));

        w.key("steps");
        w.begin_array();
        for (size_t i = 0; i < stats.stages.size(); ++i)
        {
            w.begin_object();
                w.member("name", stats.stages[i].first);
                w.member("seconds", stats.stages[i].second);
            w.end_object();
        }
        w.end_array();

        w.key("resources");
        w.begin_object();
            w.member("wall_time_s", total_s);
            w.member("peak_memory_MiB", RunStats::peak_rss_mb());
            w.member("warnings", (long long)warnings);
            w.member("errors", (long long)errors);
        w.end_object();

        // Provenance link: each output names the input files that produced it.
        w.key("outputs");
        w.begin_array();
        for (const OutputRecord& o : rep.outputs)
        {
            w.begin_object();
                w.member("name", o.name);
                w.member("path", o.path);
                w.member("observable", o.observable);
                w.member("component", o.component);
                w.member("units", o.units);
                w.member("provenance", o.provenance);
                w.key("sources");
                w.begin_array();
                for (const std::string& s : o.sources) w.str(s);
                w.end_array();
            w.end_object();
        }
        w.end_array();

        // The user-declared provenance block, echoed so the receipt is self-contained.
        w.key("manual_provenance");
        if (rep.manual.present) emit_manual_provenance(w, rep.manual);
        else                    w.null();
    w.end_object();
    os << "\n";
}

#endif
