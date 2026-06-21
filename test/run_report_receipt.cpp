// .out run-receipt test. Verifies write_run_report serializes the four sections the
// receipt promises: invocation (label/mode/input file), steps (per-stage timings),
// resources (wall time / peak memory / warning+error tally), and the outputs ledger
// — each output naming the input files that produced it (the provenance link), plus
// the echoed manual-provenance block. The serialized string is the assertion
// surface; this is a pure unit test (no model load, no fixtures).
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include "run_report.hpp"
#include "run_stats.hpp"

using namespace std;

static bool has(const string& hay, const string& needle) { return hay.find(needle) != string::npos; }

int main()
{
    RunReport rep;
    rep.label = "gr"; rep.mode = "sparse"; rep.input_file = "gr.w2s";
    rep.output_dir = "out"; rep.command_line = "wannier2sparse gr.w2s";
    rep.started_at = "2026-01-01 00:00:00"; rep.version = "1.0";
    rep.success = true;

    rep.manual.present = true;
    rep.manual.code = "QE 7.2";
    rep.manual.has_kpoint_grid = true; rep.manual.kpoint_grid = {{8, 8, 1}};

    rep.add_output("HAM", "out/gr.HAM.CSR", "hamiltonian", "", "eV",
                   "wannier90 _hr.dat", std::vector<string>{"gr_hr.dat"});
    rep.add_output("VX", "out/gr.VX.CSR", "velocity", "X", "eV*Angstrom",
                   "H_ij * (-i) dr", std::vector<string>{"gr_hr.dat", "gr.uc", "gr.xyz"});

    RunStats stats;
    stats.record("load model", 0.01);
    stats.record("write HAM", 0.02);

    ostringstream os;
    write_run_report(os, rep, stats, /*warnings*/ 1, /*errors*/ 0, /*total_s*/ 0.05);
    const string body = os.str();

    // Invocation + status.
    assert(has(body, "\"label\": \"gr\""));
    assert(has(body, "\"mode\": \"sparse\""));
    assert(has(body, "\"input_file\": \"gr.w2s\""));
    assert(has(body, "\"status\": \"success\""));

    // Steps carry the recorded stage names.
    assert(has(body, "\"name\": \"load model\""));
    assert(has(body, "\"name\": \"write HAM\""));

    // Resources: the warning tally and a memory figure are present.
    assert(has(body, "\"warnings\": 1"));
    assert(has(body, "\"errors\": 0"));
    assert(has(body, "peak_memory_MiB"));

    // Outputs ledger: each operator names its producing input files.
    assert(has(body, "\"name\": \"HAM\"") && has(body, "out/gr.HAM.CSR"));
    assert(has(body, "\"name\": \"VX\""));
    assert(has(body, "gr_hr.dat") && has(body, "gr.uc") && has(body, "gr.xyz"));

    // Manual provenance is echoed.
    assert(has(body, "\"manual_provenance\"") && has(body, "QE 7.2") && has(body, "\"kpoint_grid\""));

    // A failed run is reported as such.
    {
        RunReport bad; bad.label = "x"; bad.mode = "sparse"; bad.version = "1.0"; bad.success = false;
        RunStats s2; ostringstream o2;
        write_run_report(o2, bad, s2, 0, 1, 0.0);
        assert(has(o2.str(), "\"status\": \"failed\""));
    }

    cout << "RUN REPORT RECEIPT TEST PASSED" << endl;
    return 0;
}
