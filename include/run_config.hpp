/**
 * @file run_config.hpp
 * @brief Config-file (`run.json`) driver for the bundle output mode.
 *
 * `RunConfig` is the parsed form of a small `run.json` that describes a bundle
 * run without the positional CLI: which model to read, which operators to build,
 * and where the DFT/Wannier provenance sources live. It is read by a minimal,
 * dependency-free JSON parser (see run_config.cpp) so the format stays golden-
 * stable and adds no heavy dependency. `apply_to` projects the config onto a
 * W2SP_arguments so the existing run_bundle_mode path is reused verbatim.
 *
 * Example run.json:
 * @code
 * {
 *   "label": "graphene", "project_dir": ".", "seed": "graphene", "output_dir": "out",
 *   "mode": "bundle", "operators": ["VX", "VY", "SZ"],
 *   "exact_spin": false, "orbital_L": false, "emit_bounds": true,
 *   "truncation_threshold": 1e-8,
 *   "provenance": { "qe_xml": "scf.save/data-file-schema.xml", "win": "graphene.win" }
 * }
 * @endcode
 *
 * Every field is optional: a key that is absent leaves the corresponding
 * `has_*` flag false, so apply_to() only overrides the W2SP_arguments defaults
 * for keys the file actually set.
 */
#ifndef RUN_CONFIG
#define RUN_CONFIG

#include <string>
#include <vector>

#include "w2sp_arguments.hpp"

/**
 * @brief Parsed `run.json` for bundle mode.
 */
struct RunConfig
{
    // Each field carries a presence flag so apply_to() only overrides defaults
    // for keys the file actually set.
    bool        has_label;        std::string label;
    bool        has_project_dir;  std::string project_dir;
    bool        has_seed;         std::string seed;
    bool        has_output_dir;   std::string output_dir;
    bool        has_mode;         std::string mode;
    bool        has_operators;    std::vector<std::string> operators;
    bool        has_exact_spin;   bool        exact_spin;
    bool        has_orbital_l;    bool        orbital_l;
    bool        has_emit_bounds;  bool        emit_bounds;
    bool        has_truncation;   double      truncation_threshold;
    bool        has_qe_xml;       std::string qe_xml;
    bool        has_win;          std::string win;

    RunConfig()
        : has_label(false), has_project_dir(false), has_seed(false),
          has_output_dir(false), has_mode(false), has_operators(false),
          has_exact_spin(false), exact_spin(false), has_orbital_l(false),
          orbital_l(false), has_emit_bounds(false), emit_bounds(false),
          has_truncation(false), truncation_threshold(0.0),
          has_qe_xml(false), has_win(false) {}

    /**
     * @brief Project the config onto a W2SP_arguments, overriding only set keys.
     * @param a arguments object to populate (its defaults remain for unset keys)
     *
     * After this call `a` drives the existing run_bundle_mode exactly as an
     * equivalent positional `--mode bundle` invocation would.
     */
    void apply_to(W2SP_arguments& a) const
    {
        if (has_label)       a.label       = label;
        if (has_project_dir) a.project_dir = project_dir;
        if (has_seed)        a.seed        = seed;
        if (has_output_dir)  a.output_dir  = output_dir;
        if (has_mode)        a.mode        = mode;
        if (has_operators)   a.operators   = operators;
        if (has_exact_spin)  a.exact_spin  = exact_spin;
        if (has_orbital_l)   a.orbital_l   = orbital_l;
        if (has_emit_bounds) a.emit_descriptor = emit_bounds;
        if (has_truncation)  { a.truncation_threshold = truncation_threshold; a.has_truncation = true; }
        if (has_qe_xml)      a.qe_xml_path = qe_xml;
        if (has_win)         a.win_path    = win;
    }
};

/**
 * @brief Parse a `run.json` file into a RunConfig.
 * @param path path to the JSON config file
 * @return parsed config
 * @throws std::runtime_error if the file cannot be opened, is not valid JSON, or
 *         a known key has the wrong type / an unknown operator name.
 */
RunConfig read_run_config(const std::string& path);

#endif
