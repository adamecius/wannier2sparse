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
 * One input file drives either output mode. A sparse (supercell CSR) run:
 * @code
 * {
 *   "label": "graphene", "mode": "sparse", "output_dir": "out",
 *   "supercell": [80, 80, 1],
 *   "operators": ["VX", "VY", "SZ"],
 *   "spin_currents": [["X","Z"]],
 *   "op_files": [{"name": "OX", "path": "ox_hr.dat"}],
 *   "checks": "all", "emit_bounds": true,
 *   "log_level": "info"
 * }
 * @endcode
 * A bundle (primitive operators + provenance manifest) run:
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
#include <array>
#include <utility>
#include <ostream>

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
    bool        has_supercell;    std::array<int,3> supercell;            ///< sparse-mode N1 N2 N3
    bool        has_operators;    std::vector<std::string> operators;
    bool        has_spin_currents; std::vector<std::pair<char,char> > spin_currents; ///< derived J=1/2{V,S}
    bool        has_op_files;     std::vector<std::pair<std::string,std::string> > op_files; ///< (NAME, PATH)
    bool        has_exact_spin;   bool        exact_spin;
    bool        has_orbital_l;    bool        orbital_l;
    bool        has_emit_bounds;  bool        emit_bounds;
    bool        has_checks;       std::string checks;                     ///< --check selector
    bool        has_truncation;   double      truncation_threshold;
    bool        has_qe_xml;       std::string qe_xml;
    bool        has_win;          std::string win;
    bool        has_manual;       ManualProvenance manual;  ///< user-declared provenance block
    bool        has_log_level;    std::string log_level;
    bool        has_log_file;     std::string log_file;

    RunConfig()
        : has_label(false), has_project_dir(false), has_seed(false),
          has_output_dir(false), has_mode(false), has_supercell(false),
          supercell({{1,1,1}}), has_operators(false), has_spin_currents(false),
          has_op_files(false), has_exact_spin(false), exact_spin(false),
          has_orbital_l(false), orbital_l(false), has_emit_bounds(false),
          emit_bounds(false), has_checks(false), has_truncation(false),
          truncation_threshold(0.0), has_qe_xml(false), has_win(false),
          has_manual(false), has_log_level(false), has_log_file(false) {}

    /**
     * @brief Project the config onto a W2SP_arguments, overriding only set keys.
     * @param a arguments object to populate (its defaults remain for unset keys)
     *
     * After this call `a` drives the existing run_bundle_mode exactly as an
     * equivalent positional `--mode bundle` invocation would.
     */
    void apply_to(W2SP_arguments& a) const
    {
        if (has_label)        a.label        = label;
        if (has_project_dir)  a.project_dir  = project_dir;
        if (has_seed)         a.seed         = seed;
        if (has_output_dir)   a.output_dir   = output_dir;
        if (has_mode)         a.mode         = mode;
        if (has_supercell)    a.cellDim      = supercell;
        if (has_operators)    a.operators    = operators;
        if (has_spin_currents) a.spin_currents = spin_currents;
        if (has_op_files)     a.op_files     = op_files;
        if (has_exact_spin)   a.exact_spin   = exact_spin;
        if (has_orbital_l)    a.orbital_l    = orbital_l;
        if (has_emit_bounds)  a.emit_descriptor = emit_bounds;
        if (has_checks)       a.check        = checks;
        if (has_truncation)   { a.truncation_threshold = truncation_threshold; a.has_truncation = true; }
        if (has_qe_xml)       a.qe_xml_path  = qe_xml;
        if (has_win)          a.win_path     = win;
        if (has_manual)       a.manual       = manual;
        if (has_log_level)    a.log_level    = log_level;
        if (has_log_file)     a.log_file     = log_file;
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

/**
 * @brief Serialize a fully-resolved W2SP_arguments as a JSON input file.
 * @param a  arguments (typically from the positional CLI)
 * @param os destination stream
 *
 * Emits a `run.json`-shaped document that, fed back through read_run_config +
 * apply_to (i.e. `--run`), reproduces the same run. This backs the `--write`
 * scaffolding command: the positional CLI becomes a generator for the input file
 * rather than an execution path. Only set / non-default keys are written, in a
 * fixed deterministic key order.
 */
void write_run_config(const W2SP_arguments& a, std::ostream& os);

/**
 * @brief Write a fully-annotated, self-documenting `.w2s` input template.
 * @param os destination stream
 *
 * Emits every supported key with a sensible default value and an inline `//`
 * comment describing it (units included). The output is itself a valid `.w2s`:
 * the reader tolerates `//` and `/* * /` comments, so the scaffold round-trips
 * through read_run_config unchanged. Backs the `--create-template` command.
 */
void write_run_config_template(std::ostream& os);

#endif
