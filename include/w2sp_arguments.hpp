/**
 * @file w2sp_arguments.hpp
 * @brief Command-line interface parser for wannier2sparse.
 *
 *   wannier2sparse LABEL N1 N2 N3 [OP ... | all] [options]
 *
 * LABEL      system label; the tool reads LABEL.uc, LABEL.xyz and LABEL_hr.dat.
 * N1 N2 N3   supercell dimensions along each lattice vector (integers >= 1).
 * OP ...     operators to generate (e.g. VX SZ); 'all' generates every operator.
 *            If omitted, only the Hamiltonian (LABEL.HAM.CSR) is written.
 *
 * parse() never aborts: invalid input prints a clear message and returns a
 * non-zero status so the caller decides how to exit.  The list of valid
 * operators below reflects exactly what this build of tbmodel can generate
 * (velocity, spin density and spin current; there is no torque operator here).
 */
#ifndef W2SP_ARGUMENTS
#define W2SP_ARGUMENTS

#include <array>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cctype>

#include "system_provenance.hpp"   // ManualProvenance (user-declared provenance block)

/**
 * @brief Command-line argument state and parser.
 */
class W2SP_arguments
{
public:
    enum Status { PROCEED = 0, EXIT_OK = 1, EXIT_ERROR = 2 }; ///< parse() outcome

    std::array<int, 3>       cellDim;       ///< supercell dimensions
    std::string              label;         ///< positional system label
    std::string              output_dir;    ///< -o/--output-dir
    std::string              project_dir;   ///< directory holding the input files (default: cwd)
    std::string              seed;          ///< seedname of the input files (default: label)
    std::vector<std::string> operators;     ///< explicitly requested operators
    /**
     * @brief External operators to ingest from _hr.dat-format files: (NAME, PATH).
     *
     * Each is expanded through the same engine and written as <prefix>.NAME.CSR.
     */
    std::vector<std::pair<std::string, std::string> > op_files;
    /**
     * @brief Derived spin currents J = 1/2{V,S}: (velocity axis, spin axis).
     *
     * Each axis is in {X,Y,Z}. Written as <prefix>.J<V>S<S>.CSR.
     */
    std::vector<std::pair<char, char> > spin_currents;
    bool                     emit_descriptor; ///< --bounds: write .desc sidecars
    bool                     exact_spin;      ///< --exact-spin: use gauge transform
    bool                     orbital_l;       ///< --orbital-L: build orbital angular momentum
    std::string              check;           ///< --check selector
    std::string              mode;            ///< --mode: "sparse" (default) or "bundle"
    std::string              config_path;     ///< -x/--input-file: path to the .w2s input file
    bool                     has_truncation;  ///< truncation_threshold was set (from run.json)
    double                   truncation_threshold; ///< bundle truncation threshold (manifest)
    std::string              qe_xml_path;     ///< QE data-file-schema.xml for DFT provenance
    std::string              win_path;        ///< .win for Wannier provenance
    ManualProvenance         manual;          ///< user-declared provenance (input file only)
    std::string              log_level;       ///< console verbosity: trace|debug|info|warn|error
    std::string              log_file;        ///< explicit log-file path ("" => default location)
    bool                     no_log_file;     ///< --no-log-file: console only, no file sink
    bool                     write_config;    ///< --write: scaffold an input file from the CLI, don't run
    bool                     create_template; ///< --create-template: write an annotated .w2s template
    bool                     create_from_string; ///< --create "...": scaffold a .w2s from a CLI line
    std::string              create_string;   ///< the quoted CLI line passed to --create
    std::string              create_inp;      ///< -inp NAME: target stem for --create / --create-template
    std::string              program_name;    ///< argv[0]
    bool                     covariant_velocity; ///< legacy flag == (velocity_mode == "covariant")
    /**
     * @brief Velocity construction mode (the velocity ladder), applied to VX/VY/VZ
     *        and to the velocity inside any spin current J=1/2{V,S}:
     *   - "bare"            : v = -i(R.lat)H            (pure gradient, no positions)
     *   - "berry_connection": v = -i(R.lat + dr_ij)H    (adds Wannier centres; default)
     *   - "covariant"       : v = -i(R.lat)H - i[H,A]   (full Berry connection from _r.dat)
     */
    std::string              velocity_mode;
    std::string              r_dat_path;      ///< override path to the position matrix _r.dat
    std::vector<KpathNode>   kpoint_path;     ///< band high-symmetry path (provenance), if recorded
    std::string              kpoint_path_source; ///< where it came from: "win" | "qe_bands" | "w2s"
    bool                     make_provenance; ///< --provenance: extract provenance into a .w2s and exit
    std::string              qe_bands_path;   ///< --qe-bands: QE bands.in for the crystal_b k-path

    /**
     * @brief Default constructor.
     */
    W2SP_arguments()
        : cellDim({{1, 1, 1}}), output_dir("."), emit_descriptor(false),
          exact_spin(false), orbital_l(false), mode("sparse"),
          has_truncation(false), truncation_threshold(0.0),
          log_level("info"), no_log_file(false), write_config(false),
          create_template(false), create_from_string(false),
          program_name("wannier2sparse"),
          covariant_velocity(false), velocity_mode("berry_connection"),
          make_provenance(false) {}

    /**
     * @brief Whether a string ends with a given suffix.
     * @param s   string to test
     * @param suf suffix
     * @return true if s ends with suf
     */
    static bool ends_with(const std::string& s, const std::string& suf)
    {
        return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
    }

    /**
     * @brief Resolved input file stem: <project_dir>/<seed>.
     * @return input prefix
     *
     * seed defaults to the positional LABEL and project_dir defaults to the current
     * directory. Input files are <prefix>_hr.dat, <prefix>.uc, <prefix>.xyz.
     */
    std::string input_prefix() const
    {
        const std::string base = seed.empty() ? label : seed;
        return project_dir.empty() ? base : (project_dir + "/" + base);
    }

    /**
     * @brief Operators this tool knows how to build (HAM is always written).
     * @return reference to the operator list
     */
    static const std::vector<std::string>& available_operators()
    {
        static const std::vector<std::string> ops = {
            "VX", "VY", "VZ",
            "SX", "SY", "SZ",
            "VXSX", "VXSY", "VXSZ",
            "VYSX", "VYSY", "VYSZ",
            "VZSX", "VZSY", "VZSZ"
        };
        return ops;
    }

    /**
     * @brief Validate an operator code.
     * @param op operator string
     * @return true if op is in available_operators()
     */
    static bool is_valid_operator(const std::string& op)
    {
        for (const auto& a : available_operators())
            if (op == a) return true;
        return false;
    }

    /**
     * @brief Normalize a single-letter cartesian axis.
     * @param s input string
     * @return upper-case axis 'X'/'Y'/'Z', or '?' if invalid
     */
    static char to_axis(const std::string& s)
    {
        if (s.size() != 1) return '?';
        const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
        return (c == 'X' || c == 'Y' || c == 'Z') ? c : '?';
    }

    /**
     * @brief Whether a token is a plain (optionally signed) integer literal.
     * @param s input string
     * @return true if s parses entirely as an integer
     */
    static bool looks_like_int(const std::string& s)
    {
        if (s.empty()) return false;
        size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
        if (i == s.size()) return false;
        for (; i < s.size(); ++i)
            if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
        return true;
    }

    /**
     * @brief Parse command-line arguments into this object.
     * @param argc argument count
     * @param argv argument vector
     * @return PROCEED, EXIT_OK, or EXIT_ERROR
     */
    Status parse(int argc, char* argv[])
    {
        if (argc > 0) program_name = argv[0];

        std::vector<std::string> positional;
        bool want_all = false;

        for (int i = 1; i < argc; ++i)
        {
            const std::string a = argv[i];
            if (a == "-h" || a == "--help")        { print_help();      return EXIT_OK; }
            else if (a == "--version")             { print_version();   return EXIT_OK; }
            else if (a == "--list-operators")      { print_operators(); return EXIT_OK; }
            else if (a == "-o" || a == "--output-dir")
            {
                if (i + 1 >= argc) { error("missing directory after '" + a + "'"); return EXIT_ERROR; }
                output_dir = argv[++i];
            }
            else if (a == "-p" || a == "--project")
            {
                if (i + 1 >= argc) { error("missing directory after '" + a + "'"); return EXIT_ERROR; }
                project_dir = argv[++i];
            }
            else if (a == "--seed")
            {
                if (i + 1 >= argc) { error("missing name after '--seed'"); return EXIT_ERROR; }
                seed = argv[++i];
            }
            else if (a == "--op-file")
            {
                if (i + 2 >= argc) { error("'--op-file' needs NAME and PATH"); return EXIT_ERROR; }
                const std::string nm = argv[++i];
                const std::string pth = argv[++i];
                op_files.push_back(std::make_pair(nm, pth));
            }
            else if (a == "--spin-current")
            {
                if (i + 2 >= argc) { error("'--spin-current' needs a velocity axis and a spin axis (X|Y|Z)"); return EXIT_ERROR; }
                const char vd = to_axis(argv[++i]);
                const char sd = to_axis(argv[++i]);
                if (vd == '?' || sd == '?') { error("'--spin-current' axes must be X, Y or Z"); return EXIT_ERROR; }
                spin_currents.push_back(std::make_pair(vd, sd));
            }
            else if (a == "--mode")
            {
                if (i + 1 >= argc) { error("missing value after '--mode' (sparse|bundle)"); return EXIT_ERROR; }
                mode = argv[++i];
                if (mode != "sparse" && mode != "bundle")
                {
                    error("unknown --mode '" + mode + "' (expected 'sparse' or 'bundle')");
                    return EXIT_ERROR;
                }
            }
            else if (a == "-x" || a == "--input-file")
            {
                if (i + 1 >= argc) { error("missing path after '" + a + "'"); return EXIT_ERROR; }
                config_path = argv[++i];
            }
            else if (a == "--write")               { write_config = true; }
            else if (a == "--create-template")
            {
                create_template = true;
                // Optional NAME: consume the next token unless it is an option.
                if (i + 1 < argc && argv[i+1][0] != '-') create_inp = argv[++i];
            }
            else if (a == "--create")
            {
                if (i + 1 >= argc) { error("'--create' needs a quoted command line, e.g. --create \"graphene 50 50 1 VX SZ\""); return EXIT_ERROR; }
                create_from_string = true;
                create_string = argv[++i];
            }
            else if (a == "-inp" || a == "--inp")
            {
                if (i + 1 >= argc) { error("missing name after '" + a + "'"); return EXIT_ERROR; }
                create_inp = argv[++i];
            }
            else if (a == "--provenance")
            {
                make_provenance = true;
                // Optional SEED: consume the next token unless it is an option.
                if (i + 1 < argc && argv[i+1][0] != '-') label = argv[++i];
            }
            else if (a == "--win")
            {
                if (i + 1 >= argc) { error("missing path after '--win'"); return EXIT_ERROR; }
                win_path = argv[++i];
            }
            else if (a == "--qe-bands")
            {
                if (i + 1 >= argc) { error("missing path after '--qe-bands'"); return EXIT_ERROR; }
                qe_bands_path = argv[++i];
            }
            else if (a == "--verbose")             { log_level = "debug"; }
            else if (a == "--quiet")               { log_level = "warn"; }
            else if (a == "--log-level")
            {
                if (i + 1 >= argc) { error("missing level after '--log-level' (trace|debug|info|warn|error)"); return EXIT_ERROR; }
                log_level = argv[++i];
                if (log_level != "trace" && log_level != "debug" && log_level != "info" &&
                    log_level != "warn" && log_level != "error")
                {
                    error("unknown --log-level '" + log_level + "' (expected trace|debug|info|warn|error)");
                    return EXIT_ERROR;
                }
            }
            else if (a == "--log-file")
            {
                if (i + 1 >= argc) { error("missing path after '--log-file'"); return EXIT_ERROR; }
                log_file = argv[++i];
            }
            else if (a == "--no-log-file")         { no_log_file = true; }
            else if (a == "--bounds")              { emit_descriptor = true; }
            else if (a == "--exact-spin")          { exact_spin = true; }
            else if (a == "--covariant-velocity")  { velocity_mode = "covariant"; covariant_velocity = true; }
            else if (a == "--velocity-mode")
            {
                if (i + 1 >= argc) { error("'--velocity-mode' needs bare|berry_connection|covariant"); return EXIT_ERROR; }
                velocity_mode = argv[++i];
                if (velocity_mode != "bare" && velocity_mode != "berry_connection" && velocity_mode != "covariant")
                    { error("unknown --velocity-mode '" + velocity_mode + "' (use bare|berry_connection|covariant)"); return EXIT_ERROR; }
                covariant_velocity = (velocity_mode == "covariant");
            }
            else if (a == "--r-dat")
            {
                if (i + 1 >= argc) { error("'--r-dat' needs a path to the _r.dat position matrix"); return EXIT_ERROR; }
                r_dat_path = argv[++i];
            }
            else if (a == "--orbital-l" || a == "--orbital-L") { orbital_l = true; }
            else if (a == "--check")
            {
                // optional selector; default "all"
                static const char* known[] = {"all","hermiticity","sum_rules","algebra","aliasing","bounds"};
                check = "all";
                if (i + 1 < argc)
                {
                    const std::string nxt = argv[i+1];
                    for (const char* kw : known) if (nxt == kw) { check = nxt; ++i; break; }
                }
            }
            else if (a == "all")                   { want_all = true; }
            else if (a.size() > 1 && a[0] == '-' && !std::isdigit(static_cast<unsigned char>(a[1])))
                                                   { error("unknown option '" + a + "'"); return EXIT_ERROR; }
            else                                   { positional.push_back(a); }
        }

        // --create-template / --create / --provenance are scaffolding actions: they
        // write a `.w2s` input file and exit, so the positional LABEL N1 N2 N3
        // contract does not apply. main acts on them. (--provenance takes its SEED
        // inline, so no positional is required.)
        if (create_template || create_from_string) return PROCEED;
        if (make_provenance)
        {
            if (label.empty() && !positional.empty()) { label = positional[0]; positional.clear(); }
            if (label.empty()) { error("--provenance needs a SEED (e.g. --provenance graphene)"); return EXIT_ERROR; }
            return PROCEED;
        }

        // The `.w2s` interface is `-x FILE` (the official run flag). A bare `.w2s`
        // positional is no longer auto-run; catch it and point at -x rather than
        // failing with a confusing "expected LABEL N1 N2 N3".
        if (config_path.empty() && !positional.empty() && ends_with(positional[0], ".w2s"))
        {
            error("to run an input file, use '-x " + positional[0] + "'");
            return EXIT_ERROR;
        }

        // With --config the run.json supplies the label/operators/provenance, so
        // positional LABEL and the supercell dimensions are not required here.
        const bool config_mode = !config_path.empty();

        if (positional.empty() && !config_mode) { print_usage(); return EXIT_ERROR; }

        const bool bundle_mode = (mode == "bundle");

        // In bundle (or config) mode the operator is not expanded, so the supercell
        // dimensions are optional (and ignored). Accept LABEL [N1 N2 N3] [OP...]: if
        // the three tokens after LABEL all parse as integers, treat them as (ignored)
        // dims; otherwise the remaining positionals are operators.
        bool has_dims = !bundle_mode && !config_mode;
        if (bundle_mode || config_mode)
            has_dims = positional.size() >= 4 &&
                       looks_like_int(positional[1]) &&
                       looks_like_int(positional[2]) &&
                       looks_like_int(positional[3]);

        if (!bundle_mode && !config_mode && positional.size() < 4)
        {
            error("expected LABEL N1 N2 N3, but got " +
                  std::to_string(positional.size()) + " positional argument(s)");
            print_usage();
            return EXIT_ERROR;
        }

        if (!positional.empty()) label = positional[0];
        if (has_dims)
            for (int i = 0; i < 3; ++i)
            {
                const std::string& tok = positional[1 + i];
                try                           { cellDim[i] = std::stoi(tok); }
                catch (const std::exception&) { error("supercell dimension '" + tok + "' is not an integer"); return EXIT_ERROR; }
                if (cellDim[i] < 1)           { error("supercell dimension must be >= 1 (got '" + tok + "')"); return EXIT_ERROR; }
            }

        const size_t op_start = has_dims ? 4 : 1;
        if (want_all)
            operators = available_operators();
        else
            for (size_t i = op_start; i < positional.size(); ++i) operators.push_back(positional[i]);

        for (const auto& op : operators)
            if (!is_valid_operator(op))
            {
                error("unknown operator '" + op + "'. Run --list-operators for the valid names.");
                return EXIT_ERROR;
            }

        return PROCEED;
    }

private:
    void error(const std::string& msg) const
    {
        std::cerr << program_name << ": error: " << msg << "\n";
    }

    void print_usage() const
    {
        std::cerr << "usage: " << program_name << " -x INPUT.w2s\n"
                  << "       " << program_name << " --create-template [NAME]\n"
                  << "       " << program_name << " --create \"LABEL N1 N2 N3 [OP ...]\" [-inp NAME]\n"
                  << "       " << program_name << " --provenance SEED [--win FILE] [--qe-bands FILE]\n"
                  << "       " << program_name << " LABEL N1 N2 N3 [OP ... | all] [-o DIR]\n"
                  << "       " << program_name << " --help\n";
    }

    void print_version() const
    {
        std::cout << "wannier2sparse (LinQT utility) version 1.0\n";
    }

    void print_operators() const
    {
        std::cout << "Available operators:\n"
                  << "  velocity     : VX VY VZ\n"
                  << "  spin         : SX SY SZ\n"
                  << "  spin-current : VXSX VXSY VXSZ VYSX VYSY VYSZ VZSX VZSY VZSZ\n"
                  << "  (or 'all' to generate every operator)\n";
    }

    void print_help() const
    {
        std::cout <<
"wannier2sparse - expand a Wannier90 tight-binding model into a supercell\n"
"and export the Hamiltonian and operators as sparse (CSR) matrices.\n"
"\n"
"USAGE\n"
"  " << program_name << " -x INPUT.w2s                           run from an input file\n"
"  " << program_name << " --create-template [NAME]              scaffold an annotated NAME.w2s\n"
"  " << program_name << " --create \"LABEL N1 N2 N3 [OP...]\" [-inp NAME]   scaffold from a CLI line\n"
"  " << program_name << " --provenance SEED [--win FILE] [--qe-bands FILE]   extract provenance into SEED.w2s\n"
"  " << program_name << " LABEL N1 N2 N3 [OP ... | all] [options]   legacy positional run\n"
"\n"
"INPUT FILE (.w2s)\n"
"  The recommended interface is a JSON `.w2s` input file (// and /* */ comments\n"
"  allowed). It concentrates every option (label, mode, supercell, operators,\n"
"  checks, logging, and DFT/Wannier provenance) in one validated, traceable place.\n"
"  Create one with --create-template (blank, documented) or --create \"...\" (from a\n"
"  command line), edit it, then run it: `" << program_name << " -x input.w2s`.\n"
"  Every run also writes a `<LABEL>.out` JSON receipt: per-step timing, peak memory,\n"
"  warning/error tally, and the list of operators written with the input files that\n"
"  produced each (provenance link).\n"
"\n"
"ARGUMENTS\n"
"  LABEL        System label. The tool reads:\n"
"                 LABEL.uc      unit-cell lattice vectors\n"
"                 LABEL.xyz     orbital positions\n"
"                 LABEL_hr.dat  Wannier90 Hamiltonian\n"
"  N1 N2 N3     Supercell dimensions along each lattice vector (integers >= 1).\n"
"  OP ...       Operators to generate (e.g. VX SZ). 'all' generates every\n"
"               operator. If omitted, only LABEL.HAM.CSR is written.\n"
"\n"
"OPTIONS\n"
"  -o, --output-dir DIR   Directory for the .CSR files (default: current dir).\n"
"  -p, --project DIR      Directory holding the input files (default: current dir).\n"
"      --seed NAME        Seedname of the input files (default: LABEL). Inputs are\n"
"                         <DIR>/<NAME>_hr.dat, <DIR>/<NAME>.uc, <DIR>/<NAME>.xyz.\n"
"      --op-file NAME PATH Ingest an external operator in _hr.dat format from\n"
"                         PATH and write it as <LABEL>.NAME.CSR. Repeatable.\n"
"                         Enables hand-built tight-binding models and operators\n"
"                         produced outside this tool.\n"
"      --spin-current V S Write the derived spin current J = 1/2(V*S + S*V) for\n"
"                         velocity axis V and spin axis S (each X|Y|Z), formed by\n"
"                         sparse matrix product after expansion, as\n"
"                         <LABEL>.J<V>S<S>.CSR. Repeatable.\n"
"      --bounds           Also write a physical descriptor (.desc) next to each\n"
"                         CSR, including spectral bounds (a,b) for the Hamiltonian\n"
"                         (from <seed>.eig if present, else a Lanczos estimate).\n"
"      --exact-spin       Build the exact spin operators SXexact/SYexact/SZexact\n"
"                         from <seed>.spn + <seed>_u.mat (+ _u_dis.mat) via the\n"
"                         gauge transform (units hbar/2). Requires those files.\n"
"      --velocity-mode M  Velocity ladder for VX/VY/VZ and J=1/2{V,S}: 'bare'\n"
"                         (-i(R.lat)H), 'berry_connection' (default; adds Wannier\n"
"                         centres) or 'covariant' (-i(R.lat)H - i[H,A], full Berry\n"
"                         connection; needs <seed>_r.dat). Interband responses\n"
"                         (anomalous/spin Hall) want 'covariant'.\n"
"      --covariant-velocity  alias for --velocity-mode covariant.\n"
"      --r-dat PATH       Position matrix _r.dat (Wannier90 write_rmn) for the\n"
"                         covariant velocity. Default <seed>_r.dat.\n"
"      --orbital-L        Build orbital angular momentum LX/LY/LZ from <seed>.amn\n"
"                         + <seed>_u.mat + <seed>.win projections (units hbar).\n"
"                         Pure p/d shells only; errors on hybrids (e.g. sp3d2).\n"
"      --check [NAME]     Self-verify each operator (hermiticity, sum_rules,\n"
"                         algebra, aliasing, bounds) -> <op>.check sidecars.\n"
"                         NAME selects one check; default is all. CSR unchanged.\n"
"      --mode MODE        Output mode: 'sparse' (default; expand to supercell CSR)\n"
"                         or 'bundle' (emit the primitive operators O_ij(R) plus a\n"
"                         JSON manifest with provenance to <out>/<LABEL>.w2sp/, for\n"
"                         lsquant to build the Hamiltonian itself). In bundle mode\n"
"                         the supercell dimensions are optional and ignored.\n"
"  -x, --input-file PATH  Run from a `.w2s` input file: label, mode, supercell,\n"
"                         operators, output dir, checks and DFT/Wannier provenance\n"
"                         are all read from it. This is the recommended interface.\n"
"      --create-template [NAME]\n"
"                         Write a fully-annotated template to NAME.w2s (default\n"
"                         template.w2s) and exit. The scaffold is itself runnable.\n"
"      --create \"CMDLINE\" Parse a quoted command line (LABEL N1 N2 N3 [OP...] ...)\n"
"                         and serialize the equivalent <stem>.w2s, then exit.\n"
"  -inp, --inp NAME       Target stem for --create / --create-template (else the\n"
"                         LABEL from the command line, else 'template').\n"
"      --provenance SEED  Extract provenance from the model's side files and write\n"
"                         (or merge) it into SEED.w2s, then exit. Currently records\n"
"                         the band high-symmetry k-path: from SEED.win's kpoint_path\n"
"                         block (preferred), else a QE bands.in --qe-bands file. The\n"
"                         k-path travels with the model into the bundle manifest and\n"
"                         is read back by tools/hr_exactdiag.py bands.\n"
"      --win FILE         Wannier90 .win to read provenance from (default SEED.win).\n"
"      --qe-bands FILE    Quantum ESPRESSO bands.in (K_POINTS crystal_b) k-path\n"
"                         source, used as a fallback when no .win kpoint_path exists.\n"
"      --write            Do not run; serialize the equivalent input file to\n"
"                         <output>/<LABEL>.w2s from the positional arguments.\n"
"      --verbose          Console log level DEBUG (default INFO).\n"
"      --quiet            Console log level WARN (errors and warnings only).\n"
"      --log-level LVL    Set console level explicitly: trace|debug|info|warn|error.\n"
"      --log-file PATH    Write the full log here (default <output>/<LABEL>.run.log).\n"
"      --no-log-file      Do not write a log file (console only).\n"
"  -h, --help             Show this help and exit.\n"
"      --list-operators   List valid operator names and exit.\n"
"      --version          Show version and exit.\n"
"\n"
"NOTE\n"
"  LABEL.uc and LABEL.xyz are only required when a velocity/spin operator is\n"
"  requested; a Hamiltonian (or an --op-file operator) needs only its _hr.dat.\n"
"\n"
"OUTPUT\n"
"  LABEL.HAM.CSR and one LABEL.<OP>.CSR per requested operator, plus a LABEL.out\n"
"  JSON run receipt and (unless --no-log-file) a LABEL.run.log.\n"
"\n"
"EXAMPLES\n"
"  " << program_name << " --create-template graphene        # -> graphene.w2s (edit it)\n"
"  " << program_name << " -x graphene.w2s                    # run it\n"
"  " << program_name << " --create \"graphene 50 50 1 VX SZ\"  # -> graphene.w2s\n"
"  " << program_name << " graphene 50 50 1                    # legacy direct run\n"
"  " << program_name << " graphene 50 50 1 all -o out\n";
    }
};

#endif
