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
    std::string              config_path;     ///< --config: path to a run.json (bundle driver)
    bool                     has_truncation;  ///< truncation_threshold was set (from run.json)
    double                   truncation_threshold; ///< bundle truncation threshold (manifest)
    std::string              qe_xml_path;     ///< QE data-file-schema.xml for DFT provenance
    std::string              win_path;        ///< .win for Wannier provenance
    std::string              program_name;    ///< argv[0]
    // Input-file CLI (--create / --write / --run): see include/input_file.hpp.
    bool                     do_create;       ///< --create: write a template input file
    bool                     do_write;        ///< --write: set one option in the input file
    std::vector<std::string> write_tokens;    ///< tokens after --write (option [value...])
    std::string              run_path;        ///< --run <input-file>
    std::string              inp_path;        ///< -inp/--inp <input-file> (create/write target)
    bool                     covariant_velocity; ///< covariant-velocity option (recorded; run errors)

    /**
     * @brief Default constructor.
     */
    W2SP_arguments()
        : cellDim({{1, 1, 1}}), output_dir("."), emit_descriptor(false),
          exact_spin(false), orbital_l(false), mode("sparse"),
          has_truncation(false), truncation_threshold(0.0),
          program_name("wannier2sparse"),
          do_create(false), do_write(false), covariant_velocity(false) {}

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
            else if (a == "--config")
            {
                if (i + 1 >= argc) { error("missing path after '--config'"); return EXIT_ERROR; }
                config_path = argv[++i];
            }
            else if (a == "--create")   // optional filename follows (else use -inp / default)
            {
                do_create = true;
                if (i + 1 < argc && argv[i + 1][0] != '-') inp_path = argv[++i];
            }
            else if (a == "--run")
            {
                if (i + 1 >= argc) { error("missing input file after '--run'"); return EXIT_ERROR; }
                run_path = argv[++i];
            }
            else if (a == "-inp" || a == "--inp")
            {
                if (i + 1 >= argc) { error("missing file after '-inp'"); return EXIT_ERROR; }
                inp_path = argv[++i];
            }
            else if (a == "--write")    // collect the option + its value tokens (until next flag)
            {
                do_write = true;
                while (i + 1 < argc &&
                       !(argv[i + 1][0] == '-' && !std::isdigit(static_cast<unsigned char>(argv[i + 1][1]))))
                    write_tokens.push_back(argv[++i]);
            }
            else if (a == "--bounds")              { emit_descriptor = true; }
            else if (a == "--exact-spin")          { exact_spin = true; }
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

        // Input-file CLI: --create/--write are pure file edits (no run), and --run
        // takes all its options from the input file (applied in main). None of them
        // require positional LABEL / supercell dims, so short-circuit here.
        if (do_create || do_write || !run_path.empty()) return PROCEED;

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
        std::cerr << "usage: " << program_name
                  << " LABEL N1 N2 N3 [OP ... | all] [-o DIR]\n"
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
"  " << program_name << " LABEL N1 N2 N3 [OP ... | all] [options]\n"
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
"      --config PATH      Drive a bundle run from a run.json config file instead\n"
"                         of positional arguments (label, operators, output dir\n"
"                         and DFT/Wannier provenance sources are read from it).\n"
"\n"
"INPUT-FILE WORKFLOW (the old options recorded in an editable `key = value` file)\n"
"      --create [FILE]    Write a commented template input file (default w2sp.inp,\n"
"                         or the name given, or the -inp FILE).\n"
"      --write OPT [VAL]  Set/append one option in the input file (-inp FILE).\n"
"                         OPT may be a flag (covariant-velocity, exact-spin,\n"
"                         orbital-L, bounds), an operator name (VX, SZ, ...),\n"
"                         'spin-current X Z', 'op-file NAME PATH', or key=value\n"
"                         (e.g. supercell 50 50 1, label=foo, mode=bundle).\n"
"      -inp FILE          Input file targeted by --create / --write.\n"
"      --run FILE         Run the input file; writes <output_dir>/<label>.w2sp.out\n"
"                         (a summary of the resolved options and produced files).\n"
"                         Produces output byte-identical to the positional CLI.\n"
"  -h, --help             Show this help and exit.\n"
"      --list-operators   List valid operator names and exit.\n"
"      --version          Show version and exit.\n"
"\n"
"NOTE\n"
"  LABEL.uc and LABEL.xyz are only required when a velocity/spin operator is\n"
"  requested; a Hamiltonian (or an --op-file operator) needs only its _hr.dat.\n"
"\n"
"OUTPUT\n"
"  LABEL.HAM.CSR and one LABEL.<OP>.CSR per requested operator.\n"
"\n"
"EXAMPLES\n"
"  " << program_name << " graphene 50 50 1\n"
"  " << program_name << " graphene 50 50 1 VX VY SZ\n"
"  " << program_name << " graphene 50 50 1 all -o out\n";
    }
};

#endif
