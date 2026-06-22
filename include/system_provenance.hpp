/**
 * @file system_provenance.hpp
 * @brief Plain data structures describing the physical/DFT provenance of a model.
 *
 * These structs aggregate everything the bundle manifest records besides the
 * operators themselves: the crystal structure, the wannier sites, and the DFT
 * (Quantum ESPRESSO) and Wannier90 conditions that produced the model. They hold
 * no logic and no Eigen types; the parsers (qe_xml_parser, win_parser) and the
 * tbmodel populate them, and bundle_writer serializes them to JSON.
 *
 * Each provenance block carries a `present` flag: when its source file is absent
 * the block is left default-constructed and the manifest emits `null` for it, so
 * a bundle is always well-formed even with partial provenance.
 */
#ifndef SYSTEM_PROVENANCE
#define SYSTEM_PROVENANCE

#include <string>
#include <vector>
#include <array>

/**
 * @brief One atom from the DFT structure (Cartesian + fractional, Angstrom).
 */
struct AtomSite
{
    std::string          species;
    std::array<double,3> frac;   ///< fractional coordinates
    std::array<double,3> cart;   ///< Cartesian coordinates (Angstrom)
    AtomSite() : frac({{0,0,0}}), cart({{0,0,0}}) {}
};

/**
 * @brief One crystal symmetry operation (rotation + fractional translation).
 */
struct SymmetryOp
{
    std::array<std::array<double,3>,3> rotation;     ///< rotation matrix (crystal axes)
    std::array<double,3>               translation;  ///< fractional translation
    bool                               time_reversal;
    SymmetryOp() : rotation({{{{0,0,0}},{{0,0,0}},{{0,0,0}}}}),
                   translation({{0,0,0}}), time_reversal(false) {}
};

/**
 * @brief One pseudopotential entry.
 */
struct PseudoInfo
{
    std::string species;
    std::string file;
    double      z_valence;
    bool        has_z_valence;
    PseudoInfo() : z_valence(0.0), has_z_valence(false) {}
};

/**
 * @brief DFT (Quantum ESPRESSO) conditions, parsed from data-file-schema.xml.
 */
struct DftProvenance
{
    bool        present;            ///< false => manifest emits null
    std::string code;               ///< "Quantum ESPRESSO"
    std::string version;
    std::string source_file;        ///< path to the parsed XML
    std::string xc_functional;
    bool        spin_orbit;
    bool        noncolin;
    bool        has_ecutwfc;
    double      ecutwfc_Ry;
    bool        has_k_mesh;
    std::array<int,3> k_mesh;
    std::array<int,3> k_mesh_shift;
    std::vector<PseudoInfo> pseudopotentials;

    // DFT-level structure (supplements the .uc/.xyz structure).
    bool                               has_structure;
    std::array<std::array<double,3>,3> lattice;     ///< Angstrom, rows are lattice vectors
    std::vector<AtomSite>              atoms;
    std::vector<SymmetryOp>            symmetry;

    DftProvenance()
        : present(false), code("Quantum ESPRESSO"), spin_orbit(false), noncolin(false),
          has_ecutwfc(false), ecutwfc_Ry(0.0), has_k_mesh(false),
          k_mesh({{0,0,0}}), k_mesh_shift({{0,0,0}}), has_structure(false),
          lattice({{{{0,0,0}},{{0,0,0}},{{0,0,0}}}}) {}
};

/**
 * @brief One Wannier90 projection (e.g. "C:p").
 */
struct Projection
{
    std::string site;
    std::string orbital;
    std::string raw;     ///< the raw projection line
};

/**
 * @brief One high-symmetry node on the band k-path (label + fractional k).
 *
 * The band path is recorded as an ordered list of these nodes (e.g. G-M-K-G),
 * extracted from the Wannier90 `.win` `kpoint_path` block or a Quantum ESPRESSO
 * `K_POINTS crystal_b` block, so the Wannier bands can be drawn on the same
 * high-symmetry path the DFT used. Coordinates are fractional (crystal) and the
 * label is the symmetry-point name as written in the source (G/Gamma kept verbatim).
 */
struct KpathNode
{
    std::string          label;   ///< symmetry-point label (e.g. "G", "M", "K")
    std::array<double,3> k;       ///< fractional (crystal) coordinates
    KpathNode() : k({{0,0,0}}) {}
};

/**
 * @brief Wannier90 conditions, parsed from the .win file.
 */
struct WannierProvenance
{
    bool             present;        ///< false => manifest emits null
    std::string      source_file;
    bool             has_num_wann;
    int              num_wann;
    bool             has_num_bands;
    int              num_bands;
    std::vector<int> exclude_bands;
    bool             has_mp_grid;
    std::array<int,3> mp_grid;
    std::vector<Projection> projections;
    bool             dis_enabled;
    bool             has_dis_win;
    double           dis_win_min, dis_win_max;
    bool             has_dis_froz;
    double           dis_froz_min, dis_froz_max;
    bool             use_ws_distance;
    std::vector<KpathNode> kpoint_path;       ///< band high-symmetry path (G-M-K-...), if present
    std::string      kpoint_path_source;      ///< where the path came from ("win" or "qe_bands")

    WannierProvenance()
        : present(false), has_num_wann(false), num_wann(0), has_num_bands(false),
          num_bands(0), has_mp_grid(false), mp_grid({{0,0,0}}), dis_enabled(false),
          has_dis_win(false), dis_win_min(0), dis_win_max(0), has_dis_froz(false),
          dis_froz_min(0), dis_froz_max(0), use_ws_distance(false) {}
};

/**
 * @brief One manually-declared orbital (used when no .xyz / DFT source is given).
 *
 * The user supplies these in the input file's `provenance.manual.orbitals` block so
 * a model that arrives as a bare `_hr.dat` (no Wannier90 / QE side-files) can still
 * record which basis function each index is and where it sits. Positions are in
 * Angstrom, matching the .xyz / lattice convention used everywhere else.
 */
struct ManualOrbital
{
    int                  index;    ///< zero-based orbital index
    std::string          label;    ///< e.g. "Fe-d_xy", "pz"
    std::string          element;  ///< chemical species, e.g. "Fe"
    std::array<double,3> xyz;      ///< Cartesian position (Angstrom)
    ManualOrbital() : index(0), xyz({{0,0,0}}) {}
};

/**
 * @brief User-declared provenance, supplied verbatim in the input file.
 *
 * Whereas DftProvenance / WannierProvenance are *parsed* from QE/W90 side-files,
 * this block is what the user types into the `.w2s` input under
 * `provenance.manual`. It exists for the common case of a model that is only a
 * `_hr.dat` with no machine-readable DFT/Wannier sources: the user can still pin
 * the pseudopotentials, k-point grid, basis set, exchange-correlation functional,
 * plane-wave cutoff and the per-orbital identity/position. Every field is optional;
 * the block is echoed unchanged into the bundle manifest and the `.out` receipt so
 * the record of *how this model was made* travels with the operators. Nothing here
 * affects the numerics — it is documentation that is validated for shape only.
 */
struct ManualProvenance
{
    bool                     present;          ///< false => omitted from manifest/.out
    std::string              code;             ///< DFT/Wannier code + version, free text
    std::vector<PseudoInfo>  pseudopotentials; ///< species -> file (z_valence optional)
    bool                     has_kpoint_grid;
    std::array<int,3>        kpoint_grid;      ///< SCF k-point mesh
    std::string              basis;            ///< basis-set description (free text)
    std::string              xc_functional;    ///< e.g. "PBE"
    bool                     has_ecutwfc;
    double                   ecutwfc_Ry;       ///< plane-wave cutoff (Ry)
    std::vector<ManualOrbital> orbitals;       ///< per-index basis identity + position
    std::string              notes;            ///< any extra free-text provenance

    ManualProvenance()
        : present(false), has_kpoint_grid(false), kpoint_grid({{0,0,0}}),
          has_ecutwfc(false), ecutwfc_Ry(0.0) {}
};

/**
 * @brief One wannier site (orbital), from the .xyz file.
 */
struct WannierSite
{
    int                  index;   ///< zero-based orbital index
    std::string          label;
    std::array<double,3> cart;    ///< Cartesian position (Angstrom)
    int                  spin;    ///< +1 (_s+_), -1 (_s-_), 0 (unmarked)
    WannierSite() : index(0), cart({{0,0,0}}), spin(0) {}
};

/**
 * @brief Aggregate system + provenance for the bundle manifest.
 */
struct SystemProvenance
{
    bool                               has_lattice;   ///< lattice from .uc loaded
    std::array<std::array<double,3>,3> lattice;       ///< Angstrom, rows are lattice vectors
    std::vector<WannierSite>           wannier_sites; ///< from .xyz
    int                                num_wann;

    DftProvenance     dft;
    WannierProvenance wann;
    ManualProvenance  manual;        ///< user-declared provenance from the input file

    SystemProvenance()
        : has_lattice(false),
          lattice({{{{0,0,0}},{{0,0,0}},{{0,0,0}}}}), num_wann(0) {}
};

#endif
