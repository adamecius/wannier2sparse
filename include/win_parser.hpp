/**
 * @file win_parser.hpp
 * @brief Parse a Wannier90 `.win` input file into Wannier provenance.
 *
 * The `.win` file records the conditions under which the Wannier functions were
 * built: the number of wannier/bands, the Monkhorst-Pack grid, the excluded
 * bands, the projection block, the disentanglement windows, and whether the
 * minimum-image (`use_ws_distance`) convention was on. parse_win_provenance()
 * scans those into a WannierProvenance for the bundle manifest, reusing the same
 * block-scan idiom as parse_projection_shells().
 */
#ifndef WIN_PARSER
#define WIN_PARSER

#include <string>
#include "system_provenance.hpp"

/**
 * @brief Parse a `.win` file into a WannierProvenance.
 * @param winfile path to the Wannier90 `.win` file
 * @return populated WannierProvenance with `present == true`
 * @throws std::runtime_error if the file cannot be opened
 *
 * Keys that are absent leave their `has_*` flag false (the manifest omits or
 * nulls them); the parser never fails on a missing key, only on a missing file.
 */
WannierProvenance parse_win_provenance(const std::string& winfile);

/**
 * @brief Parse the band high-symmetry k-path from a Quantum ESPRESSO `bands.in`.
 * @param bandsfile path to a QE input with a `K_POINTS crystal_b` block
 * @return the ordered list of high-symmetry nodes (label + fractional k)
 * @throws std::runtime_error if the file cannot be opened or has no crystal_b block
 *
 * The `crystal_b` form lists one node per line as `kx ky kz nseg ! LABEL` in
 * fractional (crystal) coordinates — the same path the DFT bands were computed on.
 * This is the fallback k-path source for `--provenance` when the Wannier90 `.win`
 * carries no `kpoint_path` block.
 */
std::vector<KpathNode> parse_qe_bands_kpath(const std::string& bandsfile);

#endif
