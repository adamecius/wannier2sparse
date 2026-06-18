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

#endif
