/**
 * @file qe_xml_parser.hpp
 * @brief Parse a Quantum ESPRESSO `data-file-schema.xml` into DFT provenance.
 *
 * Pulls the structure (lattice + atoms, Bohr->Angstrom), crystal symmetry
 * operations, and DFT conditions (code/version, XC functional, spin-orbit /
 * noncollinear flags, wavefunction cutoff, k-mesh, pseudopotentials) out of the
 * QE XML and into a DftProvenance for the bundle manifest. The actual XML is
 * read with the vendored pugixml (QE's namespaced, version-varying XML is too
 * brittle to hand-parse). Handles the QE 6.x and 7.x flavors; any node that is
 * absent simply leaves its field unset rather than failing.
 */
#ifndef QE_XML_PARSER
#define QE_XML_PARSER

#include <string>
#include "system_provenance.hpp"

/**
 * @brief Parse a QE `data-file-schema.xml` into a DftProvenance.
 * @param xml_path path to the QE data-file-schema.xml
 * @return populated DftProvenance with `present == true`
 * @throws std::runtime_error if the file cannot be opened or is not valid XML
 *
 * Lattice vectors and atomic positions are converted from Bohr to Angstrom.
 * Fractional atomic coordinates are computed from the parsed lattice. Energies
 * (ecutwfc) are reported in Rydberg. Missing optional nodes leave their
 * `has_*` flags false; only a missing/invalid file throws.
 */
DftProvenance parse_qe_xml(const std::string& xml_path);

#endif
