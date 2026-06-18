/**
 * @file bundle_writer.hpp
 * @brief Write an lsquant-readable provenance bundle.
 *
 * A bundle is a directory `<label>.w2sp/` holding a JSON manifest plus one
 * `_hr.dat`-format data file per operator. Unlike the CSR export path, the bundle
 * ships the PRIMITIVE real-space operators O_ij(R) together with the crystal
 * structure, symmetry, and DFT/Wannier provenance, so a consumer (lsquant) can
 * build the Hamiltonian / supercell itself. The supercell expansion engine is not
 * involved here by design.
 *
 *   <label>.w2sp/
 *     manifest.json
 *     operators/HAM.hr.dat, VX.hr.dat, ...
 *     wsvec.dat            (copied verbatim if present)
 */
#ifndef BUNDLE_WRITER
#define BUNDLE_WRITER

#include <string>
#include <vector>
#include "hopping_list.hpp"
#include "descriptor.hpp"
#include "system_provenance.hpp"

/**
 * @brief One operator to write into the bundle, with its physical descriptor.
 */
struct BundleOperator
{
    std::string        name;   ///< file/operator name, e.g. "HAM", "VX", "SZ"
    OperatorDescriptor desc;   ///< observable/component/units/provenance (+ optional bounds)
    hopping_list       hl;     ///< primitive-cell operator O_ij(R)
};

/**
 * @brief Top-level bundle settings recorded in the manifest.
 */
struct BundleSpec
{
    std::string label;
    double      truncation_threshold;  ///< echoed into the manifest
    bool        ndegen_applied;        ///< values are post-ndegen (always true here)
    bool        wsvec_applied;         ///< whether the wsvec correction was folded in
    std::string wsvec_src;             ///< path to _wsvec.dat to copy, or empty

    BundleSpec() : truncation_threshold(0.0), ndegen_applied(true), wsvec_applied(false) {}
};

/**
 * @brief Write the bundle directory, its operator data files, and manifest.json.
 * @param spec  top-level bundle settings
 * @param prov  structure + DFT/Wannier provenance
 * @param ops   operators to emit (HAM first by convention)
 * @param out_dir parent directory in which `<label>.w2sp/` is created
 * @return path to the created bundle directory
 */
std::string write_bundle(const BundleSpec& spec, const SystemProvenance& prov,
                         const std::vector<BundleOperator>& ops, const std::string& out_dir);

#endif
