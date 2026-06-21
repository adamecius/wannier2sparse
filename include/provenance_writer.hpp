/**
 * @file provenance_writer.hpp
 * @brief Shared JSON emitter for the user-declared (manual) provenance block.
 *
 * The `provenance.manual` block of a `.w2s` input is echoed verbatim into three
 * surfaces: the input file itself (write_run_config / --create), the bundle
 * manifest (bundle_writer), and the `.out` run receipt (run_report). To keep those
 * three byte-identical and avoid drift, the serialization lives here once as a
 * small inline helper over the project's JsonWriter. It writes nothing when the
 * block is absent (present == false); the caller decides whether to emit a `null`
 * placeholder instead.
 *
 * This is documentation data only: nothing here influences the numerics. Units
 * follow the rest of the project (Angstrom for positions, Ry for the cutoff as in
 * the QE block). Keys are emitted in a fixed order so the output is deterministic.
 */
#ifndef PROVENANCE_WRITER
#define PROVENANCE_WRITER

#include "json_writer.hpp"
#include "system_provenance.hpp"

/**
 * @brief Emit a ManualProvenance as a JSON object value (no leading key).
 * @param w  the JSON writer (a value slot must be open, e.g. after w.key("manual"))
 * @param mp the manual provenance block (present == true expected)
 *
 * Only set fields are written, in a fixed key order, so absent data does not clutter
 * the record and the layout stays byte-stable for goldens.
 */
inline void emit_manual_provenance(JsonWriter& w, const ManualProvenance& mp)
{
    w.begin_object();
        if (!mp.code.empty())          w.member("code", mp.code);
        if (!mp.basis.empty())         w.member("basis", mp.basis);
        if (!mp.xc_functional.empty()) w.member("xc_functional", mp.xc_functional);
        if (mp.has_ecutwfc)            w.member("ecutwfc_Ry", mp.ecutwfc_Ry);
        if (mp.has_kpoint_grid)
        {
            w.key("kpoint_grid");
            w.begin_array(); for (int x : mp.kpoint_grid) w.inum(x); w.end_array();
        }
        if (!mp.pseudopotentials.empty())
        {
            w.key("pseudopotentials");
            w.begin_array();
            for (const PseudoInfo& p : mp.pseudopotentials)
            {
                w.begin_object();
                    w.member("species", p.species);
                    w.member("file", p.file);
                    if (p.has_z_valence) w.member("z_valence", p.z_valence);
                w.end_object();
            }
            w.end_array();
        }
        if (!mp.orbitals.empty())
        {
            w.key("orbitals");
            w.begin_array();
            for (const ManualOrbital& o : mp.orbitals)
            {
                w.begin_object();
                    w.member("index", (long long)o.index);
                    if (!o.label.empty())   w.member("label", o.label);
                    if (!o.element.empty()) w.member("element", o.element);
                    w.key("xyz"); w.array_d(o.xyz);
                w.end_object();
            }
            w.end_array();
        }
        if (!mp.notes.empty()) w.member("notes", mp.notes);
    w.end_object();
}

#endif
