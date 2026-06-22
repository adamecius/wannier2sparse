# Documentation overhaul — items for you to evaluate

This file collects every open decision, placeholder, and follow-up from the docs
overhaul (branch `docs/readme-overhaul`, docs-first per our plan). Read it, mark
each item, and I'll act on your answers. Nothing here is merged.

Status legend: **[NEEDS YOU]** = I cannot resolve it · **[CONFIRM]** = decided, want
your OK · **[CODE]** = needs the follow-up code PR · **[DOING]** = doc work still to do.

---

## A. Placeholders left in the docs (search `TODO`)

- **[NEEDS YOU] TODO[1] — LinQT / lsquant URLs.** In `README.md` the toolchain line
  links `[LinQT](TODO-LINQT-URL)` and `[lsquant](TODO-LSQUANT-URL)`. Give me the two
  canonical URLs (repo or paper) and I'll fill them everywhere.

- **[NEEDS YOU] TODO[2] — provenance fields.** See section C below; confirm the list.

- **[CODE] TODO[3] — planned-CLI notice.** The README's `-x` run flag, the
  `--provenance`/`-p` section, and `VXSZ` as the anticommutator all describe the
  *target* CLI; they need the code PR before this README is correct on master.

---

## B. Decisions you already made — recorded here so you can re-check

- **[CONFIRM] Spin current = `VXSZ` operator code**, built as the expanded
  anticommutator ½{V,S} (the old separate `spin_currents` input key is dropped). The
  docs now present spin currents only as operator codes.
- **[CONFIRM] `emit_bounds` and `checks` become ON by default** (and `log_level`
  default stays). The minimal `.w2s` no longer shows them.
- **[CONFIRM] Docs-first**, code (`-x`, `--provenance`, schema, defaults) in a second
  PR once the CLI shape here is agreed.
- **[CONFIRM] Run flag is `-x` / `--run`** (replacing the bare `wannier2sparse file.w2s`).

---

## C. Provenance: what to pull from QE + Wannier90  — RESOLVED

**Principle (your answer): record only the minimum needed to reproduce the model.**

**From Wannier90:**
- the **basis** — the projections that define the Wannier functions
- the **on-disk paths of the Wannier matrices** on this machine
  (`_hr.dat`, `_u.mat`, `_r.dat`, `.spn`, and any others used)

**From Quantum ESPRESSO (enough to reproduce the DFT run):**
- the **final structure** (lattice + atomic positions)
- the **pseudopotentials**
- the **non-default input settings** (cutoffs, k-mesh, spin-orbit, smearing, …
  — i.e. anything that differs from QE's defaults)
- the **band k-path** (the `K_POINTS crystal_b` high-symmetry path used for the
  bands) — so the Wannier bands can be drawn on the same k-points as the DFT
  (`hr_exactdiag.py bands --kpath` already reads this from a QE `bands.in`).

This replaces my earlier enumerated list. The follow-up `--provenance` code extracts
exactly these. Implementation note: "non-default" means diffing the parsed input
against QE's defaults, or recording the controlled `&system`/`&control`/`K_POINTS`
keys that were explicitly set — I'll pick the simplest robust route when coding.

---

## D. CLI shape to implement in the follow-up code PR  (so the docs are true)

- **[CODE] `-x` / `--run`** as the run flag (and decide whether the bare
  `wannier2sparse file.w2s` keeps working or is removed).
- **[CODE] `-p` / `--provenance`** to read QE/Wannier90 output and emit a `.w2s`;
  **`-ws2 NAME`** chooses the output file. Note: `-p` currently means `--project`;
  proposal is to free `-p` for provenance and keep `--project` long-form only.
- **[CODE] `VXSZ`** (and the family) to build the expanded anticommutator ½{V,S}.
- **[CODE] Defaults**: `emit_bounds` and `checks` ON by default.
- **[CODE] `operators/` folder** as the default location for user-supplied operator
  `_hr.dat` files (replacing the explicit `op_files` list).

---

## E. Documentation work remaining on this branch

- **[DOING] Rename** `examples/example_9_SHC_in_PdSe2` → `examples/05_wannier_shc_pdse2`
  and fix every cross-reference. (README already links to the new path.)
- **[DOING] Rewrite the tutorial gallery** `examples/README.md` in the new voice.
- **[DOING] Rewrite the five tutorials** in the voice: physics-first, `-x`/`.w2s`
  commands, and **one arXiv-linked reference each**. Proposed references —
  **[NEEDS YOU]** confirm or replace:
  - 01_chain1d — a textbook tight-binding reference (no arXiv); suggest Marzari RMP
    [arXiv:1112.5411] for the Wannier framing, or you give one.
  - 02_graphene — Castro Neto *et al.*, RMP 81, 109 (2009), [arXiv:0709.1163].
  - 03_cubic — none obvious; KPM review Weiße *et al.*, RMP 78, 275 (2006),
    [arXiv:cond-mat/0504627].
  - 04_haldane — Haldane, PRL 61, 2015 (1988) (no arXiv); pair with the AHC review
    Nagaosa *et al.*, RMP 82, 1539 (2010), [arXiv:0904.4154].
  - 05_wannier_shc_pdse2 — the PdSe₂ / spin-Hall reference you have in mind?
    (and the covariant velocity: Wang, Yates, Souza, Vanderbilt, PRB 74, 195118
    (2006), [arXiv:cond-mat/0608257]).

---

## E2. PdSe2 figures to regenerate (combined SHC+DOS, DFT-path bands)

The PdSe2 model is DFT-derived and not committed, so I could not render its PNGs
here. The tooling is done and **verified on the committed Haldane model** (combined
σ+DOS figure renders correctly). On a machine with the DFT data, regenerate with
`bash get_inputs.sh 50` (the figure steps are appended to it):

- **[NEEDS YOU] `img/pdse2_shc_dos.png`** — the new combined figure (FIG. 3): σ^z_xy
  and the DOS on one shared energy axis, trivial gap (grey) + topological plateau
  (hatched) + DOS. The tutorial references it; it is not yet rendered in the repo.
- **[NEEDS YOU] `img/pdse2_bands.png`** — rerun on the DFT k-path
  (`--kpath qe/bands.in`); optionally overlay the QE bands for the DFT-vs-Wannier
  comparison the caption describes.
- Once `pdse2_shc_dos.png` exists, the old **`img/pdse2_shc.png`** (the zoom, now
  superseded) can be deleted; I left it in place to avoid a dangling reference.
- The exact trivial-gap window in the plot command (`--gap -0.05 0.05`) is a
  placeholder; set it to the real Fermi-level gap edges.

## F. Housekeeping

- **[NEEDS YOU]** `PR #25` (the earlier "tutorials → .w2s" PR) is **superseded** by
  this overhaul. OK to close it once this branch's PR is up?
- **[CONFIRM]** "Using wannier2sparse as a library" moved off the landing page into
  `legacy.md`; logging/run-summary prose also moved to `legacy.md`.
- The scratch files `documentation_plan.md` and `documentation_plan.mdç` are your
  notes — I have not committed them. Delete the stray `…ç` duplicate?
