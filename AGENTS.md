# AGENTS.md — agent operating contract (wannier2sparse)

This file defines how an agent **behaves**. The behavior is project-agnostic and
is shared with the sibling repositories (it originates in `adamecius/lsquant`);
the **project facts** for `wannier2sparse` live in the docs named in §1, never
here. Auto-read by Claude Code; `CLAUDE.md` is a symlink to this file, and other
agents (Kimi, etc.) should be pointed at it explicitly. **Read it fully before
acting.**

On conflict: **this file governs behavior, the project docs govern facts.**

---

## §0 — Session preamble (run first, every session)

You operate **in place** inside an existing checkout, do not clone elsewhere.
Verify the working tree is the right one and not stale, then read context (§1).

```bash
git rev-parse --is-inside-work-tree            # must print: true
git remote get-url origin                      # must be github.com/adamecius/wannier2sparse
git rev-parse --abbrev-ref HEAD                # branch
git rev-parse HEAD                             # HEAD — record it in the PR body
git fetch --quiet origin
git rev-list --left-right --count HEAD...@{u}  # right number = commits behind upstream
```

If the remote is not `adamecius/wannier2sparse`, or you are behind upstream,
**stop and report** before doing anything.

Build and test (the gate referenced throughout this file):

```bash
cmake -S . -B build && cmake --build build -j      # C++11, Eigen header-only, CMake >= 3.14
ctest --test-dir build --output-on-failure         # tests are GLOB-discovered from test/*.cpp
```

---

## §1 — Read project context before assuming anything

Do not reconstruct state from the diff. Read the project's context docs first;
they are the source of truth for all project facts. A missing doc → say so,
don't invent state. The living docs for this repo:

| Role | Doc(s) here | Holds |
|------|-------------|-------|
| Orient | [README.md](README.md) | build/run, CLI, input/output files, repo map |
| State / Plan | [roadmap_plans.md](roadmap_plans.md) | phases, invariants, key signatures (Spanish; facts are authoritative) |
| **Reference (conventions)** | [docs/conventions.md](docs/conventions.md) | **source-verified** FT sign, `ndegen`, `.spn`/`.amn` packing, spin/`L` units, WS convention |
| Reference (operators) | [docs/operators.md](docs/operators.md) | how each operator is built, geometric vs gauge (`U`-matrix) route, what QE/W90 files each needs |
| Findings | [docs/FINDINGS.md](docs/FINDINGS.md) | point-in-time discovery notes |
| Tests | [test/CMakeLists.txt](test/CMakeLists.txt) + `test/*.cpp` | the living gate; golden fixtures under `test/golden/`, `test/fixtures/` |
| Tutorials / plotting | [examples/README.md](examples/README.md), `agent_template/` | the worked tutorial gallery and the style guides |

**Anti-rot.** `docs/conventions.md` records the versions it was verified against
(Wannier90 3.1.0, QE 7.2) and the empirical anchor fixture (Fe `example17`). If a
convention claim is not traceable to a source line there, treat it as **possibly
stale** and reconcile before trusting it. The sacred invariant of this repo: the
expansion engine (`wrap_in_supercell` / `save_supercell_as_csr` in
`src/hopping_list.cpp`) is **byte-stable and never touched** without an explicit
task to change it; its output is a golden.

---

## §2 — Language, communication, and documentation

- Write **everything in English**: PR bodies, commit messages, code comments,
  docstrings, docs. (The maintainer may converse in Spanish; agent output is
  English regardless. The legacy `roadmap_plans.md` is Spanish and stays as-is.)
- **Templates (match them; never invent a structure when one exists):**
  - Tutorials / examples → `agent_template/tutorial_*`
  - In-code comments / PR text → `agent_template/documentation_*`
  - Figures → `agent_template/plotting_style.md`
  - If a referenced template is absent, flag it rather than guessing.
- **Documentation = in-code comments by default.** They are the source Doxygen
  (C++) can later turn into rendered docs. This repo already uses Javadoc-style
  `@brief`/`@param` blocks (see [include/hopping_list.hpp](include/hopping_list.hpp));
  match that style.
  - **Do not run a doc generator or emit generated/rendered docs unless
    explicitly requested.** When requested, read both `documentation_*` templates
    first and follow them.
  - When not requested, the in-code comments must be an **extensive,
    problem-oriented description**: the problem solved, the approach, assumptions,
    **units/conventions** (this project lives or dies on units: eV, eV·Å, ħ/2 for
    spin, ħ for orbital `L`), complexity — not a line-by-line restatement.
  - Always include, where applicable: **Pitfall**, **To-do**, **Warning**.

---

## §3 — Change discipline

- **Coverage before change.** Never refactor a path that lacks a regression /
  golden test; add coverage first. Tests are auto-discovered from `test/*.cpp`.
- **Green before commit.** Run the full `ctest` (§0) before every commit; state
  the pass count in the PR body. Keep exact-match goldens **byte-identical**
  unless the change is explicitly intended to move them, and then justify it.
  The CSR export path and the `_hr.dat` bundle writer are byte-deterministic by
  design (canonical sort, fixed precision); a diff there is a red flag.
- **No generated or bulk artifacts.** No `build/`, no `*.CSR`, no `*.o`, no large
  binaries or data blobs. If you find committed cruft, propose removing and
  gitignoring it. Committed fixtures under `test/golden/`, `test/fixtures/`, and
  the small `test/*_hr.dat` models are the deliberate exception.
- **Respect placeholders.** Intentional throwing stubs / not-implemented guards
  (e.g. the orbital-`L` hybrid-shell error path) stay as they are unless the task
  is to implement them. Don't silently revive or re-stub them.

---

## §4 — PR and merge-commit contract

The durable record lives in the **squash-merge commit body**, not only the GitHub
PR UI (`git log` is forever). Every PR body **and** the squash-merge message must
contain all four (mirror to `.github/pull_request_template.md`; see
`agent_template/documentation_pr.md`):

1. **Why these files changed** — per file or logical group, the *reason*, not the
   *what*.
2. **Tests** — which `ctest` targets ran, the pass count, and confirmation that
   exact-match goldens (CSR, bundle, `.ref`) stayed byte-identical. Name any new
   gate added.
3. **Reason for acceptance** — why this is safe to merge now (behavior-preserving?
   covered? gated?), tied to a project invariant (§1) or a roadmap phase.
4. **To-do / follow-ups** — what this PR leaves open and what it unblocks. Note
   any placeholder touched.

Record the HEAD hash from §0 and the relevant commit hashes. End commit messages
and PR bodies with the project's standard co-author / generated-by trailer.

---

## §5 — Reference material (mild RAG, never committed)

Reference books/papers live in `agent_references/` — **gitignored, local only.**
Only the manifest [agent_references.md](agent_references.md) is committed; it
carries title + topics + a "use-for" tag per entry, plus the arXiv id where one
exists, so you can grep it to pick what to open. Never `git add` a reference
file; never paste its content into the repo. When a change needs a literature or
convention justification (a Fourier sign, a spin unit, a disentanglement
detail), grep the manifest by topic, open only the matching local file or its
arXiv page, and cite it in the PR body by title.

The primary references for this project are the **Quantum ESPRESSO** and
**Wannier90** user guides and their method papers; see `agent_references.md`.

`.gitignore`:  `agent_references/`

---

## §6 — Definition of done

§0 verification recorded · full `ctest` green with goldens byte-identical · PR
body + squash-merge message satisfy all four points of §4 · no banned artifacts
committed (§3, §5) · output in English (§2) · the relevant project doc
(`docs/conventions.md`, `docs/operators.md`, or `roadmap_plans.md`) updated if the
change advances a phase or alters a convention · any figure satisfies §7.

---

## §7 — Figures and plots

Behavior here; the **specifications** (column widths, fonts, palette, encoding,
caption format) live in `agent_template/plotting_style.md`, a Reference-role doc.
On conflict: this section governs *when/how you act*, the style doc governs
*the numbers*.

- **Trigger.** Any figure you generate or edit, for any surface — article,
  GitHub Pages, slides, README, an `examples/` tutorial. No exceptions for
  "quick" or "throwaway" plots.
- **Read the style doc first.** Before producing a figure, read
  `agent_template/plotting_style.md` and use its presets. Missing doc → **say so,
  don't invent defaults.** Do not reconstruct style from existing figures.
- **Figure source is code, not raster.** Commit the generating script; rendered
  figures are regenerable cache. The `examples/` demos already do this
  (`gen_models.py`, `w2s_dos.py`, `validate.py`); extend that pattern.
- **No title; the caption carries everything.** Never set a figure/panel title.
  Every figure ships **with** a self-contained caption (begins `FIG. N.`, defines
  each colour/linestyle/marker, and lists the parameters needed to reproduce the
  panel: supercell size, broadening, KPM moment count, k-mesh). **A figure
  delivered without its caption is not done.**
- **Redundant encoding is mandatory.** Traces separable by colour **and**
  linestyle **and** marker; legible in grayscale and under colour-vision
  deficiency. Colour-only encoding, `jet`/rainbow, and sole red–green contrast
  are rejected.
- **Self-check gate.** Run the style doc's §8 checklist and report PASS/FAIL in
  the PR body next to the test gate. Any failed box → not shipped.
- **No post-export rescaling.** Export at the target width and include at natural
  size; rescaling on `\includegraphics` breaks the lettering-height rule.
</content>
