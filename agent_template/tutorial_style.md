# How to write a wannier2sparse tutorial

This guide explains how to turn `tutorial_template.md` into a finished tutorial,
and it is opinionated about voice on purpose. The worked reference is the
`examples/` gallery (`examples/README.md`): each demo builds a canonical
tight-binding model, expands it with `wannier2sparse`, and plots a density of
states from the CSR. Read one, then read this.

## The voice

Write as if you were explaining the physics to a smart friend who has never met
the code, in the spirit of Sagan, Feynman, or Tyson. That is not decoration. It
is a discipline with consequences for every paragraph.

Lead with the phenomenon, not the tool. The reader should feel a small puzzle
before they ever see a command. "Graphene has two atoms per cell, yet the bands
touch at a single point; what does that do to the density of states?" earns
attention. "This tutorial runs `wannier2sparse graphene 80 80 1`" does not.

Carry one idea. A tutorial is not a tour of features. It is the unfolding of a
single physical lesson, and every step exists to advance that one lesson. If a
step does not move the idea forward, cut it. A good spine for this project: the
real-space operator `O_ij(R)` is the model, and the supercell size is the
resolution dial on its spectrum.

Earn every command. No command appears without a physical reason for running it.
The reader should always know what question the next line answers. The code is
the answer to a physics question, never the subject.

Respect the reader's expertise, and only theirs. Assume a working physicist who
programs. Do not explain Bloch's theorem, a Fourier transform, or what a Wannier
function is at the level of a textbook. Do explain anything specific to this
pipeline: the `_hr.dat` gauge, the `ndegen` normalization, the spin-label
convention, the Wigner-Seitz minimum image. Define a term the first time it
earns its place.

Trust the phenomenon to be interesting. You do not need to tell the reader that
a result is "striking" more than once, and you never need filler. Compact prose
reads as confidence. Repetition reads as padding.

Name the trap. The most useful sentence in a tutorial is often the one that
stops a smart reader from drawing the obvious wrong conclusion. For this project
that sentence is often a convention: the velocity operator carries the Wannier
center displacement, not just the cell vector; the label-spin operator is zero
for a model whose `.xyz` labels do not mark `_s+_`/`_s-_`. Find that sentence in
your topic and give it room.

## Hard rules

These are not stylistic preferences. Apply them without exception.

- All mathematics is LaTeX: `$...$` inline, `$$...$$` for displayed equations.
  Never write equations as plain ASCII in prose or in code fences.
- Everything the reader runs goes in a code block. Outputs the reader will refer
  to (a file name, a CSR header, a manifest key) also go in a plain code block.
- No em dashes anywhere. Use a comma, a colon, parentheses, or a full stop.
- No setup paragraphs. Do not write about `PATH`, build directories, CMake,
  Eigen, or threading. Installation lives in the main README, and the tutorial
  points there once, in References.
- Every figure is embedded with a relative path, and its alt text states the
  takeaway, not the file name.
  `![Graphene DOS vanishes linearly at the Dirac point](fig_dos.png)`,
  not `![figure 2](fig_dos.png)`.
- Minimal formatting. No bold scattered through prose, no nested bullets.
  Headings and the occasional short list, nothing more.

## Structure

Follow `tutorial_template.md`. The shape is fixed because it works:

1. Title as a physical question in plain words.
2. A two paragraph hook: the puzzle, then what we are after and the one lesson.
3. `## The physics`: the minimum theory, the closed form or known limit to
   compare against, the one conceptual result stated in words, and an opening
   figure that shows the problem.
4. `## Step 1: build ...`: light, the system as physical objects (a `_hr.dat`
   Hamiltonian, its band edges, its orbital positions).
5. `## Step 2: ...`: the `wannier2sparse` call that produces the supercell CSR
   (or the bundle), framed by the operator it builds.
6. `## Step k: <a physical claim>`: later steps are titled with the claim they
   demonstrate, run a command, embed a figure, read the figure as physics.
7. `## What to take away`: three to five bullets, each a physics statement.
8. One sentence forward to the next tutorial.
9. `## References and links`: the fixed footer below.

Step titles are claims, not verbs. "The spectrum is intensive in the supercell
size" tells the reader what they are about to learn. "Plot the DOS" tells them
nothing.

## Figures

Figures are the argument, so make them carry it. Read `plotting_style.md` before
producing any figure and use its presets; the figure is not done without its
caption.

Ship a small script with each tutorial that produces every embedded PNG, either
by calling `wannier2sparse` and the example KPM plotter (`examples/w2s_dos.py`)
or, where the physics has a closed form, by computing it directly. A closed-form
overlay is faithful, not a mock-up; commit the PNGs next to the README so the
page renders offline.

Name figures for their content: `fig_bands.png`, `fig_dos.png`,
`fig_spin_texture.png`. Keep one idea per figure or per panel. A two panel figure
should make a single comparison the reader can state in one sentence. Put the
physical parameters the reader needs (the supercell size, the broadening, the
KPM moment count, the k-mesh) into the legend, so the figure is readable on its
own.

## Naming, so commands stay copy-pasteable

Keep labels and file patterns consistent across tutorials, because readers paste
commands and expect the output names to match. The examples use `graphene` /
`chain1d` / `cubic` / `haldane` for system labels, `LABEL.HAM.CSR` and
`LABEL.<OP>.CSR` for outputs, and `<LABEL>.w2sp/` for a bundle. When you add a
system, extend the convention rather than inventing a new one.

## References

Every tutorial ends with the same footer, nothing more by default:

```markdown
## References and links

- wannier2sparse source and documentation: https://github.com/adamecius/wannier2sparse
- Operator and gauge conventions: docs/conventions.md and docs/operators.md.
- Wannier functions: N. Marzari et al., Rev. Mod. Phys. 84, 1419 (2012),
  arXiv:1112.5411. Wannier90: G. Pizzi et al., J. Phys. Condens. Matter 32,
  165902 (2020), arXiv:1907.09788.
- Transport methodology: Z. Fan, J. H. Garcia, A. W. Cummings et al., Linear
  scaling quantum transport methodologies, Phys. Rep. 903, 1 (2021),
  arXiv:1811.07387.
- Installation: see the main README of the repository.
```

If, and only if, the tutorial leans on a specific method or result, add a short
`## Further reading` with that one citation, grepped from `agent_references.md`.
The disentanglement step belongs to Souza, Marzari, Vanderbilt, Phys. Rev. B 65,
035109 (2001), arXiv:cond-mat/0108084; the orbital-magnetization Wannier route to
Lopez et al., Phys. Rev. B 85, 014435 (2012). Do not build a bibliography. Two or
three links is the ceiling.

## Before you commit

A finished tutorial passes all of these:

- A reader could state the one lesson in a single sentence after reading.
- Every command has a physical reason stated before it.
- Every equation is LaTeX, every runnable line is in a code block.
- Every figure is embedded, with alt text that states its takeaway, the PNG is
  committed, and the figure satisfies `plotting_style.md` (caption included).
- There are no em dashes and no setup, PATH, CMake, or build paragraphs.
- The references footer is the standard block above.
- Nothing repeats, and nothing is there that the lesson does not need.
</content>
