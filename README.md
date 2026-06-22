# wannier2sparse

wannier2sparse turns a Wannier90 tight-binding model into the sparse matrices that
large-scale quantum-transport codes need. It reads the real-space Hamiltonian
$H_{ij}(\mathbf{R})$ from a Wannier90 `_hr.dat`, replicates every hopping across an
$N_1\times N_2\times N_3$ supercell under periodic boundary conditions, and writes
the Hamiltonian and the operators you ask for — velocity, spin, spin current — as
sparse matrices in CSR format.

It is the bridge in a standard chain:

**[Quantum ESPRESSO](https://www.quantum-espresso.org/) → [Wannier90](https://wannier.org/) → wannier2sparse → KPM / Chebyshev transport ([LinQT](TODO-LINQT-URL), [lsquant](TODO-LSQUANT-URL))**

<!-- TODO[1]: replace TODO-LINQT-URL and TODO-LSQUANT-URL with the real links (see documentation_todo.md) -->


The output feeds kernel-polynomial evaluations of the density of states,
conductivity, and spin transport at sizes far beyond exact diagonalization. The
Fourier-transform sign, the `ndegen` normalization, the spin/orbital units, and the
Wigner-Seitz correction are documented in [`docs/conventions.md`](docs/conventions.md).

New here? Start with the [tutorials](examples/README.md).

---

## Requirements

- A C++11 compiler (GCC ≥ 4.8, Clang ≥ 3.4)
- CMake ≥ 3.14
- Eigen ≥ 3.3 (header-only)

| Platform | Install Eigen |
|----------|---------------|
| Debian / Ubuntu | `sudo apt install libeigen3-dev` |
| Fedora | `sudo dnf install eigen3-devel` |
| macOS (Homebrew) | `brew install eigen` |
| Conda | `conda install -c conda-forge eigen` |

## Building

```bash
mkdir build && cd build
cmake ..          # add -DEIGEN_DIR=/path/to/eigen3 if Eigen is not found
cmake --build . -j
```

This produces `build/wannier2sparse`. Run the tests with `ctest --output-on-failure`.

---

## Running

A calculation is described by a single input file, a `.w2s`, and run with:

```bash
wannier2sparse -x model.w2s          # long form: --input-file model.w2s
```

The `.w2s` is a JSON file (comments allowed). Its full set of keys is documented in
the [input-file reference](docs/input_file.md); a minimal run needs only:

```json
{
  "label": "graphene",
  "mode": "sparse",
  "output_dir": "out",
  "supercell": [80, 80, 1],
  "operators": ["VX", "VY", "VXSZ"]
}
```

- **`label`** — the model name. The tool reads `graphene_hr.dat`, plus
  `graphene.xyz` and `graphene.uc` when an operator needs orbital positions.
- **`mode`** — `sparse` writes the expanded supercell matrices; `bundle` ships the
  unexpanded $H_{ij}(\mathbf{R})$ instead, for a consumer to expand at its own size.
- **`supercell`** — the supercell size $[N_1, N_2, N_3]$.
- **`operators`** — which operators to build, by code (see [Operators](#operators)).
  The Hamiltonian is always written.

The run writes `out/graphene.HAM.CSR` and one `out/graphene.<OP>.CSR` per operator.
Spectral bounds, self-consistency checks, and a log of the run are produced
automatically — you do not configure them. To bring in an operator you built
yourself, drop its `_hr.dat` in an `operators/` folder next to the model; it is
picked up by name (details in the [input-file reference](docs/input_file.md)).

---

## Operators

Operators are requested by code in the `operators` list:

| Family | Codes | Units |
|--------|-------|-------|
| Hamiltonian | written always | eV |
| Velocity | `VX` `VY` `VZ` | eV·Å |
| Spin | `SX` `SY` `SZ` | ħ/2 |
| Spin current | `VXSX` … `VZSZ` ($J^{a}_{b}=\tfrac12\{v_b,S_a\}$) | eV·Å·ħ/2 |

`all` requests every operator. Spin operators need a model whose orbital labels mark
the spin channel (`_s+_` / `_s-_` in the `.xyz`). The exact gauge-transform spin
($S$ from `.spn` + `_u.mat`) and the orbital angular momentum $L$ are available too.
How each operator is built — geometric vs. gauge-matrix route, and which Wannier90
files it needs — is in [`docs/operators.md`](docs/operators.md).

---

## Provenance tracking

A tight-binding model is more useful when it remembers where it came from. The
crystal symmetry, the $k$-mesh, the pseudopotentials, and the Wannier centres are
what let you later exploit symmetries, reconstruct the band structure, or reproduce
the model — none of which can be recovered from the `_hr.dat` alone.

The principle is to record only what is needed to reproduce the model, nothing more.
Two surfaces carry it. A **bundle run** (`"mode": "bundle"` with
`provenance: { "win": ..., "qe_xml": ... }` in the `.w2s`) parses those side files
and writes the full record into the bundle `manifest.json`: from **Wannier90** the
*basis* (the projections that define the Wannier functions), the disentanglement
windows, `use_ws_distance`, and now the band high-symmetry $k$-path; from **Quantum
ESPRESSO** the structure, pseudopotentials, and $k$-mesh.

The **band $k$-path** can also be baked directly into the `.w2s` input so it travels
even without the side files at run time:

```bash
wannier2sparse --provenance graphene                 # reads graphene.win's kpoint_path -> graphene.w2s
wannier2sparse --provenance graphene --qe-bands qe/bands.in   # fallback: QE 'K_POINTS crystal_b'
```

`--provenance SEED` extracts the high-symmetry path — preferring the Wannier90
`.win` `kpoint_path` block, falling back to a QE `bands.in` — and writes (or merges)
it into `SEED.w2s` under `provenance.kpoint_path`, preserving any existing run keys.
The path then flows into the bundle `manifest.json` at run time, and
[`hr_exactdiag.py bands`](tools/hr_exactdiag.py) reads it straight from the `.w2s`,
so the Wannier bands are drawn on exactly the $k$-points the DFT used — with no extra
argument. (`-p` remains `--project`; the provenance command is the long-form
`--provenance`.)

<!-- TODO[3]: the broader --provenance (basis/structure/pseudopotentials/non-default
     settings written into the .w2s) is still planned; today --provenance records the
     band k-path, and bundle mode captures the rest into the manifest. -->


---

## Outputs

A `sparse` run writes one CSR file per operator (`label.HAM.CSR`,
`label.VX.CSR`, …). The CSR text format and the storage convention are in the
[input-file reference](docs/input_file.md).

A `bundle` run (`"mode": "bundle"`) instead writes the unexpanded operators
$O_{ij}(\mathbf{R})$ plus a JSON manifest to `label.w2sp/`, so a downstream code
holds the model itself and forms $H(\mathbf{k})=\sum_{\mathbf{R}}
e^{i\mathbf{k}\cdot\mathbf{R}}H(\mathbf{R})$ at any $\mathbf{k}$, choosing its own
supercell later.

---

## What you can compute

The [`examples/`](examples/README.md) directory is a set of short, self-contained
tutorials. Each builds a model, expands it, and reads a piece of physics off the
result — and each is checked against an analytic or exact reference.

| Tutorial | Physics |
|----------|---------|
| [1D chain](examples/01_chain1d/) | a band and its van Hove edges from a single hopping |
| [graphene](examples/02_graphene/) | the Dirac point, cross-checked against exact bands |
| [simple cubic](examples/03_cubic/) | linear scaling and the 3D resolution trade-off |
| [Haldane](examples/04_haldane/) | a complex hopping, a gap, and a quantized Hall plateau |
| [PdSe₂ spin Hall](examples/05_wannier_shc_pdse2/) | a real SOC material end to end, from DFT to $\sigma^{z}_{xy}$ |

---

## Conventions and references

The operator conventions are source-verified against Wannier90 3.1.0 and Quantum
ESPRESSO 7.2 and recorded in [`docs/conventions.md`](docs/conventions.md). For the
transport methodology the output feeds, see Z. Fan, J. H. Garcia, A. W. Cummings
*et al.*, *Linear scaling quantum transport methodologies*, Phys. Rep. **903**, 1
(2021), [arXiv:1811.07387](https://arxiv.org/abs/1811.07387); for Wannier functions,
N. Marzari *et al.*, Rev. Mod. Phys. **84**, 1419 (2012),
[arXiv:1112.5411](https://arxiv.org/abs/1112.5411).

## Earlier command-line interface

The positional command line of earlier versions still works and is documented in
[`legacy.md`](legacy.md).
