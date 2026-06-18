# wannier2sparse

A command-line tool (part of the **LinQT** package) that expands a Wannier90
`tight-binding model into a supercell and exports the Hamiltonian and a set of
operators as sparse matrices in CSR format, ready for KPM / Chebyshev transport
calculations.

The connectivity is taken directly from the Wannier90 `_hr.dat` file (it is not
searched geometrically): each non-zero `H_ij(R)` is replicated across the
supercell and wrapped under periodic boundary conditions. Operator conventions
— Fourier-transform sign, `ndegen` normalization, spin units, and the
Wigner-Seitz minimum-image correction — are documented in
[`docs/conventions.md`](docs/conventions.md).

---

## Requirements

- A C++11 compiler (GCC ≥ 4.8, Clang ≥ 3.4)
- CMake ≥ 3.14
- Eigen ≥ 3.3 (header-only)

Install Eigen:

| Platform | Command |
|----------|---------|
| Debian / Ubuntu | `sudo apt install libeigen3-dev` |
| Fedora | `sudo dnf install eigen3-devel` |
| macOS (Homebrew) | `brew install eigen` |
| Conda | `conda install -c conda-forge eigen` |

---

## Building

```bash
mkdir build && cd build
cmake ..
cmake --build . -j
```

The executable is produced at `build/wannier2sparse`.

If Eigen is not found automatically, point CMake at it:

```bash
cmake .. -DEIGEN_DIR=/path/to/eigen3      # directory containing the Eigen/ folder
```

Other options:

| Option | Default | Effect |
|--------|---------|--------|
| `-DW2SP_BUILD_TESTS=OFF` | `ON` | Skip building the unit tests |
| `-DCMAKE_INSTALL_PREFIX=<path>` | system default | Install location |

Run the tests:

```bash
ctest --output-on-failure
```

Install:

```bash
cmake --install . --prefix /your/prefix
```

---

## Using wannier2sparse as a library

Installing exports a CMake config-package, so a downstream project can consume
the library with `find_package`:

```cmake
find_package(wannier2sparse 1.0 CONFIG REQUIRED)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE wannier2sparse::wannierlib)
```

Point CMake at the install prefix with `-DCMAKE_PREFIX_PATH=/your/prefix`. The
imported target carries its include directory and, when Eigen was found through
its own `Eigen3Config.cmake`, propagates `Eigen3::Eigen` transitively. If Eigen
was instead located through the manual `-DEIGEN_DIR=...` fallback, the installed
package cannot re-find Eigen and the consumer must provide it.

---

## Usage

```
wannier2sparse LABEL N1 N2 N3 [OP ... | all] [options]
```

| Argument | Description |
|----------|-------------|
| `LABEL`     | System label. |
| `N1 N2 N3`  | Supercell dimensions along each lattice vector (integers ≥ 1). |
| `OP ...`    | Operators to generate, e.g. `VX SZ`. `all` generates every operator. If omitted, only the Hamiltonian is written. |

| Option | Description |
|--------|-------------|
| `-o, --output-dir DIR` | Directory for the `.CSR` output (default: current dir). |
| `-p, --project DIR`    | Directory holding the input files (default: current dir). |
| `--seed NAME`          | Seedname of the input files (default: `LABEL`). |
| `--op-file NAME PATH`  | Ingest an external operator in `_hr.dat` format from `PATH` and write it as `<LABEL>.NAME.CSR`. Repeatable. |
| `--spin-current V S`   | Write the derived spin current `J = 1/2 {V_V, S_S}` as `<LABEL>.JVSS.CSR` (`V,S ∈ {X,Y,Z}`). Repeatable. |
| `--bounds`             | Write a physical descriptor (`.desc`) next to each CSR, including spectral bounds `[a,b]` for the Hamiltonian. |
| `--exact-spin`         | Build exact spin operators from `<seed>.spn` + `<seed>_u.mat` via the gauge transform. |
| `--orbital-L`          | Build orbital angular momentum from `<seed>.amn` + `<seed>_u.mat` + `<seed>.win` (pure p/d shells only). |
| `--check [NAME]`       | Self-verify operators and write `<op>.check` sidecars (`NAME` = `hermiticity`, `sum_rules`, `algebra`, `aliasing`, `bounds`; default `all`). |
| `--mode MODE`          | Output mode: `sparse` (default; expand to the supercell CSR) or `bundle` (emit the primitive operators `O_ij(R)` plus a JSON manifest with provenance to `<out>/<LABEL>.w2sp/`). In `bundle` mode `N1 N2 N3` are optional and ignored. See [Bundle output](#bundle-output-provenance-for-lsquant). |
| `--config PATH`        | Drive a `bundle` run from a `run.json` config file (label, operators, output dir, and DFT/Wannier provenance sources) instead of positional arguments. |
| `-h, --help`           | Show help and exit. |
| `--list-operators`     | List valid operator names and exit. |
| `--version`            | Show version and exit. |

Invalid input (unknown operator, non-integer or non-positive dimension, missing
input file, too few arguments) prints a clear message and returns a non-zero
exit code — the tool never aborts with an assertion.

`LABEL.uc` and `LABEL.xyz` are only required when a velocity or spin operator is
requested; a plain Hamiltonian (or an `--op-file` operator) needs only its
`_hr.dat`.

### Examples

```bash
# Hamiltonian only, 50x50x1 supercell
wannier2sparse graphene 50 50 1

# Hamiltonian + selected operators
wannier2sparse graphene 50 50 1 VX VY SZ

# Every operator, written to ./out
wannier2sparse graphene 50 50 1 all -o out

# Exact spin operators from Wannier90 gauge data
wannier2sparse Fe 4 4 4 --exact-spin

# Bundle: primitive operators O_ij(R) + a provenance manifest, no expansion
# (the supercell dimensions are accepted but ignored in bundle mode)
wannier2sparse graphene 1 1 1 VX SZ --mode bundle -o out
# -> out/graphene.w2sp/{manifest.json, operators/HAM.hr.dat, operators/VX.hr.dat, ...}

# The same bundle run, driven entirely by a config file
wannier2sparse --config run.json
```

---

## Operators

| Family | Codes |
|--------|-------|
| Velocity (current) | `VX` `VY` `VZ` |
| Spin density | `SX` `SY` `SZ` |
| Spin current | `VXSX` `VXSY` `VXSZ` `VYSX` `VYSY` `VYSZ` `VZSX` `VZSY` `VZSZ` |

`all` expands to the full list above. Spin operators require a model whose
orbital labels mark the spin channel (`_s+_` / `_s-_` in `LABEL.xyz`); for a
spinless model they evaluate to zero. Exact gauge-transform spin operators
(`--exact-spin`) and orbital angular momentum (`--orbital-L`) are produced under
their own flags.

---

## Input files

For a given `LABEL`, the tool reads:

| File | Required | Content |
|------|----------|---------|
| `LABEL_hr.dat` | yes | Wannier90 real-space Hamiltonian: rows `R1 R2 R3 i j Re Im`. |
| `LABEL.uc`     | for V/S operators | Three lines, each three floats: the unit-cell lattice vectors (Å). |
| `LABEL.xyz`    | for V/S operators | First line: number of orbitals. Then `label x y z` per orbital. |
| `LABEL_wsvec.dat` | no | Wigner-Seitz minimum-image correction from `use_ws_distance`. |
| `LABEL.eig`    | no | Eigenvalues for exact spectral bounds with `--bounds`. |

## Output files

| File | Content |
|------|---------|
| `LABEL.HAM.CSR` | Supercell Hamiltonian, sparse CSR. |
| `LABEL.<OP>.CSR` | One file per requested operator. |
| `LABEL.<OP>.desc` | Physical descriptor (with `--bounds`). |
| `LABEL.<OP>.check` | Self-check report (with `--check`). |
| `LABEL.w2sp/` | Provenance bundle directory (with `--mode bundle`); see [Bundle output](#bundle-output-provenance-for-lsquant). |

CSR text format (one matrix per file):

```
<dim> <nnz>
<value.real value.imag> ...        # nnz complex values
<column index> ...                 # nnz column indices
<row pointer> ...                  # dim+1 row pointers
```

The FULL Hermitian matrix is stored (both triangles). The consumer reconstructs
H = ½(M + M†); for a Hermitian operator this returns M unchanged. A future
upper-triangle mode (flagged storage=upper in the .desc sidecar) would store
col ≥ row only and be reconstructed as U + U† − diag(U); it is enabled only when
the operator passes the Hermiticity check (`--check hermiticity`).

---

## Bundle output (provenance for lsquant)

The default `sparse` mode folds every `H(R)` into one supercell matrix under
periodic boundary conditions, so the CSR carries the *spectrum* of a chosen
supercell but no longer the cell index `R`, and none of the structure or DFT
history that produced the model. `--mode bundle` is the additive alternative: it
ships the **primitive** real-space operators `O_ij(R)` unexpanded, together with a
JSON manifest of full provenance, so a consumer (the twin KPM package **lsquant**)
can build `H(k) = Σ_R e^{ik·R} H(R)`, choose its own supercell, and keep the
crystal structure, symmetry, and DFT/Wannier conditions attached. The positional
CLI and the CSR path are unchanged; bundle mode bypasses the supercell engine by
design, so the supercell dimensions are accepted but ignored.

A bundle is a directory `<LABEL>.w2sp/`:

```
<LABEL>.w2sp/
  manifest.json          # structure + symmetry + DFT/Wannier provenance + operator index
  operators/
    HAM.hr.dat           # primitive O_ij(R), Wannier90 _hr.dat-shaped text (re-ingestible)
    VX.hr.dat  SZ.hr.dat ...
  wsvec.dat              # copied verbatim if a _wsvec.dat was applied
```

Each `operators/<NAME>.hr.dat` is written in the same `_hr.dat` shape the tool
reads, so it round-trips through `--op-file` and the `read_wannier_file` parser.
Values are emitted **post-`ndegen` with `ndegen = 1`** for every block, so the
consumer must not divide again; the manifest's `normalization` block records
`ndegen_applied` / `wsvec_applied` / `truncation_threshold`. Output is
byte-deterministic (canonical `R`/orbital sort, one-based orbital indices, fixed
precision). The `manifest.json` carries `structure` (lattice + reciprocal vectors,
atoms, wannier sites, `num_wann`), `symmetry` (rotation + fractional translation
per operation), `dft_provenance` (code/version, XC functional, spin-orbit /
noncollinear, `ecutwfc`, k-mesh, pseudopotentials), `wannier_provenance`
(num_wann/num_bands, mp_grid, exclude_bands, projections, disentanglement,
`use_ws_distance`), and an `operators` index. Any missing provenance source leaves
that block `null` and sets `provenance_complete: false`; the bundle is always
well-formed.

The provenance sources are parsed from a Quantum ESPRESSO `data-file-schema.xml`
(structure, symmetry, DFT conditions) and a Wannier90 `.win` (wannierisation
conditions). They are supplied through a `run.json` config:

```json
{
  "label": "graphene", "project_dir": ".", "seed": "graphene", "output_dir": "out",
  "mode": "bundle", "operators": ["VX", "VY", "SZ"],
  "exact_spin": false, "orbital_L": false, "emit_bounds": true,
  "truncation_threshold": 1e-8,
  "provenance": { "qe_xml": "scf.save/data-file-schema.xml", "win": "graphene.win" }
}
```

| `run.json` key | Meaning |
|----------------|---------|
| `label`        | System label; names the `<LABEL>.w2sp/` bundle. |
| `project_dir`, `seed` | Resolve the inputs as `<project_dir>/<seed>` (default seed = `label`). |
| `output_dir`   | Parent directory of the bundle. |
| `mode`         | `bundle` (the `--config` driver defaults to bundle if omitted). |
| `operators`    | Operators to build, e.g. `["VX","VY","SZ"]` (same codes as the CLI). |
| `exact_spin`, `orbital_L` | Build the gauge-transform spin / orbital-`L` operators. |
| `emit_bounds`  | Record spectral bounds for the Hamiltonian when a `<seed>.eig` is present. |
| `truncation_threshold` | Echoed into the manifest's `normalization` block. |
| `provenance.qe_xml` | Path to the QE `data-file-schema.xml` (DFT provenance). |
| `provenance.win`    | Path to the Wannier90 `.win` (Wannier provenance). |

Every key is optional: a `run.json` with only a `label` produces a well-formed
bundle with `null` provenance blocks. The two builds below are byte-identical:

```bash
wannier2sparse graphene 1 1 1 VX SZ --mode bundle -o out   # positional
wannier2sparse --config run.json                           # config-driven
```

---

## Physics background

`wannier2sparse` sits at the boundary between Wannier interpolation and
KPM/Chebyshev linear-response methods. The Hamiltonian `H(R)` is read in the
Wannier90 real-space gauge, operators such as velocity `V = -i[H, r]` and spin
`S` are built in the same real-space representation, and the supercell CSR
output feeds kernel-polynomial evaluations of the Kubo–Greenwood / Bastin
conductivity formulas. The Fourier-transform convention and spin/orbital-
angular-momentum definitions used here are source-verified against Wannier90
3.1.0 and Quantum ESPRESSO 7.2 and recorded in `docs/conventions.md`. For the
underlying methodology see the project articles (npj Computational Materials
`s41524…`, `TEGADPS.pdf`) and the Wannier90 user guide.

---

## Runnable examples

The [`examples/`](examples/README.md) directory has self-contained, self-validating
demos: each builds a canonical tight-binding model, expands it with
`wannier2sparse`, and plots the **density of states** from the resulting CSR via
the Kernel Polynomial Method (the same sparse-mat-vec technique LinQT uses
downstream). `examples/validate.py` checks each DOS against its analytic spectrum.

[![graphene DOS](examples/img/graphene_dos.png)](examples/README.md)

```bash
cd examples && bash run.sh graphene 80     # -> models/tb/graphene/graphene_dos.png
```

See [`examples/README.md`](examples/README.md) for the full gallery (chain1d,
graphene, cubic, haldane) and the real-Wannier (DFT) workflow.

---

## Notes on performance

The supercell expansion (`save_supercell_as_csr`) writes the CSR directly from a
single pass over the primitive hoppings — replicate + PBC-wrap into flat arrays
— instead of first building an intermediate string-keyed supercell container.
For large supercells this uses roughly **3× less memory and runs ~3× faster**
than the two-stage expansion, with byte-identical output.
