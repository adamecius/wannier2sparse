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

A run is described by a single **`.w2s` input file** — JSON that tolerates `//`
and `/* */` comments — and executed by passing it directly:

```
wannier2sparse input.w2s
```

You do not write that file by hand. Scaffold a fully-annotated one with
`--create-template`, or generate it from a short command with `--create`, edit it,
then run it:

```bash
wannier2sparse --create-template graphene            # writes graphene.w2s (every key, documented)
# or, from a command line:
wannier2sparse --create "graphene 50 50 1 VX SZ" -inp graphene   # writes graphene.w2s
# edit graphene.w2s to taste, then run it:
wannier2sparse graphene.w2s
```

Every run also leaves a **`<LABEL>.out` JSON receipt** next to the outputs: the
resolved invocation, per-stage timing, peak memory, the warning/error tally, and
the list of operators written — each naming the input files that produced it. See
[Run receipt](#run-receipt-out).

The input file is documented in [Input file](#input-file-w2s). The scaffolders and
the legacy direct run share the positional grammar `LABEL N1 N2 N3 [OP ... | all]`:

| Argument | Description |
|----------|-------------|
| `LABEL`     | System label. |
| `N1 N2 N3`  | Supercell dimensions along each lattice vector (integers ≥ 1). |
| `OP ...`    | Operators to generate, e.g. `VX SZ`. `all` generates every operator. If omitted, only the Hamiltonian is written. |

The same positional command run directly (e.g. `wannier2sparse graphene 50 50 1 VX
SZ`) still works without an input file, for quick one-offs and backward
compatibility.

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
| `INPUT.w2s` (positional) | Run from a `.w2s` input file (drives both `sparse` and `bundle` modes). Equivalent to `--run INPUT.w2s`. See [Input file](#input-file-w2s). |
| `--run PATH`           | Run from a `.w2s` (or legacy `run.json`) input file. `--input` / `--config` are synonyms. |
| `--create-template [NAME]` | Write a fully-annotated `NAME.w2s` template (default `template.w2s`) and exit. The scaffold is itself runnable. |
| `--create "CMDLINE"`   | Parse a quoted command line (`LABEL N1 N2 N3 [OP...] ...`) and serialize the equivalent `.w2s`, then exit. |
| `-inp, --inp NAME`     | Target stem for `--create` / `--create-template` (else the `LABEL` from the command line, else `template`). |
| `--write`              | Scaffold the equivalent input file to `<output>/<LABEL>.w2s` from the positional arguments and exit, instead of running. |
| `--verbose` / `--quiet` | Console log level `DEBUG` / `WARN` (default `INFO`). |
| `--log-level LVL`      | Set the console level explicitly: `trace`, `debug`, `info`, `warn`, `error`. |
| `--log-file PATH`      | Write the full (TRACE-level) log here. Default: `<output>/<LABEL>.run.log`. |
| `--no-log-file`        | Do not write a log file (console only). |
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
# Scaffold an input file, then run it (the normal workflow)
wannier2sparse --create "graphene 50 50 1 VX VY SZ" -inp graphene   # -> graphene.w2s
wannier2sparse graphene.w2s

# Start from a fully-documented blank template instead
wannier2sparse --create-template graphene                          # -> graphene.w2s
# edit graphene.w2s, then:
wannier2sparse graphene.w2s

# Scaffold a bundle run (primitive operators O_ij(R) + a provenance manifest)
wannier2sparse --create "graphene 1 1 1 VX SZ --mode bundle" -inp graphene
wannier2sparse graphene.w2s
# -> graphene.w2sp/{manifest.json, operators/HAM.hr.dat, operators/VX.hr.dat, ...}

# Scaffold a run with the exact gauge-transform spin operators
wannier2sparse --create "Fe 4 4 4 --exact-spin" -inp Fe
wannier2sparse Fe.w2s
```

An input file written by hand (see [Input file](#input-file-w2s)) is run the same
way: `wannier2sparse my_run.w2s`.

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
| `LABEL.out` | JSON run receipt: invocation, per-step timing, peak memory, warn/error tally, and the per-operator provenance ledger. See [Run receipt](#run-receipt-out). |
| `LABEL.run.log` | Full TRACE-level run log (unless `--no-log-file`). |
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

## Input file (.w2s)

A whole run is described by one `.w2s` input file and executed by passing it
directly (`wannier2sparse input.w2s`; the synonyms `--run` / `--input` / `--config`
also work, and a legacy `run.json` is still accepted). It is JSON that tolerates
`//` line and `/* */` block comments, so a hand-edited file can be self-documenting.
This is more organizable and traceable than a long command line, and it speaks the
same language as the JSON manifest the bundle writes: input and output are
symmetric. The same file drives either output mode via its `"mode"` key. Scaffold
one with `--create-template` (blank, documented) or `--create "..."` (from a command
line); an input file overrides only the keys it sets.

A `sparse` (supercell CSR) run:

```json
{
  "label": "graphene", "mode": "sparse", "output_dir": "out",
  "supercell": [80, 80, 1],
  "operators": ["VX", "VY", "SZ"],
  "spin_currents": [["X", "Z"]],
  "op_files": [{ "name": "OX", "path": "ox_hr.dat" }],
  "checks": "all", "emit_bounds": true,
  "log_level": "info"
}
```

```bash
wannier2sparse run.w2s                   # identical output to the equivalent CLI
```

| Key | Meaning |
|-----|---------|
| `label` | System label (names the outputs / the bundle). |
| `mode` | `sparse` (default for the positional CLI) or `bundle`. A `--run` input file defaults to `bundle` if `mode` is omitted. |
| `project_dir`, `seed` | Resolve inputs as `<project_dir>/<seed>` (seed defaults to `label`). |
| `output_dir` | Where outputs and the log file are written. |
| `supercell` | `[N1, N2, N3]` supercell dimensions (sparse mode). |
| `operators` | Operator codes, e.g. `["VX","VY","SZ"]`. |
| `velocity_mode` | Velocity ladder for `VX/VY/VZ` and the `V` in any spin current: `bare`, `berry_connection` (default), or `covariant` (needs `r_dat`). |
| `r_dat` | Position matrix `_r.dat` (Wannier90 `write_rmn`) for `velocity_mode: covariant`. |
| `spin_currents` | Derived `J = 1/2{V,S}` as `[[V,S], ...]`, e.g. `[["X","Z"]]` (sparse mode). |
| `op_files` | External `_hr.dat` operators as `[{"name","path"}, ...]`. |
| `exact_spin`, `orbital_L` | Build the gauge-route spin / orbital-`L` operators. |
| `checks` | Self-check selector (same values as `--check`). |
| `emit_bounds` | Write `.desc` sidecars (spectral bounds for the Hamiltonian). |
| `truncation_threshold` | Recorded in a bundle manifest's `normalization` block. |
| `provenance` | `{ "qe_xml": ..., "win": ..., "manual": {...} }` — provenance sources to auto-parse, plus a user-declared `manual` block (see below). |
| `log_level`, `log_file` | Console verbosity and an explicit log-file path. |

### Manual provenance

When a model arrives as a bare `_hr.dat` with no machine-readable DFT/Wannier
sources, the `provenance.manual` block lets you declare the conditions by hand. It
is **documentation only** — never read by the numerics — and is echoed unchanged
into the bundle manifest (`manual_provenance`) and the `.out` receipt, so the record
of *how the model was made* travels with the operators. Every field is optional:

```json
"provenance": {
  "manual": {
    "code": "Quantum ESPRESSO 7.2 + Wannier90 3.1.0",
    "xc_functional": "PBE",
    "basis": "plane waves",
    "ecutwfc_Ry": 90,
    "kpoint_grid": [12, 12, 12],
    "pseudopotentials": [{ "species": "Fe", "file": "Fe.UPF", "z_valence": 16 }],
    "orbitals": [{ "index": 0, "label": "Fe-d_xy", "element": "Fe", "xyz": [0.0, 0.0, 0.0] }],
    "notes": "..."
  }
}
```

Positions are Cartesian Å and the cutoff is in Ry, matching the rest of the project.

## Run receipt (.out)

Every run writes `<output_dir>/<LABEL>.out`, a machine-readable JSON receipt that
complements the human console / `.run.log` summary. It records, in one document:

- **invocation** — resolved `label` / `mode`, the input file that drove the run,
  the output directory, the command line, the tool version, and a start timestamp;
- **steps** — the per-stage wall-clock times, in order (a failed run still lists
  how far it got);
- **resources** — total wall time, peak resident memory (MiB), and the
  warning / error tally;
- **outputs** — the ledger: one entry per operator written, with its on-disk path,
  observable / component / units, build provenance, and the **list of input files
  that produced it**. Given an operator, the receipt names exactly which
  `_hr.dat` / `.spn` / `.amn` / op-file went into it;
- **manual_provenance** — the user-declared block above, echoed for a
  self-contained record.

The `.out` is written even with `--no-log-file` (it is a run deliverable, not a
log). Its layout is deterministic, so two runs of the same input differ only in the
timing / memory / timestamp fields.

## Logging and run summary

Every run is tracked. Leveled messages go to the console and, by default, to a
complete log file at `<output_dir>/<LABEL>.run.log` (override with `--log-file`,
disable with `--no-log-file`). Warnings and errors are always reported and
counted, regardless of console verbosity, and any error makes the tool exit
non-zero. The run ends with a summary of per-stage wall times, peak memory, and
the warning/error tally:

```
[2026-01-01 12:00:00.000] [INFO ] ==== run summary ====
[2026-01-01 12:00:00.000] [INFO ] total wall time: 1.204 s
[2026-01-01 12:00:00.000] [INFO ]   load model            0.012 s
[2026-01-01 12:00:00.000] [INFO ]   write HAM             1.100 s
[2026-01-01 12:00:00.000] [INFO ]   write operators       0.090 s
[2026-01-01 12:00:00.000] [INFO ] peak memory: 142.5 MiB
[2026-01-01 12:00:00.000] [INFO ] warnings: 0   errors: 0
```

Console verbosity is set with `--verbose` (DEBUG), `--quiet` (WARN), or
`--log-level LVL`; the file sink always records the full TRACE-level trace.

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
conditions). They are supplied through a `.w2s` config:

```json
{
  "label": "graphene", "project_dir": ".", "seed": "graphene", "output_dir": "out",
  "mode": "bundle", "operators": ["VX", "VY", "SZ"],
  "exact_spin": false, "orbital_L": false, "emit_bounds": true,
  "truncation_threshold": 1e-8,
  "provenance": { "qe_xml": "scf.save/data-file-schema.xml", "win": "graphene.win" }
}
```

| `.w2s` key | Meaning |
|----------------|---------|
| `label`        | System label; names the `<LABEL>.w2sp/` bundle. |
| `project_dir`, `seed` | Resolve the inputs as `<project_dir>/<seed>` (default seed = `label`). |
| `output_dir`   | Parent directory of the bundle. |
| `mode`         | `bundle` (a `--run` input file defaults to bundle if omitted). |
| `operators`    | Operators to build, e.g. `["VX","VY","SZ"]` (same codes as the CLI). |
| `exact_spin`, `orbital_L` | Build the gauge-transform spin / orbital-`L` operators. |
| `emit_bounds`  | Record spectral bounds for the Hamiltonian when a `<seed>.eig` is present. |
| `truncation_threshold` | Echoed into the manifest's `normalization` block. |
| `provenance.qe_xml` | Path to the QE `data-file-schema.xml` (DFT provenance). |
| `provenance.win`    | Path to the Wannier90 `.win` (Wannier provenance). |
| `provenance.manual` | User-declared provenance for a bare `_hr.dat` model; see [Manual provenance](#manual-provenance). |

Every key is optional: a `.w2s` with only a `label` produces a well-formed bundle
with `null` provenance blocks. The two builds below are byte-identical:

```bash
wannier2sparse --create "graphene 1 1 1 VX SZ --mode bundle" -inp graphene   # scaffold graphene.w2s
wannier2sparse graphene.w2s                                                  # run it
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
