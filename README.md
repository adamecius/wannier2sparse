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

CSR text format (one matrix per file):

```
<dim> <nnz>
<value.real value.imag> ...        # nnz complex values
<column index> ...                 # nnz column indices
<row pointer> ...                  # dim+1 row pointers
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

## Notes on performance

The supercell expansion (`save_supercell_as_csr`) writes the CSR directly from a
single pass over the primitive hoppings — replicate + PBC-wrap into flat arrays
— instead of first building an intermediate string-keyed supercell container.
For large supercells this uses roughly **3× less memory and runs ~3× faster**
than the two-stage expansion, with byte-identical output.
