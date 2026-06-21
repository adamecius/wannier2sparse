# Building operators from Wannier90 — carrying the full `U` matrix

How `wannier2sparse` turns a Wannier90 / Quantum ESPRESSO run into the real-space
operators `O_ij(R)` it exports. The emphasis here is the **gauge route**: the path
that keeps the *whole* disentanglement+rotation matrix `V(k) = U_dis(k)·U(k)` and
uses it to map any Bloch-space (DFT) operator into the Wannier basis exactly,
rather than reconstructing the operator from geometry or orbital labels.

This is the "how/why/what-you-need" companion to
[docs/conventions.md](conventions.md), which holds the source-verified format and
sign conventions. Where the two overlap, `conventions.md` is authoritative; this
doc explains the construction and points at the code.

---

## 1 — The single primitive object: `O_ij(R)`

Every operator this tool produces — Hamiltonian, velocity, spin, orbital `L` — is
the same data structure: a real-space, primitive-cell operator

```
O_ij(R)        i,j = Wannier orbitals in the home cell,   R = integer cell vector
```

stored as a [hopping_list](../include/hopping_list.hpp) (a flat list of
`(R, value, (i,j))` records). The physics of a tight-binding model lives entirely
in this object: the Bloch-space operator at any `k` is recovered by the Wannier
interpolation

```
O(k) = Σ_R e^{+i 2π k·R} O(R) / ndegen(R)
```

and the Hamiltonian is just the special case `H(k) = Σ_R e^{+i2πk·R} H(R)/ndegen(R)`.
Keeping `R` (instead of collapsing it into an expanded supercell CSR) is exactly
what lets a downstream consumer rebuild `O(k)` itself — see `--mode bundle` below
and [include/bundle_writer.hpp](../include/bundle_writer.hpp).

There are **two ways** an operator gets into `O_ij(R)` form, and the difference is
the whole point of this document:

| Route | Builds the operator from… | Needs the `U` matrix? | Operators |
|---|---|---|---|
| **Geometric / label** | the Hamiltonian `H_ij(R)` + Wannier-center geometry + orbital spin labels | no | velocity `V`, density, label-spin `S`, spin-current |
| **Gauge (`U`-matrix)** | the Bloch operator `O_B(k)` projected through `V(k)=U_dis U` | **yes — the full V(k)** | exact spin `S`, orbital `L` |

The geometric route is cheap and self-contained but approximate/heuristic (it
infers spin from `_s+_`/`_s-_` label tags, and velocity from finite Wannier-center
displacements). The gauge route is exact: it carries the complete unitary that
Wannier90 used to go from Bloch states to Wannier functions, so any operator known
in the DFT/Bloch basis is transported into the Wannier basis with no modeling
assumption.

---

## 2 — The gauge transform (the "whole `U` matrix" route)

### 2.1 What the `U` matrices are

Wannier90 builds Wannier functions from Bloch states in two steps, both `k`-local:

- **Disentanglement** `U_dis(k)` — a `num_bands × num_wann` (semi-)unitary that
  projects the larger entangled band manifold onto the optimal `num_wann`-dim
  subspace. Present only when `num_bands > num_wann` (an *entangled* run); written
  to `<seed>_u_dis.mat`.
- **Gauge rotation** `U(k)` — a `num_wann × num_wann` unitary (the maximal-
  localization / smoothness rotation) inside that subspace; written to
  `<seed>_u.mat`.

Their product is the **composite gauge matrix**

```
V(k) = U_dis(k) · U(k)          (num_bands × num_wann)
```

and for a non-entangled run (`num_bands == num_wann`, no `_u_dis.mat`) we take
`V(k) = U(k)`. This is exactly [gauge_data::V](../include/gauge.hpp) — note the
code keeps the *full per-`k`* list of `U`, `U_dis` matrices, not a single averaged
gauge. That is what "keeping track of the whole `U` matrix" means in practice: the
reader [read_gauge](../src/gauge.cpp) stores `std::vector<Eigen::MatrixXcd> U,
Udis` indexed by `k`, with the per-`k` k-point in `kpt[k]`.

### 2.2 Bloch → Wannier → real space

Given any operator known in the **Bloch basis** at each `k`, `O_B(k)` (a
`num_bands × num_bands` matrix), the Wannier-gauge operator is

```
O_W(k) = V(k)† · O_B(k) · V(k)            (num_wann × num_wann)
```

and the real-space operator is the Wannier90 forward Fourier transform

```
O_W(R) = (1/N_k) Σ_k e^{−i 2π k·R} O_W(k)
```

with `k` in fractional coordinates and `R` the integer Wigner-Seitz vectors taken
from `H`'s `_hr.dat` (so `S`/`L` land on the **same** `R`-grid and gauge as `H`).
Phase sign `−i`, normalization `1/N_k`, and **no `ndegen` on the forward write** —
`ndegen` is applied only on the inverse interpolation. All of this is fixed in
[conventions.md §1–§2](conventions.md) and implemented verbatim in
[exact_spin_operator / orbital_L_operator](../src/gauge.cpp).

### 2.3 Exact spin `S_α(R)`

The Bloch operator is the spin (Pauli) matrix `S_B^α(k) = ⟨ψ_m|σ_α|ψ_n⟩` read from
the `pw2wannier90` `<seed>.spn` file. Then

```
S_W^α(k) = V(k)† · S_B^α(k) · V(k)
S_α(R)   = (1/N_k) Σ_k e^{−i2πk·R} S_W^α(k)
```

Units are **ħ/2** (bare Pauli; the file carries no ħ/2 factor). Code:
[exact_spin_operator](../src/gauge.cpp); `.spn` packing/units in
[conventions.md §3](conventions.md). This is what `--exact-spin` produces
(`SXexact/SYexact/SZexact`), and it is the spin operator to trust for SOC /
noncollinear systems — the label route below cannot represent those.

### 2.4 Orbital angular momentum `L_α(R)` (projector route)

`L` is built on the *same* gauge `V(k)`, plus the projection overlaps
`A_{m,α}(k)` from `<seed>.amn` (`num_bands × num_proj`):

```
C(k)   = A(k)† · V(k)                      (num_proj × num_wann)
L_W^α(k) = C(k)† · L_local^α · C(k)        (num_wann × num_wann)
L_α(R) = (1/N_k) Σ_k e^{−i2πk·R} L_W^α(k)
```

`L_local^α` is the on-site angular-momentum generator, block-diagonal over shells,
built symbolically in the `|l,m⟩` basis and rotated into Wannier90's real-harmonic
order ([shell_L](../include/local_operators.hpp)). The shell list (which projector
column is `s`/`p`/`d`) is parsed from the `.win` `begin projections` block
([parse_projection_shells](../src/gauge.cpp)). Units **ħ**. Only complete pure
`s`/`p`/`d` shells are in scope; hybrids (`sp3`, `sp3d2`, …) and partial shells
raise a clear error — see [conventions.md §5](conventions.md) for why this is a
convention boundary, not a punt. Code:
[orbital_L_operator](../src/gauge.cpp).

> Physical caveat: integer eigenvalues hold on `L_local`, not on the emitted
> `L_W(R)` — Wannier functions are not pure atomic harmonics, so `L_W = P_W L P_W`
> need not have integer eigenvalues. The meaningful checks on `L_W` are
> Hermiticity and `L_W(R) = L_W(−R)†` ([conventions.md §5](conventions.md)).

---

## 3 — The geometric / label route (no `U` matrix)

These operators are derived from `H_ij(R)` directly and need only `.uc`, `.xyz`,
`_hr.dat`. Implemented in [tbmodel](../include/tbmodel.hpp).

- **Hamiltonian `H_ij(R)`** — read from `<seed>_hr.dat`, divided by `ndegen(R)`
  on ingest ([create_hopping_list](../include/hopping_list.hpp)); optional
  Wigner-Seitz minimum-image refinement from `_wsvec.dat`. Units eV.
- **Velocity / current `V_d = −i · dr_d · H`** — each hopping is multiplied by the
  Cartesian displacement between its initial and final Wannier centers
  (`r_j + R − r_i`), so it needs the centers (`.xyz`) and lattice (`.uc`). Units
  eV·Å. [createHoppingCurrents_list](../include/tbmodel.hpp).
  - **Warning — interband quantities need the Berry connection.** This geometric
    velocity is the gradient `∂H/∂k` (plus the diagonal Wannier centers). It gives
    the correct band group velocity and is sufficient for the DOS and the
    longitudinal `σ_xx`. It is **not** sufficient for interband responses
    (anomalous/spin Hall): the off-diagonal velocity matrix elements `⟨n|v|m⟩`
    require the inter-Wannier Berry connection `A_d(R)=⟨0i|r_d|Rj⟩`. The full
    **covariant velocity** is `v_d(R) = i(R·lat)_d H(R) − i[H, A_d](R)`, with
    `A_d(R)` from Wannier90's `<seed>_r.dat` (`write_rmn=.true.`). To-do: a
    `--covariant-velocity` route that reads `_r.dat` and adds the `−i[H,A_d]`
    commutator; until then, build it externally and ingest via `--op-file`, and
    verify the sign/normalization against Wannier90 `postw90/get_oper.F90`
    (`get_AA_R`)/`berry.F90`.
- **Density** — keep onsite (`R=0`) terms, zero translations.
  [createHoppingDensity_list](../include/tbmodel.hpp).
- **Label-spin `S_α`** — built from the *spin doubling*, not from `H`'s onsite
  graph: orbitals are paired up/down by position + base label via the `_s+_` /
  `_s-_` tags, and each spin pair gets the onsite 2×2 Pauli block
  (`S_z=diag(+1,−1)`, etc.). Pauli convention, units ħ/2. Returns zero for a
  spinless model. [createHoppingSpinDensity_list](../include/tbmodel.hpp).
  *This is a heuristic*: it requires the orbital labels to mark spin and assumes a
  collinear spin frame. For SOC/noncollinear models use `--exact-spin` (§2.3).
  If the up/down pairing is ambiguous or incomplete (common for real spinor
  Wannier models), `map_id2partner` now **throws a clear error** (it no longer
  asserts/aborts) pointing to the `--exact-spin` (`.spn`) route or to ingesting a
  precomputed spin operator via `--op-file S{X,Y,Z} <seed>_S?_hr.dat`. A
  **tight-binding** user with `*_spin_hr.dat` (or per-component `*_S?_hr.dat`)
  passes them the same way and skips the label heuristic entirely.
- **Spin current `J_d S_α = ½{V_d, S_α}`** — in CSR mode the production form is the
  **matrix anticommutator** `½(V_d S_α + S_α V_d)` assembled after supercell
  expansion (`--spin-current`), which is correct for any `α`. The primitive,
  element-wise `createHoppingSpinCurrents_list` is exact only for the **diagonal**
  `S_z` (`J^z` = `v·s_i` on the same-spin block, **zero on the opposite-spin
  block** — the latter was previously left un-zeroed, a bug for SOC, now fixed);
  for off-diagonal `S_{x,y}` it cannot be written element-wise, so use
  `--spin-current`. Units eV·Å·ħ/2.
  [createHoppingSpinCurrents_list](../include/tbmodel.hpp).

---

## 4 — What you need from Quantum ESPRESSO / Wannier90

The DFT half of the pipeline is a standard `pw.x` → `pw2wannier90.x` → `wannier90.x`
run. Per operator, the inputs are:

| File | Produced by | Carries | Needed for |
|---|---|---|---|
| `<seed>_hr.dat` | `wannier90.x` | `H_ij(R)`, `ndegen` | **everything** (sets the `R`-grid + gauge) |
| `<seed>.uc` | you / converter | lattice vectors (Å) | velocity, spin-current, manifest structure |
| `<seed>.xyz` | you / converter | Wannier-center positions + spin label tags | velocity, label-spin, spin-current |
| `<seed>_wsvec.dat` | `wannier90.x` (`use_ws_distance`) | minimum-image `T`-vectors | optional WS refinement of every operator |
| `<seed>.eig` | `pw2wannier90.x` | band energies | optional Hamiltonian spectral bounds `(a,b)` |
| **`<seed>_u.mat`** | `wannier90.x` | **`U(k)`** (num_wann²) | **gauge route** (exact spin, orbital `L`) |
| **`<seed>_u_dis.mat`** | `wannier90.x` (entangled only) | **`U_dis(k)`** (num_bands×num_wann) | **gauge route** when `num_bands > num_wann` |
| **`<seed>.spn`** | `pw2wannier90.x` (`write_spn`) | `S_B^α(k) = ⟨ψ|σ_α|ψ⟩` | **exact spin** |
| `<seed>.amn` | `pw2wannier90.x` | projection overlaps `A(k)` | **orbital `L`** |
| `<seed>.win` | you | `begin projections` block | **orbital `L`** (shell list) |

To get the gauge files written, the relevant switches are:

- `wannier90.x`: `write_u_matrices = .true.` (writes `_u.mat`, and `_u_dis.mat`
  when disentangling). The exact-spin/`L` routes need the **per-`k`** `U` matrices,
  so these must be on disk — a `.chk`-only run is not enough for this tool.
- `pw2wannier90.x`: `write_spn = .true.` for `.spn` (requires a noncollinear/SOC
  `pw.x` run for a meaningful spin texture); `.amn` is written by the standard
  projection step.

Provenance you may *also* want recorded in the bundle manifest (lattice, atoms,
symmetry ops, k-mesh, XC functional, SOC flag, ecutwfc, pseudopotentials) is read
from the QE `data-file-schema.xml` and the `.win`; see
[include/system_provenance.hpp](../include/system_provenance.hpp) and the bundle
plan. That metadata is a sidecar — it never enters the numeric kernel.

### Convention traps (read before trusting numbers)

- **Same gauge for `H` and `S`/`L`.** `S_α(R)` and `L_α(R)` are Fourier-transformed
  on `H`'s own `R`-set with the identical forward transform, so they live in the
  Hamiltonian's gauge. Mixing transforms silently rotates the operator.
- **`ndegen` once.** Applied on ingest of `_hr.dat` / on the inverse interpolation
  only — never double-divide ([conventions.md §1](conventions.md)).
- **`.spn` units are ħ/2** (bare Pauli), `L` is in ħ — the descriptor records the
  factor so downstream `½{V,S}` has unambiguous units.
- **WS convention must match** across any cross-check (`use_ws_distance`,
  TB convention I vs II) — [conventions.md §6–§7](conventions.md).
- **`_u_dis.mat` is absent for non-entangled runs**; the reader falls back to
  `V=U` automatically.

---

## 5 — Compatible commands

The CLI is `wannier2sparse LABEL N1 N2 N3 [OP ... | all] [options]`
([include/w2sp_arguments.hpp](../include/w2sp_arguments.hpp)). The geometric
operators are positional; the `U`-matrix (gauge) operators are opt-in flags.

```bash
# Hamiltonian only (50x50x1 supercell)
wannier2sparse graphene 50 50 1

# Geometric route: velocity + label-spin
wannier2sparse graphene 50 50 1 VX VY SZ

# Spin current J = 1/2{V,S}  (derived, post-expansion)
wannier2sparse graphene 50 50 1 --spin-current X Z

# Gauge route: EXACT spin from the full U matrix + .spn   (-> SXexact/SYexact/SZexact)
#   reads Fe_u.mat (+ Fe_u_dis.mat if entangled) + Fe.spn
wannier2sparse Fe 4 4 4 --exact-spin

# Gauge route: orbital angular momentum L from U matrix + .amn + .win projections
wannier2sparse Cu 4 4 4 --orbital-L

# Everything, with spectral bounds + self-checks
wannier2sparse graphene 50 50 1 all --bounds --check all -o out

# Bundle mode: emit primitive O_ij(R) + JSON manifest (lsquant rebuilds H(k) itself).
# Supercell dims are ignored here; R is preserved.
wannier2sparse graphene 1 1 1 VX SZ --exact-spin --mode bundle -o out
#   -> out/graphene.w2sp/{manifest.json, operators/*.hr.dat}
```

Useful selectors:

- `-p, --project DIR` / `--seed NAME` — locate inputs at `<DIR>/<NAME>...`.
- `--op-file NAME PATH` — ingest an external `_hr.dat`-format operator verbatim.
- `--bounds` — write `.desc` sidecars (units, provenance, `(a,b)` for `H`).
- `--check [all|hermiticity|sum_rules|algebra|aliasing|bounds]` — self-verify;
  `algebra` checks `[Lx,Ly]=iLz` / spin SU(2) without an external reference.
- `--mode bundle` — the `R`-preserving export described in §1.

### Quick verification

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
```

The gauge route is cross-checked end-to-end: the emitted `S_W(R)` inverse-transforms
back to `V(k)†S_B(k)V(k)` on the mp_grid
([test/spin_roundtrip_crosscheck.cpp](../test/spin_roundtrip_crosscheck.cpp), max
error ≈ 2.7e-9 on the Fe fixture), and matches an independent WannierBerri golden to
~1e-9 at the matrix level ([conventions.md §7](conventions.md)).

---

## 6 — Pointers

- Gauge data + exact operators: [include/gauge.hpp](../include/gauge.hpp),
  [src/gauge.cpp](../src/gauge.cpp)
- Geometric operators: [include/tbmodel.hpp](../include/tbmodel.hpp)
- Local `L` generators: [include/local_operators.hpp](../include/local_operators.hpp)
- Core data structure + FT: [include/hopping_list.hpp](../include/hopping_list.hpp)
- Format/sign conventions (authoritative): [docs/conventions.md](conventions.md)
- CLI: [include/w2sp_arguments.hpp](../include/w2sp_arguments.hpp), [src/main.cpp](../src/main.cpp)
</content>
</invoke>
