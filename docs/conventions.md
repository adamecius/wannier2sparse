# wannier2sparse — operator conventions (source-verified)

Conventions for building real-space operators `O(R)` (especially the exact spin
operator, Plan 7, and orbital angular momentum, Plan 8) from Wannier90 / QE
output. Every item is verified against the source actually used to generate the
development fixtures and cross-checked on the Fe SOC fixture (`example17`).

**Versions verified against** (built from source on this machine):
- Wannier90 **3.1.0** — `wannier-developers/wannier90`, tag `v3.1.0` (`/tmp/w90build/wannier90-3.1.0`).
- Quantum ESPRESSO **7.2** — `QEF/q-e`, tag `qe-7.2` (`/tmp/qebuild/q-e-qe-7.2`, `include/qe_version.h: version_number='7.2'`).

**Fe fixture** (the empirical anchor): `num_bands = 28`, `num_wann = 18`,
`mp_grid = 8 8 8` → `num_kpts = 512`, `exclude_bands` absent → `nexband = 0`,
disentangled (`28 > 18`).

> The user's confirmations cited the W90/QE `develop` branch; where a literal
> (e.g. a Fortran format string) differs between `develop` and the 3.1.0/7.2
> releases used here, **the release value below is authoritative for our
> fixtures** and the difference is noted.

---

## 1 — Fourier transform for `O(R)` (S_W and H share it)

```
O(R) = (1/N_k) Σ_k exp(−i · 2π · k·R) O(k)
```
- Phase sign `−i`; normalization `1/N_k`; `k` in fractional/reciprocal-lattice
  coordinates (`kpt_latt`); `R` the integer Wigner-Seitz vector (`irvec`).
- **Do NOT divide by `ndegen` on the forward write.** `ndegen` is applied only on
  the inverse interpolation `H(k) = Σ_R e^{+ik·R} H(R)/ndegen(R)` (that is Plan 1,
  on ingest). `S_W(R)` and `H(R)` must use the identical forward transform so the
  two operators live in the same gauge.

**Source** — W90 3.1.0 `src/hamiltonian.F90:299-300` (the `get_HH_R` build of
`H(R)`; identical phase/normalization also at `:320-321`, `:660-661`, `:759-760`):
```fortran
rdotk = twopi*dot_product(kpt_latt(:, loop_kpt), real(irvec(:, irpt), dp))
fac   = exp(-cmplx_i*rdotk)/real(num_kpts, dp)
```

**Empirical anchor**: `Fe_hr.dat` was produced by W90 with exactly this
transform, so reusing it for `S_W(R)` keeps S in the H gauge. Implementation-time
check: `Σ_R e^{+i2πk·R} S_W(R)` must reproduce `S_W(k)` on the WS mesh.

---

## 2 — Disentanglement bookkeeping

```
V(k)   = U_dis(k) · U(k)                    ! (num_bands×num_wann)·(num_wann×num_wann)
S_W(k) = V(k)† · S_B(k) · V(k)              ! num_wann × num_wann
```
- `U_dis(k)`: `num_bands × num_wann` (Fe: 28×18). `U(k)`: `num_wann × num_wann`
  (18×18). `S_B(k)`: `num_bands × num_bands` (the .spn matrix at k).
- The band sets of `.spn`, `.eig`, and the **rows** of `u_dis.mat` must all refer
  to the same bands after `exclude_bands` from `.win`. Assert per k:
  `S_B` is `num_bands²`, `U_dis` is `num_bands×num_wann`, `U` is `num_wann²`;
  abort on mismatch.
- **Non-entangled runs**: `u_dis.mat` is written **only** `if (have_disentangled)`
  (i.e. `num_bands > num_wann`). When absent (`num_bands == num_wann`), take
  `V(k) = U(k)`.

**Source** — W90 3.1.0 `src/plot.F90` (`plot_u_matrices`):
- `:1659` `open(... '_u.mat', form='formatted')`; `:1662` header `num_kpts, num_wann, num_wann`;
  per k: `:1665` blank line, `:1666` `kpt_latt(:,nkp)` as `(f15.10,sp,f15.10,sp,f15.10)`,
  `:1667` `((u_matrix(i,j,nkp), i=1,num_wann), j=1,num_wann)` as `(f15.10,sp,f15.10)`
  — **column-major, row index `i` fastest**.
- `:1671` `if (have_disentangled)`; `:1673` `open(... '_u_dis.mat' ...)`;
  `:1674` header `num_kpts, num_wann, num_bands`;
  `:1679` `((u_matrix_opt(i,j,nkp), i=1,num_bands), j=1,num_wann)` — `num_bands×num_wann`,
  **row index `i` (band) fastest**.

> ⚠️ Header field order ≠ data shape: `_u_dis.mat` header is
> `num_kpts num_wann num_bands` but the data block is `num_bands × num_wann`.
> Each per-k block is: one blank line, one k-point line (3 floats), then the matrix.

**Empirical (Fe)**: `Fe_u.mat` header `512 18 18`; `Fe_u_dis.mat` header
`512 18 28`; `Fe.win` has no `exclude_bands` (`num_bands = nbnd = 28`).

---

## 3 — `.spn` packing order, format, and spin units

**Header** (formatted): line 1 = comment; line 2 = `num_bands  iknum`, where
`num_bands = nbnd − nexband` (post-`exclude_bands`) and `iknum = num_kpts`.

**Packing**: upper triangle, column `m` outer, row `n ≤ m` inner; for each
`(n,m)` the three Pauli components are written σx, σy, σz, each on one line.
Lower triangle recovered by Hermiticity `S(n,m) = conj(S(m,n))`.

**Spin units**: the values are bare Pauli matrix elements `⟨ψ_n|σ_α|ψ_m⟩`, i.e.
spin in units of **ħ/2** (no ħ/2 factor in the file). Record this factor in the
descriptor so `½{V,S}` downstream has unambiguous units.

**Source** — QE 7.2 `PP/src/pw2wannier90.f90`:
- header `:2468-2470` (`spn_formatted`): `OPEN(...".spn",form='formatted')`,
  `WRITE header`, `WRITE nbnd-nexband, iknum`.
- compute/pack `:2499-2512`:
```fortran
counter=0
DO m=1,nbnd
   if(excluded_band(m)) cycle
   DO n=1,m
      if(excluded_band(n)) cycle
      sigma_x = cdum1 + cdum2          ! <n↑|m↓> + <n↓|m↑>
      sigma_y = i*(cdum2 - cdum1)
      sigma_z = <n↑|m↑> − <n↓|m↓>
      counter = counter + 1
      spn(1,counter)=sigma_x ; spn(2,counter)=sigma_y ; spn(3,counter)=sigma_z
```
- write `:2588-2596` (per k-block):
```fortran
counter=0
do m=1,num_bands
   do n=1,m
      counter=counter+1
      do s=1,3
         write(iun_spn,'(2es26.16)') spn(s,counter) + spn_aug(s,counter)
```
> Format is `(2es26.16)` in QE 7.2 (the `develop`/older `(2ES20.10)` differs);
> our fixture uses `2es26.16`.

**Empirical (Fe)**: formatted `Fe.spn` has
`3 · num_bands(num_bands+1)/2 · num_kpts + 2 = 3·406·512 + 2 = 623618` lines (exact).
First counter `(n=m=1)` → σx≈0, σy≈0, σz≈−1.0 (real diagonal ⇒ Hermiticity holds).

---

## 4 — Validation

Self-contained checks to run now:
- **Hermiticity**: `S_W(R) = S_W(−R)†` (component-wise), and each `S_α(k)` Hermitian.
- **Spin sum-rules**: e.g. `Σ_R Tr S_α(R)` against the expected per-cell value.

Defer the quantitative cross-check against an external reference (WannierBerri or
similar) until a reference is available. **Do not fabricate a golden.**

---

## Additional (resolved while reading source)

### `.amn` packing (Plan 8)
- Header: line 1 comment; line 2 `num_bands  num_kpts  num_wann`
  (`nbnd-nexband, iknum, n_wannier`).
- Body: one line per `(band, projector, k)` —
  `(3i5,2f18.12)  m  α  k  Re  Im` — so each entry is **self-indexed** (parse is
  order-independent). Internal loop order is `k` → `ipol` → projector `iw` →
  band `ibnd` (band innermost); for spinor/SOC the column is `iw + n_proj*(ipol-1)`.
- **Source** — QE 7.2 `PP/src/pw2wannier90.f90`: header `:3414-3415`; data write
  `:3499` (spinor) and `:3539` (non-spinor).
- **Empirical (Fe)**: `Fe.amn` header `28 512 18`.

### Spin units
`pw2wannier90` writes σ (Pauli) → spin in **ħ/2**, not `(ħ/2)σ`. Carry the ħ/2
factor explicitly in the operator descriptor (Plan 6 sidecar).

### `u_dis.mat` presence
Absent for non-entangled runs (`num_bands == num_wann`) ⇒ use `V = U`. Present
for Fe (disentangled). See `plot.F90:1671` (`if (have_disentangled)`).

---

## 5 — Orbital angular momentum `L` (Plan 8): real harmonics, fixed-ℓ only

### Principle — why hybrids are out of scope (not a punt)
On-site `L` is well-defined **only within a fixed-ℓ shell**: `L` is block-diagonal
in ℓ (`⟨ℓ‖L‖ℓ′⟩ = 0` for ℓ≠ℓ′) and mixes only the real harmonics inside one ℓ.
A hybrid projector (`sp3`, `sp3d2`, …) is a linear combination spanning ℓ=0,1,2
on the same atom; projecting `L` onto a hybrid basis requires choosing how the
hybrid decomposes into pure-ℓ parts and what "orbital L" even means across an
s/p/d mixture — a **modeling choice, not a convention** you can read off a table.

Therefore this first cut supports a **complete pure shell**: `s` (trivially
`L=0`), a full `p` shell (pz,px,py), or a full `d` shell
(dz2,dxz,dyz,dx2-y2,dxy) per atom. `s` is included because ℓ=0 has a
well-defined (zero) on-site `L`; it is the complete-shell case, not a hybrid.
Anything else **must raise a clear error naming the offending projection**:
- hybrids (`sp`, `sp2`, `sp3`, `sp3d`, `sp3d2`) and `f`+;
- an **incomplete** shell (e.g. only `dxy`, or `dxy;dxz;dyz`) — `L` mixes all
  members of the shell, so a partial shell cannot represent it.

### Real-harmonic order + phase (source)
- W90 3.1.0 `src/parameters.F90:5367-5383` maps labels to `(l, mr)`:
  p → mr 1,2,3 = `pz, px, py`; d → mr 1..5 = `dz2, dxz, dyz, dx2-y2, dxy`.
- `doc/user_guide/projections.tex` Table (l, mr, name, Θ_lmr) gives the explicit
  **Condon–Shortley** real harmonics: `pz∝z/r`, `px∝x/r`, `py∝y/r`;
  `dz2∝3z²−r²`, `dxz∝xz`, `dyz∝yz`, `dx2-y2∝x²−y²`, `dxy∝xy`.
- Cross-check: `QE 7.2 PP/src/pw2wannier90.f90:69` — `l_w, mr_w … as from table
  3.1,3.2`; the `.amn` columns / `.win` projection expansion are emitted in this
  `(l, mr)` order. **The `L_local` basis order must equal the `.amn` column order**
  (a mismatch silently transposes L into the wrong basis → wrong-but-plausible
  values, the worst failure mode).

### `Lx, Ly, Lz` — derived symbolically (do not transcribe literals)
Work in `ħ = 1`. In the complex `|l,m⟩` basis (`m = −l..+l`):
```
Lz|l,m⟩ = m|l,m⟩ ;  L±|l,m⟩ = sqrt(l(l+1) − m(m±1)) |l,m±1⟩ ;  Lx=(L+ +L-)/2 ;  Ly=(L+ −L-)/(2i)
```
Let `C` have as columns the W90 real harmonics expressed in `|l,m⟩` (Condon–Shortley):
- p: `pz = |1,0⟩`, `px = (|1,−1⟩ − |1,+1⟩)/√2`, `py = i(|1,−1⟩ + |1,+1⟩)/√2`.
- d: `dz2 = |2,0⟩`, `dxz = (|2,−1⟩ − |2,+1⟩)/√2`, `dyz = i(|2,−1⟩ + |2,+1⟩)/√2`,
     `dx2-y2 = (|2,−2⟩ + |2,+2⟩)/√2`, `dxy = i(|2,−2⟩ − |2,+2⟩)/√2`.

Then in the real basis (W90 column order):
```
L_α(real) = C† · L_α(complex) · C        (α = x,y,z)
```
For p this yields the standard vector (so(3)) generators `(L_a)_{bc} = −i ε_{abc}`
in the (z,x,y) ordering, eigenvalues of `Lz` = {−1,0,+1}. The implementation
**constructs `C` and `L_α(complex)` from the formulas above and forms `C† L C`**;
the literal matrices are never hand-typed. Units of the emitted operator: **ħ**.

### Projector route + FT (reuses P7)
`A_{m,α}(k)` from `.amn` (`num_bands × num_proj`); `V(k)` from P7
(`num_bands × num_wann`). Then `C(k) = A(k)† V(k)` (`num_proj × num_wann`),
`L_W^α(k) = C(k)† · L_local^α · C(k)` (`num_wann × num_wann`) with `L_local^α`
block-diagonal over the per-shell `L_α(real)` blocks (zero across shells/atoms),
and `L_W^α(R) = (1/N_k) Σ_k e^{−i 2π k·R} L_W^α(k)` — **identical FT to §1**.

### Validation (analogous to §4)
Self-contained, decisive: `L_α = L_α†`; `[Lx,Ly] = i Lz` (catches ordering/phase
errors immediately); `Tr L_α = 0` per shell; `eig(Lz) = {−1,0,1}` (p) /
`{−2,−1,0,1,2}` (d). Defer any external quantitative cross-check.

**Physical caveat — integer eigenvalues are on `L_local`, not `L_W`.** The
emitted `L_W(R)` is the operator *projected onto the Wannier subspace*
(`P_W L P_W`); Wannier functions are not pure atomic harmonics, so `L_W` does
**not** have integer eigenvalues. The integer-eigenvalue gate ({−ℓ..ℓ}) is on
the local generator `L_local`. For `L_W` the meaningful checks are Hermiticity
and `L_W(R)=L_W(−R)†`. Validated end-to-end on the real pure-`d` copper fixture
(example04, `Cu:d`): `Lz_local` eigenvalues `{−2,−1,0,1,2}`, `L_W(R)=L_W(−R)†`
to ~1e-16. Fe still hits the hybrid error (`sp3d2`), as designed.

### ⚠️ Forward note — Fe is *not* a P8 validation target
The Fe `example17` projections are `sp3d2;dxy;dxz;dyz` — a **hybrid (`sp3d2`) plus a
partial d** — so Fe **deliberately hits the hybrid/incomplete-shell error path**.
That is expected and correct for this cut. P8 must be validated against a
**pure p-shell or d-shell-only** Wannier model; **"P8 done" does NOT mean "Fe
orbital L works"** — full Fe orbital L needs the deferred hybrid handling.
