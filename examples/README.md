# Runnable examples — see wannier2sparse in action

Each example builds a canonical tight-binding model, expands it into a supercell
with `wannier2sparse`, and plots its **density of states (DOS)** from the
resulting CSR using the **Kernel Polynomial Method** — the same technique LinQT
uses downstream, and one that needs only sparse matrix–vector products, so it
scales to large supercells without dense diagonalization.

## Requirements

- `wannier2sparse` built (see the top-level `README.md`). The examples look for
  `../build/wannier2sparse`; override with `W2SP_BIN=/path/to/wannier2sparse`.
- Python 3 with `numpy`, `scipy`, `matplotlib`.

## Run one

```bash
cd examples
bash run.sh graphene 80        # model + supercell size N
# -> graphene/graphene_dos.png
```

Models: `chain1d`, `graphene`, `cubic`, `haldane`. `N` is the supercell size
(1D uses N×1×1, 2D N×N×1, 3D uses ~N^(1/3) per side).

## What each example shows (and the analytic check)

| Model | What you see in the DOS | Analytic |
|---|---|---|
| **chain1d** | a 1D band with van Hove peaks at the edges | support `[-2,2]` (`E=2t cos k`, `t=-1`) |
| **graphene** | **Dirac dip to zero at E=0**, van Hove peaks at `±|t|` | edges `±3`, peaks `±1` |
| **cubic** | 3D band, van Hove kink near center | edges `±6` (`2t(cos kx+cos ky+cos kz)`) |
| **haldane** | a **gap around E=0** (vs graphene's gapless dip) | gap `≈ 3√3·t₂ ≈ 0.78` for `t₂=0.15` |

These are produced by `gen_models.py`; edit it to change parameters
(`t`, `t₂`, `φ`, mass) and re-run.

## What you should see

The plots below are produced by `w2s_dos.py` (KPM, Jackson kernel, fixed
`seed=0`), so a correct run reproduces them. Validated by `validate.py`
(support within tolerance, gap/dip/metal class correct).

### chain1d — 1D van Hove edges
![chain1d DOS](img/chain1d_dos.png)
A single cosine band: support `[-2,2]`, with the characteristic 1D
`1/√(…)` van Hove divergences at the band edges `±2|t|`.

### graphene — Dirac point
![graphene DOS](img/graphene_dos.png)
Van Hove peaks at `E=±|t|`, DOS → 0 at the Dirac point (`E=0`), band edges `±3`.

### cubic — 3D band, van Hove kink
![cubic DOS](img/cubic_dos.png)
Simple-cubic band over `[-6,6]` (`2t(cos kₓ+cos k_y+cos k_z)`), with a van Hove
feature near the centre. A small supercell, so the KPM density is coarser.

### haldane — topological gap
![haldane DOS](img/haldane_dos.png)
The complex next-nearest hopping opens a gap (`≈ ±0.78 = 3√3·t₂` for `t₂=0.15`);
contrast with graphene's gapless Dirac dip above.

## Plot any CSR yourself

`w2s_dos.py` works on any `*.HAM.CSR` (or any operator CSR) the tool writes:

```bash
python3 w2s_dos.py graphene/graphene.HAM.CSR --title "graphene" --out g.png
python3 w2s_dos.py chain1d/chain1d.HAM.CSR --mode spectrum --out band.png   # small cells
```

- `--mode dos` (default): KPM spectral density (`--moments`, `--vectors` control
  resolution/smoothness).
- `--mode spectrum`: exact eigenvalues (dense; use only on small supercells). For
  1D the sorted spectrum traces the band `E(k)`.

## A note on bands vs. DOS

`wannier2sparse` exports the **real-space supercell Hamiltonian** (all `H(R)`
folded into one matrix under PBC), so the quantity that comes straight out of the
CSR is the **spectrum / DOS**, not `E(k)`. That is exactly what a KPM transport
calculation consumes. A full `E(k)` band path would need the per-`R` blocks
before folding; for the 1D chain, `--mode spectrum` already recovers the
dispersion because its sorted eigenvalues are the band sampled at the supercell's
allowed `k`-points.

## Real Wannier90 model — silicon (WS minimum-image **on**)

Unlike the ideal TB models above, this is a **real DFT-derived** Wannier model
(silicon, `sp3`, SOC-free) with genuine Wigner–Seitz degeneracies and a
`_wsvec.dat`. The tool **auto-detects** WS handling from the presence of
`<seed>_wsvec.dat` (no flag): it applies the minimum-image correction, exactly the
real workflow.

```bash
# fixture lives outside the repo (regenerate it — see below); point the tool at it
cd models/wannier/silicon                       # or wherever the fixture is
$W2SP_BIN silicon 8 8 8                          # WS auto-on (silicon_wsvec.dat present)
python3 /path/to/examples/w2s_dos.py silicon.HAM.CSR \
        --title "silicon (W90, WS) DOS" --out silicon_dos.png
```

![silicon DOS](img/silicon_dos.png)

- **Why `8 8 8`, not `4 4 4`?** The real Wannier `H(R)` reaches out to `|R| = 3`,
  so the minimum-image guard requires a supercell `N ≥ 2·range+1 = 7` per axis
  (smaller cells would alias distinct bonds and corrupt the operator — the tool
  refuses them). Ideal TB models have `range = 1`, hence the small `N` there.
- **Validation (against its own `.eig`, not a closed form):** the DOS support
  `[-5.83, 16.40] eV` lies inside the DFT eigenvalue window `[-5.82, 19.44] eV`
  and shows a clear **gap** — consistent with silicon being an insulator. A wrong
  WS treatment typically smears or closes that gap.
- **Inputs are not committed** (each `_hr.dat`/`_wsvec.dat` is ~0.3 MB; the repo
  keeps only the plot). Regenerate the fixture with
  `../../../test/fixtures/gen_fixture.sh` (see `../test/fixtures/README.md`). A
  HAM/DOS run needs only `silicon_hr.dat` + `silicon_wsvec.dat` (+ a `silicon.uc`
  with the three lattice vectors from `silicon.win`). A velocity operator would
  also need `silicon.xyz` — Wannier90's `_centres.xyz` is **not** this tool's
  format, so build one as: first line `num_wann`, then the first `num_wann`
  `label x y z` rows of `_centres.xyz` (the Wannier centres; drop the comment line).

## Orbital angular momentum — copper (`--orbital-L`)

Copper (`Cu:d` complete ℓ=2 shell + two interstitial `s`) is the **orbital-L**
example. Correctness is anchored to the repo's own test, not re-derived here:

```bash
W2SP_ORBITAL_FIXTURE=/path/to/copper/copper \
  ctest --test-dir ../../../build -R orbital_L_fixture --output-on-failure   # PASS
$W2SP_BIN copper 7 7 7 --orbital-L          # writes copper.{LX,LY,LZ}.CSR (range 3 -> N>=7)
```

![copper Lz ladder](img/copper_Lz.png)

- The **on-site local** `Lz` on a complete shell has the integer ladder
  `{−2,−1,0,1,2}` (d) and `0` (each s) — exactly what `orbital_L_fixture`
  verifies. The figure shows those eigenvalues for copper's 7 Wannier orbitals.
- **Honesty note:** the operator the tool *exports*, the **projected** supercell
  `L_W(R)` (`copper.LZ.CSR`), is Hermitian but **not** integer-valued — Wannier
  functions are not pure spherical harmonics. The integer ladder lives on the
  local block, not the projected one. (See `docs/conventions.md` §5.)

## Exact spin / spin–orbit — iron (`--exact-spin`)

BCC iron with SOC (`spinors`) is the **spin** example — the only fixture with a
`.spn`. Spin-operator correctness is carried by tests, not by a figure:

```bash
ctest --test-dir ../../../build -R gauge_spin --output-on-failure              # PASS (self-contained)
# with the Fe fixture + committed goldens present, the independent cross-checks:
W2SP_SPIN_FIXTURE=/path/to/Fe/Fe \
  ctest --test-dir ../../../build -R 'wberri_(matrix|texture)_crosscheck'      # PASS
$W2SP_BIN Fe 13 13 13 --exact-spin           # writes Fe.{SXexact,SYexact,SZexact}.CSR
```

![iron spin texture](img/iron_spin_texture.png)

- The plot is the **band spin texture** `⟨S_α⟩_{nk}` (ħ/2) from the
  WannierBerri-validated golden (`test/golden/Fe_S_texture.ref`; matrix agreement
  `2.6e-9`, texture `1.3e-3` — see `docs/conventions.md` §7).
- **Why no DOS plot for Fe?** Its `H(R)` reaches `|R| = 6`, so the guard requires
  `N ≥ 13`; an `18`-Wannier, range-6 supercell at `13³` produces a **>1 GB** CSR —
  impractical for the KPM demo on a laptop. The spin operator is validated by the
  tests above; the heavy Fe fixture (~80 MB) is regenerated by
  `../../../test/fixtures/gen_fixture.sh`, never committed.

For the full fixture catalogue and the WannierBerri cross-validation, see
`../test/fixtures/README.md`.
