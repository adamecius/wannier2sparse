# Tutorial 1: how does a single hopping become a band, and a band a smooth curve?

Take the simplest crystal you can imagine: identical sites on a line, each one
talking only to its two neighbours through a single hopping $t$. There is one
number in the whole model, yet out of it comes a continuous band, and out of the
band a density of states that piles up sharply at the two band edges. A finite
chain has only discrete levels, so where does the smooth curve come from, and how
do we know it is faithful to the model and not to the size of the box we drew it in?

This tutorial answers that for the 1D chain, and in doing so walks the entire
`wannier2sparse` pipeline once, slowly. We start from the primitive real-space
operator $O_{ij}(R)$ that a Wannier90 `_hr.dat` records, expand it into a sparse
supercell matrix, and reconstruct the density of states from that matrix by the
Kernel Polynomial Method. The lesson here is that $O_{ij}(R)$ is the whole model
and the supercell size is the resolution dial on its spectrum.

## The physics

The chain has one orbital per site and a single nearest-neighbour hopping $t$, so
its only non-zero real-space blocks are the on-site term and the two neighbours,

$$ H(R{=}0) = 0, \qquad H(R{=}\pm 1) = t . $$

Bloch's theorem turns these three numbers into a band by summing the real-space
blocks with their phases,

$$ H(k) = \sum_{R} e^{ikR}\, H(R) = t\,(e^{ik} + e^{-ik}) = 2t\cos k . $$

For $t=-1$ the band runs over $E \in [-2, 2]$. The density of states follows from
$\rho(E) = \tfrac{1}{\pi}\,|dk/dE|$, and because $dE/dk = -2t\sin k$ vanishes at
the zone centre and edge, $\rho(E)$ carries the hallmark 1D van Hove signature,

$$ \rho(E) = \frac{1}{\pi\sqrt{4t^2 - E^2}} , $$

an inverse-square-root divergence at each band edge $E = \pm 2|t|$. That closed
form is the thing every later step is checked against.

The conceptual result this tutorial turns on is simpler than the formula: the
three blocks $H(R)$ are the entire model, and nothing the pipeline does adds
physics to them. Expanding to a supercell only chooses how finely $k$ is sampled,
and KPM only chooses how finely $E$ is resolved.

![A single cosine band whose density of states diverges as an inverse square root at both band edges](../img/chain1d_dos.png)

FIG. 1. Density of states $\rho(E)$ of the 1D tight-binding chain. The
characteristic 1D inverse-square-root van Hove divergences sit at the band edges
$E = \pm 2|t|$ and the support is $[-2,2]$ for $t=-1$. KPM reconstruction with
$M = 2048$ Chebyshev moments, $R = 20$ stochastic vectors, Jackson kernel, on a
$400\times1\times1$ supercell.

## Step 1: build the chain as a primitive operator

```bash
bash run.sh chain1d 400       # generates models/tb/chain1d/ (and a first DOS plot)
cd models/tb/chain1d          # the rest of this tutorial runs from here
```

This writes the model in the Wannier90 real-space gauge: a `_hr.dat` holding the
three blocks $H(0)=0$ and $H(\pm 1)=-1$, a `.uc` with the lattice vector, and a
`.xyz` with the single orbital at the origin. Every degeneracy `ndegen` is $1$, so
the numbers are the hoppings themselves with no normalization to undo. The next
step reads these files.

## Step 2: describe the run in an input file, then run it

The recommended way to drive a run is the input-file workflow: record every option
in an editable `key = value` file, then execute it. Create a template, populate the
label and the supercell, and run it.

```bash
wannier2sparse --create chain1d.inp                    # 1. write a template
wannier2sparse --write label=chain1d    -inp chain1d.inp   # 2. populate it
wannier2sparse --write supercell 400 1 1 -inp chain1d.inp
wannier2sparse --run chain1d.inp                       # 3. run -> chain1d.HAM.CSR
```

The supercell engine replicates each primitive block $H(R)$ across $N$ cells and
PBC-wraps it, turning the three-number model into one $N\times N$ sparse matrix
whose eigenvalues are $2t\cos k$ sampled at the $N$ allowed momenta $k=2\pi n/N$.
`--run` writes the expanded `chain1d.HAM.CSR` and a provenance summary of the
resolved options and produced files, `chain1d.w2sp.out`. Because the input file is
self-documenting and reproducible, prefer it; the older one-line positional form
`wannier2sparse chain1d 400 1 1` produces byte-identical output.

## Step 3: the supercell size is the resolution dial on the spectrum

The expanded matrix samples the band at $N$ momenta, so $N$ sets how finely the
spectrum is resolved. A small cell gives a coarse comb of levels; a large cell
fills in the continuous band and lets the van Hove edges sharpen toward their
divergence.

```bash
python3 ../../../w2s_dos.py chain1d.HAM.CSR --out chain1d_dos.png
```

KPM expands $\rho(E)$ in $M=2048$ Chebyshev moments using only sparse
matrix-vector products, the same operation a LinQT transport calculation runs, so
the curve in FIG. 1 scales to large $N$ without any dense diagonalization. Two
knobs are physically distinct and easy to conflate: $N$ sets which energies exist
(the sampling of $k$), while $M$ and the Jackson kernel set how sharply each is
drawn (the broadening). Raising $M$ on too small an $N$ resolves a comb of spikes,
not a band; the smooth curve needs both a large cell and enough moments.

## Step 4: the closed form, checked by exact diagonalization

For a model this small the band can be recovered exactly, with no supercell and no
stochastic noise, by diagonalizing $H(k)$ densely on a $k$-grid:

```bash
python3 ../../../../tools/hr_exactdiag.py bands chain1d
python3 ../../../../tools/hr_exactdiag.py dos   chain1d
```

The `bands` subcommand traces $E = 2t\cos k$ directly, and `dos` reproduces the
inverse-square-root density of FIG. 1. This is the oracle: the KPM curve from the
supercell CSR and the dense-diagonalization curve from the same `_hr.dat` are two
routes to one spectrum, and where they agree the pipeline has added no physics of
its own.

## What to take away

- The 1D chain has a single cosine band $E = 2t\cos k$ over $[-2, 2]$ for $t=-1$,
  with inverse-square-root van Hove divergences at both edges.
- The primitive blocks $H(R)$ are the entire model; expansion and KPM only choose
  how it is sampled and drawn, never what it contains.
- The supercell size $N$ sets the energy sampling and the KPM moment count $M$
  sets the broadening; a smooth, faithful curve needs both.
- The input-file workflow (`--create` / `--write` / `--run`) records a run as an
  editable file and reproduces the positional CLI byte for byte.

Later tutorials reuse this exact pipeline, primitive $O_{ij}(R)$ to supercell CSR
to KPM density, on models where the band structure is no longer a single cosine.

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
