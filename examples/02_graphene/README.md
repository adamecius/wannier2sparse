# Tutorial 2: how do we know the smooth DOS curve is the real one?

Graphene has two carbon atoms per cell, yet its two bands touch at a single point
and the density of states drops smoothly to zero there. Run the Kernel Polynomial
Method (KPM) on a finite supercell and you get a clean curve with a dip at zero
energy and two sharp peaks. It looks right. But KPM is stochastic: it probes the
spectrum with a handful of random vectors and a truncated polynomial expansion, so
the curve is an estimate, not the spectrum itself. How do we know it is faithful to
the model and not an artifact of the random seed or the moment count?

We want the density of states of honeycomb graphene, and we want to trust it. The
way to trust an approximation is to hold it against something exact. The same
`_hr.dat` that feeds the supercell can be Fourier transformed back to $H(\mathbf{k})$
and diagonalized densely on a $\mathbf{k}$-mesh, with no supercell and no random
vectors. That dense reference is the oracle. The lesson here is that the exact
diagonalization is the truth the stochastic KPM is measured against, and graphene
is the perfect test because its features are sharp and analytic.

## The physics

Graphene is a honeycomb lattice: two sublattices $A$ and $B$, one $p_z$ orbital
each, nearest-neighbour hopping $t = -1$. The real-space model is a set of blocks
$H(\mathbf{R})$, and the Bloch Hamiltonian is their phased sum,

$$ H(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\, H(\mathbf{R}). $$

With only off-diagonal $A$-$B$ hopping the two bands are $E_\pm(\mathbf{k}) = \pm\,t\,
\lvert f(\mathbf{k})\rvert$, and the density of states follows in closed form: it
vanishes linearly at the Dirac point, diverges logarithmically at the van Hove
energies, and has hard band edges,

$$ \rho(E) \xrightarrow{E\to 0} 0, \qquad \text{van Hove peaks at } E = \pm\lvert t\rvert,
\qquad \text{edges at } E = \pm 3\lvert t\rvert. $$

The conceptual result of this tutorial is that two routes start from the same
$H(\mathbf{R})$ and must agree where both are valid. Dense diagonalization of
$H(\mathbf{k})$ is exact: no supercell, no broadening you did not choose, no
stochastic noise. KPM trades that exactness for linear scaling, since it needs only
sparse matrix-vector products and never forms or diagonalizes a dense matrix. Where
the exact curve is the oracle, the KPM curve is the candidate, and graphene's Dirac
dip and van Hove peaks are sharp enough that any disagreement would show.

![Graphene DOS vanishing linearly at the Dirac point with van Hove peaks at plus and minus the hopping](../img/graphene_dos.png)

FIG. 1. Density of states $\rho(E)$ of honeycomb graphene from the supercell CSR.
$\rho(E)$ vanishes at the Dirac point $E = 0$, with van Hove peaks at $E = \pm\lvert
t\rvert$ and band edges $\pm 3\lvert t\rvert$ for $t = -1$. KPM reconstruction with
$M = 2048$ Chebyshev moments, $R = 20$ stochastic vectors, Jackson kernel, on an
$80\times80\times1$ supercell.

## Step 1: build the honeycomb model

```bash
bash run.sh graphene 80        # generates models/tb/graphene/ (and a first DOS plot)
cd models/tb/graphene          # the rest of this tutorial runs from here
```

This writes the primitive model `graphene_hr.dat` (the two-sublattice Hamiltonian
$H(\mathbf{R})$), its lattice vectors `graphene.uc`, and its orbital positions
`graphene.xyz`. The next steps read those same three files, once as the source of
the supercell and once as the source of the exact reference.

## Step 2: describe the run in an input file, then run it

The operator here is the Hamiltonian itself: `wannier2sparse` replicates every
non-zero $H_{ij}(\mathbf{R})$ across the supercell and PBC-wraps it into one sparse
matrix. The recommended way to drive that is the input-file workflow, which records
the options in an editable `key = value` file and then executes it.

```bash
wannier2sparse --create graphene.inp                     # 1. template
wannier2sparse --write label=graphene    -inp graphene.inp   # 2. populate
wannier2sparse --write supercell 80 80 1  -inp graphene.inp
wannier2sparse --run graphene.inp                        # 3. run -> graphene.HAM.CSR
```

`--run` writes `graphene.HAM.CSR` and a provenance summary `graphene.w2sp.out`.
Prefer the input file because it is self-documenting and reproducible; the older
one-line positional form `wannier2sparse graphene 80 80 1` gives byte-identical
output. (The `--write` forms follow `wannier2sparse --help`: a flag or operator is
a bare token, a setting like the supercell is `key value`.)

## Step 3: the exact bands have no supercell and no noise

`tools/hr_exactdiag.py` takes the opposite route on the same `graphene_hr.dat`: it
rebuilds $H(\mathbf{k})$ and diagonalizes it densely along a $\mathbf{k}$-path.

```bash
python3 ../../../../tools/hr_exactdiag.py bands graphene
```

The two bands meet at the Dirac point and split to the band edges $\pm 3\lvert
t\rvert$. There is no supercell folding and no random seed in this curve, so it is
the analytic dispersion sampled as finely as the path, the reference the supercell
spectrum must reproduce.

## Step 4: the KPM density matches the exact mesh density

The same tool computes the density of states on a dense $\mathbf{k}$-mesh, with only
a Gaussian broadening you set explicitly.

```bash
python3 ../../../../tools/hr_exactdiag.py dos graphene
```

Overlay this exact mesh density on the KPM curve of FIG. 1 and they agree: the
linear vanishing at $E = 0$, the van Hove peaks at $E = \pm\lvert t\rvert$, and the
edges at $\pm 3\lvert t\rvert$ sit at the same energies. The exact density has no
stochastic vectors and no Chebyshev truncation, so any KPM peak that landed at the
wrong energy, or a Dirac dip that failed to reach zero, would stand out against it.
The trap to avoid is reading the KPM curve as exact in its own right: its peak
heights are set by the moment count $M$ and its smoothness by the vector count $R$,
and only the comparison to the dense mesh tells you the broadening is resolution,
not physics.

## What to take away

- Graphene's DOS vanishes linearly at the Dirac point $E = 0$, with van Hove peaks
  at $E = \pm\lvert t\rvert$ and band edges $\pm 3\lvert t\rvert$ for $t = -1$.
- The dense diagonalization of $H(\mathbf{k})$ is the oracle: no supercell, no
  stochastic noise, exact up to the $\mathbf{k}$-sampling you choose.
- KPM trades that exactness for linear scaling, using only sparse matrix-vector
  products, and its broadening is a resolution dial set by the moment count, not a
  physical width.
- Both routes start from the same $H(\mathbf{R})$, so where they overlap they must
  agree, and graphene's sharp features make the agreement a real test.

The next tutorial reuses this exact-versus-KPM habit on a 3D band, where the
supercell is small and the mesh comparison is what keeps the coarse density honest.

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
