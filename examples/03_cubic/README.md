# Tutorial 3: how big a crystal do you need before its energy levels blur into a band?

A single atom has sharp levels. Bring a few hundred together and the levels crowd
into a band; bring Avogadro's number and the band looks perfectly smooth. Somewhere
between one site and $10^{23}$ the discrete spectrum turns into a continuous density
of states. How many sites is enough, and what sets the price of asking?

We will compute the density of states $\rho(E)$ of the simplest three-dimensional
solid, the simple-cubic tight-binding band, and watch the answer come out coarser
than its two-dimensional cousins did. The reason is not physics but arithmetic: the
Kernel Polynomial Method (KPM) reconstructs $\rho(E)$ from sparse matrix-vector
products alone, so its cost grows linearly with the number of sites, and in three
dimensions the number of sites is $N^3$. The lesson here is that linear scaling
buys you large systems cheaply, but in 3D the site count grows so fast that the
supercell stays modest and its size becomes the resolution dial you trade against.

## The physics

One orbital per site on a cubic lattice, nearest-neighbour hopping $t$ along each
of the three axes. The Bloch Hamiltonian is a single number,

$$ H(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\, H(\mathbf{R}) = 2t\,(\cos k_x + \cos k_y + \cos k_z), $$

so the band runs over $[-6\lvert t\rvert, 6\lvert t\rvert]$, that is $[-6, 6]$ for
$t = -1$. The density of states is the closed form

$$ \rho(E) = \int \frac{d^3k}{(2\pi)^3}\, \delta\!\big(E - 2t(\cos k_x + \cos k_y + \cos k_z)\big). $$

Unlike the inverse-square-root divergences of the 1D chain, the 3D van Hove
singularities are mild: $\rho(E)$ stays finite, and the only visible feature is a
kink near the band centre $E = 0$, where a saddle point of $E(\mathbf{k})$ changes
the topology of the constant-energy surface. That faint kink is the thing we are
trying to resolve, and it is exactly the kind of feature a coarse spectrum smears.

The conceptual result of this tutorial is about cost, not the band. KPM never
diagonalizes anything: it applies the sparse Hamiltonian to random vectors $R$
times and reads off $M$ Chebyshev moments, so the work is $O(M R \cdot \text{nnz})$
and the number of non-zeros is proportional to the number of sites. Doubling the
sites roughly doubles the cost, with no cubic blow-up from diagonalization. But the
site count itself is $N^3$ in 3D, so a supercell that is trivial in 2D becomes
expensive, and we are forced to keep $N$ small. The moment count $M$ sets the
energy broadening, while the supercell size sets how finely $\mathbf{k}$-space is
sampled. With few sites the sampling is coarse, and the 3D DOS comes out blurrier
than the 2D panels at the same $M$.

![Simple-cubic 3D band over minus six to plus six with a van Hove kink near the centre, coarser than the 2D panels because only a modest supercell fits in 3D](../img/cubic_dos.png)

FIG. 1. Density of states $\rho(E)$ of the simple-cubic band,
$E = 2t(\cos k_x + \cos k_y + \cos k_z)$ with $t = -1$, over the support $[-6, 6]$,
showing the mild van Hove kink near $E = 0$. KPM reconstruction with $M = 2048$
Chebyshev moments, $R = 20$ stochastic vectors, Jackson kernel, on a modest
$18\times18\times18$ supercell. The curve is coarser than the 1D and 2D panels of
the gallery because three dimensions force a small per-side $N$: at $18^3$ there are
$5832$ sites, where an 80-side 2D cell already holds $6400$.

## Step 1: build the cubic model

```bash
bash run.sh cubic 18           # generates models/tb/cubic/ (and a first DOS plot)
cd models/tb/cubic             # the rest of this tutorial runs from here
```

This writes a one-orbital Wannier90 model: an on-site term plus six hoppings of
strength $t = -1$ to the nearest neighbours along $\pm x$, $\pm y$, $\pm z$. The
next step reads these as the primitive operator $H_{ij}(\mathbf{R})$.

## Step 2: run it from an input file

The Hamiltonian is the only operator we need for a density of states, so we expand
it alone. Put

```json
{ "label": "cubic", "mode": "sparse", "supercell": [18, 18, 18] }
```

in `cubic.w2s` and run it:

```bash
wannier2sparse -x cubic.w2s        # -> cubic.HAM.CSR
```

This replicates each primitive $H_{ij}(\mathbf{R})$ across the supercell, PBC-wraps
it, and folds everything into one sparse `cubic.HAM.CSR`.

## Step 3: the cost scales with the sites, but 3D forces a small cell

```bash
python3 ../../../w2s_dos.py cubic.HAM.CSR --title "cubic DOS (KPM)" --out cubic_dos.png
```

![Simple-cubic DOS over minus six to plus six with a van Hove kink near the centre](../img/cubic_dos.png)

The CSR has one non-zero on the diagonal and six off-diagonal per site, so its size
grows linearly with the $N^3$ sites, and so does the KPM work: going from $14^3$ to
$18^3$ roughly doubles the site count and roughly doubles the runtime, with no extra
penalty from a dense solve. That linear scaling is the whole reason KPM is used for
transport. The catch is the $N^3$ itself. An $80\times80\times1$ graphene cell and
an $18\times18\times18$ cubic cell hold comparable numbers of sites, yet one samples
$\mathbf{k}$-space finely in two directions and the other coarsely in three. With
fewer allowed $\mathbf{k}$-points per axis, the 3D spectrum is sampled sparsely, and
at the same $M = 2048$ moments the curve in FIG. 1 looks visibly grainier than the
2D panels: the broadening from $M$ no longer hides the discreteness of the
supercell. To check the shape against an exact, supercell-free reference, diagonalize
$H(\mathbf{k})$ densely on a $\mathbf{k}$-mesh with `python3 ../../../../tools/hr_exactdiag.py dos cubic`,
which reconstructs the same closed-form $\rho(E)$ with no stochastic noise.

## Step 4: ship the model unexpanded with a bundle

Expanding to a CSR commits to one supercell size now. The alternative is to ship the
primitive operator and let a downstream consumer choose its own resolution later.
Set `"mode": "bundle"` in `cubic.w2s` and run it again:

```bash
wannier2sparse -x cubic.w2s
```

In `bundle` mode the supercell dimensions are ignored, because nothing
is expanded. Instead of a folded CSR, the run writes the primitive $H_{ij}(\mathbf{R})$
plus a JSON manifest to `cubic.w2sp/`. The consumer (the KPM package lsquant) holds the
model itself, forms $H(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}} H(\mathbf{R})$
for any $\mathbf{k}$, and rebuilds the supercell at whatever $N$ its own budget allows.
For a 3D model this is the natural division of labour: the resolution dial, which we
just saw is the binding constraint, is set where the calculation actually runs, not
frozen into the file. The trap to avoid is treating the bundle's numbers as raw
Wannier90 values: they are written post-`ndegen` with `ndegen = 1`, so a consumer must
not divide by the Wigner-Seitz degeneracy a second time.

## What to take away

- The simple-cubic band fills $[-6, 6]$ for $t = -1$ and shows only a mild van Hove
  kink near $E = 0$, the finite 3D saddle-point feature, not a divergence.
- KPM costs $O(M R \cdot \text{nnz})$ from sparse matrix-vector products alone, so
  doubling the sites roughly doubles the work, with no cubic cost from diagonalization.
- In 3D the site count is $N^3$, so the supercell stays modest and its size, not the
  moment count, is the binding resolution dial: the 3D DOS is coarser than the 2D one.
- `"mode": "sparse"` expands to a supercell CSR now and fixes that resolution;
  `"mode": "bundle"` keeps the primitive $O_{ij}(\mathbf{R})$ so the consumer sets it later.

The next tutorial reuses this same expand-versus-bundle choice on a model where the
unexpanded $H(\mathbf{R})$ carries topology the folded spectrum alone would hide.

## References and links

- Kernel polynomial method: A. Weiße, G. Wellein, A. Alvermann, H. Fehske,
  Rev. Mod. Phys. 78, 275 (2006),
  [arXiv:cond-mat/0504627](https://arxiv.org/abs/cond-mat/0504627).
- Transport methodology: Z. Fan, J. H. Garcia, A. W. Cummings et al., Linear
  scaling quantum transport methodologies, Phys. Rep. 903, 1 (2021),
  [arXiv:1811.07387](https://arxiv.org/abs/1811.07387).
- wannier2sparse source and documentation: https://github.com/adamecius/wannier2sparse
