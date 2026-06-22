#!/usr/bin/env python3
"""Load a wannier2sparse CSR file and plot its spectral density (DOS) via KPM,
or the sorted eigenvalue spectrum. Only needs numpy/scipy/matplotlib.

Usage:
    python3 w2s_dos.py PATH.HAM.CSR [--title T] [--out fig.png]
                       [--mode dos|spectrum] [--moments 2048] [--vectors 20]

The CSR text format written by wannier2sparse is:
    line 1: dim nnz
    line 2: 2*nnz reals  -> (real, imag) of each stored value, row-major
    line 3: nnz ints     -> column index of each value
    line 4: dim+1 ints   -> row pointers
"""
import sys, argparse, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh

def load_csr(path):
    with open(path) as f:
        dim, nnz = map(int, f.readline().split())
        v   = np.fromstring(f.readline(), sep=' '); val = v[0::2] + 1j*v[1::2]
        col = np.fromstring(f.readline(), sep=' ', dtype=np.int64)
        rp  = np.fromstring(f.readline(), sep=' ', dtype=np.int64)
    H = sp.csr_matrix((val, col, rp), shape=(dim, dim))
    return H

def kpm_dos(H, M=2048, R=20, npts=2000, pad=1.05, seed=0):
    r"""Estimate the density of states of a sparse Hamiltonian by KPM.

    Reconstructs :math:`\rho(E)` from Chebyshev moments
    :math:`\mu_m = \mathrm{Tr}\,T_m(\tilde H)` estimated by a stochastic trace
    over ``R`` random-phase vectors, damped by the Jackson kernel to suppress
    Gibbs oscillations. The Hamiltonian is first rescaled into :math:`[-1,1]`
    from its extremal eigenvalues (with a small ``pad`` margin) so the Chebyshev
    recursion stays bounded. The cost is ``R*M`` sparse matrix-vector products
    and no dense factorization, which is why it scales to the large supercells
    ``wannier2sparse`` emits.

    Parameters
    ----------
    H : scipy.sparse.csr_matrix, shape (n, n)
        Hermitian Hamiltonian in eV; only its mat-vec action is used.
    M : int
        Chebyshev moment count; the resolution dial (broadening ~ bandwidth/M).
    R : int
        Number of stochastic vectors; the trace variance falls as 1/R.
    pad : float
        Multiplicative margin on the half-bandwidth so the spectrum stays
        strictly inside [-1, 1] after rescaling.
    seed : int
        Seed for the random-phase vectors, fixed so the figure is reproducible.

    Returns
    -------
    E : ndarray, shape (npts,)
        Energies in eV.
    rho : ndarray, shape (npts,)
        Density of states per unit energy per site.

    Warnings
    --------
    Returns a meaningless curve if ``H`` is not Hermitian; the routine assumes
    real moments and does not check. The endpoints :math:`E=\pm` band edge are
    excluded because the KPM weight :math:`1/\sqrt{1-x^2}` diverges there.
    """
    n = H.shape[0]
    emax = float(eigsh(H, k=1, which='LA', return_eigenvectors=False)[0])
    emin = float(eigsh(H, k=1, which='SA', return_eigenvectors=False)[0])
    a = (emax - emin) / 2 * pad
    b = (emax + emin) / 2
    Hs = (H - b * sp.identity(n, format='csr')) / a       # rescaled to [-1,1]
    mu = np.zeros(M)
    rng = np.random.default_rng(seed)
    for _ in range(R):
        v = np.exp(2j*np.pi*rng.random(n)); v /= np.linalg.norm(v)
        t0 = v; t1 = Hs @ t0
        mu[0] += np.vdot(v, t0).real
        mu[1] += np.vdot(v, t1).real
        for m in range(2, M):
            t2 = 2*(Hs @ t1) - t0
            mu[m] += np.vdot(v, t2).real
            t0, t1 = t1, t2
    mu /= R
    m = np.arange(M)                                       # Jackson kernel
    g = ((M-m+1)*np.cos(np.pi*m/(M+1)) + np.sin(np.pi*m/(M+1))/np.tan(np.pi/(M+1)))/(M+1)
    mu *= g
    x = np.linspace(-1, 1, npts + 2)[1:-1]  # open interval: avoid 1/sqrt(1-x^2) blow-up at +-1
    T = np.cos(np.outer(np.arccos(np.clip(x, -1, 1)), m))
    rho = (mu[0] + 2*np.sum(mu[1:]*T[:, 1:], axis=1)) / (np.pi*np.sqrt(1 - x**2))
    return a*x + b, rho/(a*n)

def exact_dos(H, eta=0.03, npts=2000, emin=None, emax=None):
    r"""Exact density of states of the sparse Hamiltonian by *dense* diagonalization.

    Diagonalizes the full supercell matrix (``H.toarray()``) and Gaussian-broadens
    its eigenvalues, :math:`\rho(E)=\tfrac1N\sum_n g_\eta(E-E_n)`. This is the
    no-approximation route the KPM is checked against: feasible only for the small
    supercells used in a tutorial (it is O(N^3)), but exact up to the broadening.
    """
    w = np.linalg.eigvalsh(H.toarray())
    if emin is None:
        emin = float(w.min()) - 10 * eta
    if emax is None:
        emax = float(w.max()) + 10 * eta
    E = np.linspace(emin, emax, npts)
    nrm = 1.0 / (np.sqrt(2 * np.pi) * eta)
    rho = nrm * np.exp(-((E[:, None] - w[None, :]) ** 2) / (2 * eta ** 2)).sum(1)
    return E, rho / H.shape[0]                              # per site


def chain1d_analytic_dos(E, t=1.0):
    r"""Closed-form 1D-chain DOS :math:`\rho(E)=1/(\pi\sqrt{4t^2-E^2})`, per site.

    Zero outside the band :math:`|E|<2|t|`; the inverse-square-root van Hove
    divergences at the band edges are the analytic fingerprint the tutorial checks.
    """
    rho = np.zeros_like(E)
    m = np.abs(E) < 2 * abs(t)
    rho[m] = 1.0 / (np.pi * np.sqrt(4 * t * t - E[m] ** 2))
    return rho


OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]


def _style(usetex):
    """agent_template/plotting_style.md single-column preset (matches plot_hall.py)."""
    import matplotlib
    matplotlib.rcParams.update({
        "text.usetex": usetex, "font.family": "serif",
        "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm", "axes.unicode_minus": False,
        "font.size": 9, "axes.labelsize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "legend.frameon": False,
        "lines.linewidth": 1.4, "lines.markersize": 5, "axes.linewidth": 0.6,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "figure.dpi": 600, "savefig.dpi": 600,
    })


def cmd_compare(args):
    """Overlay the DOS from three routes that must agree, then KPM, on one figure:
       closed-form analytic, exact diagonalization of the sparse CSR, and KPM."""
    import shutil, json as _json
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    H = load_csr(args.csr)
    Ek, rk = kpm_dos(H, M=args.moments, R=args.vectors)            # KPM (stochastic)
    rk = rk * H.shape[0]   # kpm_dos returns per-state (1/N); rescale to per-site to match exact/analytic
    Ee, re = exact_dos(H, eta=args.eta)                             # exact diag of the CSR
    usetex = bool(shutil.which("latex") and shutil.which("dvipng"))
    _style(usetex)
    print(f"[w2s_dos] LaTeX path: {'usetex' if usetex else 'mathtext/cm'}")
    fig, ax = plt.subplots(figsize=(3.375, 2.6))
    if args.analytic == "chain1d":                                 # closed form, on its own grid
        Ea = np.linspace(Ee.min(), Ee.max(), 1500)
        ax.plot(Ea, chain1d_analytic_dos(Ea), color="0.55", ls="-", lw=2.6, zorder=1,
                label="analytic")
    if args.overlay:                                               # hr_exactdiag.py dos JSON (on-the-fly H(k))
        d = _json.load(open(args.overlay))
        ax.plot(np.array(d["energy_eV"]), np.array(d["dos"]) / args.norm, color=OKABE_ITO[2],
                ls="--", zorder=2, label="exact diag, $H(k)$")
    me = max(1, len(Ee) // 24)
    ax.plot(Ee, re, color=OKABE_ITO[0], ls="none", marker="o", markevery=me,
            markerfacecolor="white", markeredgecolor=OKABE_ITO[0], markeredgewidth=1.0,
            zorder=4, label="exact diag, sparse CSR")
    ax.plot(Ek, rk, color=OKABE_ITO[1], ls="-.", zorder=3, label="KPM")
    ax.set_xlabel(r"$E\;[\mathrm{eV}]$"); ax.set_ylabel(r"$\rho(E)\;[\mathrm{eV}^{-1}\,\mathrm{site}^{-1}]$")
    if args.ymax: ax.set_ylim(0, args.ymax)
    ax.set_xlim(Ee.min(), Ee.max()); ax.legend(loc="upper center")
    fig.savefig(args.out, bbox_inches="tight", pad_inches=0.02)
    onthefly = ("dashed green, exact diagonalization of H(k) on a dense k-mesh "
                "(on-the-fly, no supercell); " if args.overlay else "")
    cap = args.out.rsplit(".", 1)[0] + ".caption.txt"
    open(cap, "w").write(
        "FIG. N. Density of states rho(E) per site of the 1D tight-binding chain "
        "from routes that agree: solid grey, the closed-form "
        "rho=1/(pi sqrt(4t^2-E^2)); " + onthefly +
        "open blue circles, exact (dense) diagonalization of the expanded sparse "
        "supercell matrix; dash-dotted orange, the kernel polynomial method "
        f"(M={args.moments} moments, R={args.vectors} vectors, Jackson kernel) on "
        f"that same matrix. The exact curves use a Gaussian broadening eta={args.eta} "
        "eV. Chain of t=-1, band [-2,2]; the van Hove edges are resolved as the "
        "broadening allows.\n")
    print(f"[w2s_dos] wrote {args.out} and {cap}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('csr'); p.add_argument('--title', default='')
    p.add_argument('--out', default='dos.png'); p.add_argument('--mode', default='dos',
                   choices=['dos', 'spectrum', 'dos-exact', 'compare'])
    p.add_argument('--moments', type=int, default=2048)
    p.add_argument('--vectors', type=int, default=20)
    p.add_argument('--eta', type=float, default=0.03, help="Gaussian broadening for the exact DOS (eV)")
    p.add_argument('--analytic', default=None, choices=[None, 'chain1d'],
                   help="overlay a known closed-form DOS (compare mode)")
    p.add_argument('--overlay', default=None, help="overlay an hr_exactdiag.py dos JSON (compare mode)")
    p.add_argument('--norm', type=float, default=1.0, help="divide the --overlay dos by this (per-site)")
    p.add_argument('--ymax', type=float, default=None)
    args = p.parse_args()

    if args.mode == 'compare':
        cmd_compare(args); return

    H = load_csr(args.csr)
    herr = abs(H - H.getH()).max() if H.nnz else 0.0
    print("loaded %s  dim=%d  nnz=%d  hermiticity_err=%.2e" % (args.csr, H.shape[0], H.nnz, herr))

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    if args.mode == 'spectrum':                            # exact, for SMALL supercells
        E = np.linalg.eigvalsh(H.toarray())
        plt.plot(np.linspace(0, 1, E.size), E, '.', ms=2)
        plt.xlabel('state index (normalized)'); plt.ylabel('E')
    elif args.mode == 'dos-exact':                         # exact DOS by dense diagonalization
        E, rho = exact_dos(H, eta=args.eta)
        plt.fill_between(E, rho, alpha=.4); plt.plot(E, rho)
        plt.xlabel('E'); plt.ylabel('DOS (exact diag)')
    else:
        E, rho = kpm_dos(H, M=args.moments, R=args.vectors)
        plt.fill_between(E, rho, alpha=.4); plt.plot(E, rho)
        plt.xlabel('E'); plt.ylabel('DOS')
    plt.title(args.title or args.csr); plt.tight_layout(); plt.savefig(args.out, dpi=120)
    print("saved", args.out)

if __name__ == '__main__':
    main()
