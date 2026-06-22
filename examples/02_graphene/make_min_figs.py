#!/usr/bin/env python3
"""Generate the two density-of-states figures of tutorial 02 from the committed
inputs, contrasting the *minimal Wannier* file set against the *minimal
tight-binding* file set.

The point of tutorial 02 is the minimal input for a real Wannier calculation. A
Wannier model carries genuine Wigner-Seitz degeneracies (``ndegen > 1`` in
``graphene_hr.dat``) and ships a ``graphene_wsvec.dat`` that the tool auto-detects
for the minimum-image correction; that file is part of the minimal set, not an
extra. This script makes both halves of that statement visible:

  * Panel A (``graphene_wannier_dos.png``) -- the DOS of the real DFT-derived
    Wannier graphene, computed once with the ``_wsvec.dat`` correction applied
    (the correct, minimal-Wannier result) and once with it withheld. For this
    well-localized model the two curves nearly coincide (the WS effect on the
    integrated DOS is small, ~0.004/eV/site, shown in the zoom inset); the file
    still belongs in the minimal set because the tool auto-applies it and the
    corrected operator is what downstream code expects.

  * Panel B (``graphene_tb_dos.png``) -- the DOS of the idealized two-band
    tight-binding graphene (``graphene_tb_hr.dat`` only, ``ndegen = 1``, no
    ``_wsvec.dat``): the minimal tight-binding file set of tutorial 01 applied to
    the honeycomb lattice. Clean Dirac dip at E=0, van Hove peaks at +-|t|, edges
    at +-3|t| for t=-1.

Each figure obeys agent_template/plotting_style.md (single-column width, Okabe-Ito
palette with redundant linestyle+marker encoding, no title, 600 dpi) and is shipped
with a self-contained ``FIG.`` caption sidecar.

Usage:  python3 make_min_figs.py          # writes ../img/graphene_{wannier,tb}_dos.png
Set W2SP_BIN=/path/to/wannier2sparse if the binary is not at ../../build/wannier2sparse.
The tool is invoked as `wannier2sparse -x <file>.w2s`.

Pitfall: the Wannier _hr.dat reaches |R| = 8 in-plane, so the minimum-image guard
requires a supercell N >= 2*range+1 = 17 per axis; smaller cells alias distinct
bonds and the tool refuses them. The idealized TB model has range 1, hence any N
works there.
"""
import os, sys, shutil, subprocess, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))                 # examples/ -> w2s_dos
import w2s_dos as W

BIN = os.environ.get("W2SP_BIN", os.path.join(HERE, "..", "..", "build", "wannier2sparse"))
IMG = os.path.join(HERE, "..", "img")

# Both panels use dense diagonalization (exact_dos): on a few-thousand-state cell it
# is fast and gives a smooth, broadening-controlled curve, whereas KPM at the fine
# resolution needed to expose the ~0.06 eV Wigner-Seitz shift would instead resolve a
# comb of individual levels. N >= 17 satisfies the in-plane minimum-image guard.
N_WANNIER = 60            # 2*60*60 = 7200 states
ETA_W     = 0.06          # eV, Gaussian broadening (resolves the WS shift, stays smooth)
EMIN_W, EMAX_W = -9.8, 4.9   # common energy window so the on/off curves share a grid
N_TB      = 60            # 2*60*60 = 7200 states
ETA_TB    = 0.04          # |t| (eV), broadening for the idealized TB DOS


def run_w2s(workdir, seed, n1, n2, n3):
    """Run wannier2sparse in ``workdir`` on ``seed`` and return the CSR path."""
    w2s = os.path.join(workdir, seed + ".w2s")
    with open(w2s, "w") as f:
        f.write('{ "label": "%s", "mode": "sparse", "supercell": [%d, %d, %d] }\n'
                % (seed, n1, n2, n3))
    subprocess.run([BIN, "-x", w2s], cwd=workdir, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return os.path.join(workdir, seed + ".HAM.CSR")


def stage(workdir, *names):
    for n in names:
        shutil.copy(os.path.join(HERE, n), os.path.join(workdir, n))


def panel_wannier(ax):
    """DOS of the real Wannier graphene, Wigner-Seitz correction on vs off, with a
    zoom inset auto-placed where the two curves separate most."""
    with tempfile.TemporaryDirectory() as on, tempfile.TemporaryDirectory() as off:
        # ON: ship the _wsvec.dat so the tool auto-applies the minimum-image correction
        stage(on, "graphene_hr.dat", "graphene_wsvec.dat", "graphene.uc", "graphene.xyz")
        Eon, ron = W.exact_dos(W.load_csr(run_w2s(on, "graphene", N_WANNIER, N_WANNIER, 1)),
                               eta=ETA_W, emin=EMIN_W, emax=EMAX_W)
        # OFF: identical inputs minus the _wsvec.dat -> no correction
        stage(off, "graphene_hr.dat", "graphene.uc", "graphene.xyz")
        Eoff, roff = W.exact_dos(W.load_csr(run_w2s(off, "graphene", N_WANNIER, N_WANNIER, 1)),
                                 eta=ETA_W, emin=EMIN_W, emax=EMAX_W)
    roff_i = roff                                           # same grid (shared emin/emax)

    def draw(a, lw=1.4, me_on=None, me_off=None):
        a.plot(Eon, ron, color=W.OKABE_ITO[0], ls="-", lw=lw, marker="o", markevery=me_on,
                markerfacecolor="white", markeredgecolor=W.OKABE_ITO[0], markeredgewidth=1.0,
                label="WS on")
        a.plot(Eon, roff_i, color=W.OKABE_ITO[1], ls="--", lw=lw, marker="s", markevery=me_off,
                markerfacecolor="white", markeredgecolor=W.OKABE_ITO[1], markeredgewidth=1.0,
                label="WS off")

    me = max(1, len(Eon) // 26)
    draw(ax, me_on=me, me_off=(me // 2, me))
    ax.set_xlabel(r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(r"$\rho(E)\;[\mathrm{eV}^{-1}\,\mathrm{site}^{-1}]$")
    ax.set_xlim(Eon.min(), Eon.max())
    ax.legend(loc="upper left")

    # auto-locate the 1.5 eV window of largest sustained on/off difference for the inset
    half = np.searchsorted(Eon - Eon[0], 0.75)
    d = np.convolve(np.abs(ron - roff_i), np.ones(2 * half + 1) / (2 * half + 1), mode="same")
    c = Eon[d.argmax()]
    axin = ax.inset_axes([0.50, 0.46, 0.46, 0.46])
    draw(axin, lw=1.1, me_on=None, me_off=None)
    sel = (Eon > c - 0.9) & (Eon < c + 0.9)
    axin.set_xlim(c - 0.9, c + 0.9)
    axin.set_ylim(0, 1.15 * max(ron[sel].max(), roff_i[sel].max()))
    axin.tick_params(labelsize=8)            # >= the 8 pt rule (plotting_style.md §4)
    ax.indicate_inset_zoom(axin, edgecolor="0.4")
    print("[make_min_figs] WS on/off max|dDOS|=%.4f (per site/eV) near E=%.2f eV" % (d.max(), c))


def panel_tb(ax):
    """DOS of the idealized two-band tight-binding graphene (KPM)."""
    with tempfile.TemporaryDirectory() as wd:
        stage(wd, "graphene_tb_hr.dat", "graphene_tb.uc", "graphene_tb.xyz")
        Ek, rk = W.exact_dos(W.load_csr(run_w2s(wd, "graphene_tb", N_TB, N_TB, 1)), eta=ETA_TB)
    ax.plot(Ek, rk, color=W.OKABE_ITO[2], ls="-.", marker="^",
            markevery=max(1, len(Ek) // 24), markerfacecolor="white",
            markeredgecolor=W.OKABE_ITO[2], markeredgewidth=1.0, label="exact diagonalization")
    ax.set_xlabel(r"$E\;[\mathrm{eV}]$")
    ax.set_ylabel(r"$\rho(E)\;[\mathrm{eV}^{-1}\,\mathrm{site}^{-1}]$")
    ax.set_xlim(-3.3, 3.3)
    ax.legend(loc="upper center")


def main():
    if not (shutil.which(BIN) or os.path.isfile(BIN)):
        sys.exit("wannier2sparse not found at %s; build it or set W2SP_BIN=..." % BIN)
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    usetex = bool(shutil.which("latex") and shutil.which("dvipng"))
    W._style(usetex)
    print("[make_min_figs] LaTeX path:", "usetex" if usetex else "mathtext/cm")
    os.makedirs(IMG, exist_ok=True)

    fig, ax = plt.subplots(figsize=(3.375, 2.6)); panel_wannier(ax)
    out_a = os.path.join(IMG, "graphene_wannier_dos.png")
    fig.savefig(out_a, bbox_inches="tight", pad_inches=0.02); plt.close(fig)
    open(out_a.rsplit(".", 1)[0] + ".caption.txt", "w").write(
        "FIG. N. Density of states rho(E) per site of the real DFT-derived Wannier "
        "graphene model (two Wannier functions, 149 R-points, genuine Wigner-Seitz "
        "degeneracies ndegen up to 2). Solid blue with open circles (WS on): the "
        "minimum-image correction from graphene_wsvec.dat applied, the correct "
        "minimal-Wannier result. Dashed orange with open squares (WS off): the same "
        "graphene_hr.dat with the _wsvec.dat withheld. The two curves nearly "
        "coincide; the inset zooms on the window of largest difference (near "
        "E=-3.5 eV), where they part by at most ~0.004 eV^-1 site^-1, within the "
        "finite-size ripple of the cell. For this well-localized model the "
        "Wigner-Seitz correction is thus a small perturbation on the integrated "
        "DOS, while it matters far more for k-resolved quantities (the band "
        "structure); it is auto-applied whenever _wsvec.dat is present and is part "
        "of the minimal Wannier input because the corrected operator is the one "
        "downstream code expects. Both curves are dense diagonalizations of the "
        "expanded supercell Hamiltonian, Gaussian-broadened with eta=0.06 eV, on a "
        "60x60x1 supercell (7200 states). Energies are absolute DFT eigenvalues in "
        "eV.\n")
    print("[make_min_figs] wrote", out_a)

    fig, ax = plt.subplots(figsize=(3.375, 2.6)); panel_tb(ax)
    out_b = os.path.join(IMG, "graphene_tb_dos.png")
    fig.savefig(out_b, bbox_inches="tight", pad_inches=0.02); plt.close(fig)
    open(out_b.rsplit(".", 1)[0] + ".caption.txt", "w").write(
        "FIG. N. Density of states rho(E) per site of the idealized two-band "
        "tight-binding graphene (graphene_tb_hr.dat only: nearest-neighbour t=-1, "
        "every ndegen=1, no _wsvec.dat -- the minimal tight-binding file set of "
        "tutorial 01 on the honeycomb lattice). Dash-dotted green with open "
        "triangles: dense diagonalization of the expanded supercell Hamiltonian, "
        "Gaussian-broadened with eta=0.04. rho(E) vanishes linearly at the Dirac "
        "point E=0, with van Hove peaks at E=+-|t| and band edges +-3|t|. Supercell "
        "60x60x1 (7200 states). Energies in units of |t| (eV).\n")
    print("[make_min_figs] wrote", out_b)


if __name__ == "__main__":
    main()
