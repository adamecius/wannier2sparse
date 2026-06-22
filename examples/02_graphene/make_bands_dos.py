#!/usr/bin/env python3
"""Tutorial 02 figure: the band structure and density of states of the real
DFT-derived Wannier graphene model, side by side on a shared energy axis.

The two panels are the same Wannier Hamiltonian seen two ways. Left: the band
structure E(k) along the high-symmetry path G-M-K-G, the two pi (p_z) bands that
touch at the Dirac point K. Right: the density of states rho(E), with the Dirac dip
lining up with the band touching and van Hove peaks at the M-point saddle. Energies
are referenced to the Fermi level E_F read from graphene.win.

Both quantities use the *Wigner-Seitz-corrected* operator, which is the default and
the only physically correct choice — the tool applies it automatically from
graphene_wsvec.dat, so there is nothing to switch on. The bands come from the
unfolded H(R) (a bundle run writes the WS-corrected, ndegen=1 operator) Fourier
transformed along the path that `wannier2sparse --provenance` records from the .win;
the DOS comes from the expanded supercell CSR. Both are drawn as solid lines.

Usage:  python3 make_bands_dos.py        # writes ../img/graphene_bands_dos.png
Set W2SP_BIN=/path/to/wannier2sparse if the binary is not at ../../build/wannier2sparse.

Pitfall: the Wannier H(R) reaches |R| = 8 in-plane, so the supercell DOS needs
N >= 2*range+1 = 17 per axis (the minimum-image guard); the band panel has no such
limit. The band route must use the WS-corrected operator (the bundle's HAM.hr.dat),
not the raw graphene_hr.dat, or the minimum-image bonds are aliased.
"""
import os, sys, re, json, shutil, subprocess, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))                 # examples/ -> w2s_dos
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "..", "tools"))  # -> hr_exactdiag
import w2s_dos as W
import hr_exactdiag as HX

BIN = os.environ.get("W2SP_BIN", os.path.join(HERE, "..", "..", "build", "wannier2sparse"))
IMG = os.path.join(HERE, "..", "img")
N_DOS = 60                     # 2*60*60 = 7200 states; >= 17 satisfies the WS guard
ETA   = 0.05                   # eV, Gaussian broadening of the supercell DOS
PATH  = [("G", (0, 0, 0)), ("M", (0.5, 0, 0)), ("K", (1/3, 1/3, 0)), ("G", (0, 0, 0))]
NSEG  = 150                    # k-points per path segment


def fermi_energy():
    for line in open(os.path.join(HERE, "graphene.win")):
        m = re.match(r"\s*fermi_energy\s*[=:]\s*([-\d.eE+]+)", line, re.I)
        if m:
            return float(m.group(1))
    return 0.0


def run_w2s(workdir, w2s_text, label):
    w2s = os.path.join(workdir, label + ".w2s")
    open(w2s, "w").write(w2s_text)
    subprocess.run([BIN, "-x", w2s], cwd=workdir, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return workdir


def compute_bands():
    """E(k) along G-M-K-G from the WS-corrected (bundle) Hamiltonian."""
    with tempfile.TemporaryDirectory() as wd:
        for f in ("graphene_hr.dat", "graphene_wsvec.dat", "graphene.uc",
                  "graphene.xyz", "graphene.win"):
            shutil.copy(os.path.join(HERE, f), os.path.join(wd, f))
        run_w2s(wd, '{ "label":"graphene", "mode":"bundle", '
                    '"provenance":{ "win":"graphene.win" } }\n', "graphene")
        ham = os.path.join(wd, "graphene.w2sp", "operators", "HAM.hr.dat")
        iR, H = HX.read_hr(ham)
    lat = HX.read_uc(os.path.join(HERE, "graphene.uc"))
    nodes = [(lbl, k, NSEG) for (lbl, k) in PATH]
    kf, xd, tx, lbl = HX.build_kpath(lat, nodes=nodes)
    Hk = HX.opk(iR, H, kf); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
    return xd, np.linalg.eigvalsh(Hk), tx, lbl


def compute_dos():
    """rho(E) per site from the expanded (WS-corrected) supercell CSR."""
    with tempfile.TemporaryDirectory() as wd:
        for f in ("graphene_hr.dat", "graphene_wsvec.dat", "graphene.uc", "graphene.xyz"):
            shutil.copy(os.path.join(HERE, f), os.path.join(wd, f))
        run_w2s(wd, '{ "label":"graphene", "mode":"sparse", '
                    '"supercell":[%d,%d,1] }\n' % (N_DOS, N_DOS), "graphene")
        return W.exact_dos(W.load_csr(os.path.join(wd, "graphene.HAM.CSR")), eta=ETA)


def main():
    if not (shutil.which(BIN) or os.path.isfile(BIN)):
        sys.exit("wannier2sparse not found at %s; build it or set W2SP_BIN=..." % BIN)
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    usetex = bool(shutil.which("latex") and shutil.which("dvipng"))
    W._style(usetex)
    print("[make_bands_dos] LaTeX path:", "usetex" if usetex else "mathtext/cm")

    EF = fermi_energy()
    xd, w, tx, lbl = compute_bands()
    Edos, rho = compute_dos()

    fig, (axb, axd) = plt.subplots(
        1, 2, sharey=True, figsize=(3.375, 2.7),
        gridspec_kw={"width_ratios": [2.3, 1], "wspace": 0.06})

    # Left: band structure (both bands, solid blue).
    for b in range(w.shape[1]):
        axb.plot(xd, w[:, b] - EF, color=W.OKABE_ITO[0], ls="-", lw=1.4)
    for x in tx:
        axb.axvline(x, color="0.85", lw=0.6)
    axb.axhline(0.0, color="0.5", lw=0.6, ls=":")
    gam = lambda s: r"$\Gamma$" if s.upper() in ("G", "GAMMA") else f"${s}$"
    axb.set_xticks(tx); axb.set_xticklabels([gam(s) for s in lbl])
    axb.set_xlim(xd[0], xd[-1]); axb.set_ylabel(r"$E - E_F\;[\mathrm{eV}]$")

    # Right: DOS (solid orange), energy on the shared vertical axis.
    axd.plot(rho, Edos - EF, color=W.OKABE_ITO[1], ls="-", lw=1.4)
    axd.axhline(0.0, color="0.5", lw=0.6, ls=":")
    axd.set_xlabel(r"$\rho\;[\mathrm{eV}^{-1}\,\mathrm{site}^{-1}]$")
    axd.set_xlim(0, None)
    axb.set_ylim(-9.8, 4.9)

    os.makedirs(IMG, exist_ok=True)
    out = os.path.join(IMG, "graphene_bands_dos.png")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02); plt.close(fig)
    open(out.rsplit(".", 1)[0] + ".caption.txt", "w").write(
        "FIG. N. Band structure and density of states of the DFT-derived Wannier "
        "graphene model (two p_z Wannier functions), referenced to the Fermi level "
        "E_F = %.4f eV. Left (solid blue): the two pi bands E(k) along G-M-K-G, "
        "touching at the Dirac point K on E_F. Right (solid orange): the density of "
        "states rho(E) per site; the Dirac dip lines up with the band touching and "
        "the van Hove peaks sit at the M-point saddles. Both use the Wigner-Seitz "
        "minimum-image correction (applied automatically): the bands from the "
        "unfolded WS-corrected H(R) Fourier-transformed along the path "
        "wannier2sparse --provenance records from the .win, the DOS from the "
        "expanded %dx%dx1 supercell CSR (7200 states) with Gaussian broadening "
        "eta=%.2f eV.\n" % (EF, N_DOS, N_DOS, ETA))
    print("[make_bands_dos] wrote", out)


if __name__ == "__main__":
    main()
