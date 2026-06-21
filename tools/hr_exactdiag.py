#!/usr/bin/env python3
"""
hr_exactdiag.py -- exact-diagonalization reference for wannier2sparse operators.

Problem solved
--------------
wannier2sparse expands a primitive Wannier90 model `O_ij(R)` into a supercell CSR
that lsquant feeds to the (linear-scaling, stochastic) KPM. This tool is the
*opposite, exact* route on the same operators: it reconstructs `O(k)` from the
`_hr.dat` files and diagonalizes `H(k)` densely on a k-path or k-mesh, giving
band structures, density of states, and Kubo conductivities with no supercell and
no stochastic noise. It is the oracle the KPM is validated against, and the
"exact" curve used in the SHC example (examples/example_9_SHC_in_PdSe2).

Conventions (docs/conventions.md)
---------------------------------
  O(k) = sum_R e^{+i 2 pi k.R} O(R) / ndegen(R)      (k fractional, R integer)
  velocity operators v_a(R) are stored in eV*Angstrom (= dH/dk_a, i.e. hbar*v)
  spin operators in hbar/2 (Pauli, eigenvalues +-1)

Subcommands
-----------
  bands <seed>            band structure along a k-path (default 2D Gamma-X-S-Y-Gamma)
  dos   <seed>            density of states (Gaussian broadening eta)
  sigma <seed> [--comp xx|xy]   Kubo-Greenwood conductivity (e^2/h), needs v_a
  shc   <seed> --jop J --vop V  intrinsic Kubo (Fermi-sea) spin Hall (a.u.), needs J,v

Inputs: <seed>_hr.dat (H), <seed>.uc (lattice, for k-distance/units), and the
operator files <seed>_<op>_hr.dat (e.g. _vx_hr.dat, _vy_hr.dat, _Sz_hr.dat,
_JXSZ_hr.dat). A tight-binding user who supplies these _hr.dat files directly
gets every quantity below with no DFT.

Pitfall: this tool does the PLAIN Fourier transform; if the model still needs the
Wigner-Seitz `_wsvec.dat` correction it must be folded into `_hr.dat` first (see
docs/conventions.md). Re-applying a stray `_wsvec.dat` double-counts it.
"""
import argparse, sys, json
import numpy as np


# ----------------------------------------------------------------------------- I/O
def read_hr(path):
    """Wannier90 _hr.dat -> (iR (nR,3) int, O (nR,nw,nw) complex), ndegen folded."""
    L = open(path).read().split('\n')
    nw = int(L[1]); nr = int(L[2])
    i = 3; deg = []
    while len(deg) < nr:
        deg += [int(x) for x in L[i].split()]; i += 1
    deg = np.array(deg, float)
    rec = {}; iR = []; O = np.zeros((nr, nw, nw), complex)
    for line in L[i:]:
        s = line.split()
        if len(s) < 7:
            continue
        R = (int(s[0]), int(s[1]), int(s[2]))
        if R not in rec:
            rec[R] = len(rec); iR.append(R)
        r = rec[R]
        O[r, int(s[3]) - 1, int(s[4]) - 1] = (float(s[5]) + 1j * float(s[6])) / deg[r]
    return np.array(iR), O


def read_uc(path):
    """<seed>.uc -> 3x3 lattice (Angstrom). Falls back to identity if absent."""
    try:
        rows = [list(map(float, l.split())) for l in open(path) if l.split()]
        return np.array(rows[:3])
    except OSError:
        return np.eye(3)


def opk(iR, OR, kf):
    """O(k) = sum_R e^{+i2pi k.R} O(R) for k-points kf (nk,3)."""
    return np.einsum('kr,rij->kij', np.exp(2j * np.pi * (np.atleast_2d(kf) @ iR.T)), OR)


# ------------------------------------------------------------------- k-path helper
DEFAULT_PATH = [("G", (0, 0, 0)), ("X", (.5, 0, 0)), ("S", (.5, .5, 0)),
                ("Y", (0, .5, 0)), ("G", (0, 0, 0)), ("S", (.5, .5, 0))]


def build_kpath(lattice, nodes=DEFAULT_PATH, npts=400):
    """Return (kfrac (Nk,3), xdist (Nk,), tick_x, tick_lbl) on the reciprocal metric."""
    recip = 2 * np.pi * np.linalg.inv(lattice).T          # rows = b1,b2,b3 (1/Ang)
    labels = [n[0] for n in nodes]
    fr = [np.array(n[1], float) for n in nodes]
    kf, xd, tick_x = [], [], [0.0]
    x = 0.0
    for a in range(len(fr) - 1):
        seg = np.linspace(0, 1, npts, endpoint=(a == len(fr) - 2))
        d = (fr[a + 1] - fr[a])
        dcart = d @ recip
        seglen = np.linalg.norm(dcart)
        for t in seg:
            kf.append(fr[a] + t * d); xd.append(x + t * seglen)
        x += seglen; tick_x.append(x)
    return np.array(kf), np.array(xd), np.array(tick_x), labels


# -------------------------------------------------------------------- subcommands
def cmd_bands(a):
    iR, H = read_hr(f"{a.seed}_hr.dat"); lat = read_uc(f"{a.seed}.uc")
    kf, xd, tx, lbl = build_kpath(lat, npts=a.npts)
    Hk = opk(iR, H, kf); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
    w = np.linalg.eigvalsh(Hk)                            # (Nk, nw)
    out = {"kdist": xd.tolist(), "bands": w.tolist(), "tick_x": tx.tolist(),
           "tick_label": lbl, "EF": a.ef}
    json.dump(out, open(a.out + ".json", "w"))
    if not a.no_plot:
        _plot_bands(xd, w, tx, lbl, a.ef, a.out + ".png", a.seed)
    print(f"bands: {w.shape[1]} bands x {len(xd)} k -> {a.out}.json"
          + ("" if a.no_plot else f", {a.out}.png"))


def cmd_dos(a):
    iR, H = read_hr(f"{a.seed}_hr.dat")
    E, dos, _ = _spectral(iR, H, a.nk, a.eta, a.ngrid, a.emin, a.emax)
    json.dump({"energy_eV": E.tolist(), "dos": dos.tolist(),
               "meta": {"nk": a.nk, "eta_eV": a.eta, "integral": float(np.trapezoid(dos, E))}},
              open(a.out + ".json", "w"))
    print(f"dos: nk={a.nk}^2, eta={a.eta*1e3:.0f} meV, integral={np.trapezoid(dos,E):.3f} (~num_wann)")


def cmd_sigma(a):
    iR, H = read_hr(f"{a.seed}_hr.dat"); lat = read_uc(f"{a.seed}.uc")
    A_cell = abs(lat[0, 0] * lat[1, 1])                   # 2D cell area (Ang^2)
    iRa, va = read_hr(f"{a.seed}_v{a.comp[0]}_hr.dat")
    iRb, vb = read_hr(f"{a.seed}_v{a.comp[1]}_hr.dat")
    E = np.linspace(a.emin, a.emax, a.ngrid); sig = np.zeros(a.ngrid)
    nrm = 1 / (np.sqrt(2 * np.pi) * a.eta); inv2 = 1 / (2 * a.eta ** 2)
    kx = (np.arange(a.nk) + .5) / a.nk
    KX, KY = np.meshgrid(kx, kx, indexing='ij')
    kf = np.stack([KX.ravel(), KY.ravel(), np.zeros(a.nk ** 2)], 1)
    for s in range(0, len(kf), 1024):
        kb = kf[s:s + 1024]
        Hk = opk(iR, H, kb); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
        w, U = np.linalg.eigh(Hk)
        Va = np.einsum('kmi,kmn,knj->kij', np.conj(U), opk(iRa, va, kb), U)
        Vb = np.einsum('kmi,kmn,knj->kij', np.conj(U), opk(iRb, vb, kb), U)
        g = nrm * np.exp(-((E[None, None, :] - w[:, :, None]) ** 2) * inv2)
        if a.comp == "xx":
            sig += np.einsum('knm,kne,kme->e', np.abs(Va) ** 2, g, g)
        else:                                              # xy: Re[Va_nm Vb_mn]
            sig += np.einsum('knm,kne,kme->e', np.real(Va * np.transpose(Vb, (0, 2, 1))), g, g)
    sig *= (2 * np.pi ** 2 / A_cell) / len(kf)
    json.dump({"energy_eV": E.tolist(), f"sigma_{a.comp}_e2h": sig.tolist()}, open(a.out + ".json", "w"))
    print(f"sigma_{a.comp}: nk={a.nk}^2, eta={a.eta*1e3:.0f} meV -> {a.out}.json")


def cmd_shc(a):
    """Intrinsic (Fermi-sea / spin-Berry) SHC: sigma(E_F)=cumulative over occupied of
       Omega_n(k) = sum_{m!=n} 2 Im[<n|J|m><m|v|n>]/(E_n-E_m)^2."""
    iRH, H = read_hr(f"{a.seed}_hr.dat")
    iRJ, J = read_hr(a.jop); iRV, V = read_hr(a.vop)   # each operator keeps its own R-grid
    A_cell = abs(read_uc(f"{a.seed}.uc")[0, 0] * read_uc(f"{a.seed}.uc")[1, 1])  # per-cell area (2D)
    E = np.linspace(a.emin, a.emax, a.ngrid); omE = np.zeros(a.ngrid)
    nrm = 1 / (np.sqrt(2 * np.pi) * a.eta); inv2 = 1 / (2 * a.eta ** 2); tol = 1e-4
    kx = (np.arange(a.nk) + .5) / a.nk
    KX, KY = np.meshgrid(kx, kx, indexing='ij')
    kf = np.stack([KX.ravel(), KY.ravel(), np.zeros(a.nk ** 2)], 1)
    for s in range(0, len(kf), 512):
        kb = kf[s:s + 512]
        Hk = opk(iRH, H, kb); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
        w, U = np.linalg.eigh(Hk)
        Je = np.einsum('kmi,kmn,knj->kij', np.conj(U), opk(iRJ, J, kb), U)
        Ve = np.einsum('kmi,kmn,knj->kij', np.conj(U), opk(iRV, V, kb), U)
        dE = w[:, :, None] - w[:, None, :]; mask = np.abs(dE) > tol
        term = np.where(mask, 2 * np.imag(Je * np.transpose(Ve, (0, 2, 1))) / np.where(mask, dE ** 2, 1), 0)
        Om = term.sum(2)                                   # Omega_n(k)  (nk,nw)
        g = nrm * np.exp(-((E[None, None, :] - w[:, :, None]) ** 2) * inv2)   # (nk,nw,ngrid)
        omE += np.einsum('kn,kne->e', Om, g)
    omE /= len(kf)
    sig = np.cumsum(omE) * (E[1] - E[0]) / A_cell      # per-cell (matches sigma + the exact reference)
    json.dump({"energy_eV": E.tolist(), "shc": sig.tolist(),
               "meta": {"nk": a.nk, "eta_eV": a.eta, "jop": a.jop, "vop": a.vop}}, open(a.out + ".json", "w"))
    print(f"shc: nk={a.nk}^2, eta={a.eta*1e3:.0f} meV, J={a.jop} v={a.vop} -> {a.out}.json")


# ---------------------------------------------------------------------- internals
def _spectral(iR, H, nk, eta, ngrid, emin, emax):
    E = np.linspace(emin, emax, ngrid); dos = np.zeros(ngrid)
    nrm = 1 / (np.sqrt(2 * np.pi) * eta); inv2 = 1 / (2 * eta ** 2)
    kx = (np.arange(nk) + .5) / nk
    KX, KY = np.meshgrid(kx, kx, indexing='ij')
    kf = np.stack([KX.ravel(), KY.ravel(), np.zeros(nk ** 2)], 1)
    allw = []
    for s in range(0, len(kf), 2048):
        Hk = opk(iR, H, kf[s:s + 2048]); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
        w = np.linalg.eigvalsh(Hk); allw.append(w)
        dos += (nrm * np.exp(-((E[None, None, :] - w[:, :, None]) ** 2) * inv2)).sum((0, 1))
    return E, dos / len(kf), np.concatenate(allw)


def _plot_bands(xd, w, tx, lbl, ef, out, title):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    for b in range(w.shape[1]):
        ax.plot(xd, w[:, b] - ef, color="#0072B2", lw=0.9)
    for x in tx:
        ax.axvline(x, color="0.85", lw=0.6)
    ax.axhline(0, color="0.6", lw=0.6, ls=":")
    ax.set_xticks(tx); ax.set_xticklabels([r"$\Gamma$" if l == "G" else l for l in lbl])
    ax.set_xlim(xd[0], xd[-1]); ax.set_ylabel(r"$E-E_F$ (eV)")
    fig.tight_layout(); fig.savefig(out, dpi=200)


def main(argv=None):
    p = argparse.ArgumentParser(description="exact-diagonalization reference for wannier2sparse _hr.dat operators")
    sub = p.add_subparsers(dest="cmd", required=True)
    def common(q):
        q.add_argument("seed"); q.add_argument("--out", default=None)
        q.add_argument("--nk", type=int, default=200); q.add_argument("--eta", type=float, default=0.02)
        q.add_argument("--ngrid", type=int, default=600)
        q.add_argument("--emin", type=float, default=-3.5); q.add_argument("--emax", type=float, default=1.5)
    b = sub.add_parser("bands"); b.add_argument("seed"); b.add_argument("--out", default=None)
    b.add_argument("--npts", type=int, default=400); b.add_argument("--ef", type=float, default=0.0)
    b.add_argument("--no-plot", action="store_true"); b.set_defaults(fn=cmd_bands)
    d = sub.add_parser("dos"); common(d); d.set_defaults(fn=cmd_dos)
    s = sub.add_parser("sigma"); common(s); s.add_argument("--comp", default="xx", choices=["xx", "xy"]); s.set_defaults(fn=cmd_sigma)
    h = sub.add_parser("shc"); common(h); h.add_argument("--jop", required=True); h.add_argument("--vop", required=True); h.set_defaults(fn=cmd_shc)
    a = p.parse_args(argv)
    if getattr(a, "out", None) is None:
        a.out = f"{a.seed}_{a.cmd}"
    a.fn(a)


if __name__ == "__main__":
    main()
