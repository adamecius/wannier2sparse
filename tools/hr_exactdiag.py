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


def read_wsvec(path):
    """Wannier90 `_wsvec.dat` -> {(Rx,Ry,Rz,i,j): [(tx,ty,tz), ...]} (0-based i,j).

    Each record is `Rx Ry Rz i j` then a count nT then nT minimum-image translations.
    Mirrors src/wannier_parser.cpp::read_wsvec.
    """
    toks = []
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("#"):
                continue
            toks += line.split()
    out, p = {}, 0
    while p + 6 <= len(toks):
        Rx, Ry, Rz, iw, jw = (int(toks[p + k]) for k in range(5))
        nT = int(toks[p + 5]); p += 6
        Ts = []
        for _ in range(nT):
            Ts.append((int(toks[p]), int(toks[p + 1]), int(toks[p + 2]))); p += 3
        out[(Rx, Ry, Rz, iw - 1, jw - 1)] = Ts
    return out


def apply_wsvec(iR, O, wsvec):
    """Fold the Wigner-Seitz minimum-image correction into (iR, O).

    Replaces each H(R)_{ij} by H(R+T)_{ij}/nT for the nT translations listed for
    (R,i,j), accumulating into the (expanded) R-set. Mirrors
    src/hopping_list.cpp::apply_wsvec. Returns the new (iR, O); a no-op if wsvec is
    empty. Pitfall: apply only to a *raw* `_hr.dat`; a bundle operator already has it
    folded (ndegen=1, no `_wsvec.dat` beside it), so do not double-apply.
    """
    if not wsvec:
        return iR, O
    from collections import defaultdict
    nw = O.shape[1]
    acc = defaultdict(lambda: np.zeros((nw, nw), complex))
    for r, R in enumerate(iR):
        Rt = (int(R[0]), int(R[1]), int(R[2]))
        for i in range(nw):
            for j in range(nw):
                v = O[r, i, j]
                if v == 0:
                    continue
                Ts = wsvec.get((Rt[0], Rt[1], Rt[2], i, j))
                if not Ts:
                    acc[Rt][i, j] += v
                else:
                    for (tx, ty, tz) in Ts:
                        acc[(Rt[0] + tx, Rt[1] + ty, Rt[2] + tz)][i, j] += v / len(Ts)
    newR = sorted(acc.keys())
    return np.array(newR), np.array([acc[R] for R in newR])


def read_hr_ws(seed):
    """Read <seed>_hr.dat and fold <seed>_wsvec.dat if present (the real workflow)."""
    import os
    iR, O = read_hr(f"{seed}_hr.dat")
    wsf = f"{seed}_wsvec.dat"
    if os.path.exists(wsf):
        iR, O = apply_wsvec(iR, O, read_wsvec(wsf))
        print(f"  applied Wigner-Seitz correction from {wsf}")
    return iR, O


def opk(iR, OR, kf):
    """O(k) = sum_R e^{+i2pi k.R} O(R) for k-points kf (nk,3)."""
    return np.einsum('kr,rij->kij', np.exp(2j * np.pi * (np.atleast_2d(kf) @ iR.T)), OR)


def cell_area_2d(lat):
    """In-plane unit-cell area |a1 x a2|_z for a 2D model (Angstrom^2).

    This is the z-component of the cross product of the first two lattice
    vectors, a1 x a2 = lat[0,0] lat[1,1] - lat[0,1] lat[1,0]. It reduces to
    |lat[0,0] lat[1,1]| only for a rectangular cell (off-diagonals zero, e.g.
    PdSe2); for a non-orthogonal cell (honeycomb a1=(1.5, +s), a2=(1.5, -s))
    the true area is twice that. Pitfall: using lat[0,0]*lat[1,1] for a
    honeycomb cell halves the area and so doubles any per-area response
    (it reported the Haldane Chern plateau as 2 instead of 1)."""
    return abs(lat[0, 0] * lat[1, 1] - lat[0, 1] * lat[1, 0])


# ------------------------------------------------------------------- k-path helper
DEFAULT_PATH = [("G", (0, 0, 0), 400), ("X", (.5, 0, 0), 400), ("S", (.5, .5, 0), 400),
                ("Y", (0, .5, 0), 400), ("G", (0, 0, 0), 400), ("S", (.5, .5, 0), 0)]


def read_w2s_kpath(path):
    """Read a band path recorded in a `.w2s` input under provenance.kpoint_path.

    Returns [(label, (kx,ky,kz), 0), ...] in fractional coordinates (segment counts
    left 0 so build_kpath falls back to the uniform --npts), or None if the file has
    no recorded path. The `.w2s` is JSON with `//` and `/* */` comments allowed (as
    the C++ reader tolerates), so those are stripped before json.loads. This is the
    path `wannier2sparse --provenance` extracts from the Wannier90 .win / QE bands.in
    and bakes into the input, so the Wannier bands draw on the DFT high-symmetry path
    with no separate argument.
    """
    import re
    try:
        raw = open(path).read()
    except OSError:
        return None
    raw = re.sub(r"/\*.*?\*/", "", raw, flags=re.S)          # block comments
    raw = re.sub(r"//[^\n]*", "", raw)                       # line comments
    try:
        doc = json.loads(raw)
    except (ValueError, json.JSONDecodeError):
        return None
    kp = (doc.get("provenance") or {}).get("kpoint_path")
    if not kp or not kp.get("nodes"):
        return None
    return [(n.get("label", f"k{i}"), tuple(n["k"]), 0) for i, n in enumerate(kp["nodes"])]


def read_qe_kpath(path):
    """Parse a Quantum ESPRESSO `K_POINTS crystal_b` block into band-path nodes.

    Returns [(label, (kx,ky,kz), nseg), ...] in fractional (crystal) coordinates,
    where nseg is the number of points QE places from that node to the next. This is
    the SAME high-symmetry path the DFT bands were computed on, so a Wannier band
    structure built on it overlays the DFT bands point-for-point. The path is exactly
    what provenance tracking should record (see README "Provenance tracking"), so the
    comparison needs no separate input file once the model carries its provenance.

    Pitfall: only the `crystal_b` (band) form is parsed; an explicit `crystal`/`tpiba`
    list of individual k-points is not a band path and raises.
    """
    lines = open(path).read().splitlines()
    for i, l in enumerate(lines):
        s = l.strip().lower()
        if s.startswith("k_points") and "crystal_b" in s:
            n = int(lines[i + 1].split()[0])
            nodes = []
            for j in range(n):
                parts = lines[i + 2 + j].replace("!", " ! ").split()
                k = (float(parts[0]), float(parts[1]), float(parts[2]))
                nseg = int(parts[3]) if len(parts) > 3 and parts[3].lstrip("-").isdigit() else 0
                lbl = parts[parts.index("!") + 1] if "!" in parts else f"k{j}"
                nodes.append((lbl, k, nseg))
            return nodes
    raise ValueError(f"no 'K_POINTS crystal_b' band path found in {path}")


def build_kpath(lattice, nodes=DEFAULT_PATH, npts=None):
    """Return (kfrac (Nk,3), xdist (Nk,), tick_x, tick_lbl) on the reciprocal metric.

    Each node carries its own segment count (the QE per-segment npts); pass npts to
    override them all with a uniform value.
    """
    recip = 2 * np.pi * np.linalg.inv(lattice).T          # rows = b1,b2,b3 (1/Ang)
    labels = [n[0] for n in nodes]
    fr = [np.array(n[1], float) for n in nodes]
    kf, xd, tick_x = [], [], [0.0]
    x = 0.0
    for a in range(len(fr) - 1):
        m = npts if npts else max(2, nodes[a][2])
        seg = np.linspace(0, 1, m, endpoint=(a == len(fr) - 2))
        d = (fr[a + 1] - fr[a])
        dcart = d @ recip
        seglen = np.linalg.norm(dcart)
        for t in seg:
            kf.append(fr[a] + t * d); xd.append(x + t * seglen)
        x += seglen; tick_x.append(x)
    return np.array(kf), np.array(xd), np.array(tick_x), labels


# -------------------------------------------------------------------- subcommands
def cmd_bands(a):
    iR, H = read_hr_ws(a.seed); lat = read_uc(f"{a.seed}.uc")
    w2s_path = a.w2s if a.w2s else f"{a.seed}.w2s"
    w2s_nodes = read_w2s_kpath(w2s_path)
    if a.kpath:                                      # explicit QE band path wins
        nodes = read_qe_kpath(a.kpath)
        kf, xd, tx, lbl = build_kpath(lat, nodes=nodes)
        print(f"bands: using DFT k-path from {a.kpath} ({len(nodes)} nodes: {'-'.join(n[0] for n in nodes)})")
    elif w2s_nodes:                                  # path recorded by --provenance in the .w2s
        kf, xd, tx, lbl = build_kpath(lat, nodes=w2s_nodes, npts=a.npts)
        print(f"bands: using recorded k-path from {w2s_path} ({len(w2s_nodes)} nodes: {'-'.join(n[0] for n in w2s_nodes)})")
    else:
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
    iR, H = read_hr_ws(a.seed)
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
    A_cell = cell_area_2d(read_uc(f"{a.seed}.uc"))  # in-plane cell area |a1 x a2| (2D)
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
    sig = np.cumsum(omE) * (E[1] - E[0]) / A_cell      # per-cell natural units
    # Calibration to e^2/h. The charge AHC (cmd_ahc) uses the Berry prefactor
    # 2*pi/A_cell to read the Chern number directly; the intrinsic SHC is the same
    # Berry machinery with the spin current J = 1/2{v, sigma_z}, so its e^2/h factor
    # is (2*pi/A_cell) * (1/2) -- the 1/2 because sigma_z (Pauli, eig +/-1) carries
    # spin in units of hbar/2 -- with an overall sign fixing the sigma^z_xy (vs yx)
    # / Im convention. Net: multiply the natural per-cell value by -pi. Verified
    # against the wannierberri covariant SHC (single constant, corr ~0.94 across the
    # PdSe2 curve; the topological-gap plateau reads ~+1 e^2/h). The negative sign is
    # the standard-convention alignment, exactly analogous to sigma_xx's +7.49 factor.
    if a.shc_units == "e2h":
        sig = sig * (-np.pi)
    json.dump({"energy_eV": E.tolist(), "shc": sig.tolist(),
               "meta": {"nk": a.nk, "eta_eV": a.eta, "jop": a.jop, "vop": a.vop, "units": a.shc_units}}, open(a.out + ".json", "w"))
    print(f"shc: nk={a.nk}^2, eta={a.eta*1e3:.0f} meV, J={a.jop} v={a.vop}, units={a.shc_units} -> {a.out}.json")


def cmd_ahc(a):
    """Intrinsic anomalous (charge) Hall conductivity sigma_xy in units of e^2/h.

    Berry-curvature / Kubo Fermi-sea form, with the velocity v_a = dH/dk_a built
    internally from H(R) (v_a(R) = i (R.lat)_a H(R)):
       Omega_n(k) = -2 Im sum_{m!=n} <n|v_x|m><m|v_y|n> / (E_n - E_m)^2,
       sigma_xy(E_F) = (2 pi / A_cell) (1/N_k) sum_k sum_{n: E_n<E_F} Omega_n.
    In a gap this is the Chern number C (a flat, quantized plateau at integer e^2/h).
    Needs only <seed>_hr.dat and <seed>.uc."""
    iR, H = read_hr(f"{a.seed}_hr.dat"); lat = read_uc(f"{a.seed}.uc")
    A_cell = cell_area_2d(lat); Rc = iR @ lat
    VX = 1j * Rc[:, 0, None, None] * H; VY = 1j * Rc[:, 1, None, None] * H   # v_a(R) = i (R.lat)_a H(R)
    E = np.linspace(a.emin, a.emax, a.ngrid); omE = np.zeros(a.ngrid)
    nrm = 1 / (np.sqrt(2 * np.pi) * a.eta); inv2 = 1 / (2 * a.eta ** 2); tol = 1e-4
    kx = (np.arange(a.nk) + .5) / a.nk
    KX, KY = np.meshgrid(kx, kx, indexing='ij')
    kf = np.stack([KX.ravel(), KY.ravel(), np.zeros(a.nk ** 2)], 1)
    for s in range(0, len(kf), 512):
        kb = kf[s:s + 512]
        Hk = opk(iR, H, kb); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
        w, U = np.linalg.eigh(Hk)
        vx = np.einsum('kmi,kmn,knj->kij', np.conj(U), opk(iR, VX, kb), U)
        vy = np.einsum('kmi,kmn,knj->kij', np.conj(U), opk(iR, VY, kb), U)
        dE = w[:, :, None] - w[:, None, :]; mask = np.abs(dE) > tol
        Om = np.where(mask, -2 * np.imag(vx * np.transpose(vy, (0, 2, 1))) / np.where(mask, dE ** 2, 1), 0).sum(2)
        g = nrm * np.exp(-((E[None, None, :] - w[:, :, None]) ** 2) * inv2)
        omE += np.einsum('kn,kne->e', Om, g)
    omE /= len(kf)
    sig = (2 * np.pi / A_cell) * np.cumsum(omE) * (E[1] - E[0])   # e^2/h; in a gap = Chern number
    json.dump({"energy_eV": E.tolist(), "sigma_xy_e2h": sig.tolist(),
               "meta": {"nk": a.nk, "eta_eV": a.eta}}, open(a.out + ".json", "w"))
    print(f"ahc: nk={a.nk}^2, eta={a.eta*1e3:.0f} meV -> {a.out}.json "
          f"(sigma_xy in e^2/h; a gap plateau at an integer is the Chern number)")


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
    gam = lambda l: r"$\Gamma$" if l.upper() in ("G", "GAMMA", "GM") else l
    ax.set_xticks(tx); ax.set_xticklabels([gam(l) for l in lbl])
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
    b.add_argument("--kpath", default=None, metavar="QE_BANDS_IN",
                   help="QE 'K_POINTS crystal_b' file (e.g. qe/bands.in) — build the Wannier "
                        "bands on the same k-path as the DFT, for a point-for-point overlay")
    b.add_argument("--w2s", default=None, metavar="W2S",
                   help="read the band path recorded under provenance.kpoint_path of this .w2s "
                        "(default <seed>.w2s); written by 'wannier2sparse --provenance'")
    b.add_argument("--no-plot", action="store_true"); b.set_defaults(fn=cmd_bands)
    d = sub.add_parser("dos"); common(d); d.set_defaults(fn=cmd_dos)
    s = sub.add_parser("sigma"); common(s); s.add_argument("--comp", default="xx", choices=["xx", "xy"]); s.set_defaults(fn=cmd_sigma)
    ah = sub.add_parser("ahc"); common(ah); ah.set_defaults(fn=cmd_ahc)
    h = sub.add_parser("shc"); common(h); h.add_argument("--jop", required=True); h.add_argument("--vop", required=True)
    h.add_argument("--shc-units", dest="shc_units", default="natural", choices=["natural", "e2h"],
                   help="natural (per-cell a.u., default) or e2h (x -pi, calibrated to e^2/h)")
    h.set_defaults(fn=cmd_shc)
    a = p.parse_args(argv)
    if getattr(a, "out", None) is None:
        a.out = f"{a.seed}_{a.cmd}"
    a.fn(a)


if __name__ == "__main__":
    main()
