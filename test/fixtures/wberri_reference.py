#!/usr/bin/env python3
"""Plan 11 -- generate WannierBerri committed goldens for the cross-check tests.

DECOUPLED COMMITTED-GOLDEN architecture (docs/conventions.md sec 7): this script
is run ONCE, manually, by someone with WannierBerri installed. It writes two
reference files per (fixture, operator) that are versioned in the repo. The C++
tests (test/wberri_matrix_crosscheck.cpp, test/wberri_texture_crosscheck.cpp)
only READ those .ref files and never import WannierBerri.

Confirmed against WannierBerri 26.4.6 (v26 pipeline: Wannier90data.from_w90_files
-> System_w90; see docs/conventions.md sec 7 for the exact calls).

Usage:
    pip install wannierberri fortio      # fortio reads the binary .chk
    # <prefix> must point at a dir holding <seed>.chk + <seed>.mmn + .eig + .spn
    # (+ _u.mat for the k-points). The .chk/.mmn live in the gen working dir
    # (e.g. /tmp/fixrun.*), not the slim /tmp/fix/<seed>. The .spn here is FORMATTED.
    W2SP_WS_CONVENTION=II python3 wberri_reference.py <prefix> spin <out_dir>
    #   -> writes <out_dir>/<seed>_S_matrix.ref  and  <seed>_S_texture.ref
    # Optional: W2SP_MATRIX_KSTRIDE (default 8) subsets the matrix golden's k;
    #           W2SP_DEGENERACY_TOL (default 1e-5; use 1e-3 -- see sec 7) for texture.

This is normally invoked via `make regen-wberri-golden`. Orbital ('L') is NOT
emitted (the WannierBerri<->projector-route correspondence is unconfirmed; the
script raises rather than fabricate -- sec 7 Deferred).

------------------------------------------------------------------------------
.ref formats (must stay in lockstep with the two C++ tests)

  <seed>_<op>_matrix.ref   (Test 1 -- Wannier-gauge matrix element)
    line 1 (header):  seed operator n_k num_wann WS_convention wberri_version
    then one row per (ik, m, n, alpha):  ik m n alpha Re Im
      m,n are 0-based Wannier indices; alpha in {0,1,2} = x,y,z. n_k is the FULL
      grid size; the rows may be a stride subset of k (ik says which).

  <seed>_<op>_texture.ref  (Test 2 -- band-resolved expectation)
    line 1 (header):  seed operator n_k num_wann WS_convention wberri_version degeneracy_tol
    then one row per (ik, ibnd, alpha):  ik ibnd alpha block value
      block = degenerate-block id within (ik): bands with |E_i-E_j|<degeneracy_tol
      share a block id so the comparator sums <O_alpha> over the block (subspace
      trace -- band-by-band is ill-defined inside a multiplet, e.g. Fe+SOC).

------------------------------------------------------------------------------
WS CONVENTION (docs/conventions.md sec 6/7): WannierBerri and wannier2sparse must
use the SAME use_ws_distance convention or O_W(k) differs by trivial per-orbital
phases. This script writes W2SP_WS_CONVENTION into every header; the C++ tests
assert it matches their own. Empirically (Fe) the plain inverse FT of WannierBerri's
Ham_R reproduces the w2s H(k) eigenvalues to ~4e-5 (the _hr.dat text precision),
confirming the conventions agree -- but the token is still a build-time fact: set
it to match how <seed>_hr.dat was made.
"""
import os
import sys


def _kpts_from_umat(prefix):
    """Read the mp_grid k-points (fractional kpt_latt) from <prefix>_u.mat so the
    golden lives on EXACTLY the k-points the C++ side reconstructs on."""
    import numpy as np
    path = prefix + "_u.mat"
    with open(path) as f:
        f.readline()                                   # comment
        nk, nw, _ = (int(x) for x in f.readline().split())
        kpts = []
        for _ik in range(nk):
            line = f.readline()
            while line.strip() == "":                  # skip the blank separator
                line = f.readline()
            kpts.append([float(x) for x in line.split()[:3]])
            for _ in range(nw * nw):                    # skip the matrix block
                f.readline()
    return np.array(kpts), nw


def _build_system(prefix, which):
    """Build a WannierBerri System from <prefix>.{chk,eig,spn} via the v26 pipeline
    (Wannier90data.from_w90_files -> System_w90). Returns (iRvec, Ham_R, Op_R, nw)
    with Ham_R shaped (nR, nw, nw) and Op_R (nR, nw, nw, 3); Op_R is the spin SS_R.
    Verified end-to-end against wannierberri 26.4.6 + fortio 0.4 (.chk reader) on
    the Fe SOC fixture: the plain inverse FT below reproduces the w2s H(k)
    eigenvalues to ~4e-5 (the _hr.dat text precision), i.e. same R-set / ndegen /
    convention. .chk is binary (needs `pip install fortio`); .spn here is FORMATTED
    (the fixture was written with spn_formatted=.true.), so it is flagged as such.

    Orbital L is deliberately NOT emitted: WannierBerri's orbital operator has no
    guaranteed one-to-one with the w2s projector route (C(k)=A(k)^dag V(k),
    L_local) -- raises rather than emit an unverified golden (conventions sec 7)."""
    import numpy as np
    import wannierberri as wb
    import wannierberri.w90files as w90f

    if which != "spin":
        raise SystemExit(
            "wberri_reference: orbital-L in the Wannier gauge has no guaranteed "
            "one-to-one with the wannier2sparse projector route (C(k)=A(k)^dag V(k), "
            "L_local). Confirm the WannierBerri orbital operator matches that "
            "definition before trusting an L golden (docs/conventions.md sec 7, "
            "Deferred). Refusing to emit an unverified orbital golden.")

    try:
        wandata = w90f.Wannier90data.from_w90_files(
            prefix, files=("chk", "eig", "spn"), formatted=("spn",))
    except ModuleNotFoundError as e:
        if "fortio" in str(e):
            raise SystemExit("wberri_reference: reading the binary .chk needs fortio "
                             "-> pip install fortio")
        raise
    # spin=True makes System_w90 build SS_R; symmetrize=False (no symmetrizer here).
    system = wb.system.System_w90(wandata, symmetrize=False, spin=True)
    Ham_R = np.asarray(system.get_R_mat("Ham"))         # (nR, nw, nw)
    Op_R = np.asarray(system.get_R_mat("SS"))           # (nR, nw, nw, 3)
    iRvec = np.asarray(system.rvec.iRvec)               # (nR, 3)
    nw = int(system.num_wann)
    if Ham_R.shape != (iRvec.shape[0], nw, nw):
        raise SystemExit("wberri_reference: unexpected Ham_R shape %s (expected "
                         "(%d,%d,%d)); the API may differ in this WannierBerri "
                         "version -- confirm before trusting the golden."
                         % (Ham_R.shape, iRvec.shape[0], nw, nw))
    return iRvec, Ham_R, Op_R, nw


def _ft(kpt, iRvec, mat):
    """O(k) = sum_R e^{+i 2pi k.R} O(R)  (inverse interpolation, docs sec 1).

    No 1/ndegen weight: WannierBerri folds the WS treatment into its stored
    real-space matrices, so the plain sum is the correct inverse (confirmed by the
    H(k) eigenvalue match above). mat is (nR, nw, nw) or (nR, nw, nw, 3)."""
    import numpy as np
    phase = np.exp(2j * np.pi * (iRvec @ np.asarray(kpt)))   # (nR,)
    return np.tensordot(phase, mat, axes=([0], [0]))


def main():
    if len(sys.argv) < 4:
        sys.exit("usage: wberri_reference.py <prefix> <spin|orbital> <out_dir>")
    prefix, which, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    if which not in ("spin", "orbital"):
        sys.exit("operator must be 'spin' or 'orbital'")
    op_tag = "S" if which == "spin" else "L"
    seed = os.path.basename(prefix)
    ws_conv = os.environ.get("W2SP_WS_CONVENTION", "")
    if not ws_conv:
        sys.exit("set W2SP_WS_CONVENTION to the use_ws_distance convention used to "
                 "build <seed>_hr.dat (must match the C++ test); see conventions sec 7")

    try:
        import numpy as np
        import wannierberri as wb
    except ImportError:
        sys.exit("wannierberri / numpy not installed: pip install wannierberri")
    wb_version = getattr(wb, "__version__", "unknown")

    iRvec, Ham_R, Op_R, nw = _build_system(prefix, which)
    kpts, nw_umat = _kpts_from_umat(prefix)
    if nw != nw_umat:
        sys.exit("num_wann mismatch: u.mat says %d, System says %d" % (nw_umat, nw))
    nk = kpts.shape[0]
    os.makedirs(out_dir, exist_ok=True)

    # ---- Test 1 golden: Wannier-gauge matrix element O_W(k) ------------------
    # The full mp_grid matrix is large (nk*3*nw^2 rows); store a BZ-spread stride
    # subset (W2SP_MATRIX_KSTRIDE, default 8) -- element-by-element on a sample of k
    # is already the strict check. The k-index is the first column, so the C++ test
    # reconstructs exactly the k present. n_k in the header is the FULL grid size.
    kstride = int(os.environ.get("W2SP_MATRIX_KSTRIDE", "8"))
    mpath = os.path.join(out_dir, "%s_%s_matrix.ref" % (seed, op_tag))
    with open(mpath, "w") as f:
        f.write("%s %s %d %d %s %s\n" % (seed, op_tag, nk, nw, ws_conv, wb_version))
        for ik in range(0, nk, kstride):
            Ok = _ft(kpts[ik], iRvec, Op_R)                  # (nw, nw, 3)
            for a in range(3):
                for m in range(nw):
                    for n in range(nw):
                        v = Ok[m, n, a]
                        f.write("%d %d %d %d %.16e %.16e\n"
                                % (ik, m, n, a, v.real, v.imag))

    # ---- Test 2 golden: band texture <O_alpha>_{nk} with degeneracy blocks ----
    deg_tol = float(os.environ.get("W2SP_DEGENERACY_TOL", "1e-5"))
    tpath = os.path.join(out_dir, "%s_%s_texture.ref" % (seed, op_tag))
    with open(tpath, "w") as f:
        f.write("%s %s %d %d %s %s %.3e\n"
                % (seed, op_tag, nk, nw, ws_conv, wb_version, deg_tol))
        for ik in range(nk):
            Hk = _ft(kpts[ik], iRvec, Ham_R)
            Hk = 0.5 * (Hk + Hk.conj().T)
            evals, evecs = np.linalg.eigh(Hk)
            Ok = _ft(kpts[ik], iRvec, Op_R)                  # (nw, nw, 3)
            block, b = [0] * nw, 0
            for i in range(1, nw):
                if evals[i] - evals[i - 1] >= deg_tol:
                    b += 1
                block[i] = b
            for ibnd in range(nw):
                psi = evecs[:, ibnd]
                for a in range(3):
                    val = np.vdot(psi, Ok[:, :, a] @ psi).real
                    f.write("%d %d %d %d %.16e\n" % (ik, ibnd, a, block[ibnd], val))

    sys.stderr.write("wrote %s\nwrote %s\n" % (mpath, tpath))


if __name__ == "__main__":
    main()
