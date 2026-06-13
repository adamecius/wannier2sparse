#!/usr/bin/env python3
"""Plan 11 -- generate WannierBerri committed goldens for the cross-check tests.

DECOUPLED COMMITTED-GOLDEN architecture (docs/conventions.md sec 7): this script
is run ONCE, manually, by someone with WannierBerri installed. It writes two
reference files per (fixture, operator) that are versioned in the repo. The C++
tests (test/wberri_matrix_crosscheck.cpp, test/wberri_texture_crosscheck.cpp)
only READ those .ref files and never import WannierBerri.

Usage:
    pip install wannierberri
    # the fixture must include <seed>.chk (collect it in gen_fixture.sh) plus the
    # usual _hr.dat / _u.mat / _u_dis.mat / .eig and .spn (spin) / .amn (orbital).
    python3 wberri_reference.py <prefix> <spin|orbital> <out_dir>
    #   -> writes <out_dir>/<seed>_<S|L>_matrix.ref  and  <seed>_<S|L>_texture.ref

This is normally invoked via `make regen-wberri-golden`.

------------------------------------------------------------------------------
.ref formats (must stay in lockstep with the two C++ tests)

  <seed>_<op>_matrix.ref   (Test 1 -- Wannier-gauge matrix element)
    line 1 (header):  seed operator n_k num_wann WS_convention wberri_version
    then one row per (ik, m, n, alpha):  ik m n alpha Re Im
      m,n are 0-based Wannier indices; alpha in {0,1,2} = x,y,z.

  <seed>_<op>_texture.ref  (Test 2 -- band-resolved expectation)
    line 1 (header):  seed operator n_k num_wann WS_convention wberri_version degeneracy_tol
    then one row per (ik, ibnd, alpha):  ik ibnd alpha block value
      block = degenerate-block id within (ik): bands with |E_i-E_j|<degeneracy_tol
      share a block id so the comparator sums <O_alpha> over the block (subspace
      trace -- band-by-band is ill-defined inside a multiplet, e.g. Fe+SOC).

------------------------------------------------------------------------------
TWO things that are deliberately NOT assumed here (flag, do not fabricate):

1. WS CONVENTION (docs/conventions.md sec 6/7). WannierBerri and wannier2sparse
   must use the SAME use_ws_distance convention or O_W(k) differs by trivial
   per-orbital phases. This script writes the convention token it used into every
   header; the C++ tests assert it matches their own. The matching value is a
   build-time fact of how <seed>_hr.dat was made -- set W2SP_WS_CONVENTION to
   match and VERIFY against the installed WannierBerri, not from memory.

2. The VERSION-SPECIFIC real-space-matrix attribute names on the WannierBerri
   System object (Ham_R / SS_R / iRvec / degeneracies). These have changed across
   releases. We introspect and raise a clear, actionable error if the expected
   attribute is absent, rather than emit unverified numbers.
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


def _real_space_matrices(system, which):
    """Return (iRvec, weights, Ham_R, Op_R) from a WannierBerri System, introspecting
    the version-specific attribute names. Op_R is the spin (SS_R) or orbital matrix
    in the Wannier gauge, shape (num_wann, num_wann, nRvec, 3). Raises with guidance
    if an expected attribute is missing -- we never guess silently."""
    import numpy as np

    def first_attr(obj, names, label):
        for n in names:
            if hasattr(obj, n):
                return getattr(obj, n)
            getter = "get_R_mat"
            if hasattr(obj, getter):
                try:
                    return getattr(obj, getter)(n)
                except Exception:
                    pass
        raise SystemExit(
            "wberri_reference: could not find the '%s' real-space matrix on this "
            "WannierBerri System (tried %s). The attribute/API is version-specific; "
            "confirm it against the installed WannierBerri source and update "
            "_real_space_matrices(). Do NOT guess -- a wrong matrix yields a "
            "wrong-but-plausible golden." % (label, names))

    iRvec = np.asarray(first_attr(system, ["iRvec"], "iRvec"))
    Ham_R = np.asarray(first_attr(system, ["Ham_R", "HH_R"], "Ham_R"))

    # WannierBerri may fold the WS degeneracy into the stored matrices already, or
    # expose it separately. If a degeneracy vector exists, divide it out so our
    # explicit FT below does not double-count.
    weights = np.ones(iRvec.shape[0])
    for n in ("Ndegen", "ndegen", "_NKFFT_recommended"):
        if hasattr(system, n):
            cand = np.asarray(getattr(system, n))
            if cand.ndim == 1 and cand.shape[0] == iRvec.shape[0]:
                weights = 1.0 / cand
                break

    if which == "spin":
        Op_R = np.asarray(first_attr(system, ["SS_R"], "SS_R (spin)"))
    else:
        Op_R = np.asarray(first_attr(
            system, ["OAM_R", "AA_R"], "orbital-L (Wannier gauge)"))
        raise SystemExit(
            "wberri_reference: orbital-L in the Wannier gauge has no guaranteed "
            "one-to-one with the wannier2sparse projector route (C(k)=A(k)^dag V(k), "
            "L_local). Confirm the WannierBerri orbital operator matches that "
            "definition before trusting an L golden (docs/conventions.md sec 7, "
            "Deferred). Refusing to emit an unverified orbital golden.")

    # Normalize Ham_R/Op_R orientation to (nw, nw, nR[, 3]).
    nw = Ham_R.shape[0] if Ham_R.shape[0] == Ham_R.shape[1] else Ham_R.shape[1]
    if Ham_R.shape[0] != nw:                            # stored as (nR, nw, nw)
        Ham_R = np.moveaxis(Ham_R, 0, -1)
    if Op_R.shape[0] != nw:
        Op_R = np.moveaxis(Op_R, 0, -3 if Op_R.ndim == 4 else -1)
    return iRvec, weights, Ham_R, Op_R, nw


def _ft(kpt, iRvec, weights, mat):
    """O(k) = sum_R w_R e^{+i 2pi k.R} O(R)  (inverse interpolation, docs sec 1)."""
    import numpy as np
    phase = np.exp(2j * np.pi * (iRvec @ np.asarray(kpt)))   # (nR,)
    wphase = weights * phase
    # mat: (nw, nw, nR[, 3]) -> contract over R.
    return np.tensordot(mat, wphase, axes=([2], [0]))


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

    # System from the W90 .chk (+ .eig, and .spn for spin). Match the WS convention
    # used to build <seed>_hr.dat -- verify the keyword against the installed version.
    sys_w90 = wb.system.System_w90(prefix, spin=(which == "spin"),
                                   transl_inv=False)
    kpts, nw_umat = _kpts_from_umat(prefix)
    iRvec, weights, Ham_R, Op_R, nw = _real_space_matrices(sys_w90, which)
    if nw != nw_umat:
        sys.exit("num_wann mismatch: u.mat says %d, System says %d" % (nw_umat, nw))
    nk = kpts.shape[0]
    os.makedirs(out_dir, exist_ok=True)

    # ---- Test 1 golden: Wannier-gauge matrix element O_W(k) ------------------
    mpath = os.path.join(out_dir, "%s_%s_matrix.ref" % (seed, op_tag))
    with open(mpath, "w") as f:
        f.write("%s %s %d %d %s %s\n" % (seed, op_tag, nk, nw, ws_conv, wb_version))
        for ik in range(nk):
            Ok = np.stack([_ft(kpts[ik], iRvec, weights, Op_R[..., a])
                           for a in range(3)], axis=0)        # (3, nw, nw)
            for a in range(3):
                for m in range(nw):
                    for n in range(nw):
                        v = Ok[a, m, n]
                        f.write("%d %d %d %d %.16e %.16e\n"
                                % (ik, m, n, a, v.real, v.imag))

    # ---- Test 2 golden: band texture <O_alpha>_{nk} with degeneracy blocks ----
    deg_tol = float(os.environ.get("W2SP_DEGENERACY_TOL", "1e-5"))
    tpath = os.path.join(out_dir, "%s_%s_texture.ref" % (seed, op_tag))
    with open(tpath, "w") as f:
        f.write("%s %s %d %d %s %s %.3e\n"
                % (seed, op_tag, nk, nw, ws_conv, wb_version, deg_tol))
        for ik in range(nk):
            Hk = _ft(kpts[ik], iRvec, weights, Ham_R)
            Hk = 0.5 * (Hk + Hk.conj().T)
            evals, evecs = np.linalg.eigh(Hk)
            Ok = [_ft(kpts[ik], iRvec, weights, Op_R[..., a]) for a in range(3)]
            block, b = [0] * nw, 0
            for i in range(1, nw):
                if evals[i] - evals[i - 1] >= deg_tol:
                    b += 1
                block[i] = b
            for ibnd in range(nw):
                psi = evecs[:, ibnd]
                for a in range(3):
                    val = np.vdot(psi, Ok[a] @ psi).real
                    f.write("%d %d %d %d %.16e\n" % (ik, ibnd, a, block[ibnd], val))

    sys.stderr.write("wrote %s\nwrote %s\n" % (mpath, tpath))


if __name__ == "__main__":
    main()
