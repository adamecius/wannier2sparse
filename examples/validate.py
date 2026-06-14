#!/usr/bin/env python3
"""Automated acceptance check for the tight-binding examples (no human needed).
For each model: load its HAM.CSR, assert Hermiticity, compute the KPM DOS, and
assert the spectral support / gap / Dirac dip match the analytic prediction.
Exit code 0 = all pass, 1 = any failure. Run after run.sh has produced the CSRs."""
import sys, os, numpy as np
import importlib.util as IU
_here = os.path.dirname(os.path.abspath(__file__))
spec = IU.spec_from_file_location("w2s", os.path.join(_here, "w2s_dos.py"))
w2s = IU.module_from_spec(spec); spec.loader.exec_module(w2s)

# model -> (csr_relpath, analytic_edge, kind)   kind: 'metal'|'dirac'|'gapped'
CASES = {
    "chain1d":  ("tb/chain1d/chain1d.HAM.CSR",   2.0, "metal"),
    "graphene": ("tb/graphene/graphene.HAM.CSR", 3.0, "dirac"),
    "cubic":    ("tb/cubic/cubic.HAM.CSR",       6.0, "metal"),
    "haldane":  ("tb/haldane/haldane.HAM.CSR",   3.0, "gapped"),
}
PAD = 1.05            # KPM rescaling pad used by w2s_dos.kpm_dos
EDGE_TOL = 0.15       # relative tolerance on the band edge

def support(E, rho, frac=1e-2):
    f = np.isfinite(rho); E, rho = E[f], rho[f]
    m = rho > frac*rho.max(); return E[m].min(), E[m].max()

def check(name, rel, edge, kind, root):
    path = os.path.join(root, rel)
    if not os.path.exists(path):
        print(f"  [SKIP] {name}: {rel} not found (run run.sh first)"); return None
    H = w2s.load_csr(path)
    herr = abs(H - H.getH()).max()
    assert herr < 1e-9, f"{name}: not Hermitian ({herr:.1e})"
    E, rho = w2s.kpm_dos(H, M=1024, R=8)
    lo, hi = support(E, rho)
    exp = edge*PAD
    ok_edges = abs(abs(lo)-exp) < EDGE_TOL*exp and abs(hi-exp) < EDGE_TOL*exp
    # near-zero behaviour
    mid = rho[np.abs(E) < 0.4]; midmax = mid.max() if mid.size else 0.0
    peak = rho.max()
    detail = f"support=[{lo:.2f},{hi:.2f}] exp~±{exp:.2f}  DOS(|E|<0.4)/peak={midmax/peak:.2e}"
    if kind == "gapped":
        ok = ok_edges and midmax < 0.05*peak                      # a real gap
    elif kind == "dirac":
        ok = ok_edges and 0.0 < midmax < 0.5*peak                 # dip, not a gap, not flat
    else:
        ok = ok_edges and midmax > 0.05*peak                      # states at E=0
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:9s} {detail}")
    return ok

def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_here, "models")
    results = [check(n, *v, root) for n, v in CASES.items()]
    ran = [r for r in results if r is not None]
    if not ran: print("no CSRs found; run run.sh first"); sys.exit(2)
    if all(ran): print("ALL EXAMPLES VALIDATED"); sys.exit(0)
    print("VALIDATION FAILED"); sys.exit(1)

if __name__ == "__main__":
    main()
