#!/usr/bin/env python3
"""Plan 10C (independent path) -- generate a reference O_W(k) with WannierBerri.

This is the *independent implementation* cross-check (beyond the self-contained
round-trip test, which already gates the FT/gauge conventions). WannierBerri reads
the same Wannier90 files and fixes the same gauge, so a match confirms the .spn
packing / real-harmonic ordering against a second codebase.

Status: scaffolding. NOT run in CI. To use:
    pip install wannierberri
    # regenerate the fixture so it includes <seed>.chk (and .mmn/.eig/.spn):
    #   wannier90.x writes <seed>.chk; ensure gen_fixture collects it.
    python3 wberri_reference.py <prefix> <operator: spin|orbital>  > ref_Ok.dat

ref_Ok.dat layout (read by test/wberri_crosscheck.cpp):
    line 1:  num_wann  num_kpts  ncomp(=3)
    then per (k, component, i, j):  Re  Im        (i fastest, then j, then comp, then k)

IMPORTANT -- WS convention trap: WannierBerri and wannier2sparse must use the
SAME minimum-image (use_ws_distance) convention, or O_W(k) differs by trivial
phases. Build the System with the convention matching how <seed>_hr.dat was made
(here: use_ws_distance applied -> set transl_inv / ws_dist accordingly; verify
against the WannierBerri version in use). See docs/conventions.md sec 6.
"""
import sys

def main():
    if len(sys.argv) < 3:
        sys.exit("usage: wberri_reference.py <prefix> <spin|orbital>")
    prefix, which = sys.argv[1], sys.argv[2]
    try:
        import numpy as np
        import wannierberri as wb
    except ImportError:
        sys.exit("wannierberri / numpy not installed: pip install wannierberri")

    # System from the W90 .chk (+ .eig, and .spn for spin). berry=False; we only
    # need the gauge matrices to reconstruct O_W(k). Match the WS convention used
    # to build <seed>_hr.dat.
    sys_w90 = wb.system.System_w90(prefix, spin=(which == "spin"),
                                   transl_inv=False)   # <-- verify vs how hr.dat was made
    # NOTE: extracting O_W(k) on the mp_grid from a given WannierBerri version may
    # require its internal API (Ham_R / SS_R Fourier-summed to k). Left explicit
    # here because the exact call path is version-dependent and must be confirmed
    # against the installed WannierBerri source, not assumed.
    raise SystemExit("wberri_reference: fill in the version-specific O_W(k) export "
                     "for the installed WannierBerri, then emit ref_Ok.dat")

if __name__ == "__main__":
    main()
