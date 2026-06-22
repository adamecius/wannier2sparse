#!/usr/bin/env bash
# Regenerate the PdSe2 Wannier model and its operators from the committed inputs.
# Nothing bulky is shipped with this example (repo policy): this script reproduces
# everything from qe/*.in and w90/*.{win,pw2wan*.in}.
#
# Requires on PATH: Quantum ESPRESSO (pw.x, pw2wannier90.x) and Wannier90
# (wannier90.x), plus the fully-relativistic PBE PAW pseudopotentials for Pd and
# Se. The pseudos are NOT committed (~9 MB); fetch them (PSL 1.0.0 / pseudo-dojo):
#     Pd.rel-pbe-n-kjpaw_psl.1.0.0.UPF
#     Se.rel-pbe-dn-kjpaw_psl.1.0.0.UPF
# and put them where scf.in's pseudo_dir points.
#
# Outputs (none committed): pdse2_proj_hr.dat, .chk, .xyz, .spn, .mmn, .amn, .eig,
# and the operator _hr.dat / CSR files the SHC workflow consumes.
set -euo pipefail
SEED=pdse2_proj
W2SP="${W2SP_BIN:-../../build/wannier2sparse}"

# 1. DFT: SCF then NSCF on the FULL uniform 4x4x1 grid (noncollinear + SOC).
pw.x < qe/scf.in  > scf.out
pw.x < qe/nscf.in > nscf.out

# 2. Wannier90 pre-processing -> <seed>.nnkp
cp w90/${SEED}.win .
wannier90.x -pp ${SEED}

# 3. overlaps/projections/eigenvalues (.mmn .amn .eig) and the spin matrices (.spn)
pw2wannier90.x < w90/pw2wan.in     > pw2wan.out      # writes .mmn .amn .eig
pw2wannier90.x < w90/pw2wan_spn.in > pw2wan_spn.out  # writes .spn (needs SOC run)

# 4. Wannierize -> ${SEED}_hr.dat (write_hr), ${SEED}.xyz (write_xyz),
#    ${SEED}.chk, ${SEED}_u.mat (need write_u_matrices=.true. for --exact-spin),
#    ${SEED}_band.dat (bands_plot). For the covariant velocity (Berry connection)
#    also set write_rmn=.true. to emit ${SEED}_r.dat.
wannier90.x ${SEED}

# 5. Operators in CSR for a finite supercell (NxNx1). Velocity is geometric here;
#    spin is the exact gauge transform from .spn + _u.mat; J^z_x is the supercell
#    anticommutator 1/2{V_x,S_z} (authoritative for the off-diagonal spin).
N="${1:-50}"
"$W2SP" ${SEED} "$N" "$N" 1 VX VY --exact-spin --spin-current X Z -o out
echo "done -> out/${SEED}.{HAM,VX,VY,SXexact,SYexact,SZexact,JXSZ}.CSR"

# 6. Figures (exact-diagonalization reference; needs the operator _hr.dat files).
#    Bands on the SAME k-path the DFT used (from qe/bands.in); SHC + DOS; then one
#    combined figure (FIG. 1 and FIG. 3 of the tutorial).
ED=../../tools/hr_exactdiag.py
EF=-1.3162
python3 $ED bands ${SEED} --kpath qe/bands.in --ef $EF --out img/pdse2_bands
python3 $ED shc   ${SEED} --jop ${SEED}_JXSZ_hr.dat --vop ${SEED}_vy_hr.dat \
        --shc-units e2h --nk 300 --eta 0.006 --ngrid 1200 --emin -3.0 --emax 1.0 --out pdse2_shc
python3 $ED dos   ${SEED} --nk 300 --eta 0.006 --ngrid 1200 --emin -3.0 --emax 1.0 --out pdse2_dos
python3 ../plot_hall.py pdse2_shc.json --dos pdse2_dos.json --ef $EF \
        --gap -0.05 0.05 --gap2 1.175 1.230 --plateau \
        --ylabel '$\sigma^z_{xy}\;[e^2/h]$' --out img/pdse2_shc_dos.png
echo "figures -> img/pdse2_bands.png, img/pdse2_shc_dos.png"
