#!/usr/bin/env bash
# Plan 4b - generate a wannier2sparse fixture from a Wannier90 (+QE) tutorial.
#
# Usage:
#   gen_fixture.sh wsdist  SEED TUTORIAL_DIR OUTDIR     # Level 1: hr + wsvec + eig
#   gen_fixture.sh soc     SEED TUTORIAL_DIR OUTDIR     # Level 2: + .spn + .amn (QE)
#
#   SEED          wannier90 seedname (e.g. gaas, Fe, GaAs)
#   TUTORIAL_DIR  directory with the tutorial inputs (SEED.win, SCF/NSCF for soc,
#                 and precomputed SEED.amn/.mmn/.eig for the Level-1 wsdist mode)
#   OUTDIR        destination seedname folder that Plan 4's --project resolves
#
# Requires W90_BIN (from build_w90.sh) and, for `soc`, QE_BIN (from build_qe.sh).
# This is dev scaffolding: exact tutorial file names differ across W90/QE
# versions, so adjust SEED/paths to your tutorial set (see README.md).
set -euo pipefail

MODE="${1:?mode: wsdist|soc}"
SEED="${2:?seedname}"
TUT="${3:?tutorial dir}"
OUT="${4:?output dir}"

WORK_W90="${W90_WORK:-/tmp/w90build}"
WORK_QE="${QE_WORK:-/tmp/qebuild}"
W90_BIN="${W90_BIN:-$(cat "$WORK_W90/.w90_bin_dir" 2>/dev/null || true)}"
QE_BIN="${QE_BIN:-$(cat "$WORK_QE/.qe_bin_dir" 2>/dev/null || true)}"

# Resolve binaries: prefer the provisioned dirs, else fall back to PATH (e.g. a
# system-installed quantum-espresso).
W90X="$W90_BIN/wannier90.x"; [ -x "$W90X" ] || W90X="$(command -v wannier90.x || true)"
PWX="$QE_BIN/pw.x";          [ -x "$PWX" ]  || PWX="$(command -v pw.x || true)"
P2WX="$QE_BIN/pw2wannier90.x"; [ -x "$P2WX" ] || P2WX="$(command -v pw2wannier90.x || true)"
[ -x "$W90X" ] || { echo "wannier90.x not found; run build_w90.sh first"; exit 1; }

# First existing of a set of candidate input names (W90 examples use SEED.scf;
# some setups use SEED.scf.in).
pick() { for f in "$@"; do [ -f "$f" ] && { echo "$f"; return; }; done; echo "$1"; }

mkdir -p "$OUT"
run="$(mktemp -d /tmp/fixrun.XXXXXX)"
cp -r "$TUT"/. "$run"/
cd "$run"

# Ensure wannier2sparse-relevant outputs are written by appending flags idempotently.
ensure_flag() { grep -qiE "^\s*$1\s*[:=]" "$SEED.win" || echo "$1 = $2" >> "$SEED.win"; }
ensure_flag write_hr           .true.
ensure_flag write_xyz          .true.

if [ "$MODE" = "soc" ]; then
    [ -x "$PWX" ] && [ -x "$P2WX" ] || { echo "pw.x/pw2wannier90.x not found; run build_qe.sh or apt-install quantum-espresso"; exit 1; }
    ensure_flag write_u_matrices .true.
    # Pseudopotentials ship next to the W90 tutorials (TUT/../../pseudo). The
    # example's relative pseudo_dir breaks once inputs are copied here, so pin it
    # to the absolute location (and ESPRESSO_PSEUDO as a fallback).
    PSEUDO_DIR="$(cd "$TUT/../../pseudo" 2>/dev/null && pwd || true)"
    if [ -n "$PSEUDO_DIR" ]; then
        export ESPRESSO_PSEUDO="$PSEUDO_DIR"
        for inp in "$SEED.scf" "$SEED.scf.in" "$SEED.nscf" "$SEED.nscf.in"; do
            [ -f "$inp" ] && sed -i "s|pseudo_dir *=.*|pseudo_dir = '$PSEUDO_DIR/',|" "$inp"
        done
    fi
    # scf -> nscf -> overlaps/projections/spin -> wannierise
    "$PWX"  -in "$(pick "$SEED.scf"  "$SEED.scf.in")"  > scf.out
    "$PWX"  -in "$(pick "$SEED.nscf" "$SEED.nscf.in")" > nscf.out
    "$W90X" -pp "$SEED"                                            # writes SEED.nnkp
    "$P2WX" -in "$(pick "$SEED.pw2wan" "$SEED.pw2wan.in")" > pw2wan.out  # amn/mmn/eig/spn
fi

# use_ws_distance is the point of the Level-1 wsdist fixture (produces SEED_wsvec.dat)
[ "$MODE" = "wsdist" ] && ensure_flag use_ws_distance .true.

"$W90X" "$SEED" > "wannier90.$SEED.out"

# Collect what wannier2sparse consumes (plus provenance).
for f in "${SEED}_hr.dat" "${SEED}_wsvec.dat" "${SEED}.eig" "${SEED}_centres.xyz" \
         "${SEED}.spn" "${SEED}.amn" "${SEED}_u.mat" "${SEED}_u_dis.mat" \
         "${SEED}.win" "${SEED}.wout"; do
    [ -f "$f" ] && cp "$f" "$OUT"/ && echo "  collected $f"
done
echo "[gen_fixture] fixture for '$SEED' ($MODE) -> $OUT"
ls -la "$OUT"
