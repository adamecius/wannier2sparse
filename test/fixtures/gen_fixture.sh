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
[ -x "$W90_BIN/wannier90.x" ] || { echo "wannier90.x not found; run build_w90.sh first"; exit 1; }

mkdir -p "$OUT"
run="$(mktemp -d /tmp/fixrun.XXXXXX)"
cp -r "$TUT"/. "$run"/
cd "$run"

# Ensure wannier2sparse-relevant outputs are written by appending flags idempotently.
ensure_flag() { grep -qiE "^\s*$1\s*[:=]" "$SEED.win" || echo "$1 = $2" >> "$SEED.win"; }
ensure_flag write_hr           .true.
ensure_flag write_xyz          .true.

if [ "$MODE" = "soc" ]; then
    [ -x "$QE_BIN/pw.x" ] || { echo "pw.x not found; run build_qe.sh first"; exit 1; }
    ensure_flag write_u_matrices .true.
    # scf -> nscf -> overlaps/projections/spin -> wannierise
    "$QE_BIN/pw.x"            -in "$SEED.scf.in"      > scf.out
    "$QE_BIN/pw.x"            -in "$SEED.nscf.in"     > nscf.out
    "$W90_BIN/wannier90.x" -pp "$SEED"                            # writes SEED.nnkp
    "$QE_BIN/pw2wannier90.x" -in "$SEED.pw2wan.in"   > pw2wan.out # writes amn/mmn/eig/spn
fi

# use_ws_distance is the point of the Level-1 wsdist fixture (produces SEED_wsvec.dat)
[ "$MODE" = "wsdist" ] && ensure_flag use_ws_distance .true.

"$W90_BIN/wannier90.x" "$SEED" > "wannier90.$SEED.out"

# Collect what wannier2sparse consumes (plus provenance).
for f in "${SEED}_hr.dat" "${SEED}_wsvec.dat" "${SEED}.eig" "${SEED}_centres.xyz" \
         "${SEED}.spn" "${SEED}.amn" "${SEED}_u.mat" "${SEED}_u_dis.mat" \
         "${SEED}.win" "${SEED}.wout"; do
    [ -f "$f" ] && cp "$f" "$OUT"/ && echo "  collected $f"
done
echo "[gen_fixture] fixture for '$SEED' ($MODE) -> $OUT"
ls -la "$OUT"
