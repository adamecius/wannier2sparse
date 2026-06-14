#!/usr/bin/env bash
# Usage:  ./run.sh <model> [N]      e.g.  ./run.sh graphene 80
# Generates the model, expands it to a supercell with wannier2sparse, and
# plots the DOS (KPM) to <model>/<model>_dos.png.
set -euo pipefail
MODEL="${1:?model: chain1d|graphene|cubic|haldane}"; N="${2:-80}"
BIN="${W2SP_BIN:-../build/wannier2sparse}"          # override with W2SP_BIN=...
command -v "$BIN" >/dev/null 2>&1 || [ -x "$BIN" ] || { echo "Build the tool first (see ../README.md), or set W2SP_BIN=/path/to/wannier2sparse"; exit 1; }
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"   # -> absolute path
HERE="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$HERE/models/tb"
( cd "$HERE/models/tb" && python3 "$HERE/gen_models.py" )  # writes all models under models/tb (idempotent)
cd "$HERE/models/tb/$MODEL"
case "$MODEL" in
  chain1d|chain1d_mag) DIMS="$N 1 1" ;;                                      # 1D (incl. spin-doubled)
  cubic)   M=$(python3 -c "print(max(4,int($N**0.5)))"); DIMS="$M $M $M" ;;  # 3D, modest
  *)       DIMS="$N $N 1" ;;                                                 # 2D
esac
echo ">> wannier2sparse $MODEL $DIMS"
"$BIN" "$MODEL" $DIMS
python3 "$HERE/w2s_dos.py" "${MODEL}.HAM.CSR" --title "$MODEL DOS (KPM)" --out "${MODEL}_dos.png"
echo ">> open ${MODEL}/${MODEL}_dos.png"
