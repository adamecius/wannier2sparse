#!/usr/bin/env bash
# Plan 4b - provision Wannier90 (wannier90.x, postw90.x) locally, no sudo, in /tmp.
#
# Two routes are supported; the script picks whichever the machine allows:
#   (A) source build  - if gfortran + BLAS/LAPACK are present (e.g. after
#                        `sudo apt install gfortran liblapack-dev libblas-dev`).
#   (B) conda-forge    - otherwise, fetch a standalone micromamba binary and
#                        create an env with the prebuilt `wannier90` package.
#
# Outputs (sourced by gen_fixture.sh):
#   $W90_BIN   -> directory containing wannier90.x / postw90.x
# Everything lives under /tmp (ephemeral); re-run after a reboot.
set -euo pipefail

W90_VERSION="${W90_VERSION:-3.1.0}"
WORK="${W90_WORK:-/tmp/w90build}"
ENVDIR="${W90_ENV:-/tmp/w90env}"
MAMBA_ROOT="${MAMBA_ROOT:-/tmp/mamba}"
mkdir -p "$WORK"

have() { command -v "$1" >/dev/null 2>&1; }

build_from_source() {
    echo "[build_w90] gfortran found - building Wannier90 $W90_VERSION from source"
    local tgz="$WORK/v${W90_VERSION}.tar.gz"
    [ -f "$tgz" ] || curl -L -o "$tgz" \
        "https://github.com/wannier-developers/wannier90/archive/refs/tags/v${W90_VERSION}.tar.gz"
    tar -xzf "$tgz" -C "$WORK"
    local src="$WORK/wannier90-${W90_VERSION}"
    cp "$src/config/make.inc.gfort" "$src/make.inc"
    make -C "$src" -j"$(nproc)" wannier90.x postw90.x
    echo "$src" > "$WORK/.w90_bin_dir"
}

provision_with_micromamba() {
    echo "[build_w90] no gfortran - provisioning wannier90 via conda-forge (micromamba)"
    local mm="$MAMBA_ROOT/bin/micromamba"
    if [ ! -x "$mm" ]; then
        mkdir -p "$MAMBA_ROOT"
        # NOTE: downloads a standalone micromamba binary. In a locked-down sandbox
        # this fetch may require an explicit Bash permission rule.
        curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
            | tar -xj -C "$MAMBA_ROOT" bin/micromamba
    fi
    "$mm" create -y -p "$ENVDIR" -c conda-forge "wannier90>=3.1"
    echo "$ENVDIR/bin" > "$WORK/.w90_bin_dir"
}

if have gfortran; then build_from_source; else provision_with_micromamba; fi

W90_BIN="$(cat "$WORK/.w90_bin_dir")"
echo "[build_w90] wannier90.x -> $W90_BIN"
"$W90_BIN/wannier90.x" --version 2>/dev/null || ls -la "$W90_BIN"/wannier90.x
