#!/usr/bin/env bash
# Plan 4b - provision Quantum ESPRESSO (pw.x, pw2wannier90.x) locally, no sudo.
# Needed only for the SOC fixtures (.spn for Plan 7, .amn from scratch for Plan 8);
# the Level-1 fixtures (hr/wsvec/eig from precomputed overlaps) do NOT need QE.
#
# Routes, picked automatically:
#   (A) source build  - gfortran + BLAS/LAPACK + FFTW present:
#                        configure --disable-parallel && make pw pp w90
#   (B) conda-forge    - otherwise micromamba env with the prebuilt `qe` package.
#
# Outputs:
#   $QE_BIN -> directory containing pw.x and pw2wannier90.x
set -euo pipefail

QE_VERSION="${QE_VERSION:-7.2}"
WORK="${QE_WORK:-/tmp/qebuild}"
ENVDIR="${QE_ENV:-/tmp/qeenv}"
MAMBA_ROOT="${MAMBA_ROOT:-/tmp/mamba}"
mkdir -p "$WORK"

have() { command -v "$1" >/dev/null 2>&1; }

build_from_source() {
    echo "[build_qe] toolchain found - building QE $QE_VERSION (serial) from source"
    local tgz="$WORK/qe-${QE_VERSION}.tar.gz"
    [ -f "$tgz" ] || curl -L -o "$tgz" \
        "https://github.com/QEF/q-e/archive/refs/tags/qe-${QE_VERSION}.tar.gz"
    tar -xzf "$tgz" -C "$WORK"
    local src="$WORK/q-e-qe-${QE_VERSION}"
    ( cd "$src" && ./configure --disable-parallel && make -j"$(nproc)" pw pp )
    echo "$src/bin" > "$WORK/.qe_bin_dir"
}

provision_with_micromamba() {
    echo "[build_qe] no FFTW/gfortran toolchain - provisioning qe via conda-forge"
    local mm="$MAMBA_ROOT/bin/micromamba"
    if [ ! -x "$mm" ]; then
        mkdir -p "$MAMBA_ROOT"
        curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
            | tar -xj -C "$MAMBA_ROOT" bin/micromamba
    fi
    "$mm" create -y -p "$ENVDIR" -c conda-forge qe
    echo "$ENVDIR/bin" > "$WORK/.qe_bin_dir"
}

if have gfortran && [ -e /usr/include/fftw3.f03 -o -e /usr/include/fftw3.h ]; then
    build_from_source
else
    provision_with_micromamba
fi

QE_BIN="$(cat "$WORK/.qe_bin_dir")"
echo "[build_qe] pw.x -> $QE_BIN"
ls -la "$QE_BIN/pw.x" "$QE_BIN/pw2wannier90.x"
