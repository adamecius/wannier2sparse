# Plan 4b — local W90/QE fixture bootstrap (development scaffolding)

These scripts provision Wannier90 (and, for SOC, Quantum ESPRESSO) **locally,
without sudo, under `/tmp`**, run a tutorial, and collect the real files that
Plans 5/7/8 consume. This is a **developer tool**, not the final golden set
(curated, versioned fixtures with pinned tags + ctest integration are deferred,
see the roadmap).

## What each plan needs
| Plan | Files | Level |
|------|-------|-------|
| P1, P2, P6 | `seedname_hr.dat`, `.eig` | 1 (W90 only) |
| **P5** `wsvec` | `seedname_wsvec.dat`, `_hr.dat` (run with `use_ws_distance`) | 1 (W90 only) |
| **P7** exact spin | `.spn`, `seedname_u.mat`, `seedname_u_dis.mat`, `.eig`, `.win` (`mp_grid`) | 2 (W90 + QE, SOC) |
| **P8** orbital L | `.amn`, `.win` projection block | 2 (W90 + QE) |

## Routes (auto-selected)
- **Source build** when a Fortran toolchain is present:
  `sudo apt install gfortran cmake liblapack-dev libblas-dev` (+ `libfftw3-dev`
  for QE). `build_w90.sh`/`build_qe.sh` then download and build into `/tmp`.
- **conda-forge (no sudo)** otherwise: the scripts fetch a standalone
  `micromamba` binary and create envs with the prebuilt `wannier90` / `qe`
  packages. In a locked-down sandbox the `micromamba` download may need an
  explicit permission rule for `curl https://micro.mamba.pm/...`.

## Usage
```bash
# 1. provision (once per machine/boot; binaries cached in /tmp)
bash test/fixtures/build_w90.sh            # wannier90.x  (Level 1 + 2)
bash test/fixtures/build_qe.sh             # pw.x, pw2wannier90.x (Level 2 only)

# 2. generate fixtures into a seedname folder that `--project` resolves
#    Level 1: use_ws_distance fixture (Plan 5) from precomputed overlaps
bash test/fixtures/gen_fixture.sh wsdist gaas  /path/to/w90/tutorials/gaas  /tmp/fix/gaas
#    Level 2: SOC fixture (Plans 7/8) — pseudos ship with the W90 SOC tutorials
bash test/fixtures/gen_fixture.sh soc    Fe    /path/to/w90/tutorials/Fe    /tmp/fix/Fe

# 3. drive wannier2sparse against a generated fixture
build/wannier2sparse run 4 4 1 --project /tmp/fix/gaas --seed gaas --bounds
```

## Suggested tutorials (per the roadmap)
- Level 1: a tutorial that ships precomputed `.amn/.mmn/.eig` (e.g. GaAs), run
  with `use_ws_distance=.true.` → exercises Plan 5.
- Level 2 SOC (pseudos included, `pseudo_dir='../../pseudo/'`): `tutorial17` Fe
  (`write_spn`, smallest, for Plan 7), `tutorial30` GaAs (ac spin Hall),
  `tutorial29` Pt (heavier spin-Hall reference).

## Validation (catches FT/convention bugs a regression would miss)
Reconstruct `H(k) = Σ_R e^{ik·R} H(R)/ndegen` from the ingested model and compare
against `seedname_band.dat` / `.eig`. For spin: Hermiticity and sum rules of
`S_W(R)`.

## Notes / caveats
- `/tmp` is ephemeral; re-run `build_*.sh` after a reboot (binaries are cached
  within a session).
- Exact tutorial file names vary across W90/QE versions — adjust `SEED` and the
  tutorial path; `gen_fixture.sh` appends the `write_hr/write_xyz/use_ws_distance`
  flags it needs idempotently.
- **Always verify the `.spn` packing order and `wsvec`/`u.mat` layouts against the
  W90 / `pw2wannier90` source**, not from memory, before trusting P5/P7/P8 output.
