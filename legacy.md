# Legacy and advanced usage

This collects features that earlier versions of wannier2sparse documented on the
front page. They still work, but the [`.w2s` input file](docs/input_file.md) is now
the recommended way to drive a run.

## Positional command line

Before the input file, a run was given entirely on the command line:

```bash
wannier2sparse LABEL N1 N2 N3 [OP ... | all] [options]
```

For example, `wannier2sparse graphene 80 80 1 VX VY VXSZ -o out` is equivalent to,
and produces byte-identical output to, the `.w2s` run with the same `label`,
`supercell`, and `operators`. The flags map one-to-one onto the input-file keys
(`--exact-spin` → `exact_spin`, `--velocity-mode covariant` → `velocity_mode`,
`-o DIR` → `output_dir`, and so on); see `wannier2sparse --help`.

## Run log and summary

Every run also writes a log next to its outputs (`label.run.log`) and prints a
short summary at the end — per-stage wall time, peak memory, and the count of
warnings and errors. A machine-readable copy of the same information, together with
the list of operators written and the input files behind each, is left in
`label.out`. Console verbosity is set with `--verbose` / `--quiet`; the log file
can be redirected with `--log-file` or suppressed with `--no-log-file`.

## Using wannier2sparse as a C++ library

Installing exports a CMake config-package, so a C++ project can link the expansion
engine directly:

```cmake
find_package(wannier2sparse 1.0 CONFIG REQUIRED)
target_link_libraries(app PRIVATE wannier2sparse::wannierlib)
```

Point CMake at the install prefix with `-DCMAKE_PREFIX_PATH=/your/prefix`. The
imported target carries its include directory and propagates `Eigen3::Eigen` when
Eigen was found through its own `Eigen3Config.cmake`; if Eigen was located through
the manual `-DEIGEN_DIR=...` fallback, the consumer must provide Eigen itself.
