#ifndef WANNIER2SPARSE_INPUT_FILE_HPP
#define WANNIER2SPARSE_INPUT_FILE_HPP

// Human-editable input file for wannier2sparse (the --create / --write / --run CLI).
//
// Problem solved
// --------------
// The positional CLI (`wannier2sparse LABEL N1 N2 N3 VX VY --exact-spin ...`) is
// fine for one-off runs but is not self-documenting and is awkward to build up or
// reproduce. This module adds an lsquant-style workflow that *preserves* those old
// options by recording them, one per line, in a plain `key = value` text file:
//
//   wannier2sparse --create run.inp                  # write a commented template
//   wannier2sparse --write covariant-velocity -inp run.inp   # set/append one option
//   wannier2sparse --write VX -inp run.inp                   # append an operator
//   wannier2sparse --write supercell 50 50 1 -inp run.inp
//   wannier2sparse --run run.inp                     # execute; writes <label>.w2sp.out
//
// Format
// ------
// `key = value`, `#` starts a comment, blank lines ignored. Repeatable keys
// (`spin_current`, `op_file`) may appear on several lines. The keys are exactly the
// old CLI options (label, supercell, operators, exact_spin, ...), so nothing about
// the established formats is lost; the file is just a serialization of them.
//
// Units/conventions: unchanged from the rest of the tool (see docs/operators.md,
// docs/conventions.md). Pitfall: `covariant_velocity` is recognized and recorded
// but the operator route is not implemented yet (docs/operators.md §3); `--run`
// errors clearly if it is set, pointing to the `--op-file` workaround.

#include <string>
#include "w2sp_arguments.hpp"

namespace w2sp {

// Write a commented template input file at `path` (overwrites). Lists every key
// with its default and a one-line explanation.
void input_file_create(const std::string& path);

// Set `key` to `value` in `path` (creating the file if absent). If `append` and the
// key already has a line, the value is appended to it (used for the operator list);
// otherwise the existing line is replaced, or a new line added. Repeatable keys
// (spin_current, op_file) always add a new line.
void input_file_set(const std::string& path, const std::string& key,
                    const std::string& value, bool append);

// Parse `path` and populate `args` (label, supercell, operators, flags, ...). Throws
// std::runtime_error on a malformed file or an unknown key.
void input_file_apply(const std::string& path, W2SP_arguments& args);

// Write the run summary (the "output file", lsquant-style): the resolved options
// and the operator files produced (name + byte size). `outpath` is the file to write.
void input_file_write_output(const std::string& outpath, const W2SP_arguments& args);

} // namespace w2sp

#endif // WANNIER2SPARSE_INPUT_FILE_HPP
