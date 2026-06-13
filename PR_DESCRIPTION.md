# Replace string-keyed hopping map with flat append-only storage

## Summary

This PR removes the string-keyed `map` from the hot path used to parse, wrap, and export Wannier hoppings. The hopping list is now stored as a flat append-only vector, with duplicate edges intentionally preserved until Eigen assembles the CSR matrix from triplets.

The change follows the performance analysis: `hr.dat` already gives explicit connectivity, so supercell generation is not a geometric neighbor-search problem. The expensive part was constructing a textual tag for every replicated hopping and inserting it into a tree, even though production code never used that tag for lookup.

## What changed

- Replaced production hopping storage with `hopping_storage`, a vector-backed container.
- Changed Wannier parsing and supercell wrapping to append hoppings with `push_back`.
- Pre-reserved the generated supercell hopping capacity as `num_hoppings * cellDim[0] * cellDim[1] * cellDim[2]`.
- Updated CSR export and `tbmodel` operator builders to iterate over direct hopping tuples instead of map entries.
- Kept a small compatibility path for legacy tests/diagnostics that still call `hoppings.insert({tag, hop})` or `hoppings[tag]`.
- Documented why tag construction is no longer part of the hot path and why duplicate sparse entries are safe.

## Why this is safe

Eigen's `setFromTriplets` already combines duplicate `(row, col)` entries during sparse matrix assembly. The previous map-based implementation partially summed duplicates while wrapping; the new implementation defers that summation to the CSR assembly step.

The output semantics are therefore preserved for the sparse matrices, while avoiding per-hopping string construction and tree insertion during supercell replication.

## Verification

- `cmake -S . -B /tmp/wannier2sparse-build` succeeds.
- `git diff --check` passes.
- `cmake --build /tmp/wannier2sparse-build` could not complete in this environment because Eigen is not installed at `/usr/include/eigen3/Eigen`.

## Follow-up

A larger second step could bypass replicated `hopping_t` storage entirely and assemble CSR triplets directly during wrapping. That would also remove the remaining tuple/vector copy cost, but this PR intentionally keeps the behavioral change small.
