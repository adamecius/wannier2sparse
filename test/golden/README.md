# WannierBerri committed goldens (Plan 11)

Reference operator data generated **once** by WannierBerri and versioned here, so
the cross-check tests can run without the package. See `docs/conventions.md` §7.

The two C++ tests read these files and **never import WannierBerri**:

| golden                          | test                              | what it checks                          |
|---------------------------------|-----------------------------------|-----------------------------------------|
| `<seed>_<op>_matrix.ref`        | `wberri_matrix_crosscheck`        | Wannier-gauge `O_W(k)` element-by-element (tight, ~1e-6) |
| `<seed>_<op>_texture.ref`       | `wberri_texture_crosscheck`       | band texture `⟨O_α⟩_{nk}` (looser, ~1e-5, degeneracy subspace-trace) |

`<op>` is `S` (spin, e.g. Fe) or `L` (orbital, e.g. copper / BaTiO3).

## Regenerating (needs WannierBerri + a fixture with `<seed>.chk`)

```sh
pip install wannierberri
export W2SP_WS_CONVENTION=...   # MUST match how <seed>_hr.dat was built (see §6/§7)
make regen-wberri-golden        # CMake custom target -> writes the .ref files here
```

After regenerating, **hand-validate** the numbers against the WannierBerri run
before committing the `.ref` files (the architecture trades a one-time manual check
for a fast, dependency-free CI).

## File formats

`matrix.ref` — header `seed operator n_k num_wann WS_convention wberri_version`,
then rows `ik m n alpha Re Im` (`m,n` 0-based Wannier indices, `alpha`∈{0,1,2}=x,y,z).

`texture.ref` — header adds a trailing `degeneracy_tol`, then rows
`ik ibnd alpha block value`; bands sharing a `block` id are degenerate
(`|E_i−E_j| < degeneracy_tol`) and the comparator sums `⟨O_α⟩` over the block.

## Skip-clean

If a `.ref` (or the underlying fixture) is absent, the corresponding test returns
the ctest `SKIP_RETURN_CODE` → **SKIPPED**, not FAIL. The two tests skip
independently. A default `ctest` (no `-L wberri`) is unaffected.
