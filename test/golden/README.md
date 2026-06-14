# WannierBerri committed goldens (Plan 11)

Reference operator data generated **once** by WannierBerri and versioned here, so
the cross-check tests can run without the package. See `docs/conventions.md` §7.

The two C++ tests read these files and **never import WannierBerri**:

| golden                          | test                              | what it checks                          |
|---------------------------------|-----------------------------------|-----------------------------------------|
| `<seed>_<op>_matrix.ref`        | `wberri_matrix_crosscheck`        | Wannier-gauge `O_W(k)` element-by-element (tight, ~1e-6) |
| `<seed>_<op>_texture.ref`       | `wberri_texture_crosscheck`       | band texture `⟨O_α⟩_{nk}` (looser, ~1e-5, degeneracy subspace-trace) |

`<op>` is `S` (spin) or `L` (orbital). **Committed so far:** `Fe_S_matrix.ref` +
`Fe_S_texture.ref`, generated with WannierBerri 26.4.6 and hand-validated against
the w2s reconstruction (matrix `2.6e-9`, texture `1.3e-3`; see §7). The orbital
`L` goldens are **not** committed yet — the WannierBerri↔projector-route
correspondence is unconfirmed, so the generator refuses to emit them (§7 Deferred);
the C++ `L` tests run the moment an `L` golden appears.

`Fe_S_matrix.ref` stores a **stride-8 k-subset** (`W2SP_MATRIX_KSTRIDE`) to keep
the file small; the `ik` column says which k each row is, and the test reconstructs
exactly those. `Fe_S_texture.ref` covers the full grid.

## Regenerating (needs WannierBerri + a fixture with `<seed>.chk`)

```sh
pip install wannierberri fortio   # fortio reads the binary .chk
export W2SP_WS_CONVENTION=II       # MUST match how <seed>_hr.dat was built (see §6/§7)
export W2SP_DEGENERACY_TOL=1e-3    # texture grouping >> the ~4e-5 _hr.dat noise (§7)
make regen-wberri-golden          # CMake custom target -> writes the .ref files here
```

The fixture prefix must point at a dir holding `<seed>.chk` + `.mmn` + `.eig` +
`.spn` + `_u.mat` (the `.chk`/`.mmn` live in the gen working dir, not the slim
`/tmp/fix/<seed>`). After regenerating, **hand-validate** the numbers against the
WannierBerri run before committing (the architecture trades a one-time manual check
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
