# The `.w2s` input file

A run is described by one `.w2s` file and executed with `wannier2sparse -x file.w2s`
(or `--run file.w2s`). It is JSON, with `//` line and `/* */` block comments allowed,
so the file can be self-documenting. Only the keys you set take effect; everything
else takes its default.

A minimal sparse run:

```json
{
  "label": "graphene",
  "mode": "sparse",
  "output_dir": "out",
  "supercell": [80, 80, 1],
  "operators": ["VX", "VY", "VXSZ"]
}
```

## Keys

| Key | Meaning | Default |
|-----|---------|---------|
| `label` | Model name; the tool reads `label_hr.dat`, and `label.xyz` / `label.uc` when an operator needs positions. | required |
| `mode` | `sparse` (expand to supercell CSR) or `bundle` (ship the unexpanded $O_{ij}(\mathbf{R})$ + manifest). | `sparse` |
| `project_dir` | Directory holding the input files. | `.` |
| `seed` | Seedname of the input files, if different from `label`. | `label` |
| `output_dir` | Where the outputs are written. | `.` |
| `supercell` | `[N1, N2, N3]` supercell size (sparse mode). | `[1,1,1]` |
| `operators` | Operator codes to build (see [Operators](operators.md)); `"all"` requests every one. The Hamiltonian is always written. | none |
| `exact_spin` | Build $S$ from the exact gauge transform (`.spn` + `_u.mat`). | `false` |
| `orbital_L` | Build orbital angular momentum $L$ (`.amn` + `_u.mat` + `.win`). | `false` |
| `velocity_mode` | Velocity ladder for `VX/VY/VZ` and the velocity inside a spin current: `bare`, `berry_connection`, or `covariant`. | `berry_connection` |
| `r_dat` | Position matrix `_r.dat` (Wannier90 `write_rmn`), needed for `velocity_mode: covariant`. | `<seed>_r.dat` |
| `provenance` | DFT/Wannier provenance, usually filled by `--provenance` (see [README](../README.md#provenance-tracking)). | none |

Spectral bounds, self-consistency checks, and the run log are produced
automatically; there are no keys to turn them on.

## Operators you provide yourself

To use an operator you built outside the tool (its own `_hr.dat`), put it in an
`operators/` folder beside the model and name it `<NAME>_hr.dat`. It is read by name
and expanded through the same engine, and written as `label.<NAME>.CSR`.

## CSR output format

Each `sparse` operator is one CSR file:

```
<dim> <nnz>
<value.real value.imag> ...   # nnz complex values
<column index> ...            # nnz column indices
<row pointer> ...             # dim+1 row pointers
```

The full Hermitian matrix is stored (both triangles); a consumer reconstructs
$H = \tfrac12(M + M^\dagger)$, which returns $M$ unchanged for a Hermitian operator.

## Bundle output

With `"mode": "bundle"` the supercell size is ignored and the run writes
`label.w2sp/` holding `manifest.json` and one `operators/<NAME>.hr.dat` per operator
— the primitive $O_{ij}(\mathbf{R})$ together with the crystal structure and the
DFT/Wannier provenance. A consumer forms $H(\mathbf{k})$ and its own supercell from
these.
