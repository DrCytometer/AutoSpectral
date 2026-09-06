# Finalize an auto-generated control file for bead/cell comparison

Post-processes a control file created by
[`create.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/create.control.file.md)
so it is usable by either the legacy or automated AutoSpectral pipeline
in the context of a bead-vs-cell comparison run. Sets `control.type` for
every row – including the single unstained row – according to the true
particle type of the folder, and back-fills `universal.negative` for the
remaining rows with that row's filename.

## Usage

``` r
finalize.control.file(
  control.file.path,
  particle.type,
  reference.type = "Cells"
)
```

## Arguments

- control.file.path:

  Character scalar. Path to the control file CSV to finalize in place.

- particle.type:

  Character scalar. The particle type label for the folder this control
  file belongs to (e.g. `"Cells"`, `"Beads"`).

- reference.type:

  Character scalar. The particle type treated as the reference (cell)
  population; any `particle.type` not identical to this is treated as
  beads. Defaults to `"Cells"`.

## Value

Invisibly, the finalized control file as a data frame. The same data is
also written back to `control.file.path`.

## Details

Exactly one row must resolve to `fluorophore %in% c("AF", "Negative")`.
If zero or more than one such row is found, the function leaves
`fluorophore`/`universal.negative` unset for that file and emits a
warning instructing the user to edit the control file manually rather
than guessing which row is the background control.

The reference (cell) folder's unstained row is labelled `"AF"`, per the
package's usual convention. A bead folder's unstained row is left
labelled `"Negative"` and typed `control.type == "beads"`, like every
other row in that folder – it is never relabelled to impersonate a
cell-based AF control.
[`run.bead.cell.comparison()`](https://drcytometer.github.io/AutoSpectral/reference/run.bead.cell.comparison.md)
resolves that row via
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)'s
own `unstained.sample` argument instead of requiring a literal `"AF"`
row.
