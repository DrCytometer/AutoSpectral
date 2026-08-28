# Finalize an auto-generated control file for bead/cell comparison

Post-processes a control file created by
[`create.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/create.control.file.md)
so it is usable by either the legacy or automated AutoSpectral pipeline
in the context of a bead-vs-cell comparison run. Sets `control.type` for
every row according to the particle type of the folder, forces the
single unstained/AF row to `control.type == "cells"` (required by both
pipelines' validation), and back-fills `universal.negative` for the
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
