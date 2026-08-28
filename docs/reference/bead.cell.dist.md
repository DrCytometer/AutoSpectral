# Per-Detector Spectral Distance Between Particle Types

For each fluorophore present in `cell.variants`, computes the raw
per-detector difference (cell reference spectrum minus bead reference
spectrum) between two particle types' extracted reference spectra.

`cell.variants` and `bead.variants` are expected to come from
independent
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
calls (e.g. one per particle-type folder) and are not guaranteed to
share the same detector set or column order, so the comparison is
aligned explicitly by detector name rather than position. A warning is
issued if the two particle types were extracted with different detector
sets.

## Usage

``` r
bead.cell.dist(cell.variants, bead.variants)
```

## Arguments

- cell.variants:

  Named list of spectral variant matrices (one per fluorophore), as
  returned in the `variants` element of
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md),
  for the reference particle type (typically Cells). Row 1 of each
  matrix is treated as the reference spectrum.

- bead.variants:

  Named list of spectral variant matrices in the same format as
  `cell.variants`, for the particle type being compared against the
  reference (e.g. a bead type).

## Value

A numeric matrix of per-detector differences (cell minus bead), one row
per fluorophore common to both lists, columns are the shared detector
channels.
