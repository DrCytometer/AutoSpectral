# Spectral-Location Alignment Between Variability and Mismatch

For each fluorophore common to both `variability.mad` and
`mismatch.dist`, computes the cosine similarity between that
fluorophore's per-detector variability profile (e.g. the denoised MAD
profile from
[`assess.variability.mad()`](https://drcytometer.github.io/AutoSpectral/reference/assess.variability.mad.md))
and the per-detector magnitude of its bead-vs-cell mismatch
([`abs()`](https://rdrr.io/r/base/MathFun.html) of, e.g.,
[`bead.cell.dist()`](https://drcytometer.github.io/AutoSpectral/reference/bead.cell.dist.md)).
Both profiles are non-negative by construction, so the resulting cosine
similarity is bounded in `[0, 1]` and reflects purely where in detector
space each quantity is concentrated – a value near 1 means the
fluorophore's variant-to-variant variability and its bead-vs-cell
mismatch peak at the same detectors, regardless of either quantity's
overall magnitude; a value near 0 means they are concentrated in
different, non- overlapping parts of the spectrum.

`variability.mad` and `mismatch.dist` are expected to come from
independent computations and are not guaranteed to share the same
detector set or column order, so the comparison is aligned explicitly by
detector name rather than position. A warning is issued if the two
detector sets differ.

## Usage

``` r
assess.variability.alignment(variability.mad, mismatch.dist)
```

## Arguments

- variability.mad:

  Numeric matrix, fluorophores in rows and detectors in columns, as
  returned by
  [`assess.variability.mad()`](https://drcytometer.github.io/AutoSpectral/reference/assess.variability.mad.md).

- mismatch.dist:

  Numeric matrix, fluorophores in rows and detectors in columns, as
  returned by
  [`bead.cell.dist()`](https://drcytometer.github.io/AutoSpectral/reference/bead.cell.dist.md).
  Values are used as `abs(mismatch.dist)` so that alignment reflects
  spectral location rather than the sign of the mismatch.

## Value

A one-column numeric matrix of cosine similarity values (column name
`"Alignment"`), one row per fluorophore common to both inputs.
