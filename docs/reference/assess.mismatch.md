# Assess Spectral Mismatch Between Reference and Test Spectra

For each fluorophore common to both `reference.variants` and
`test.variants`, computes the cosine similarity between the reference
spectrum (row 1 of each fluorophore's variant matrix) and the
corresponding test spectrum. Typically used to compare a spectral
reference library built from one particle type (e.g. Cells) against one
built from another (e.g. a bead type), by passing the reference particle
type as `reference.variants`.

`reference.variants` and `test.variants` are expected to come from
independent
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
calls and are not guaranteed to share the same detector set or column
order, so the comparison is aligned explicitly by detector name rather
than position. A warning is issued if the two detector sets differ.

## Usage

``` r
assess.mismatch(reference.variants, test.variants)
```

## Arguments

- reference.variants:

  Named list of spectral variant matrices (one per fluorophore), as
  returned in the `variants` element of
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md).
  Row 1 of each matrix is treated as the reference spectrum.

- test.variants:

  Named list of spectral variant matrices in the same format as
  `reference.variants`, to be compared against it.

## Value

A one-column numeric matrix of cosine similarity values, one row per
fluorophore common to both lists.
