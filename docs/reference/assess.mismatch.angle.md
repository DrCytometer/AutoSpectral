# Assess Spectral Angle Between Reference and Test Spectra

As
[`assess.mismatch()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.md),
but reports the spectral angle (in degrees) between each fluorophore's
reference and test spectra rather than their cosine similarity. The
angle is [`acos()`](https://rdrr.io/r/base/Trig.html) of the same cosine
value
[`assess.mismatch()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.md)
returns, so the two metrics are a strictly monotonic transform of one
another – but unlike cosine similarity, spectral angle increases with
divergence, so it can be plotted directly against other
divergence-increasing metrics (e.g. mismatch distance) without an
inverted axis.

## Usage

``` r
assess.mismatch.angle(reference.variants, test.variants)
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

A one-column numeric matrix of spectral angle values in degrees (column
name `"Angle"`), one row per fluorophore common to both lists.
