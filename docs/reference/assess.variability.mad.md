# Denoised Per-Detector Variability (MAD)

As
[`assess.variability()`](https://drcytometer.github.io/AutoSpectral/reference/assess.variability.md),
but returns the full per-detector robust SD (MAD) profile across a
fluorophore's variant spectra, rather than a single summary value.
Detectors where the reference spectrum falls below `peak.height` are
zeroed out before computing MAD, so that off-peak noise does not
contribute to the profile.

## Usage

``` r
assess.variability.mad(variant.list, peak.height = 0.1)
```

## Arguments

- variant.list:

  Named list of spectral variant matrices (one per fluorophore), as
  returned in the `variants` element of
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md).
  Row 1 of each matrix is treated as the reference spectrum; rows 2+ are
  treated as variants.

- peak.height:

  Numeric in `[0, 1]`, default `0.1`. Detectors where the reference
  spectrum falls below this fraction of that fluorophore's own peak
  value (`peak.height * max(reference.spectrum)`) are set to zero before
  MAD is computed. Relative to the peak rather than an absolute cutoff,
  so behavior is unchanged whether `variant.list` is L-infinity or L2
  normalized.

## Value

A numeric matrix with one row per fluorophore and one column per
detector, giving the denoised per-detector MAD of the variant spectra.
