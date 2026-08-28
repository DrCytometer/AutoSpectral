# Assess Within-Fluorophore Spectral Variability

Summarises how much a fluorophore's spectral variants (rows 2+ of its
variant matrix, as produced by
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md))
disagree with its reference spectrum (row 1) and with each other,
returning several complementary summary metrics per fluorophore.

## Usage

``` r
assess.variability(variant.list, peak.height = 0.1)
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
  value (`peak.height * max(reference.spectrum)`) are excluded from the
  denoised variability metric (`Denoised_rSD`). Relative to the peak
  rather than an absolute cutoff, so behavior is unchanged whether
  `variant.list` is L-infinity or L2 normalized.

## Value

A numeric matrix with one row per fluorophore and four columns:

- `rSD`:

  Sum, across detectors, of the robust SD (MAD) of the variant spectra.

- `Scaled`:

  `rSD` weighted by the reference spectrum, summed across detectors, to
  downweight low-signal channels.

- `Sim_rSD`:

  Robust SD of the pairwise cosine similarities among all variants
  (including the reference).

- `Denoised_rSD`:

  As `rSD`, but restricted to detectors where the reference spectrum
  exceeds `peak.height`.

Returns `NULL` (with a message) if no fluorophore has usable variability
data.
