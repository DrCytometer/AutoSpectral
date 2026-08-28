# Correct Reference Spectra for Systematic Variant Offset

Shifts each fluorophore's reference spectrum by the mean delta captured
in `variant.basis` (see
[`build.variant.basis()`](https://drcytometer.github.io/AutoSpectral/reference/build.variant.basis.md)),
then re-normalises to L-infinity = 1. A non-zero mean delta means the
control-derived reference is systematically off-centre relative to the
variant population – a bias in `spectra`, not a source of variance – so
folding it into the covariance would inflate uncertainty while leaving
the bias itself uncorrected. This applies the correction once, up front,
so
[`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md)
sees a centred reference and `variant.basis` describes only the residual
spread around it.

The `mean.delta` entries in the returned `variant.basis` are zeroed for
every fluorophore corrected, since the mean is now baked into `spectra`;
calling this function again on its own output is therefore a no-op.

## Usage

``` r
correct.spectra.mean.delta(spectra, variant.basis, fluorophores = NULL)
```

## Arguments

- spectra:

  Numeric matrix (fluorophores x detectors).

- variant.basis:

  List from
  [`build.variant.basis()`](https://drcytometer.github.io/AutoSpectral/reference/build.variant.basis.md).

- fluorophores:

  Optional character vector restricting the correction to named rows of
  `spectra`. Default `NULL` (all rows present in both `spectra` and
  `variant.basis`).

## Value

A list with `spectra` (corrected, L-infinity re-normalised) and
`variant.basis` (mean.delta zeroed for corrected fluorophores).
