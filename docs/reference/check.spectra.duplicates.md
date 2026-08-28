# Check Spectra for Duplicate Fluorophores

Verifies that a spectral reference matrix intended for unmixing contains
exactly one row per fluorophore. Multiple controls per fluorophore are
permitted upstream (during control file validation, spectrum extraction,
and QC comparison) via `sample`-based disambiguation, but two or more
rows representing the *same* fluorophore must never reach the unmixing
solvers: near-identical reference columns are highly collinear and will
badly degrade unmixing conditioning. This function is the single
enforcement point that catches that before any unmixing math runs, and
is called automatically from
[`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md),
[`unmix.folder()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.folder.md),
[`unmix.autospectral()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.autospectral.md),
[`unmix.autospectral.joint()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.autospectral.joint.md),
`unmix.autospectral.rcpp()`,
[`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md)
and
[`unmix.af.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.af.gls.md).

If `spectra` carries a `"fluorophore"` attribute (set automatically by
[`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md)
and
[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)),
true fluorophore identity per row is checked exactly, and a duplicate
stops execution. If that attribute is absent – notably, after `spectra`
has been round-tripped through a CSV via
[`read.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/read.spectra.md),
which cannot carry the attribute, and where a `sample`-disambiguated
rowname is unique by construction and so cannot reveal a duplicate
either – this falls back to a pairwise cosine-similarity screen as a
cheap, imperfect surrogate. That screen cannot distinguish "two controls
for the same fluorophore" from "two genuinely similar fluorophores", so
it warns rather than stopping.

## Usage

``` r
check.spectra.duplicates(spectra, collinearity.threshold = 0.95)
```

## Arguments

- spectra:

  A matrix of spectral reference data (fluorophores x detectors), as
  used by
  [`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
  and related functions.

- collinearity.threshold:

  Numeric, default `0.95`. Cosine similarity above which two rows
  trigger the fallback warning, used only when `spectra` has no
  `"fluorophore"` attribute. Matches the package default for
  `asp$similarity.warning.n`.

## Value

Invisibly returns `TRUE` if no duplicates (or, in fallback mode, no
high-similarity pairs) are found. Stops with an informative error when
an exact duplicate is found via the attribute.
