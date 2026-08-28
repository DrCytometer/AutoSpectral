# Screen a Spectral Matrix for Suspiciously Similar Rows

Internal helper shared by
[`check.spectra.duplicates()`](https://drcytometer.github.io/AutoSpectral/reference/check.spectra.duplicates.md)
(as its fallback when fluorophore identity is unavailable) and
[`read.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/read.spectra.md)
(as an early, load-time check). Not exported.

## Usage

``` r
.check.spectra.collinearity(
  spectra,
  exclude = "AF",
  threshold = 0.95,
  context = ""
)
```
