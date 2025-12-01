# Unmix Using Poisson Regression

This function performs unmixing of raw data using Poisson regression,
with iterative reweighted least squares (IRLS) and fallback methods for
cells that fail to converge.

## Usage

``` r
unmix.poisson(
  raw.data,
  spectra,
  asp,
  initial.weights = NULL,
  parallel = TRUE,
  threads = 1
)
```

## Arguments

- raw.data:

  Matrix containing raw data to be unmixed.

- spectra:

  Matrix containing spectra information.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- initial.weights:

  Optional numeric vector of weights, one per fluorescent detector.
  Default is `NULL`, in which case weighting will be done by channel
  means.

- parallel:

  Logical, default is `TRUE`, which enables parallel processing for
  per-cell unmixing.

- threads:

  Numeric, default is `1`. Multi-threading is only used if `parallel` is
  `TRUE` and `threads` is greater than 1.

## Value

A matrix containing the unmixed data.
