# Build Variant Covariance Basis

Converts the per-fluorophore variant deltas from
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
into a low-rank covariance basis for
[`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md).
Each fluorophore's contribution to the detector-space covariance is
\\x_f^2 \sum_j \lambda\_{fj} b\_{fj} b\_{fj}^\top\\.

The *mean* delta is returned separately rather than folded into the
covariance: a systematic offset from the reference spectrum is a
correction to `spectra`, not a source of variance, and treating it as
variance both inflates the noise model and leaves the bias uncorrected.

When `spectra` is supplied, the centred deltas are projected onto the
orthogonal complement of the panel span before the SVD. Only the
panel-oblique component of spectral variation is visible to unmixing;
the in-span component is collinear with the fluorophore columns, and
carrying it in the covariance inflates the variance of the abundance
estimates along exactly the spillover directions. Supplying `spectra` is
strongly recommended whenever the basis will be used by
[`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md).

## Usage

``` r
build.variant.basis(
  spectra.variants,
  spectra = NULL,
  rank = 6L,
  var.explained = 0.99,
  min.lambda = 1e-04,
  pooled.fallback = TRUE,
  verbose = TRUE
)
```

## Arguments

- spectra.variants:

  Either the full list returned by
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
  (containing `$delta.list`), or a named list of delta matrices
  (variants x detectors).

- spectra:

  Optional numeric matrix (fluorophores x detectors), the panel
  reference spectra. When supplied, variant deltas are projected onto
  the orthogonal complement of the panel span before fitting, so the
  basis captures only panel-oblique variation. Default `NULL`.

- rank:

  Integer, maximum number of components retained per fluorophore.
  Default `6`.

- var.explained:

  Numeric in (0, 1\]. Components are dropped once this fraction of the
  delta variance is captured. Default `0.99`.

- min.lambda:

  Numeric, eigenvalues below this fraction of the leading eigenvalue are
  dropped. Default `1e-4`.

- pooled.fallback:

  Logical, if `TRUE`, uses other fluorophores' variability to infer
  likely variability for a fluorophore where variants where not
  sucessfully measured.

- verbose:

  Logical, whether to output messages to the console.

## Value

A named list, one entry per fluorophore with at least two variants, each
containing `basis` (components x detectors), `lambda` (numeric, one per
component), and `mean.delta` (numeric, length detectors). Fluorophores
with insufficient variants are omitted.
