# Estimate Detector Noise Model

Estimates the two parameters of the detector noise model used by
[`unmix.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.gls.md):
a per-detector additive variance floor (read noise plus any residual
dark-offset contribution) and a per-detector photon conversion factor
`counts.per.unit`, such that \$\$Var(y_d) = read.var_d + \mu_d /
\kappa_d\$\$

Estimation uses the mean-variance relationship of the *residual* after
least-squares projection onto `spectra` (optionally augmented with
`af.spectra`). Residuals are orthogonal to the fitted values by
construction, so binning on the fitted value does not bias the variance
estimate. A degrees-of-freedom correction accounts for the variance
removed by the projection.

Best estimated on bead controls (no autofluorescence, no biological
spread) or on an unstained cell control together with its AF dictionary.
On fully stained samples, unmodelled spectral variation inflates both
parameters; that is not necessarily wrong (it reflects true predictive
uncertainty) but it is no longer a pure instrument measurement.

## Usage

``` r
estimate.noise.model(
  raw.data,
  spectra,
  af.spectra = NULL,
  n.bins = 40L,
  min.bin.n = 50L,
  trim.quantile = 0.999,
  verbose = TRUE,
  af.pc.n = 5L,
  af.raw.data = NULL,
  af.basis.n.cells = 20000L,
  read.var.floor = NULL,
  file.id = NULL,
  unstained.data = NULL,
  dark.quantile = 0.95,
  n.tranche = 10L
)
```

## Arguments

- raw.data:

  Numeric matrix (events x detectors), or a path to an FCS file.
  Detector columns must include all columns of `spectra`.

- spectra:

  Numeric matrix (fluorophores x detectors).

- af.spectra:

  Optional AF dictionary (AF components x detectors), as returned by
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md).
  Strongly recommended for cell controls.

- n.bins:

  Integer, number of quantile bins per detector. Default `40`.

- min.bin.n:

  Integer, minimum events per bin. Default `50`.

- trim.quantile:

  Numeric, upper quantile of each detector discarded before binning, to
  exclude saturation. Default `0.999`.

- verbose:

  Logical. Default `TRUE`.

- af.pc.n:

  Integer, number of principal components to use to describe
  `af.spectra`. Default is `5`.

- af.raw.data:

  Optional numeric matrix, or a path to an FCS file, of individual
  autofluorescence-only events – typically the same unstained sample
  `af.spectra` was built from. When supplied, the AF shape basis is fit
  by PCA on these per-cell events (L-infinity normalised per event,
  matching
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)'s
  own convention) instead of on the rows of `af.spectra`. `af.spectra`
  rows are SOM node centroids, already compressed from the underlying
  cells; if within-node AF spread is substantially real rather than
  noise (check `shape.frac` in `attr(af.spectra, "af.model")`), the
  centroids alone understate true per-cell AF shape diversity, and the
  basis fit here from raw events will need more components to reach the
  same variance-explained target. Ignored when `af.spectra` is `NULL`.
  Default `NULL`, which fits the basis on `af.spectra` as before.

- af.basis.n.cells:

  Integer, default `20000`. Maximum events from `af.raw.data` used for
  the PCA; downsampled if more are supplied. Only used when
  `af.raw.data` is supplied.

- read.var.floor:

  Optional named numeric vector, per-detector additive

- file.id:

  Optional vector (character or factor), length `nrow(raw.data)`,
  identifying which source file/acquisition each row came from. When
  supplied, per-file residual means are subtracted before the
  mean-variance regression, so between-file offsets (gain drift, voltage
  differences across acquisition days) are not counted as photon-driven
  variance. Recommended whenever `raw.data` was pooled across multiple
  control files, e.g. from `pool.scc.for.noise.model()`. Default `NULL`.

- unstained.data:

  Optional numeric matrix (events x detectors) of an unstained/AF-only
  control, used only to anchor the lower tranche of the linearity check
  at `dark.quantile` of the true-dark distribution. When `NULL`, falls
  back to `dark.quantile` of `raw.data` itself, which is a weaker proxy
  since it is not a genuinely dark reference. Default `NULL`.

- dark.quantile:

  Numeric in (0, 1), the quantile of `unstained.data` (or of `raw.data`
  when `unstained.data` is `NULL`) used as the upper bound of the lowest
  linearity tranche. Default `0.95`.

- n.tranche:

  Integer, number of intensity tranches for the linearity check. Default
  `10L`.

## Value

A named list:

- `read.var`:

  Named numeric vector, additive variance floor per detector, in squared
  instrument units.

- `counts.per.unit`:

  Named numeric vector, \\\kappa_d\\.

- `kappa.pooled`:

  Numeric scalar, median of `counts.per.unit`.

- `fit.table`:

  Data frame of per-bin means and variances, for QC plotting.

- `r.squared`:

  Named numeric vector, per-detector regression fit.

- `curvature.coef`:

  Named numeric vector, the quadratic coefficient of
  `variance ~ mean + mean^2` fitted across `n.tranche` intensity
  tranches per detector. Negative indicates variance growing more slowly
  than linear at high signal (compression, e.g. approaching saturation);
  positive indicates faster-than-linear growth.

- `curvature.p`:

  Named numeric vector, the p-value of `curvature.coef`. Values below
  roughly `0.01` indicate the linear mean-variance model is a poor fit
  for that detector's full range.

- `curvature.table`:

  Data frame of per-tranche means, variances, and counts, one block per
  detector, for QC plotting – the same shape as `fit.table` but at
  tranche rather than bin resolution.

- `noise.floor.signal`:

  `sqrt(read.var)`, in *signal* units, for passing to the existing
  `noise.floor` argument of the C++ pipeline. Note the unit difference:
  `read.var` is a variance, the C++ `noise.floor` is a signal level.

- `read.var.source`:

  Named character vector, `"external"` or `"regression"` per detector,
  indicating whether `read.var` came from `read.var.floor` or the fitted
  intercept.

- `kappa.source`:

  Named character vector, `"regression"` or `"pooled.median"` per
  detector, flagging detectors where the fitted slope was non-positive
  and `counts.per.unit` was filled from the pooled median instead of
  measured directly. Treat these detectors' `counts.per.unit`, and any
  ratio computed from it, with caution.
