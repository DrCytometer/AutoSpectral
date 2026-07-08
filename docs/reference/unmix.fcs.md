# Unmix FCS Data

This function performs spectral unmixing on FCS data using various
methods.

## Usage

``` r
unmix.fcs(
  fcs.file,
  spectra,
  asp,
  flow.control,
  method = c("AutoSpectral", "OLS", "WLS", "Poisson", "FastPoisson"),
  weighted = FALSE,
  weights = NULL,
  af.spectra = NULL,
  spectra.variants = NULL,
  output.dir = NULL,
  file.suffix = NULL,
  include.raw = FALSE,
  include.imaging = TRUE,
  use.dist0 = TRUE,
  divergence.threshold = 10000,
  divergence.handling = "Balance",
  balance.weight = 0.5,
  speed = c("fast", "medium", "slow"),
  parallel = TRUE,
  threads = if (parallel) 0 else 1,
  verbose = TRUE,
  n.variants = NULL,
  chunk.size = 2e+06,
  pipeline = c("legacy", "joint"),
  n.passes = 1L,
  n.af.passes = 1L,
  cell.weight = if (asp$cytometer == "ID7000") TRUE else FALSE,
  noise.floor = 125,
  alpha = 0.5,
  collinear.threshold = 0.5,
  joint.pair.resolution = TRUE,
  refine.af.quantile = 0.5,
  ...
)
```

## Arguments

- fcs.file:

  A character string specifying the path to the FCS file.

- spectra:

  A matrix containing the spectral data. Fluorophores in rows, detectors
  in columns.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`.

- flow.control:

  A list containing flow cytometry control parameters.

- method:

  A character string specifying the unmixing method to use. The default
  as of version 1.0.0 is now `AutoSpectral` to avoid confusion. To use
  AutoSpectral unmixing, you must provide at least `af.spectra` to
  perform autofluorescence extraction (on a per-cell basis). To also
  optimize fluorophore spectra, provide `spectra.variants`. To perform
  other types of unmixing, select from the options: `OLS`, `WLS`,
  `Poisson` or `FastPoisson`. `FastPoisson` requires installation of
  `AutoSpectralRcpp`.

- weighted:

  Logical, whether to use ordinary or weighted least squares unmixing as
  the base algorithm in AutoSpectral legacy pipeline unmixing. Default
  is `FALSE` and will use OLS.

- weights:

  Optional numeric vector of weights (one per fluorescent detector).
  Default is `NULL`, in which case weighting will be done by channel
  means (Poisson variance). Only used for `WLS`.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`. Required for `AutoSpectral` unmixing. Default is
  `NULL` and will thus provoke failure if no spectra are provided and
  `AutoSpectral` is selected.

- spectra.variants:

  Named list (names are fluorophores) carrying matrices of spectral
  signature variations for each fluorophore. Prepare using
  `get.spectral.variants`. Default is `NULL`. Used for AutoSpectral
  unmixing. Required for per-cell fluorophore optimization.

- output.dir:

  A character string specifying the directory to save the unmixed FCS
  file. Default is `NULL`, which will use `./AutoSpectral_unmixed`.

- file.suffix:

  A character string to append to the output file name. Default is
  `NULL`.

- include.raw:

  A logical value indicating whether to include raw expression data in
  the written FCS file. Default is `FALSE` to provide smaller output
  files.

- include.imaging:

  A logical value indicating whether to include imaging parameters in
  the written FCS file. Default is `TRUE`.

- use.dist0:

  Legacy pipeline argument. Logical, controls whether the selection of
  the optimal AF signature for each cell is determined by the
  minimization of potential AF spillover into the fluorophore channels
  (`use.dist0` = `TRUE`) or by which unmixing minimizes the per-cell
  residual (`use.dist0` = `FALSE`). Default is `TRUE`. Used for legacy
  AutoSpectral unmixing.

- divergence.threshold:

  Numeric. Used for `FastPoisson` only. Threshold to trigger reversion
  towards WLS unmixing when Poisson result diverges for a given point.
  To be deprecated.

- divergence.handling:

  String. How to handle divergent cells from Poisson IRLS. Options are
  `NonNeg` (non-negativity will be enforced), `WLS` (revert to WLS
  initial unmix) or `Balance` (WLS and NonNeg will be averaged). Default
  is `Balance`. To be deprecated.

- balance.weight:

  Numeric. Weighting to average non-convergent cells. Used for `Balance`
  option under `divergence.handling`. Default is `0.5`. To be
  deprecated.

- speed:

  Selector for the precision-speed trade-off in AutoSpectral per-cell
  fluorophore optimization. Options are `fast`, `medium` and `slow`,
  with the default being `fast`. As of version 1.0.0, the backend for
  how this works has changed. Spectral variants and AF signatures are
  now pre-screened per cell to identify likely candidates, so brute
  force testing of all variants is no longer required. So, `speed`
  controls the number of variants to be tested per cell, with `fast`
  testing a single variant, `medium` testing 3 variants, and `slow`
  testing 10 variants. While this is now implemented in pure R in
  `AutoSpectral`, installation of `AutoSpectralRcpp` is strongly
  encouraged for faster processing.

- parallel:

  Logical, default is `TRUE`, which enables parallel processing for
  per-cell unmixing methods.

- threads:

  Numeric, defaults to a single thread for sequential processing
  (`parallel=FALSE`) or all available cores if `parallel=TRUE`.

- verbose:

  Logical, controls messaging. Default is `TRUE`. Set to `FALSE` to have
  it shut up.

- n.variants:

  Numeric, used for legacy AutoSpectral pipeline unmixing. Number of
  variants to test per cell. Allows explicit control over the number
  used, as opposed to `speed`, which selects from pre-defined choices.
  Providing a numeric value to `n.variants` will override `speed`,
  allowing up to `n.variants` (or the max available) variants to be
  tested. The default is `NULL`, in which case `n.variants` will be
  ignored.

- chunk.size:

  Numeric, number of events to use per chunk of unmixing. Used to manage
  memory when processing large FCS files. As a rough guide, you will
  need approximately 10x the size of the raw FCS file on disk as
  available memory. Default is set at `2e6` events, assuming ~20GB
  memory available.

- pipeline:

  Character, one of `"legacy"` (default) or `"joint"`. Passed to
  `unmix.autospectral.rcpp()`. `"joint"` uses the new
  covariance-weighted joint per-cell pipeline; `"legacy"` reproduces the
  behaviour of AutoSpectral prior to version 1.6.0.

- n.passes:

  Integer, default `1L`. Number of joint optimisation passes per cell.
  Only used when `pipeline = "joint"`. Set higher for some improvement
  in unmixing with high spillover datasets.

- n.af.passes:

  Integer, default `1L`. Number of autofluorescence extraction passes
  per cell. Only used when `pipeline = "joint"`. Passed to
  `unmix.autospectral.rcpp()`.

- cell.weight:

  Logical, default `FALSE`. Applies per-cell detector weighting to the
  joint unmixing solve. Only used when `pipeline = "joint"`. Passed to
  `unmix.autospectral.rcpp()`.

- noise.floor:

  Numeric, default `125`. Lower clamp on the denominator of the per-cell
  detector weights when `cell.weight = TRUE`. Only used when
  `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`.

- alpha:

  Numeric, default `0.5`. Weighting for balancing residual and
  covariance spillover minimization. Only used when
  `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`.

- collinear.threshold:

  Numeric, default `0.5`. Cosine similarity value to trigger conflict
  assessment for collinear fluorophore variants. Only used when
  `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`.

- joint.pair.resolution:

  Logical, default `TRUE`. Whether to perform conflict-resolution for
  collinear fluorophore pairs. Only used when `pipeline = "joint"`.
  Passed to `unmix.autospectral.rcpp()`.

- refine.af.quantile:

  Numeric, default `0.5`. Fraction of cells taken forward for additional
  AF passes (see `n.af.passes`). Only used when `pipeline = "joint"`.
  Passed to `unmix.autospectral.rcpp()`.

- ...:

  Ignored. Used to catch deprecated arguments.

## Value

None. The function writes the unmixed FCS data to a file.
