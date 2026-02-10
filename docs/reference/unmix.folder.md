# Unmix All FCS Files in a Directory

This function unmixes all FCS files in a specified directory using the
provided spectra and method, and saves the unmixed FCS files to an
output directory of the user's choice.

## Usage

``` r
unmix.folder(
  fcs.dir,
  spectra,
  asp,
  flow.control,
  method = "AutoSpectral",
  weighted = FALSE,
  weights = NULL,
  af.spectra = NULL,
  spectra.variants = NULL,
  output.dir = NULL,
  file.suffix = NULL,
  include.raw = FALSE,
  include.imaging = FALSE,
  use.dist0 = TRUE,
  divergence.threshold = 10000,
  divergence.handling = "Balance",
  balance.weight = 0.5,
  speed = c("slow", "medium", "fast"),
  parallel = FALSE,
  threads = NULL,
  verbose = TRUE,
  k = NULL,
  ...
)
```

## Arguments

- fcs.dir:

  Directory (file path) containing FCS files to be unmixed.

- spectra:

  A matrix containing the spectral data. Fluorophores in rows, detectors
  in columns.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

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
  `AutoSpectralRcpp`.There is also `Automatic`, which switches depending
  on the inputs provided: it uses `AutoSpectral` for AF extraction if
  `af.spectra` are provided, and automatically selects `OLS` or `WLS`
  depending on which is normal for the given cytometer in
  `asp$cytometer`. This means that files from the ID7000, A8 and S8 will
  be unmixed using `WLS` while others will be unmixed with `OLS`, if
  AutoSpectral unmixing is not activated.

- weighted:

  Logical, whether to use ordinary or weighted least squares unmixing as
  the base algorithm in AutoSpectral unmixing. Default is `FALSE` and
  will use OLS.

- weights:

  Optional numeric vector of weights: one per fluorescent detector.
  Default is `NULL`, in which case weighting will be done by channel
  means. Only used for `WLS`

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

  Directory to save the unmixed FCS files (default is
  `asp$unmixed.fcs.dir`, which is `./AutoSpectral_unmixed`).

- file.suffix:

  A character string to append to the output file name. Default is
  `NULL`

- include.raw:

  A logical value indicating whether to include raw expression data in
  the written FCS file. Default is `FALSE` to provide smaller output
  files.

- include.imaging:

  A logical value indicating whether to include imaging parameters in
  the written FCS file. Default is `FALSE` to provide smaller output
  files.

- use.dist0:

  Logical, controls whether the selection of the optimal AF signature
  for each cell is determined by which unmixing brings the fluorophore
  signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing
  minimizes the per-cell residual (`use.dist0` = `FALSE`). Default is
  `TRUE`. Used for AutoSpectral unmixing. The minimization of
  fluorophore signals can be thought of as a "worst-case" scenario, but
  it provides more accurate assignments, particularly with large panels.

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
  with the default being `slow`. As of version 1.0.0, the backend for
  how this works has changed. Spectral variants and AF signatures are
  now pre-screened per cell to identify likely candidates, so brute
  force testing of all variants is no longer required. So, `speed`
  controls the number of variants to be tested per cell, with `fast`
  testing a single variant, `medium` testing 3 variants, and `slow`
  testing 10 variants. While this is now implemented in pure R in
  `AutoSpectral`, installation of `AutoSpectralRcpp` is strongly
  encouraged for faster processing.

- parallel:

  Logical, default is `FALSE`. Set to `TRUE` to activate parallel
  processing for multiple FCS files.

- threads:

  Numeric, default is `NULL`, in which case `asp$worker.process.n` will
  be used. `asp$worker.process.n` is set by default to be one less than
  the available cores on the machine. Multi-threading is only used if
  `parallel` is `TRUE`. If working on a computing cluster, try
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html).

- verbose:

  Logical, controls messaging. Default is `TRUE`. Set to `FALSE` to have
  it shut up.

- k:

  Number of variants (and autofluorescence spectra) to test per cell.
  Allows explicit control over the number used, as opposed to `speed`,
  which selects from pre-defined choices. Providing a numeric value to
  `k` will override `speed`, allowing up to `k` (or the max available)
  variants to be tested. The default is `NULL`, in which case `k` will
  be ignored.

- ...:

  Ignored. Previously used for deprecated arguments such as
  `calculate.error`.

## Value

None. Saves the unmixed FCS files to the specified output directory.
