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
  method = "Automatic",
  weighted = FALSE,
  weights = NULL,
  af.spectra = NULL,
  spectra.variants = NULL,
  output.dir = NULL,
  file.suffix = NULL,
  include.raw = FALSE,
  include.imaging = FALSE,
  calculate.error = FALSE,
  use.dist0 = TRUE,
  divergence.threshold = 10000,
  divergence.handling = "Balance",
  balance.weight = 0.5,
  speed = "fast",
  parallel = FALSE,
  threads = NULL,
  verbose = TRUE
)
```

## Arguments

- fcs.dir:

  Directory containing FCS files to be unmixed.

- spectra:

  Matrix containing spectra information.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- flow.control:

  A list containing flow cytometry control parameters.

- method:

  A character string specifying the unmixing method to use. The default
  is `Automatic`, which uses `AutoSpectral` for AF extraction if
  af.spectra are provided and automatically selects `OLS` or `WLS`
  depending on which is normal for the given cytometer in
  `asp$cytometer`. This means that files from the ID7000, A8 and S8 will
  be unmixed using `WLS` while others will be unmixed with `OLS`. Any
  option can be set manually. Manual options are `OLS`, `WLS`,
  `AutoSpectral`, `Poisson` and `FastPoisson`. Default is `OLS`.
  `FastPoisson` requires installation of `AutoSpectralRcpp`.

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
  asp\$unmixed.fcs.dir).

- file.suffix:

  A character string to append to the output file name. Default is
  `NULL`

- include.raw:

  Logical indicating whether to include raw data in the written FCS
  file. Default is `FALSE`

- include.imaging:

  Logical indicating whether to include imaging data in the written FCS
  file: relevant for S8 and A8. Default is `FALSE`

- calculate.error:

  Logical, whether to calculate the RMSE unmixing model accuracy and
  include it as a keyword in the FCS file.

- use.dist0:

  Logical, controls whether the selection of the optimal AF signature
  for each cell is determined by which unmixing brings the fluorophore
  signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing
  minimizes the per-cell residual (`use.dist0` = `FALSE`). Default is
  `TRUE`. Used for AutoSpectral unmixing.

- divergence.threshold:

  Numeric. Used for `FastPoisson` only. Threshold to trigger reversion
  towards WLS unmixing when Poisson result diverges. Default is `1e4`

- divergence.handling:

  String. How to handle divergent cells from Poisson IRLS. Options are
  `NonNeg`, in which case non-negativity will be enforced, `WLS`, where
  values will revert to the WLS initial unmix or `Balance`, where `WLS`
  and `NonNeg` results will be averaged. Default is `Balance`

- balance.weight:

  Numeric. Weighting to average non-convergent cells. Used for `Balance`
  option under `divergence.handling`. Default is `0.5`

- speed:

  Selector for the precision-speed trade-off in AutoSpectral per-cell
  fluorophore optimization. Options are the default `fast`, which
  selects the best spectral fit per cell by updating the predicted
  values for each fluorophore independently without repeating the
  unnmixing, `medium` which uses a Woodbury-Sherman-Morrison rank-one
  updating of the unnmixing matrix for better results and a moderate
  slow-down, or `slow`, which explicitly recomputes the unmixing matrix
  for each variant for maximum precision. The `fast` method is only one
  available in the `AutoSpectral` package and will be slow in the pure R
  implementation. Installation of `AutoSpectralRcpp` is strongly
  encouraged.

- parallel:

  Logical, default is `FALSE`. Set to `TRUE` to activate parallel
  processing for multiple FCS files.

- threads:

  Numeric, default is `NULL`, in which case `asp$worker.process.n` will
  be used. `asp$worker.process.n` is set by default to be one less than
  the available cores on the machine. Multi-threading is only used if
  `parallel` is `TRUE`.

- verbose:

  Logical, controls messaging. Default is `TRUE`.

## Value

None. Saves the unmixed FCS files to the specified output directory.
