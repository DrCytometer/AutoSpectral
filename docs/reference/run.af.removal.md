# Run Autofluorescence Removal

This function runs the autofluorescence removal process on a list of
samples, using the specified parameters and settings.

## Usage

``` r
run.af.removal(
  clean.expr,
  af.removal.sample,
  spectral.channel,
  peak.channel,
  universal.negative,
  asp,
  scatter.param,
  negative.n = 500,
  positive.n = 1000,
  scatter.match = TRUE,
  k.neighbors = 3L,
  intermediate.figures = FALSE,
  main.figures = TRUE,
  parallel = FALSE,
  threads = 1,
  verbose = TRUE,
  diagnostics.env = NULL
)
```

## Arguments

- clean.expr:

  List containing cleaned expression data.

- af.removal.sample:

  Vector of sample names for which autofluorescence removal is to be
  performed.

- spectral.channel:

  Vector of spectral channel names.

- peak.channel:

  Vector of peak detection channels for fluorophores.

- universal.negative:

  Name of the universal negative control.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- scatter.param:

  Vector of scatter parameters.

- negative.n:

  Integer. Number of events to include in the downsampled negative
  population. Default is `500`.

- positive.n:

  Integer. Number of events to include in the downsampled positive
  population. Default is `1000`.

- scatter.match:

  Logical, default is `TRUE`. Whether to select negative events based on
  scatter profiles matching the positive events.

- k.neighbors:

  Numeric, number of scatter-matched unstained events to pair with every
  positive event for background determination. Default is `3`.

- intermediate.figures:

  Logical, if `TRUE` returns additional figures to show the inner
  workings of the cleaning, including definition of low-AF cell gates on
  the PCA-unmixed unstained and spectral ribbon plots of the AF
  exclusion from the unstained.

- main.figures:

  Logical, if `TRUE` creates the main figures to show the impact of
  intrusive autofluorescent event removal and scatter-matching for the
  negatives.

- parallel:

  Logical, default is `FALSE`, in which case parallel processing will
  not be used. Set to `TRUE` to run in parallel.

- threads:

  Number of cores to use for parallel processing, default is `1`.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

- diagnostics.env:

  Optional environment, default `NULL`. If supplied,
  [`remove.af()`](https://drcytometer.github.io/AutoSpectral/reference/remove.af.md)
  populates it (keyed by sample name) with the objects used to identify
  and exclude intrusive autofluorescence for each cell-based AF-removal
  sample: `af.peak.channel`, `fluor.peak`, `af.boundaries`,
  `expr.data.pos`/`expr.data.neg` (spectral channels only), and the
  resulting gate indices. Intended for diagnostic/manuscript figures
  (see `plot.spectra.legacy.steps()`); has no effect on the cleaning
  result. Capture is unreliable when `parallel = TRUE`.

## Value

A list containing the expression data with autofluorescent events
removed for each sample.
