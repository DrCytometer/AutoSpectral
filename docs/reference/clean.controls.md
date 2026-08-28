# Clean Controls

A multi-part function to clean single-color controls in order to extract
fluorophore signatures. Any part can be run independently:

- **Stage 1**: Autofluorescence noise removal using PCA unmixing on
  matching unstained (cells only).

- **Stage 2**: Brightest event selection from positive, universal
  negative from matching negative, and downsampling to speed up RLM
  spectra optimization.

## Usage

``` r
clean.controls(
  flow.control,
  asp,
  af.remove = TRUE,
  universal.negative = TRUE,
  downsample = TRUE,
  drop.raw = FALSE,
  negative.n = asp$negative.n,
  positive.n = asp$positive.n,
  scatter.match = TRUE,
  k.neighbors = 3L,
  intermediate.figures = FALSE,
  main.figures = TRUE,
  parallel = FALSE,
  verbose = TRUE,
  threads = NULL,
  diagnostics.env = NULL,
  ...
)
```

## Arguments

- flow.control:

  A list prepared using `define.flow.control`, containing the data and
  essential information about the cytometer and data structure.

- asp:

  The AutoSpectral parameter list, prepared using
  `get.autospectral.param`.

- af.remove:

  Logical, default is `TRUE`. Whether to remove intrusive
  autofluorescence contamination from cell controls using PCA-based
  identification and gating. Requires universal negatives to be defined
  in the control file and in `flow.control`.

- universal.negative:

  Logical, default is `TRUE`. Whether to use a universal negative sample
  as the negative for spectral extraction. Requires universal negatives
  to be defined in the control file and in `flow.control`.

- downsample:

  Logical, default is `TRUE`. Whether to reduce cell and bead control
  events to speed up processing.

- drop.raw:

  Logical, default `FALSE`. If `TRUE`, removes `flow.control$expr.data`
  (the full, undownsampled per-event data from every control) before
  returning, freeing memory once cleaning is done. Only set this if you
  don't need `get.fluorophore.spectra(use.clean.expr = FALSE)` or other
  functions that read the raw expression data afterward.

- negative.n:

  Integer. Number of events to include in the downsampled negative
  population. Default is `asp$negative.n`.

- positive.n:

  Integer. Number of events to include in the downsampled positive
  population. Default is `asp$positive.n`.

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
  exclusion from the unstained. Default is `FALSE` to speed up
  processing.

- main.figures:

  Logical, if `TRUE` creates the main figures to show the impact of
  intrusive autofluorescent event removal and scatter-matching for the
  negatives.

- parallel:

  Logical, default is `FALSE`, in which case parallel processing will
  not be used. Parallel processing will likely be faster when many small
  files are read in. If the data is larger, parallel processing may not
  accelerate the process much.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

- threads:

  Numeric, number of threads to use for parallel processing. Default is
  `NULL` which will revert to `asp$worker.process.n` if `parallel=TRUE`.

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

- ...:

  Ignored. Used to catch deprecated arguments.

## Value

Returns a modified `flow.control` with the original data intact. New,
cleaned data and corresponding factors are stored in new slots.
