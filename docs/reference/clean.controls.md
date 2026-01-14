# Clean Controls

A four-part function to clean single-color controls in order to extract
fluorophore signatures. Any part can be run independently:

- **Stage 1**: PeacoQC to eliminate flow artefacts. Not required in most
  cases.

- **Stage 2**: Trimming to eliminate extreme events. Not recommended for
  most use cases.

- **Stage 3**: Autofluorescence noise removal using PCA unmixing on
  matching unstained (cells only).

- **Stage 4**: Brightest event selection from positive, universal
  negative from matching negative, and downsampling to speed up RLM
  spectra optimization.

## Usage

``` r
clean.controls(
  flow.control,
  asp,
  time.clean = FALSE,
  trim = FALSE,
  trim.factor = NULL,
  af.remove = TRUE,
  universal.negative = TRUE,
  downsample = TRUE,
  negative.n = asp$negative.n,
  positive.n = asp$positive.n,
  scatter.match = TRUE,
  scrub = FALSE,
  intermediate.figures = FALSE,
  main.figures = TRUE,
  parallel = FALSE,
  verbose = TRUE,
  threads = NULL
)
```

## Arguments

- flow.control:

  A list prepared using `define.flow.control`, containing the data and
  essential information about the cytometer and data structure.

- asp:

  The AutoSpectral parameter list, prepared using
  `get.autospectral.param`.

- time.clean:

  Logical, default is `FALSE`. Whether to run PeacoQC to remove
  time-based inconsistencies in the controls.

- trim:

  Logical, default is `FALSE`. Whether to remove extreme events
  (positive and negative) from controls.

- trim.factor:

  Numeric. Default is `asp$rlm.trim.factor`. Required if `trim = TRUE`.

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

- negative.n:

  Integer. Number of events to include in the downsampled negative
  population. Default is `asp$negative.n`.

- positive.n:

  Integer. Number of events to include in the downsampled positive
  population. Default is `asp$positive.n`.

- scatter.match:

  Logical, default is `TRUE`. Whether to select negative events based on
  scatter profiles matching the positive events. Defines a region of FSC
  and SSC based on the distribution of selected positive events.

- scrub:

  Logical, if `TRUE` allows for re-cleaning of already cleaned data,
  provided there are clean data in `flow.control`.

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

## Value

Returns a modified `flow.control` with the original data intact. New,
cleaned data and corresponding factors are stored in new slots.
