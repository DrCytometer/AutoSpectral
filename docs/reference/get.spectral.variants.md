# Get Spectral Variations for Fluorophores

This function cycles through all the fluorophores defined in
`control.def.file`, identifying variations in spectral profiles. It does
this by performing SOM clustering on the positive events in the cleaned
control data. The output is saved as an .rds file, and figures
summarizing the variation are saved, if desired. Note that the .rds file
contains all the needed information for downstream processing (per-cell
unmixing), so you can just load that using the `readRDS` function)
rather than re-running this process.

## Usage

``` r
get.spectral.variants(
  control.dir,
  control.def.file,
  asp,
  spectra,
  af.spectra,
  n.cells = 2000,
  pos.quantile = 0.995,
  som.dim = 10,
  sim.threshold = 0.98,
  figures = TRUE,
  output.dir = NULL,
  parallel = FALSE,
  verbose = TRUE
)
```

## Arguments

- control.dir:

  File path to the single-stained control FCS files.

- control.def.file:

  CSV file defining the single-color control file names, fluorophores
  they represent, marker names, peak channels, and gating requirements.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- spectra:

  A matrix containing the spectral data. Fluorophores in rows, detectors
  in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`.

- n.cells:

  Numeric, default `2000`. Number of cells to use for defining the
  variation in spectra. Up to `n.cells` cells will be selected as
  positive events in the peak channel for each fluorophore, above the
  `pos.quantile` in the unstained sample.

- pos.quantile:

  Numeric, default `0.995`. Threshold for positivity. This quantile will
  be used to define the maximum extent of the unstained sample in the
  unmixed space. Anything above that will be considered positive in the
  single-stained controls.

- som.dim:

  Numeric, default `10`. Number of x and y dimensions to use in the SOM
  for clustering the spectral variation.

- sim.threshold:

  Numeric, default `0.98`. Threshold for cosine similarity- based
  exclusion of spectral variants. Any variant less than `sim.threshold`
  by `cosine.similarity` from the optimized spectrum for that
  fluorophore (from `spectra`) will be excluded from output. This helps
  to exclude autofluorescence contamination.

- figures:

  Logical, controls whether the variation in spectra for each
  fluorophore is plotted in `output.dir`. Default is `TRUE`.

- output.dir:

  File path to whether the figures and .rds data file will be saved.
  Default is `NULL`, in which case `asp$variant.dir` will be used.

- parallel:

  Logical, default is `FALSE`, in which case sequential processing will
  be used. Parallel processing will likely be faster when many small
  files are read in. If the data is larger, parallel processing may not
  accelerate the process much or may fail outright.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A vector with the indexes of events inside the initial gate.
