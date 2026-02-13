# Get Fluorophore Variants

Assesses variation in the spectral signature of a single-stained flow
cytometry control sample. Uses SOM-based clustering on the brightest
positive events in the file.

## Usage

``` r
get.fluor.variants(
  fluor,
  file.name,
  control.dir,
  asp,
  spectra,
  af.spectra,
  n.cells,
  som.dim,
  figures,
  output.dir,
  verbose,
  spectral.channel,
  universal.negative,
  control.type,
  raw.thresholds,
  unmixed.thresholds,
  flow.channel,
  refine = TRUE,
  problem.quantile = 0.95
)
```

## Arguments

- fluor:

  The name of the fluorophore.

- file.name:

  A named vector of file names for the samples.

- control.dir:

  The directory containing the control files.

- asp:

  The AutoSpectral parameter list.

- spectra:

  A matrix containing the spectral data. Fluorophores in rows, detectors
  in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`.

- n.cells:

  Numeric. Number of cells to use for defining the variation in spectra.
  Up to `n.cells` cells will be selected as positive events in the peak
  channel for each fluorophore, above the 99.5th percentile level in the
  unstained sample.

- som.dim:

  Numeric. Number of x and y dimensions to use in the SOM for clustering
  the spectral variation.

- figures:

  Logical, controls whether the variation in spectra for each
  fluorophore is plotted in `output.dir`. Default is `TRUE`.

- output.dir:

  File path to whether the figures and .rds data file will be saved.
  Default is `NULL`, in which case `asp$variant.dir` will be used.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

- spectral.channel:

  A vector of spectral channels.

- universal.negative:

  A named vector of unstained negative samples, with names corresponding
  to the fluorophores.

- control.type:

  Character, either "beads" or "cells". Determines the type of control
  sample being used and the subsequent processing steps.

- raw.thresholds:

  A named vector of numerical values corresponding to the threshold for
  positivity in each raw detector channel. Determined by the 99.5th
  percentile on the unstained sample, typically.

- unmixed.thresholds:

  A named vector of numerical values corresponding to the threshold for
  positivity in each unmixed channel. Determined by the 99.5th
  percentile on the unstained sample, typically after single-cell AF
  unmixing.

- flow.channel:

  A named vector of peak raw channels, one per fluorophore.

- refine:

  Logical, default is `TRUE`. Controls whether to perform a second round
  of variation measurement on "problem cells", which are those with the
  highest spillover, as defined by `problem.quantile`.

- problem.quantile:

  Numeric, default `0.95`. The quantile for determining which cells will
  be considered "problematic" after unmixing with per-cell AF
  extraction. Cells in the `problem.quantile` or above with respect to
  total signal in the fluorophore (non-AF) channels after per-cell AF
  extraction will be used to determine additional autofluorescence
  spectra, using a second round of clustering and modulation of the
  previously selected autofluorescence spectra. A value of `0.95` means
  the top 5% of cells, those farthest from zero, will be selected for
  further investigation.

## Value

A matrix with the flow expression data.
