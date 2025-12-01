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
  spectra,
  af.spectra,
  spectral.channel,
  universal.negative,
  control.type,
  raw.thresholds,
  unmixed.thresholds,
  flow.channel,
  som.dim,
  n.cells,
  asp,
  verbose,
  output.dir,
  sim.threshold,
  figures
)
```

## Arguments

- fluor:

  The name of the fluorophore.

- file.name:

  A named vector of file names for the samples.

- control.dir:

  The directory containing the control files.

- spectra:

  A matrix containing the spectral data. Fluorophores in rows, detectors
  in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns. Prepare using
  `get.af.spectra`.

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

- som.dim:

  Numeric, default `10`. Number of x and y dimensions to use in the SOM
  for clustering the spectral variation.

- n.cells:

  Numeric, default `2000`. Number of cells to use for defining the
  variation in spectra. Up to `n.cells` cells will be selected as
  positive events in the peak channel for each fluorophore, above the
  `pos.quantile` in the unstained sample.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

- output.dir:

  File path to whether the figures and .rds data file will be saved.
  Default is `NULL`, in which case `asp$variant.dir` will be used.

- sim.threshold:

  Numeric, default `0.98`. Threshold for cosine similarity- based
  exclusion of spectral variants. Any variant less than `sim.threshold`
  by `cosine.similarity` from the optimized spectrum for that
  fluorophore (from `spectra`) will be excluded from output. This helps
  to exclude autofluorescence contamination.

- figures:

  Logical, controls whether the variation in spectra for each
  fluorophore is plotted in `output.dir`. Default is `TRUE`.

## Value

A matrix with the flow expression data.
