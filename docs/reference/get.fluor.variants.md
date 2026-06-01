# Get Fluorophore Variants

Assesses variation in the spectral signature of a single-stained flow
cytometry control sample using SOM clustering on scatter-matched,
background-corrected positive events.

Autofluorescence is characterised **in situ** from the paired
universal-negative file specified in the control table (or from the
lower 25\\ available). The AF mean vector is unit-normalised and
projected out of each event to identify the empirical peak detector,
mirroring the approach in `get.spectra.automated`. All positive events
above the raw threshold in the (empirical) peak channel are selected, up
to `n.cells` (randomly downsampled when more are present). For each
selected event, \\k\\ scatter-space nearest neighbours are found in the
unstained pool and their spectral values are averaged to form a
per-event background estimate, which is then subtracted. SOM clustering
on the resulting background-corrected matrix recovers the
population-level distribution of spectral shapes. A cosine-similarity QC
step retains only SOM centroids sufficiently similar to the reference
spectrum, followed by off-peak smoothing.

## Usage

``` r
get.fluor.variants(
  fluor,
  file.name,
  control.dir,
  asp,
  spectra,
  figures,
  output.dir,
  verbose,
  spectral.channel,
  scatter.channel,
  universal.negative,
  control.type,
  raw.thresholds,
  unmixed.thresholds,
  flow.channel,
  af.pcs,
  n.cells = 10000L,
  som.dim = 10L,
  k.neighbors = 3L,
  sim.threshold = 0.985,
  variant.fill.color = "red",
  variant.fill.alpha = 0.7,
  median.line.color = "black",
  median.linewidth = 1
)
```

## Arguments

- fluor:

  Character. Name of the fluorophore.

- file.name:

  Named character vector of control FCS filenames, named by fluorophore.

- control.dir:

  Character. Directory containing the control FCS files.

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- spectra:

  Numeric matrix. Reference spectra; fluorophores in rows, detectors in
  columns.

- figures:

  Logical. Whether to save a spectral-variant plot. Default `TRUE`.

- output.dir:

  Character. Directory for figures.

- verbose:

  Logical. Whether to print progress messages. Default `TRUE`.

- spectral.channel:

  Character vector of spectral detector channel names.

- scatter.channel:

  Character vector of scatter parameter names (e.g. `"FSC-A"`,
  `"SSC-A"`) used for KNN scatter matching against the unstained pool.

- universal.negative:

  Named character vector mapping fluorophore names to their paired
  unstained FCS filename, or `"FALSE"` / `NA` when none is available.

- control.type:

  Character, either "beads" or "cells". Determines the type of control
  sample being used and the subsequent processing steps.

- raw.thresholds:

  Named numeric vector of per-channel positivity thresholds (typically
  the 99.5th percentile of the unstained sample).

- unmixed.thresholds:

  A named vector of numerical values corresponding to the threshold for
  positivity in each unmixed channel. Determined by the 99.5th
  percentile on the unstained sample, typically after single-cell AF
  unmixing.

- flow.channel:

  Named character vector of expected peak raw channels, one per
  fluorophore.

- af.pcs:

  Matrix of autofluorescence-defining principal components.

- n.cells:

  Integer, default `10000`. Maximum number of positive events used for
  SOM clustering. Files with more events above threshold are randomly
  downsampled to this number.

- som.dim:

  Integer, default `10`. Side length of the square SOM grid. Produces up
  to `som.dim^2` candidate variant spectra before cosine QC.

- k.neighbors:

  Integer, default `3`. Number of scatter-space nearest neighbours from
  the unstained pool used to form the per-event background estimate.

- sim.threshold:

  Numeric, default `0.985`. Minimum cosine similarity between a SOM
  centroid and the reference spectrum for the centroid to be retained as
  a variant.

- variant.fill.color:

  Color for the shaded ribbon in the variant plot. Default `"red"`.

- variant.fill.alpha:

  Alpha for `variant.fill.color`. Default `0.7`.

- median.line.color:

  Color for the reference-spectrum line. Default `"black"`.

- median.linewidth:

  Width of the reference-spectrum line. Default `1`.

## Value

A numeric matrix; variants in rows, detectors in columns, values
normalised to \\\[0, 1\]\\. When no centroids survive cosine QC the
single reference spectrum is returned (one row).

## References

Van Gassen S et al. (2015). FlowSOM. *Cytometry Part A*, 87(7), 636-645.
[doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625)
