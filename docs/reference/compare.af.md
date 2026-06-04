# Compare Autofluorescence Spectra Sets

Evaluates the quality of multiple AF spectra matrices – for example,
produced by different parameter settings of
[`get.af.spectra`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
– against the same unstained FCS file. A single AF assignment function
is applied consistently across every entry in `af.spectra.list`, so that
differences in cosine similarity reflect the AF spectra themselves
rather than the assignment strategy.

## Usage

``` r
compare.af(
  unstained.fcs,
  spectra,
  af.spectra.list,
  assign.fn = "assign.af.fluorophores",
  n.downsample = 1000L,
  plot.dir = "figure_af_accuracy",
  title = "compare_af"
)
```

## Arguments

- unstained.fcs:

  Character scalar. Path to the unstained FCS file.

- spectra:

  Numeric matrix of fluorophore spectra (fluorophores x detectors). Row
  names must be fluorophore names; column names must match the detector
  channels in the FCS file. Any row named `"AF"` is removed
  automatically.

- af.spectra.list:

  Named list of AF spectra matrices. Each element must be a numeric
  matrix with columns matching `colnames(spectra)`. The first row of
  each matrix is treated as that candidate's mean AF spectrum. Names are
  used as labels throughout; if the list is unnamed, entries are
  labelled `"af1"`, `"af2"`, etc.

- assign.fn:

  Character scalar. Name of the assign-type function used to map each
  cell to an AF spectrum row. Must be available in the current search
  path and follow the assign-type calling convention:
  `fn(raw.data, spectra, af.spectra)` returning an integer vector of row
  indices into `af.spectra`. Default: `"assign.af.fluorophores"`.

- n.downsample:

  Integer scalar. Maximum number of events used from the FCS file. A
  random subsample is drawn when the file contains more events. Set to
  `Inf` to use all events. Default: `1000L`.

- plot.dir:

  Character scalar. Directory in which to save the summary plot PDF.
  Created recursively if it does not exist. Set to `NULL` to skip
  saving. Default: `"figure_af_accuracy"`.

- title:

  Character scalar. Stem used to name the output PDF
  (`<plot.dir>/<title>.pdf`). Default: `"compare_af"`.

## Value

A named list with one entry per candidate (plus `"baseline"`). Each
entry contains:

- `Assignments`:

  Integer vector of per-cell AF-spectrum row indices.

- `Similarity`:

  Numeric vector of per-cell cosine similarities between the raw
  detector signal and the assigned AF spectrum.

- `Mean_Sim`:

  Mean of `Similarity` (NAs excluded).

- `rSD_Sim`:

  Standard deviation of `Similarity` (NAs excluded).

- `n.variants`:

  Number of AF spectrum rows (variants) in this candidate's matrix.

## Details

For each `af.spectra` matrix the function:

1.  Assigns each cell to its best-matching AF spectrum row using
    `assign.fn`.

2.  Computes the cosine similarity between the cell's raw detector
    signal and its assigned AF spectrum.

3.  Summarises per-cell similarities into `Mean_Sim` and `rSD_Sim`.

A grand baseline using the first row of the first list entry (i.e. the
mean AF spectrum of the first candidate) is always prepended so every
plot has a common anchor.

Results are returned as a list and, optionally, as a summary bar chart
saved to `plot.dir`.

## See also

[`get.af.spectra`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md),
[`test.af.accuracy`](https://drcytometer.github.io/AutoSpectral/reference/test.af.accuracy.md),
[`assign.af.fluorophores`](https://drcytometer.github.io/AutoSpectral/reference/assign.af.fluorophores.md)
