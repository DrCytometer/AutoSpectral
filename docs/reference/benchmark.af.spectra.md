# Benchmark AF Assignment Accuracy Against Spectral Panel Size

Repeatedly subsamples the fluorophore spectra matrix to a range of panel
sizes, runs
[`test.af.accuracy`](https://drcytometer.github.io/AutoSpectral/reference/test.af.accuracy.md)
for each subsample, and summarises how cosine similarity between raw
detector signals and the assigned AF spectra changes with the number of
fluorophores. Results are saved as a line-plot PDF and returned as a
list.

## Usage

``` r
benchmark.af.spectra(
  unstained.fcs,
  spectra,
  af.spectra,
  asp,
  functions = c("assign.af.fluorophores", "assign.af.residuals", "assign.af.joint.cov"),
  n.fluors = c(5, 10, 15, 20, 25, 30, 35, 40),
  n.draws = 5L,
  seed = 42L,
  n.downsample = 1000L,
  plot.dir = "figure_af_accuracy",
  filename = "af_accuracy_spectra_benchmark"
)
```

## Arguments

- unstained.fcs:

  Character scalar. Path to the unstained FCS file passed through to
  [`test.af.accuracy`](https://drcytometer.github.io/AutoSpectral/reference/test.af.accuracy.md).

- spectra:

  Numeric matrix of fluorophore spectra (fluorophores x detectors). Row
  names must be fluorophore names; column names must match the detector
  channels in the FCS file. Any row named `"AF"` is removed
  automatically before processing. Column names must match those of
  `spectra`.

- af.spectra:

  Numeric matrix of AF spectra (AF variants x detectors). The first row
  is treated as the mean AF spectrum and is used for the baseline
  comparison. Column names must match those of `spectra`.

- asp:

  Aspect-ratio value passed through to
  [`test.af.accuracy`](https://drcytometer.github.io/AutoSpectral/reference/test.af.accuracy.md).

- functions:

  Character vector of AF-function names to benchmark. See
  [`test.af.accuracy`](https://drcytometer.github.io/AutoSpectral/reference/test.af.accuracy.md)
  for the required calling conventions.

- n.fluors:

  Integer vector of panel sizes (number of fluorophores) to evaluate.
  Values outside `[1, nrow(spectra)]` are silently dropped. Default:
  `c(5, 10, 15, 20, 25, 30, 35, 40)`.

- n.draws:

  Integer scalar. Number of random draws (unique fluorophore subsets)
  per panel size. A unique seed is derived for each `(n.fluors, draw)`
  combination from `seed`, so results are fully reproducible. Default:
  `5L`.

- seed:

  Integer scalar. Base seed for reproducible subsampling. Default:
  `42L`.

- n.downsample:

  Integer scalar. Maximum number of events passed to each
  [`test.af.accuracy`](https://drcytometer.github.io/AutoSpectral/reference/test.af.accuracy.md)
  call. A random subsample is drawn when the FCS file contains more
  events. Set to `Inf` to use all events. Default: `1000L`.

- plot.dir:

  Character scalar. Directory for output files. Created recursively if
  it does not exist. Default: `"figure_af_accuracy"`.

- filename:

  Character scalar. Stem used to name the summary PDF (the file will be
  `<plot.dir>/<filename>.pdf`). Default:
  `"af_accuracy_spectra_benchmark"`.

## Value

A named list with three elements:

- `raw`:

  Data frame of per-draw results with columns `n.fluors`, `draw`,
  `af_method`, `Mean_Sim`, and `SD_Sim`.

- `summary`:

  Data frame summarising `raw` across draws for each
  `(n.fluors, method)` combination, with columns `n.fluors`, `method`,
  `mean_Mean_Sim`, `sd_Mean_Sim` (variability across draws), and
  `mean_SD_Sim` (average within-draw spread).

- `plot`:

  The `ggplot2` object for the summary line plot.

## See also

[`test.af.accuracy`](https://drcytometer.github.io/AutoSpectral/reference/test.af.accuracy.md)
