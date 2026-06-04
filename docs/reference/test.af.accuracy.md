# Test Autofluorescence Assignment Accuracy

Benchmarks autofluorescence (AF) assignment and unmixing functions
against an unstained FCS file. For each function supplied, cells are
assigned to an AF spectrum, unmixed, and evaluated by cosine similarity
between the raw detector signal and the assigned AF spectrum. A mean-AF
baseline (every cell assigned to the mean AF spectrum) is always
prepended to the results. Biplot panels for each method are saved as a
single PDF.

## Usage

``` r
test.af.accuracy(
  unstained.fcs,
  spectra,
  af.spectra,
  asp,
  functions = c("assign.af.fluorophores", "assign.af.residuals", "assign.af.joint.cov"),
  n.downsample = 1000L,
  plot.dir = "figure_af_accuracy",
  title = "af_accuracy"
)
```

## Arguments

- unstained.fcs:

  Character scalar. Path to the unstained FCS file used as the reference
  data set.

- spectra:

  Numeric matrix of fluorophore spectra (fluorophores x detectors). Row
  names must be fluorophore names; column names must match the detector
  channels in the FCS file. Any row named `"AF"` is removed
  automatically before processing.

- af.spectra:

  Numeric matrix of AF spectra (AF variants x detectors). The first row
  is treated as the mean AF spectrum and is used for the baseline
  comparison. Column names must match those of `spectra`.

- asp:

  Aspect-ratio value passed to
  [`create.biplot`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).

- functions:

  Character vector of AF-function names to benchmark. Each name must
  resolve to a function in the current search path. Functions whose
  names start with `"fit."` are called with the signature
  `fn(raw.data, unmixed, unmixing.matrix, spectra, af.spectra)` and must
  return a list with elements `$unmixed` (cells x fluorophores, no AF
  column) and `$af.idx` (integer vector of per-cell AF-spectrum
  indices). All other functions are treated as assign-type and called
  with `fn(raw.data, spectra, af.spectra)`, returning an integer vector
  of AF-spectrum indices.

- n.downsample:

  Integer scalar. Maximum number of events read from the FCS file. A
  random subsample of this size is drawn when the file contains more
  events. Set to `Inf` to use all events. Default: `1000L`.

- plot.dir:

  Character scalar. Directory in which to save the biplot PDF. Created
  recursively if it does not exist. Default: `"figure_af_accuracy"`.

- title:

  Character scalar. Stem used to name the output PDF (the file will be
  `<plot.dir>/<title>_biplots.pdf`). Default: `"af_accuracy"`.

## Value

A named list with one entry per tested method (including the `"mean.af"`
baseline). Each entry is itself a list with elements:

- `Assignments`:

  Integer vector of per-cell AF-spectrum indices (all `1L` for
  `"mean.af"`).

- `Unmixed`:

  Numeric matrix of unmixed fluorophore values (cells x fluorophores, no
  AF column).

- `Similarity`:

  Numeric vector of per-cell cosine similarities between the raw
  detector signal and the assigned AF spectrum.

- `Mean_Sim`:

  Mean of `Similarity` (NAs excluded).

- `SD_Sim`:

  Standard deviation of `Similarity` (NAs excluded).

## See also

[`benchmark.af.spectra`](https://drcytometer.github.io/AutoSpectral/reference/benchmark.af.spectra.md),
[`unmix.ols`](https://drcytometer.github.io/AutoSpectral/reference/unmix.ols.md),
[`unmix.ols.fast`](https://drcytometer.github.io/AutoSpectral/reference/unmix.ols.fast.md),
[`create.biplot`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)
