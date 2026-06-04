# Assign AF Spectrum By Scatter-Matched Reference Averaging

An experimental AF assignment approach for benchmarking against the
existing covariance/residual methods (`assign.af.fluorophores`,
`assign.af.residuals`, `assign.af.joint.cov`).

For each cell in `test.data`, the `k` nearest neighbours in `ref.data`
are found by Euclidean distance in scatter space (using
[`FNN::knnx.index`](https://rdrr.io/pkg/FNN/man/knn.index.html)). Their
spectral channels are averaged to produce a cell-specific reference
spectrum. The cosine similarity between each test cell's spectral
signature and its scatter-matched reference average is returned
alongside the best-fitting AF variant index, allowing direct comparison
against the existing assignment approaches.

The rationale: cells with identical scatter profiles (size, granularity)
should have near-identical autofluorescence. If the scatter-matched
average is a better proxy for the true AF than the library-based
assignment, the returned cosine similarities will be systematically
higher than those obtained from the covariance/residual pipeline.

## Usage

``` r
assign.af.scatter.match(
  test.data,
  ref.data,
  scatter.param,
  spectra,
  k = 5L,
  af.spectra = NULL,
  scale.scatter = TRUE,
  verbose = TRUE
)
```

## Arguments

- test.data:

  Numeric matrix. Expression data from the **test** unstained FCS file.
  Cells in rows, channels in columns. Must contain both `scatter.param`
  channels and the spectral detector channels. Can also be the path to
  an FCS file (character scalar), in which case the file is read with
  [`readFCS()`](https://drcytometer.github.io/AutoSpectral/reference/readFCS.md).

- ref.data:

  Numeric matrix. Expression data from the **reference** unstained FCS
  file. Same column layout as `test.data`. Can also be an FCS file path.

- scatter.param:

  Character vector of length \>= 1. Column names that identify the
  scatter channels to use for kNN matching (e.g. `c("FSC-A", "SSC-A")`).
  These are excluded from the spectral similarity calculation.

- k:

  Integer. Number of nearest reference neighbours to average. Default
  `5`. Larger values stabilise the reference estimate but may blur
  genuine AF heterogeneity.

- af.spectra:

  Optional numeric matrix. Spectral signatures of AF variants,
  normalised `[0, 1]`, with variants in rows and detectors in columns.
  When supplied, the function also assigns each test cell to the closest
  AF variant (by cosine similarity to the scatter-matched average) and
  returns that index. Omit to skip variant assignment.

- scale.scatter:

  Logical. Whether to z-score-standardise the scatter channels before
  computing kNN distances (recommended when FSC and SSC have very
  different dynamic ranges). Default `TRUE`.

- verbose:

  Logical. Whether to emit progress messages. Default `TRUE`.

## Value

A list with the following elements:

- `cosine.similarity`:

  Numeric vector, length `nrow(test.data)`. Per-cell cosine similarity
  between the test cell's spectral signature and its scatter-matched
  reference average. Values close to 1 indicate a strong match.

- `ref.average`:

  Numeric matrix, same dimensions as the spectral portion of
  `test.data`. The scatter-matched averaged reference spectrum for each
  test cell (i.e., the mean of the `k` nearest reference neighbours).

- `af.assignment`:

  Integer vector (or `NULL` if `af.spectra` was not supplied). Row index
  into `af.spectra` of the best-fitting AF variant for each test cell,
  determined by cosine similarity between the scatter-matched reference
  average and each AF library spectrum.

- `nn.index`:

  Integer matrix, `nrow(test.data)` x `k`. The row indices in `ref.data`
  of the `k` nearest scatter neighbours for each test cell. Useful for
  diagnostics.

- `summary`:

  A one-row data frame with aggregate statistics (mean, median, sd of
  cosine similarities) for quick comparison against competing methods.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using FCS file paths:
result <- assign.af.scatter.match(
  test.data    = "path/to/test_unstained.fcs",
  ref.data     = "path/to/ref_unstained.fcs",
  scatter.param = c( "FSC-A", "SSC-A" ),
  k            = 5,
  af.spectra   = my.af.spectra   # optional
)
hist( result$cosine.similarity, main = "Scatter-match cosine similarity" )

# Or pass matrices directly:
result <- assign.af.scatter.match(
  test.data    = test.mat,
  ref.data     = ref.mat,
  scatter.param = c( "FSC-A", "SSC-A" ),
  k            = 5
)
} # }
```
