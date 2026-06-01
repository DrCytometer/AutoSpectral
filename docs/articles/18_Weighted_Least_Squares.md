# 18 Weighted Least Squares Unmixing

## Why weighting matters

Ordinary least squares (OLS) unmixing treats every detector as equally
reliable. In practice, detectors are not equally reliable. Under a
Poisson noise model, absolute noise scales with signal level: a channel
carrying high photon counts has high absolute variance and therefore a
less precise measurement. Weighted least squares (WLS) unmixing accounts
for this by down-weighting the noisier high-signal channels and
up-weighting the quieter low-signal channels, where measurements are
more precise.

### The Poisson variance model

Flow cytometry data approximate a Poisson process: the number of photons
detected in a given time window is Poisson distributed. A central
property of the Poisson distribution is that its variance equals its
mean. This means a channel with higher average signal also has higher
variance, and therefore higher absolute noise. The standard WLS weight
for channel $`j`$ is therefore:

``` math
w_j = \frac{1}{\text{mean signal in channel } j}
```

High mean signal means high absolute variance and therefore a low
weight. Low mean signal means low absolute variance: the channel
reliably reads near zero, and the solver treats that reading as precise.
Note that this is opposite to what the coefficient of variation (CV)
would suggest – a dim channel has a large CV but small absolute noise.
WLS, at least in AutoSpectral, weights by absolute Poisson noise, not
relative noise.

The lower bound of `1e-6` is applied before taking the reciprocal to
avoid division-by-zero for channels that record essentially nothing.

In practice, WLS is particularly recommended for cytometers such as the
Sony ID7000, FACSDiscover A8, and FACSDiscover S8, which have large
numbers of detectors and are therefore more over-determined.

------------------------------------------------------------------------

## `calculate.weights()` – extracting weights from a file

[`calculate.weights()`](https://drcytometer.github.io/AutoSpectral/reference/calculate.weights.md)
reads an FCS file, extracts the spectral detector columns, and computes
the inverse-mean weights described above.

``` r

library( AutoSpectral )

# weights for a single representative file
weights <- calculate.weights(
  fcs.file        = "./samples/representative_sample.fcs",
  spectral.channels = flow.control$spectral.channel
)
```

The return value is a named numeric vector, one element per spectral
detector. Names match the detector channel names in the FCS file so that
they can be matched unambiguously against the spectra matrix later.

You can optionally save the weights to a CSV for inspection or later
reuse:

``` r

weights <- calculate.weights(
  fcs.file          = "./samples/representative_sample.fcs",
  spectral.channels = flow.control$spectral.channel,
  save              = TRUE,
  output.dir        = "./table_spectra",
  filename          = "weights.csv"
)
```

The saved CSV can be read back with
[`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html) and the
weights column extracted as a named vector.

------------------------------------------------------------------------

## `unmix.wls()` – performing WLS unmixing

[`unmix.wls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.wls.md)
takes a matrix of raw spectral measurements, a spectra matrix, and an
optional weight vector, and returns the unmixed fluorophore abundances.

``` r

# raw.data: cells x detectors
# spectra:  fluorophores x detectors (normalised 0-1)
unmixed <- unmix.wls( raw.data, spectra, weights = weights )
```

If `weights` is `NULL`,
[`unmix.wls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.wls.md)
falls back to computing the inverse-column-mean weights internally from
`raw.data`. This is convenient for one-off calls but means the weights
are derived from whatever data you pass in at that moment. For
consistent results across multiple files, compute and reuse a shared
weight vector (see below).

Internally,
[`unmix.wls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.wls.md)
constructs the diagonal weight matrix $`W`$, applies it to the spectra
matrix $`A`$ to form $`A^* = W^{1/2} A`$, and solves via SVD for
numerical stability:

``` math
\hat{x} = (A^T W A)^{-1} A^T W \, y
```

The SVD path is preferred over direct matrix inversion because the
weighting can make the system poorly conditioned when some channels have
very low signal.

------------------------------------------------------------------------

## `save.unmixing.matrix()` – exporting the unmixing matrix

It is sometimes useful to inspect or archive the unmixing matrix itself
– for example, to compare OLS and WLS coefficients, or to apply the same
matrix outside of AutoSpectral.
[`save.unmixing.matrix()`](https://drcytometer.github.io/AutoSpectral/reference/save.unmixing.matrix.md)
computes and writes the matrix to a CSV, with an optional spectral
heatmap.

``` r

# OLS unmixing matrix (no weights)
save.unmixing.matrix(
  spectra    = spectra,
  filename   = "unmixing_matrix_OLS",
  output.dir = "./table_spectra",
  figures    = TRUE
)

# WLS unmixing matrix
save.unmixing.matrix(
  spectra    = spectra,
  weights    = weights,
  filename   = "unmixing_matrix_WLS",
  output.dir = "./table_spectra",
  figures    = TRUE
)
```

The heatmap (when `figures = TRUE`) shows the contribution coefficient
of each fluorophore in each detector. In a well-separated panel you will
see a strong diagonal pattern. WLS coefficients will differ from OLS in
the rows corresponding to detectors with lower signal, where WLS places
less weight.

------------------------------------------------------------------------

## Integration with `unmix.fcs()` and `unmix.folder()`

Both
[`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
and
[`unmix.folder()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.folder.md)
accept a `method` argument that selects the unmixing algorithm. Use
`method = "WLS"` to activate weighted least squares:

``` r

unmix.fcs(
  fcs.file     = "./samples/experiment.fcs",
  spectra      = spectra,
  asp          = asp,
  flow.control = flow.control,
  method       = "WLS"
)
```

When `weights = NULL` (the default),
[`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
automatically computes inverse-mean weights from the first 100,000
events of the file being unmixed. This is reasonable for a single file,
but can produce slightly inconsistent results across a batch because
each file contributes its own weight vector.

For batch processing with
[`unmix.folder()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.folder.md),
the same on-the-fly approach is used by default:

``` r

unmix.folder(
  sample.dir   = "./samples",
  spectra      = spectra,
  asp          = asp,
  flow.control = flow.control,
  method       = "WLS"
)
```

To ensure a consistent weight vector across all files in the batch, pass
pre-computed weights explicitly:

``` r

unmix.folder(
  sample.dir   = "./samples",
  spectra      = spectra,
  asp          = asp,
  flow.control = flow.control,
  method       = "WLS",
  weights      = weights    # computed once from a representative file
)
```

------------------------------------------------------------------------

## Reusing weights across a batch

For large experiments, computing weights from a single representative
file and reusing them across all samples has two advantages:

1.  **Consistency** – the unmixing matrix is identical for every file,
    so comparisons between files are not confounded by per-file weight
    variation.
2.  **Speed** –
    [`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
    skips the weight-computation read when weights are provided, saving
    one FCS read per file.

A practical strategy for choosing the representative file:

- **Concatenated file.** Use
  [`concatenateFCS()`](https://drcytometer.github.io/AutoSpectral/reference/concatenateFCS.md)
  to combine all samples into one file, then compute weights from that.

``` r
# concatenate a random subset of all sample files
concatenateFCS(
  fcs.dir    = list.files(path = "./samples", pattern = ".fcs", full.names = TRUE),
  output.name = "My merged.fcs"
  output.dir =  "./concatenated_fcs"
)

concat.file <- "./concatenated_fcs/My merged.fcs"

weights <- calculate.weights(
  fcs.file          = concat.file,
  spectral.channels = flow.control$spectral.channel,
  save              = TRUE,
  output.dir        = "./table_spectra"
)

# now unmix all files with those weights
unmix.folder(
  sample.dir   = "./samples",
  spectra      = spectra,
  asp          = asp,
  flow.control = flow.control,
  method       = "WLS",
  weights      = weights
)
```

- **A single well-stained sample.** If you know one file is
  representative of the full experiment (e.g., a fully stained PBMC
  sample), use that. Avoid files that are missing markers or have
  unusually low cell counts, as these will produce atypical per-channel
  means and therefore atypical weights.

Saved weights can also be reloaded in a later R session:

``` r

w <- utils::read.csv( "./table_spectra/weights.csv", row.names = 1 )
weights <- setNames( w[, 1], rownames( w ) )
```

------------------------------------------------------------------------

## `unmix.poisson()` – iterative per-cell refinement

WLS with global weights is still an approximation: it uses *average*
signal levels to characterise noise, whereas the actual noise in a
Poisson process is cell-specific – a brighter cell has more
photon-counting noise in every channel than a dim cell.
[`unmix.poisson()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.poisson.md)
addresses this by replacing the global weight vector with per-cell
weights derived from the cell’s own predicted signal, using iteratively
reweighted least squares (IRLS).

The algorithm:

1.  **Initialise** with a standard WLS solution (using the global
    `initial.weights` if provided, otherwise inverse column means).
2.  **Per cell:** fit a Poisson GLM with an identity link, using the
    current predicted fluorophore abundances as starting values. The
    Poisson GLM implicitly weights by the predicted mean at each step of
    the IRLS, so the weights are updated to reflect the *cell’s* photon
    count rather than the population average.
3.  If the GLM converges, accept the new coefficients. If it fails to
    converge, retain the WLS initialisation for that cell.

``` r

unmixed <- unmix.poisson(
  raw.data        = raw.data,
  spectra         = spectra,
  asp             = asp,
  initial.weights = weights,    # optional; from calculate.weights()
  parallel        = TRUE,
  threads         = asp$worker.process.n
)
```

To use Poisson unmixing via
[`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md),
set `method = "Poisson"` or `method = "FastPoisson"`. `FastPoisson`
requires the `AutoSpectralRcpp` companion package and implements the
IRLS loop in C++ for substantially faster processing:

``` r

unmix.fcs(
  fcs.file     = "./samples/experiment.fcs",
  spectra      = spectra,
  asp          = asp,
  flow.control = flow.control,
  method       = "FastPoisson",
  parallel     = TRUE
)
```

`FastPoisson` will fall back to standard `Poisson` if `AutoSpectralRcpp`
is not available, and both will fall back to their WLS initialisation
for any cell that does not converge.

### When to use Poisson unmixing

Poisson unmixing is most beneficial for panels with a wide dynamic range
– where some cells are dim in a given channel and others are very
bright, and you want the unmixing to be equally accurate for both. It
tends to help most in panels that are small relative to the number of
detectors (highly over-determined), and generally when those panels are
not particularly well designed. For most panels or when speed is
critical, WLS with a representative weight vector is often sufficient
and considerably faster. The per-cell Poisson approach is more
computationally intensive than WLS because it runs a separate
optimisation for every event; see the [Speed It
Up](https://drcytometer.github.io/AutoSpectral/articles/14_Speed_It_Up.html)
article for benchmarks and guidance on accelerating it with
`AutoSpectralRcpp`. Also, per-cell autofluorescence extraction generally
outperforms WLS and Poisson.

The reference for the Poisson unmixing model is:

> Novo D. et al. (2014). Generalized Unmixing Model for Multispectral
> Flow Cytometry Utilizing Nonsquare Compensation Matrices. *Cytometry
> Part A*, 83(5):508-520. <https://doi.org/10.1002/cyto.a.22272>
