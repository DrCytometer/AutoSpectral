# Save Unmixing Matrix

This function calculates and saves the unmixing matrix to a CSV file. By
default, it will produce an ordinary least-squares (OLS) matrix. To save
a weighted least-squares (WLS) unmixing matrix, supply detector weights
to argument `weights`. Weights should be a numeric vector the same
length as the detectors in `spectra`. Good options for this are to
use 1) detector noise (reliability) measurements, for instance, CVs from
instrument QC or 2) empirical data from an FCS file. For the latter,
take the raw data in the spectral detector channels and apply function
`colMeans`. This provides the average signal level in the detecors,
which we can assume is proportional to the noise level. In both cases,
the weights should be the inverse (1/noise). Use function
[`calculate.weights()`](https://drcytometer.github.io/AutoSpectral/reference/calculate.weights.md)
to extract `colMeans` weighting from an FCS file.

## Usage

``` r
save.unmixing.matrix(
  spectra,
  weights = NULL,
  filename = "unmixing_matrix",
  output.dir = "./table_spectra",
  figures = TRUE,
  color.palette = "mako"
)
```

## Arguments

- spectra:

  Matrix or dataframe containing spectral data. Note that the unmixing
  matrix depends on the combination of fluorophore signals present in
  the spectral mixing matrix.

- weights:

  Optional numeric vector of weights, one per fluorescent detector.
  Default is `NULL`, in which case OLS unmixing will be calculated.

- filename:

  Character string defining the output file name. Default is
  `unmixing_matrix`.

- output.dir:

  Optional output directory. Default is `./table_spectra`.

- figures:

  Logical, if `TRUE`, produces a spectral heatmap of the unmixing matrix
  in `output.dir`.

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `mako`. Options are the
  viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
  `cividis`, `rocket`, `mako` and `turbo`.
