# Unmix WLS Fast

Faster solver for per-cell optimization workflow. Performs spectral
unmixing using weighted least squares.

## Usage

``` r
unmix.wls.fast(raw.data, spectra, weights = NULL)
```

## Arguments

- raw.data:

  Expression data from raw FCS files. Cells in rows and detectors in
  columns. Columns should be fluorescent data only and must match the
  columns in spectra.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- weights:

  The weighting values for the detectors. No checks are performed in
  this function; `unmix.wls` should be used for most cases.

## Value

Unmixed data with cells in rows and fluorophores in columns.
