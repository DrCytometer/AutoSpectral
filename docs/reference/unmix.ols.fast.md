# Unmix OLS Fast

Faster solver for per-cell optimization workflow. Performs spectral
unmixing using ordinary least squares.

## Usage

``` r
unmix.ols.fast(raw.data, spectra, weights = NULL)
```

## Arguments

- raw.data:

  Expression data from raw fcs files. Cells in rows and detectors in
  columns. Columns should be fluorescent data only and must match the
  columns in spectra.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- weights:

  Dummy argument to allow dynamic switching between OLS and WLS. Default
  is `NULL`. Values passed to `weights` will be ignored.

## Value

Unmixed data with cells in rows and fluorophores in columns.
