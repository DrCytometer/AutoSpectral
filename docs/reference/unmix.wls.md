# Unmix Using Weighted Least Squares

This function performs unmixing of raw data using weighted least
squares, AKA WLS, based on the provided spectra. Weighting is by channel
power.

## Usage

``` r
unmix.wls(raw.data, spectra, weights = NULL)
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

  Optional numeric vector of weights, one per fluorescent detector.
  Default is `NULL`, in which case weighting will be done by channel
  means.

## Value

A matrix containing unnmixed data with cells in rows and fluorophores in
columns.
