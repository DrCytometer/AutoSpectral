# Unmix WLS Per Cell

Faster solver for per-cell optimization workflow. Performs spectral
unmixing using weighted least squares.

## Usage

``` r
unmix.wls.per.cell(cell.raw, spectra, weights)
```

## Arguments

- cell.raw:

  Expression data from raw FCS files. A single row of cells with
  detectors in columns. Columns should be fluorescent data only and must
  match the columns in spectra.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- weights:

  The weighting values for the detectors. No checks are performed in
  this function; `unmix.wls` should be used for standard cases.

## Value

Unmixed data (one row), fluorophores in columns.
