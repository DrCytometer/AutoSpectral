# Unmix WLS Per Cell

Faster solver for per-cell optimization workflow. Performs spectral
unmixing using weighted least squares.

## Usage

``` r
unmix.wls.per.cell(cell.raw, spectra, weights = NULL, noise.floor = 125)
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

  The weighting values for the detectors. No checks are performed on
  supplied weights; `unmix.wls` should be used for standard cases.
  Default is `NULL`, in which case weighting is computed from this
  cell's own signal, floored at `noise.floor`.

- noise.floor:

  Numeric, default `125`. Lower clamp on this cell's signal before
  inversion when `weights = NULL`. Signal units, same convention as
  `noise.floor` elsewhere in the package. Ignored if `weights` is
  supplied.

## Value

Unmixed data (one row), fluorophores in columns.
