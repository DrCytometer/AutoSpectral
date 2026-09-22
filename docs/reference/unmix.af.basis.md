# Unmix With A Continuous Autofluorescence Basis

Solves `y = t(spectra) %*% f + basis %*% k` for every cell, fitting
autofluorescence as a free combination of basis directions rather than
selecting one row from a discrete library. Because the basis is built
from the panel-oblique part of the AF library, the normal equations for
`k` are diagonal and the per-cell solve is a single projection.

Fits are unconstrained, so a cell whose reconstructed AF spectrum goes
meaningfully negative is flagged in `negative` and is a candidate for
falling back to the discrete library.

## Usage

``` r
unmix.af.basis(
  raw.data,
  spectra,
  af.basis,
  af.spectra = NULL,
  return.fitted.af = FALSE,
  negative.tol = 0.05
)
```

## Arguments

- raw.data:

  Expression data from raw FCS files. Cells in rows and detectors in
  columns.

- spectra:

  Fluorophore spectral signatures, fluorophores in rows and detectors in
  columns.

- af.basis:

  An AF basis from `get.af.basis`.

- af.spectra:

  Optional AF library used only to report a nearest-library index for
  each cell, for compatibility with the discrete workflow. Default
  `NULL`.

- return.fitted.af:

  Logical, default `FALSE`. Whether to return the fitted
  autofluorescence in detector space.

- negative.tol:

  Numeric, default `0.05`. A cell is flagged when the most negative
  detector of its reconstructed AF spectrum falls below `-negative.tol`
  times its largest.

## Value

A list with `fluorophores` (cells x fluorophores), `k` (cells x
components), `af` (total fitted AF scale per cell), `negative` (logical
vector) and, when requested, `fitted.af` and `af.index`.
