# Build A Continuous Autofluorescence Basis From The Discrete Library

Constructs a low-dimensional basis for autofluorescence that replaces
the discrete AF library in the unmixing solve. The basis is derived from
the part of the AF library that the fluorophore panel cannot explain,
because only that component is identifiable and only that component
reaches the fluorophore abundances. Ranking by out-of-span energy rather
than total energy means the leading directions are the ones that
actually matter, and it makes the resulting normal equations diagonal,
so the per-cell solve is a single projection with no matrix inverse.

## Usage

``` r
get.af.basis.library(
  af.spectra,
  spectra,
  n.components = NULL,
  var.explained = 0.99
)
```

## Arguments

- af.spectra:

  Autofluorescence spectra, variants in rows and detectors in columns.
  Prepare using `get.af.spectra`.

- spectra:

  Fluorophore spectral signatures, fluorophores in rows and detectors in
  columns.

- n.components:

  Integer, number of basis directions to retain. Default `NULL`, in
  which case the count is chosen from `var.explained`.

- var.explained:

  Numeric in (0, 1\], default `0.99`. Fraction of the out-of-span
  variance of the AF library the retained basis must capture. Ignored
  when `n.components` is supplied.

## Value

A list with `basis` (detectors x components), `sigma` (singular values
of the projected library), `directions` (the orthonormal out-of-span
directions, detectors x components), `var.explained` (cumulative
fraction captured) and `n.components`.
