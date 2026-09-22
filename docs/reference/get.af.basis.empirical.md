# Build A Continuous Autofluorescence Basis From Raw Events

As
[`get.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.basis.md),
but derives the basis directly from the panel-oblique residual of real
unstained single-cell events rather than from a SOM-derived spectral
library. A SOM library is itself a quantisation of the unstained
sample's cell-to-cell diversity into a fixed set of node shapes, so
building the basis one step further upstream, from the events the SOM
was trained on, removes that quantisation before it happens rather than
fitting around it.

Structurally identical to
[`get.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.basis.md):
both take the panel-oblique residual of a set of reference spectra
(there, SOM nodes; here, raw events), take its SVD, and keep the leading
out-of-span directions. Ranking by out-of-span variance rather than
total variance matters here more than in the library case, because raw
single-cell data also carries shot noise and any residual spillover
contamination, both of which add variance that has nothing to do with
autofluorescence shape.

## Usage

``` r
get.af.basis.empirical(
  unstained.exprs,
  spectra,
  n.components = NULL,
  var.explained = 0.99,
  trim.quantile = 0.99,
  max.events = 50000,
  seed = 1
)
```

## Arguments

- unstained.exprs:

  Numeric matrix of raw unstained expression data, cells in rows and
  detectors in columns. Columns must match `spectra`.

- spectra:

  Fluorophore spectral signatures, fluorophores in rows and detectors in
  columns.

- n.components:

  Integer, number of basis directions to retain. Default `NULL`, in
  which case the count is chosen from `var.explained`.

- var.explained:

  Numeric in (0, 1\], default `0.99`. Fraction of the out-of-span
  variance of the trimmed, subsampled event set the retained basis must
  capture. Ignored when `n.components` is supplied.

- trim.quantile:

  Numeric in (0, 1\], default `0.99`. Events whose out-of-span residual
  norm exceeds this quantile are excluded before the basis is built.
  Real unstained samples can carry debris, doublets, or residual
  spillover-contaminated events; the SVD has no built-in robustness to
  them, so they are trimmed explicitly rather than left to dominate the
  leading components. Set to `1` to disable.

- max.events:

  Integer, default `5e4`. If more events remain after trimming, a random
  subsample of this size is used to build the basis. The basis is then
  applied to every event via
  [`unmix.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.af.basis.md),
  so this only bounds the cost of basis *construction*, not of applying
  it.

- seed:

  Integer, default `1`. Seed for the subsampling step.

## Value

A list with `basis` (detectors x components, real combinations of the
sampled events), `directions` (detectors x components, the orthonormal
out-of-span directions), `sigma` (singular values), `var.explained`
(cumulative fraction captured), `n.components`, and `n.trimmed` (events
excluded by `trim.quantile`). Compatible with
[`unmix.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.af.basis.md).

## See also

`get.af.basis`, `unmix.af.basis`
