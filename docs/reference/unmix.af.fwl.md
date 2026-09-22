# Unmix With Per-Cell Autofluorescence By Frisch-Waugh

Solves the joint least-squares problem `y = t(spectra) %*% f + k * af`
for every cell at once, where each cell may use a different
autofluorescence spectrum. The Frisch-Waugh-Lovell decomposition splits
the joint solve into two precomputable library projections plus one pass
over the data, so no per-cell or per-group design matrix is ever formed.
Results are identical to solving each AF group separately with
[`unmix.ols.fast()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.ols.fast.md)
on the stacked spectra, to floating-point precision.

## Usage

``` r
unmix.af.fwl(
  raw.data,
  spectra,
  af.spectra,
  af.index,
  unmixed.no.af = NULL,
  return.fitted.af = FALSE,
  denominator.floor = 0,
  chunk.size = 100000L
)
```

## Arguments

- raw.data:

  Expression data from raw FCS files. Cells in rows and detectors in
  columns. Columns must match the columns in `spectra`.

- spectra:

  Spectral signatures of fluorophores, with fluorophores in rows and
  detectors in columns.

- af.spectra:

  Spectral signatures of autofluorescences, with variants in rows and
  detectors in columns. Prepare using `get.af.spectra`.

- af.index:

  Integer vector, one entry per cell, giving the row of `af.spectra`
  assigned to that cell.

- unmixed.no.af:

  Optional numeric matrix (cells x fluorophores) holding the AF-free
  unmixing `raw.data %*% t(U)`. Supply it when the caller has already
  computed it to avoid repeating the largest matrix product. Default
  `NULL`, in which case it is computed here.

- return.fitted.af:

  Logical, default `FALSE`. Whether to also return the fitted
  autofluorescence in detector space, `k * af.spectra[af.index, ]`.

- denominator.floor:

  Numeric, default `0`. When positive, floors each AF candidate's
  out-of-span self-dot at this fraction of the largest. An AF variant
  lying almost inside the fluorophore span has a vanishing out-of-span
  direction and an unidentifiable abundance; flooring caps the
  amplification. `0` reproduces the unregularised least-squares solution
  exactly.

- chunk.size:

  Integer, default `1e5`. Number of cells processed per block, bounding
  the size of the intermediate detector-wide gather.

## Value

A list with elements `fluorophores` (cells x fluorophores), `af`
(numeric vector of per-cell AF abundance) and, when requested,
`fitted.af` (cells x detectors).
