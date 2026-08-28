# Deconvolve Autofluorescence Background

Removes autofluorescence from raw detector data by fitting it jointly
with the panel rather than subtracting it beforehand. Each event is
unmixed against a design of the autofluorescence basis stacked on the
fluorophore spectra and only the fitted autofluorescence part is
subtracted, so a fluorophore's abundance and the background amount are
estimated together instead of competing for the same signal.

Fluorophores whose spectra lie close to the span of the autofluorescence
basis cannot be separated from it, and projecting them out removes real
signal. The returned `hotspot` matrix reports that coupling. When
`target` is supplied, components coupling to that fluorophore above
`max.hotspot` are dropped from its fit; the leading component is always
retained, because every event carries mean background.

## Usage

``` r
deconvolve.af.background(
  raw.data,
  spectra,
  af.basis,
  af.name = "AF",
  target = NULL,
  max.hotspot = 5,
  nonneg = TRUE
)
```

## Arguments

- raw.data:

  Numeric matrix (events x detectors), raw detector data.

- spectra:

  Numeric matrix (fluorophores x detectors), reference spectra. Any
  autofluorescence row named by `af.name` is removed before fitting,
  since the basis supersedes it.

- af.basis:

  Numeric matrix (components x detectors) from
  [`get.af.basis()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.basis.md).

- af.name:

  Character or `NULL`, the name of an autofluorescence row in `spectra`
  to drop. Default `"AF"`.

- target:

  Character or `NULL`. When supplied, autofluorescence components that
  this fluorophore cannot be separated from are dropped before fitting.
  Default `NULL`.

- max.hotspot:

  Numeric, the hotspot scale above which a component is considered
  inseparable from `target`. Default `5`.

- nonneg:

  Logical, whether to clamp fitted autofluorescence coefficients at zero
  before subtraction, so background can only be removed and never added.
  Default `TRUE`.

## Value

A named list:

- `residual`:

  `raw.data` with the fitted autofluorescence removed.

- `background`:

  The fitted autofluorescence contribution.

- `abundance`:

  Fluorophore abundances from the joint fit.

- `af.coef`:

  Per-event autofluorescence component coefficients, clamped at zero
  when `nonneg` is `TRUE`.

- `af.coef.raw`:

  The same coefficients before clamping, for use as regressors.

- `hotspot`:

  Hotspot matrix of the joint design.

- `af.used`, `af.dropped`:

  Component names kept and dropped.
