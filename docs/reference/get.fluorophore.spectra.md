# Get Fluorophore Spectra

This function retrieves the fluorophore spectra for flow cytometry data,
optionally using cleaned expression data. It also plots and saves the
spectra, and performs cosine similarity QC for controls.

## Usage

``` r
get.fluorophore.spectra(
  flow.control,
  asp,
  use.clean.expr = TRUE,
  af.spectra = NULL,
  title = NULL,
  figures = TRUE,
  refine = FALSE,
  control.dir = NULL,
  control.def.file = NULL,
  refine.args = list()
)
```

## Arguments

- flow.control:

  A list containing flow cytometry control data.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- use.clean.expr:

  Logical indicating whether to use cleaned expression data, default is
  `TRUE`

- af.spectra:

  Optional autofluorescence spectra to include.

- title:

  Optional prefix for plot titles, default is `NULL`, which gives
  "Initial" when `use.clean.expr` is `FALSE` and "Clean" when
  `use.clean.expr` is `TRUE`.

- figures:

  Logical, default is `TRUE`. Whether to produce plots of the
  fluorophore spectra and cosine similarity.

- refine:

  Logical, default `FALSE`. Whether to re-measure each row directly on
  its own single-color control after the first pass, via
  [`refine.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/refine.fluorophore.spectra.md).
  Unlike the first pass, this re-reads the controls fresh from
  `control.dir`/`control.def.file`, one file at a time, rather than
  using `flow.control`'s already-gated and possibly downsampled data –
  checking a spectrum against the same selection it was fit from cannot
  surface bias that selection introduced. Requires `control.dir` and
  `control.def.file`.

- control.dir, control.def.file:

  As in
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md).
  Ignored unless `refine = TRUE`, in which case both are required and
  should normally be the same files `flow.control` was built from.

- refine.args:

  Named list of further arguments passed to
  [`refine.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/refine.fluorophore.spectra.md)
  when `refine = TRUE`.

## Value

A matrix with the fluorophore spectra. When `refine = TRUE`, also
carries `attr(., "refine.log")` and `attr(., "refine.crosstalk")` with
the per-iteration diagnostics.
