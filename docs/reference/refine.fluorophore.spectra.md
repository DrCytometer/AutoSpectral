# Refine Fluorophore Spectra

Optional second pass over
[`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md)'s
or
[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)'s
output, re-measuring each row directly from its own single-color control
with
[`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md)
– the same raw-space engine
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)
uses for its own phase two.

Checking a spectrum against the same gated, downsampled population it
was fit from cannot surface bias that population selection itself
introduced – exactly the effect documented for
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md).
A fresh read is the only way this check is independent of whatever
selection produced the first-pass spectra, whether that was
[`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)'s
scatter gate,
[`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)'s
downsampling, or
[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)'s
cosine-similarity filter.

That fresh read is done one control at a time, not by building a second,
whole-panel `flow.control` in memory: `.read.fcs.clean()` (the private
FCS reader
[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)
already uses) reads one file, drops saturating and doublet events, and
returns a plain matrix; nothing else is retained. Where a control needs
autofluorescence removed, its matching negative is read the same way and
the two are handed to
[`remove.af()`](https://drcytometer.github.io/AutoSpectral/reference/remove.af.md)
in a throwaway two-entry list –
[`remove.af()`](https://drcytometer.github.io/AutoSpectral/reference/remove.af.md)
only ever touches `clean.expr[[samp]]` and
`clean.expr[[matching.negative]]`, so it does not need
[`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)'s
whole-panel structure to run correctly on a single pair. Every control's
data is discarded once that control's fit is done. The trade is disk
I/O: every one of `n.iter` passes re-reads every control file, and its
matching negative, from scratch. Nothing is cached across controls or
across iterations, so peak memory is bounded by roughly one control file
and one negative file, regardless of panel size.

Every fluorophore's control is unmixed against its own row alone
(`active = target` in
[`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md)),
not the whole panel: a single-color control has only one fluorophore
truly present, so the population every other row should read zero in is
known, not inferred, and no nuisance subtraction is needed. A companion
diagnostic,
[`estimate.residual.spillover()`](https://drcytometer.github.io/AutoSpectral/reference/estimate.residual.spillover.md),
checks the same read's abundance in every *other* fluorophore's column,
with its negative mask set to "every event", unconditionally, for the
same reason, reusing the same read rather than a second pass over the
file. That check can only ever see the component of a row's error that
lies in the span of the current library – unmixing removes a spectral
error by projecting it through the compensation operator, so only the
projection is visible to any abundance-space residual – so it is logged
and optionally warned on, not used to update a row;
[`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md),
which regresses the raw detector trace directly rather than inverting
anything, is what recovers the full error and is what actually updates
`spectra`. Because it reuses this iteration's read rather than the
just-updated matrix, it reports cross-talk under the spectra this
iteration started from, not the spectra it ends with.

Every row update passes the same acceptance stack
[`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)
phase two uses (`max.angle`, `min.explained`/`max.explained`,
`max.resid`, `max.intercept`, `min.bg.align`); a row that fails any gate
keeps its starting spectrum for that iteration. Every fluorophore is
refit against the same starting `spectra` within one iteration and all
accepted updates are applied together at the end of it, so results do
not depend on the order controls happen to be listed in.

## Usage

``` r
refine.fluorophore.spectra(
  marker.spectra,
  control.dir,
  control.def.file,
  asp,
  af.remove = TRUE,
  af.figures = FALSE,
  singlet.quantiles = c(0.85, 0.975),
  remove.doublets = TRUE,
  allow.duplicate.controls = FALSE,
  n.iter = 3L,
  intercept = TRUE,
  multivariate = TRUE,
  ridge = 1e-06,
  n.levels = 60L,
  min.bin.events = 50L,
  min.events = 200L,
  max.angle = 5,
  min.explained = 0.8,
  max.explained = 1.2,
  max.resid = 0.05,
  max.intercept = 0.05,
  min.bg.align = -0.9,
  unstained.threshold = 0.99,
  unstained.margin = 1.3,
  crosstalk.check = TRUE,
  max.crosstalk = 0.1,
  n.levels.pair = 10L,
  convergence.threshold = 0.5,
  step = 1,
  n.threads = 1L,
  verbose = TRUE
)
```

## Arguments

- marker.spectra:

  Numeric matrix (samples x detectors), the first-pass spectra. Rownames
  must be the control `sample` identifiers (not necessarily the dye
  identity – a dye can have more than one control);
  `attr(marker.spectra, "fluorophore")`, if present, is used to keep
  replicate controls of the same dye out of each other's cross-talk
  check.

- control.dir, control.def.file:

  As in
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md).
  Should normally be the same files used to build the spectra
  `marker.spectra` was fit from.

- asp:

  The AutoSpectral parameter list.

- af.remove:

  Logical, default `TRUE`. Whether to read each cell-type control's
  matching negative and run
  [`remove.af()`](https://drcytometer.github.io/AutoSpectral/reference/remove.af.md)'s
  intrusive-AF gate on it. Bead controls are never AF-removed, matching
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).
  Set `FALSE` for a bead-only panel.

- af.figures:

  Logical, default `FALSE`. Whether
  [`remove.af()`](https://drcytometer.github.io/AutoSpectral/reference/remove.af.md)
  writes its own AF-removal diagnostic figures for each control.

- singlet.quantiles, remove.doublets:

  As in
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md),
  passed to `.read.fcs.clean()` for every file read here.

- allow.duplicate.controls:

  Logical, default `FALSE`. As in
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)/[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md);
  set `TRUE` if the control file the first pass used permits multiple
  controls per dye.

- n.iter:

  Integer, maximum refine iterations. Default `3L`.

- intercept, multivariate, ridge, n.levels, min.bin.events:

  As in
  [`extract.raw.signature()`](https://drcytometer.github.io/AutoSpectral/reference/extract.raw.signature.md).

- min.events:

  Integer, minimum events for both the signature fit and the cross-talk
  pair estimator. Default `200L`.

- max.angle:

  Numeric, degrees. A candidate row is rejected if its cosine distance
  from the current row exceeds this. Default `5`.

- min.explained, max.explained, max.resid, max.intercept, min.bg.align:

  As in
  [`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md)'s
  phase-two acceptance gates.

- unstained.threshold, unstained.margin:

  Numeric, used to set the cross-talk check's source threshold from the
  control's matching negative, when one is defined. Same convention as
  [`fix.my.unmix()`](https://drcytometer.github.io/AutoSpectral/reference/fix.my.unmix.md).

- crosstalk.check:

  Logical, whether to run the diagnostic
  [`estimate.residual.spillover()`](https://drcytometer.github.io/AutoSpectral/reference/estimate.residual.spillover.md)
  pass each iteration. Default `TRUE`.

- max.crosstalk:

  Numeric, the abundance-space coefficient above which a console warning
  is printed. Diagnostic only; does not gate a row. Default `0.1`.

- n.levels.pair:

  Integer, abundance bins for the cross-talk pair estimator. Default
  `10L`.

- convergence.threshold:

  Numeric, degrees. Iteration stops early once the largest accepted
  `deg.change` in a pass falls below this. Default `0.5`.

- step:

  Numeric in (0, 1\], the fraction of each accepted change applied.
  Default `1`.

- n.threads:

  Integer, threads for the cross-talk pair estimator. Default `1L`.

- verbose:

  Logical, default `TRUE`.

## Value

A named list:

- `spectra`:

  The refined matrix, same shape and rownames as `marker.spectra`, with
  `attr(spectra, "fluorophore")` preserved.

- `log`:

  Data frame, one row per sample per iteration, with the acceptance
  decision and fit diagnostics.

- `crosstalk`:

  Data frame of
  [`estimate.residual.spillover()`](https://drcytometer.github.io/AutoSpectral/reference/estimate.residual.spillover.md)
  output for every sample and iteration, or `NULL` if
  `crosstalk.check = FALSE`.
