# Plot Standard (Manual-Gating) Spectra Extraction Workflow

Builds a manuscript-ready, multi-panel figure illustrating a "standard"
single-stained-control workflow of the kind used in vendor software
(e.g. SpectroFlo): an octagon gate on FSC-A vs SSC-A positioned on the
dominant population, selection of the brightest events on the
fluorophore's peak channel as the positive population and the dimmest
fraction as an internal negative, biplots of both fractions (negative in
black, positive coloured by cosine similarity), and a population- level
background subtraction. Intended as a direct comparison against
[`spectra.automated.steps.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.automated.steps.plot.md)
for the same fluorophore(s). Uses the same building blocks as the rest
of AutoSpectral
([`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)-style
pseudocolour density,
[`spectral.trace()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.trace.md))
so panel styling matches the package's other figures.

Requires `plot_spectra_automated_steps.R` to be loaded in the same
package namespace (several private helpers are shared).

## Usage

``` r
spectra.standard.workflow.plot(
  control.dir,
  control.def.file,
  asp,
  fluorophores = NULL,
  octagon.width.factor = 3,
  n.bright.events = 2000L,
  negative.quantile = 0.25,
  negative.quantile.min = 0.01,
  x.min.quantile = 0.005,
  singlet.quantiles = c(0.85, 0.975),
  gate.color = "darkgoldenrod1",
  density.palette = "rainbow",
  selection.fill.color = "steelblue",
  selection.line.color = "black",
  negative.bracket.color = "#377EB8",
  positive.bracket.color = "#E41A1C",
  negative.point.color = "black",
  event.point.size = NULL,
  n.highlight = 200L,
  clean.positive.color = "red",
  clean.positive.point.size = NULL,
  cells.trace.color = "#D95F02",
  beads.trace.color = "#377EB8",
  af.trace.color = "grey40",
  max.points = 50000,
  panel.width = 4,
  panel.height = 4,
  composite.width = NULL,
  composite.height = NULL,
  output.dir = NULL,
  save = TRUE,
  file.type = "jpg",
  verbose = TRUE
)
```

## Arguments

- control.dir:

  Character. Path to the directory containing the single-stained control
  FCS files.

- control.def.file:

  Character. Path to (or filename of) the control definition CSV – for
  the plotted fluorophore itself, only the `fluorophore` / `filename`
  columns are used, since this workflow always derives an internal
  negative from the same file. If the file has a `control.type` column
  (as used by
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  /
  [`check.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/check.control.file.md)),
  a row with `control.type == "beads"` sharing a fluorophore name with a
  "cells" row is treated as a paired bead control: its `filename` is the
  positive bead sample, and its `universal.negative` must point at the
  matching unstained bead file. When present, this pair supplies the "
  (Beads)" trace in panel D instead of the static spectral reference
  library.

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- fluorophores:

  Character vector of fluorophore name(s) to illustrate. Default `NULL`
  illustrates the first fluorophore found.

- octagon.width.factor:

  Numeric, default `3`. Controls how large the octagon gate is drawn, as
  a multiple of the local median-absolute- deviation spread of events
  around the density peak.

- n.bright.events:

  Integer, default `2000`. Number of top-expressing events on the peak
  channel (within the octagon gate) selected as the "positive"
  population.

- negative.quantile:

  Numeric in `(0, 1)`, default `0.25`. Upper quantile bound (on the peak
  channel, within octagon-gated events) of the internal "negative"
  population.

- negative.quantile.min:

  Numeric in `[0, negative.quantile)`, default `0.01`. Lower quantile
  bound of the internal "negative" population – the negative gate spans
  events between the `negative.quantile.min` and `negative.quantile`
  quantiles of the peak channel, excluding the extreme dim tail
  (typically debris) below `negative.quantile.min`.

- x.min.quantile:

  Numeric in `[0, 1)`, default `0.005`. Quantile of the peak-channel
  data used as the x-axis lower limit in panel B. Flow data can have a
  long, sparse negative tail, so this is set independently of
  `negative.quantile.min` to keep the histogram readable.

- singlet.quantiles:

  Numeric, default `c(0.85, 0.975)`. Quantile thresholds used only when
  cleaning a paired bead control (see `control.def.file`), matching
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md).

- gate.color:

  Colour of the octagon gate boundary in panel A, and of the gate box
  drawn around the positive fraction in panel C. Default
  `"darkgoldenrod1"` (matching
  [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)'s
  default).

- density.palette:

  Fill palette for the pseudocolour density in panel A: one of the
  viridis options (`"plasma"`, `"viridis"`, etc.) or any other value to
  use `asp$density.palette.base.color`. Default `"rainbow"` (matching
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)'s
  default).

- selection.fill.color, selection.line.color:

  Fill and line colours for the KDE histogram in panel B. Defaults
  `"steelblue"` / `"black"`.

- negative.bracket.color, positive.bracket.color:

  Colours for the negative- and positive-selection brackets in panel B.
  Defaults `"#377EB8"` (blue) / `"#E41A1C"` (red).

- negative.point.color:

  Colour for the negative-fraction events in panel C, which are plotted
  as plain points with no colour mapping (matching the "Unstained"
  convention in
  [`spectra.automated.steps.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.automated.steps.plot.md)'s
  panel D). Default `"black"`. Only the positive-fraction events are
  coloured by cosine similarity.

- event.point.size:

  Numeric or `NULL` (default). Point size for the positive-fraction
  events in panel C. If `NULL`, defaults to
  `asp$figure.gate.point.size * 1.3`.

- n.highlight:

  Integer, default `200`. Number of positive-fraction events (smallest
  cosine similarity to the negative fraction, i.e. least AF-like)
  highlighted in red in panel A.

- clean.positive.color:

  Colour for the panel A highlight. Default `"red"`.

- clean.positive.point.size:

  Numeric or `NULL` (default). Point size for the panel A highlight. If
  `NULL`, defaults to `asp$figure.gate.point.size * 1.5`.

- cells.trace.color, beads.trace.color, af.trace.color:

  Colours for the three traces in panel D: the population-level
  background-subtracted profile ("Cells"), the reference profile
  ("Beads"), and the negative- fraction trace. Defaults `"#D95F02"` /
  `"#377EB8"` / `"grey40"`.

- max.points:

  Integer. Maximum events plotted per panel (randomly downsampled beyond
  this for speed). Default `5e4`.

- panel.width, panel.height:

  Numeric. Width/height (inches) used per sub-panel when sizing the
  saved composite figure. Defaults `4` and `4`.

- composite.width, composite.height:

  Numeric or `NULL` (default). Override the overall saved figure
  dimensions (inches); if `NULL`, these are computed from `panel.width`
  / `panel.height`.

- output.dir:

  Character or `NULL` (default). Directory to save the composite
  figure(s). Defaults to the current working directory.

- save:

  Logical, default `TRUE`. Whether to save the composite figure for each
  fluorophore to `output.dir`.

- file.type:

  Character string, one of `"jpg"` (default), `"tiff"`, `"png"`, or
  `"pdf"`.

- verbose:

  Logical, default `TRUE`. Print progress messages.

## Value

Invisibly, a named list (one entry per fluorophore) each containing the
individual panel ggplot objects, the assembled `composite` cowplot
object, and the peak / peak-AF channels used.

## See also

[`spectra.automated.steps.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.automated.steps.plot.md),
[`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md),
[`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md),
[`spectral.trace()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.trace.md)
