# Plot Legacy Spectra Extraction Pipeline Steps

Builds a manuscript-ready, multi-panel figure illustrating the internal
steps of the legacy workflow
([`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md),
[`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md),
[`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md))
for one or more single-stained control samples: the automated scatter
gate, autofluorescence-exclusion gating, scatter-matched
universal-negative selection, and the robust linear model fit used to
extract the fluorophore's spectral signature. Uses the same building
blocks as the rest of AutoSpectral
([`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)-style
biexponential biplots,
[`scatter.match.plot()`](https://drcytometer.github.io/AutoSpectral/reference/scatter.match.plot.md))
so panel styling matches the package's other figures.

This function runs the full legacy pipeline itself – it calls
[`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
and
[`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)
internally for the entire control set in `control.def.file`, not just
the illustrated fluorophore(s) – so it can be slow for large panels.
`parallel = FALSE` is used internally throughout, since diagnostic
capture is unreliable under forked parallel processing.

Requires `plot_spectra_automated_steps.R` and
`plot_spectra_standard_workflow.R` to be loaded in the same package
namespace (several private helpers are shared), and requires the
`diagnostics.env`-aware versions of
[`remove.af()`](https://drcytometer.github.io/AutoSpectral/reference/remove.af.md)
/
[`run.af.removal()`](https://drcytometer.github.io/AutoSpectral/reference/run.af.removal.md)
/
[`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).

## Usage

``` r
spectra.legacy.steps.plot(
  control.dir,
  control.def.file,
  asp,
  fluorophores = NULL,
  gating.system = c("density", "landmarks"),
  af.remove = TRUE,
  universal.negative = TRUE,
  downsample = TRUE,
  scatter.match = TRUE,
  k.neighbors = 3L,
  negative.n = asp$negative.n,
  positive.n = asp$positive.n,
  gate.color = "darkgoldenrod1",
  density.palette = "rainbow",
  unstained.point.color = "black",
  cosine.point.size = NULL,
  af.gate.color = "black",
  clean.positive.color = "red",
  clean.positive.point.size = NULL,
  rlm.line.color = "blue",
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

  Character. Path to (or filename of) the control definition CSV, in the
  full legacy format required by
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  (including `control.type`, `gate.name`, `gate.define`, etc. – see
  [`check.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/check.control.file.md)).

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- fluorophores:

  Character vector of fluorophore name(s) to illustrate. Default `NULL`
  illustrates the first cell-based fluorophore with a paired universal
  negative.

- gating.system:

  Character, one of `"density"` (default) or `"landmarks"`, matching
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)'s
  argument of the same name.

- af.remove:

  Logical, default `TRUE`. Passed to
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).
  Panels B and D require this to be `TRUE` and require the illustrated
  fluorophore to have a paired universal negative; otherwise those
  panels show a placeholder.

- universal.negative, downsample, scatter.match, k.neighbors,
  negative.n, positive.n:

  Passed through to
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md).
  See that function's documentation.

- gate.color:

  Colour of the panel A gate boundary. Default `"darkgoldenrod1"`
  (matching
  [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)'s
  default).

- density.palette:

  Fill palette for the panel A pseudocolour density. Default
  `"rainbow"`.

- unstained.point.color:

  Colour for the unstained/AF events in panel B. Default `"black"`.

- cosine.point.size:

  Numeric or `NULL` (default). Point size for the single-stained control
  events in panel B. If `NULL`, defaults to
  `asp$figure.gate.point.size * 1.3`.

- af.gate.color:

  Colour of the AF-exclusion gate boundary drawn on both panel B
  biplots. Default `"black"`.

- clean.positive.color:

  Colour for the highlighted "clean" (AF-gate- excluded) events in
  panel D. Default `"red"`.

- clean.positive.point.size:

  Numeric or `NULL` (default). Point size for the panel D highlight. If
  `NULL`, defaults to `asp$figure.gate.point.size * 1.5`.

- rlm.line.color:

  Colour of the robust-linear-model fit line in panel D. Default
  `"blue"`.

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

  Logical, default `TRUE`. Print progress messages (also controls
  verbosity of the internal
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  /
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)
  calls).

## Value

Invisibly, a named list (one entry per fluorophore) each containing the
individual panel ggplot objects, the assembled `composite` cowplot
object, and the resolved gate/channel names.

## See also

[`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md),
[`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md),
[`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md),
[`spectra.automated.steps.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.automated.steps.plot.md),
[`spectra.standard.workflow.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.standard.workflow.plot.md)
