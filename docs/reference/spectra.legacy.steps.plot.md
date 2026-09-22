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
  gate.list = NULL,
  af.remove = TRUE,
  universal.negative = TRUE,
  downsample = TRUE,
  scatter.match = TRUE,
  k.neighbors = 3L,
  negative.n = asp$negative.n,
  positive.n = asp$positive.n,
  singlet.quantiles = c(0.85, 0.975),
  color.palette = NULL,
  gate.color = "darkgoldenrod1",
  density.palette = "rainbow",
  unstained.point.color = "black",
  cosine.point.size = NULL,
  af.gate.color = "black",
  clean.positive.color = "red",
  clean.positive.point.size = NULL,
  n.true.positive = 50L,
  rlm.line.color = "blue",
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
  verbose = TRUE,
  allow.duplicate.controls = TRUE
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

- gate.list:

  Optional named list of gates. To use this, pre-define the gates using
  [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)
  and/or
  [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md),
  ensure that the names of the gates correspond to the names in the
  `control.def.file`, and ensure that the `gate.name` column has been
  filled in for the `control.def.file`. Default `NULL` will revert to
  creating new gates. Passed through to
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  for the real pipeline run, and also reused directly for panel A's
  re-derived gate boundary (rather than recomputing it), so the figure
  shows the same gate the real run used.

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

- singlet.quantiles:

  Numeric, default `c(0.85, 0.975)`. Quantile thresholds for the
  two-stage FSC/SSC singlet discrimination used only when cleaning a
  paired bead control (see `control.def.file`), matching
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md).

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Use `rainbow` to be similar to FlowJo
  or SpectroFlo. Other options are the viridis color options: `magma`,
  `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako` and
  `turbo`.

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

  Colour for the highlighted "true positive" events in panels A and D.
  Default `"red"`.

- clean.positive.point.size:

  Numeric or `NULL` (default). Point size for the panels A/D highlight.
  If `NULL`, defaults to `asp$figure.gate.point.size * 1.5`.

- n.true.positive:

  Integer, default `50L`. Number of "true positive" events highlighted
  in red in panels A and D: the brightest `n.true.positive` events among
  the AF-gate-excluded population, ranked by projection onto the fitted
  RLM trend direction in (peak channel, intrusive-AF channel) space,
  rather than every AF-gate-excluded event (which is simply "not AF",
  not "positively stained").

- rlm.line.color:

  Colour of the robust-linear-model fit line in panel D. Default
  `"blue"`.

- cells.trace.color, beads.trace.color, af.trace.color:

  Colours for the three traces in panel E: the RLM-based per-channel
  signature ("Cells"), the reference profile ("Beads"), and the
  matched-negative AF trace. Defaults `"#D95F02"` / `"#377EB8"` /
  `"grey40"`, matching
  [`spectra.automated.steps.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.automated.steps.plot.md)'s
  panel F and
  [`spectra.standard.workflow.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.standard.workflow.plot.md)'s
  panel D.

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

- allow.duplicate.controls:

  Logical, default `TRUE`. Set `TRUE` to permit multiple single-stained
  controls for the same fluorophore (diagnostic/QC use only). Each is
  tracked internally under a unique `sample` identifier. The resulting
  spectral reference library still needs to be reduced to one row per
  fluorophore before unmixing – see
  [`check.spectra.duplicates()`](https://drcytometer.github.io/AutoSpectral/reference/check.spectra.duplicates.md).

## Value

Invisibly, a named list (one entry per fluorophore), each containing:

- `gate.panel`:

  Panel A, the automated scatter gate (or a placeholder if gate
  definition failed), with the `n.true.positive`
  brightest-along-the-RLM-trend events highlighted larger in red when
  AF-removal diagnostics were available.

- `af.panel`:

  Panel B, the AF-exclusion cosine-similarity biplot (or a placeholder
  if `af.remove = FALSE` or no paired universal negative was available
  for this fluorophore).

- `scatter.match.panel`:

  Panel C, the
  [`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md)
  kNN scatter-match figure embedded from its saved JPEG.

- `rlm.panel`:

  Panel D, the robust-linear-model diagnostic, fit to and displaying
  `flow.control$clean.expr` for this sample (the events as
  clean.controls() actually finalises them, not just
  `gate.population.idx`), or a placeholder alongside `af.panel` when
  AF-removal diagnostics were unavailable or too few clean.controls()
  events remained.

- `subtraction.plot`:

  Panel E, the final spectral profile comparison
  ([`spectral.trace()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.trace.md)
  of Cells / Beads / AF), fit to the same `flow.control$clean.expr`
  population as `rlm.panel`, or a placeholder when AF-removal
  diagnostics were unavailable, too few clean.controls() events
  remained, or RLM extraction failed.

- `composite`:

  The assembled five-panel cowplot object saved to `output.dir` when
  `save = TRUE`.

- `gate.name`:

  Character. The `gate.name` resolved for this fluorophore's sample, or
  `NA` if none was assigned.

- `af.peak.channel`:

  Character. The intrusive-AF channel used as panels B/D's y-axis, or
  `NA_character_` if AF-removal diagnostics were unavailable.

- `fluor.peak`:

  Character. The fluorophore's peak channel used as panels B/D's x-axis,
  or `NA_character_` if AF-removal diagnostics were unavailable.

- `reference.profile`:

  Named numeric vector (over the panel-wide spectral channels) used as
  the "Beads" trace in panel E, or `NULL` if neither a paired bead
  control nor the spectral reference library had data for this
  fluorophore.

## See also

[`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md),
[`clean.controls()`](https://drcytometer.github.io/AutoSpectral/reference/clean.controls.md),
[`get.fluorophore.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.fluorophore.spectra.md),
[`spectra.automated.steps.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.automated.steps.plot.md),
[`spectra.standard.workflow.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectra.standard.workflow.plot.md)
