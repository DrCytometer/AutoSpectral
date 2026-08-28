# Plot Automated Spectra Extraction Pipeline Steps

Builds a manuscript-ready, multi-panel figure illustrating the internal
steps of
[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)
for one or more single-stained control samples: singlet gating (linear
scale), AF-orthogonalisation peak-finding, brightest-event candidate
selection (KDE-smoothed histogram), the cosine-similarity filter
(biexponential biplot, continuous colour gradient), and the kNN
scatter-matched AF subtraction. Uses the same building blocks as the
rest of AutoSpectral
([`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md)-style
biexponential biplots,
[`spectral.trace()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.trace.md),
[`scatter.match.plot()`](https://drcytometer.github.io/AutoSpectral/reference/scatter.match.plot.md))
so panel styling matches the package's other figures.

## Usage

``` r
spectra.automated.steps.plot(
  control.dir,
  control.def.file,
  asp,
  fluorophores = NULL,
  n.candidates = 1000L,
  n.spectral = 200L,
  k.neighbors = 2L,
  singlet.quantiles = c(0.85, 0.975),
  trace.colors = c(AF = "grey40", `Pre-orthogonalization` = "#1B9E77",
    `Post-orthogonalization` = "#D95F02"),
  cells.trace.color = "#D95F02",
  beads.trace.color = "#377EB8",
  af.trace.color = "grey40",
  singlet.retained.color = "black",
  singlet.excluded.color = "grey75",
  clean.positive.color = "red",
  clean.positive.point.size = NULL,
  brightest.fill.color = "steelblue",
  brightest.line.color = "black",
  brightest.panel.width = NULL,
  brightest.panel.height = NULL,
  unstained.point.color = "black",
  cosine.point.size = NULL,
  selected.size.mult = 5,
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

  Character. Path to (or filename of) the control definition CSV, as
  used by
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md).
  If the file has a `control.type` column (as used by
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  /
  [`check.control.file()`](https://drcytometer.github.io/AutoSpectral/reference/check.control.file.md)),
  a row with `control.type == "beads"` sharing a fluorophore name with a
  "cells" row is treated as a paired bead control: its `filename` is the
  positive bead sample, and its `universal.negative` must point at the
  matching unstained bead file. When present, this pair supplies the "
  (Beads)" trace in panel F instead of the static spectral reference
  library.

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- fluorophores:

  Character vector of fluorophore name(s) (as they appear in the
  `fluorophore` column of the control file) to illustrate. Default
  `NULL` illustrates the first external-negative fluorophore found.
  Fluorophores using internal-negative mode (no `universal.negative`
  entry) are skipped with a warning, since orthogonalisation /
  cosine-filter / kNN steps do not apply to them.

- n.candidates:

  Integer, default `1000`. Number of top-expressing candidate events
  shown as the "brightest events" selection (panel C) and carried into
  the cosine-filter biplot (panel D), matching the `n.candidates`
  argument of
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md).

- n.spectral:

  Integer, default `200`. Number of the above candidates that would be
  retained after cosine filtering; these are drawn at full opacity, the
  remainder at reduced opacity, matching
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)'s
  `n.spectral` argument.

- k.neighbors:

  Integer, default `2`. Number of nearest neighbours used for the kNN
  scatter-match panel (panel E), matching
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)'s
  `k.neighbors` argument.

- singlet.quantiles:

  Numeric, default `c(0.85, 0.975)`. Quantile thresholds for the
  two-stage FSC/SSC singlet discrimination shown in panel A, matching
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md).

- trace.colors:

  Named character vector of colours for the three traces in panel B.
  Names must be `"AF"`, `"Pre-orthogonalization"` and
  `"Post-orthogonalization"`.

- cells.trace.color, beads.trace.color, af.trace.color:

  Colours for the three traces in panel F: the calculated per-cell
  profile ("Cells"), the reference profile ("Beads"), and the
  matched-unstained AF trace. Defaults `"#D95F02"` / `"#377EB8"` /
  `"grey40"`.

- singlet.retained.color, singlet.excluded.color:

  Colours for retained and excluded events in panel A. Defaults
  `"black"` / `"grey75"`.

- clean.positive.color:

  Colour for the highlighted "clean positive" events overlaid on
  panel A. Default `"red"`.

- clean.positive.point.size:

  Numeric or `NULL` (default). Point size for the panel A highlight. If
  `NULL`, defaults to `asp$figure.gate.point.size * 1.5`.

- brightest.fill.color, brightest.line.color:

  Fill and line colours for the KDE-smoothed histogram in panel C.
  Defaults `"steelblue"` / `"black"`.

- brightest.panel.width, brightest.panel.height:

  Numeric or `NULL` (default). Width/height (inches) of panel C. If
  `NULL`, defaults to `panel.width * 2` / `panel.height`. When narrower
  than the overall composite width, the panel is centred with blank
  padding on either side.

- unstained.point.color:

  Colour for unstained events in panel D. Default `"black"`.

- cosine.point.size:

  Numeric or `NULL` (default). Point size for the non-selected
  single-stained control events in panel D (these are often a small,
  pre-selected "brightest events" pool, so the default point size can be
  too small to see clearly). If `NULL`, defaults to
  `asp$figure.gate.point.size * 1.3`.

- selected.size.mult:

  Numeric, default `5`. Events selected for the final profile are drawn
  in panel D at `cosine.point.size` times this multiple, in place of a
  gate polygon (nothing is actually excluded from the plot – the
  selection feeds forward as an average).

- max.points:

  Integer. Maximum events plotted per panel (randomly downsampled beyond
  this for speed). Default `5e4`.

- panel.width, panel.height:

  Numeric. Width/height (inches) used per sub-panel when sizing the
  saved composite figure. Defaults `4` and `4`.

- composite.width, composite.height:

  Numeric or `NULL` (default). Override the overall saved figure
  dimensions (inches); if `NULL`, these are computed from `panel.width`
  / `panel.height` (and `brightest.panel.height`, if given).

- output.dir:

  Character or `NULL` (default). Directory to save the composite
  figure(s) and the kNN scatter-match panel. Defaults to the current
  working directory.

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
object, and the empirical peak / AF-collision channels used.

## See also

[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md),
[`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md),
[`spectral.trace()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.trace.md),
[`scatter.match.plot()`](https://drcytometer.github.io/AutoSpectral/reference/scatter.match.plot.md)
