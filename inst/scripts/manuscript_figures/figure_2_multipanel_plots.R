# figure_2_multipanel_plots.R
#
# Figure 2: representative biplots comparing classical unmixing against
# AutoSpectral, across an unaffected and an affected channel combination for
# each of four sources of unmixing error.
#
#   A - Unmixing-dependent spread    (unstained spleen)
#   B - Skewing                      (BV605 single-stained control)
#   C - Spillover spread             (CD19 BUV661 single-stained control)
#   D - Autofluorescence             (unstained spleen / lung / brain)
#
# Each row has three panels: Classical Unmixing on an unaffected channel
# combination, Classical Unmixing on an affected channel combination, and
# AutoSpectral on that same affected channel combination. B and C additionally
# overlay gate boxes with a within-gate metric; A and D annotate a whole-plot
# metric for the channels shown. Each row may optionally be pre-gated on
# scatter parameters; when a gate source is supplied, the same gate boundary
# is computed once and applied identically to all three panels in that row.

library( AutoSpectral )
library( cowplot )
library( ggplot2 )


# ---- setup ----------------------------------------------------------------

asp <- get.autospectral.param()

# PLACEHOLDER: directories holding the classically-unmixed and
# AutoSpectral-unmixed FCS exports. Each is expected to contain one export
# per sample below, with fluorophore-named columns (not raw detector names).
classical.dir <- "PLACEHOLDER/Classical unmixed"
autospectral.dir <- "PLACEHOLDER/AutoSpectral unmixed"

output.dir <- "./figure_2"
if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )


# ---- shared row-building helpers -------------------------------------------

# Computes one scatter gate boundary from gate.source and applies it
# identically to every dataset in data.list. Pass gate.source = NULL to skip
# gating for a row.
apply.row.scatter.gate <- function(
    data.list,
    gate.source,
    asp,
    scatter.param = asp$default.scatter.parameter,
    large.gate = FALSE,
    samp = "row"
) {

  if ( is.null( gate.source ) ) return( data.list )

  gate.boundary <- compute.scatter.gate(
    flow.data = gate.source,
    asp = asp,
    scatter.param = scatter.param,
    large.gate = large.gate,
    samp = samp,
    output.dir = file.path( output.dir, "gate_definitions" )
  )

  lapply( data.list, function( d ) {
    apply.gate( d, gate.boundary, scatter.param = scatter.param, asp = asp )
  } )
}

# Builds the three panels (unaffected, affected, AutoSpectral) for one row,
# either with a whole-plot metric (metric.mode = "whole", for A and D) or a
# gate-box metric (metric.mode = "gated", for B and C).
build.figure2.row <- function(
    data.1,
    data.2,
    data.3,
    x.dim.1,
    y.dim.1,
    x.dim.2,
    y.dim.2,
    asp,
    x.lab.1 = NULL,
    y.lab.1 = NULL,
    x.lab.2 = NULL,
    y.lab.2 = NULL,
    metric.mode = c( "whole", "gated" ),
    metric.method = c( "rsd", "mfi", "sd" ),
    gates = NULL,
    metric.channel.gated = NULL,
    x.min = -5000,
    x.max = asp$expr.data.max,
    y.min = -5000,
    y.max = asp$expr.data.max,
    panel.label = NULL
) {

  metric.mode <- match.arg( metric.mode )
  metric.method <- match.arg( metric.method )

  plot.1 <- create.biplot(
    plot.data = data.1, x.dim = x.dim.1, y.dim = y.dim.1, asp = asp,
    x.lab = x.lab.1, y.lab = y.lab.1,
    x.min = x.min, x.max = x.max, y.min = y.min, y.max = y.max,
    save = FALSE
  )
  plot.2 <- create.biplot(
    plot.data = data.2, x.dim = x.dim.2, y.dim = y.dim.2, asp = asp,
    x.lab = x.lab.2, y.lab = y.lab.2,
    x.min = x.min, x.max = x.max, y.min = y.min, y.max = y.max,
    save = FALSE
  )
  plot.3 <- create.biplot(
    plot.data = data.3, x.dim = x.dim.2, y.dim = y.dim.2, asp = asp,
    x.lab = x.lab.2, y.lab = y.lab.2,
    x.min = x.min, x.max = x.max, y.min = y.min, y.max = y.max,
    save = FALSE
  )

  if ( metric.mode == "whole" ) {

    plot.1 <- annotate.biplot.metrics(
      plot.1, data.1, c( x.dim.1, y.dim.1 ), asp, method = metric.method
    )
    plot.2 <- annotate.biplot.metrics(
      plot.2, data.2, c( x.dim.2, y.dim.2 ), asp, method = metric.method
    )
    plot.3 <- annotate.biplot.metrics(
      plot.3, data.3, c( x.dim.2, y.dim.2 ), asp, method = metric.method
    )

  } else {

    for ( gate in gates ) {
      plot.2 <- add.biplot.gated.metric(
        plot = plot.2, data = data.2, x.dim = x.dim.2, y.dim = y.dim.2,
        x.range = gate$x.range, y.range = gate$y.range,
        metric.channel = metric.channel.gated, asp = asp,
        method = metric.method, box.color = gate$color,
        x.max = x.max, y.max = y.max
      )
      plot.3 <- add.biplot.gated.metric(
        plot = plot.3, data = data.3, x.dim = x.dim.2, y.dim = y.dim.2,
        x.range = gate$x.range, y.range = gate$y.range,
        metric.channel = metric.channel.gated, asp = asp,
        method = metric.method, box.color = gate$color,
        x.max = x.max, y.max = y.max
      )
    }
  }

  if ( !is.null( panel.label ) ) {
    subtitle.theme <- theme( plot.subtitle = element_text( size = asp$figure.axis.title.size ) )
    plot.1 <- plot.1 + labs( subtitle = panel.label ) + subtitle.theme
    plot.2 <- plot.2 + labs( subtitle = panel.label ) + subtitle.theme
    plot.3 <- plot.3 + labs( subtitle = panel.label ) + subtitle.theme
  }

  list( unaffected = plot.1, affected = plot.2, autospectral = plot.3 )
}

# Assembles one row's three panels with a rotated row title on the left.
assemble.figure2.row <- function( row.plots, row.title ) {

  row.grid <- plot_grid(
    row.plots$unaffected, row.plots$affected, row.plots$autospectral,
    nrow = 1
  )

  title.grob <- ggdraw() +
    draw_label( row.title, angle = 90, size = 12, fontface = "bold" )

  plot_grid( title.grob, row.grid, ncol = 2, rel_widths = c( 0.06, 1 ) )
}


# ---- Figure 2A: Unmixing-dependent spread ----------------------------------

# PLACEHOLDER: unstained spleen file names in each unmixed directory.
fig2a.classical.file <- "PLACEHOLDER_unstained_spleen_classical.fcs"
fig2a.autospectral.file <- "PLACEHOLDER_unstained_spleen_autospectral.fcs"

fig2a.classical.data <- readFCS( file.path( classical.dir, fig2a.classical.file ) )
fig2a.autospectral.data <- readFCS( file.path( autospectral.dir, fig2a.autospectral.file ) )

# Optional scatter gate for row A. Set to NULL to skip gating this row.
fig2a.gate.source <- NULL

fig2a.gated <- apply.row.scatter.gate(
  data.list = list( fig2a.classical.data, fig2a.classical.data, fig2a.autospectral.data ),
  gate.source = fig2a.gate.source,
  asp = asp,
  samp = "Figure_2A_unstained_spleen"
)

fig2a.row <- build.figure2.row(
  data.1 = fig2a.gated[[ 1 ]], data.2 = fig2a.gated[[ 2 ]], data.3 = fig2a.gated[[ 3 ]],
  x.dim.1 = "BUV805", y.dim.1 = "APC-Fire 810",
  x.dim.2 = "BV510", y.dim.2 = "GFP",
  asp = asp,
  metric.mode = "whole",
  metric.method = "rsd"
)

fig2a <- assemble.figure2.row( fig2a.row, "Unmixing-dependent spread" )


# ---- Figure 2B: Skewing -----------------------------------------------------

# PLACEHOLDER: BV605 single-stained control file names.
fig2b.classical.file <- "PLACEHOLDER_BV605_classical.fcs"
fig2b.autospectral.file <- "PLACEHOLDER_BV605_autospectral.fcs"

fig2b.classical.data <- readFCS( file.path( classical.dir, fig2b.classical.file ) )
fig2b.autospectral.data <- readFCS( file.path( autospectral.dir, fig2b.autospectral.file ) )

fig2b.gate.source <- NULL

fig2b.gated <- apply.row.scatter.gate(
  data.list = list( fig2b.classical.data, fig2b.classical.data, fig2b.autospectral.data ),
  gate.source = fig2b.gate.source,
  asp = asp,
  samp = "Figure_2B_BV605"
)

# Negative (blue) and positive (red) gate boxes. Coordinates are modifiable.
fig2b.gates <- list(
  list( x.range = c( -1e4, 1e4 ), y.range = c( -3e3, 3e3 ), color = "blue" ),
  list( x.range = c( -1e4, 1e4 ), y.range = c( 2e4, 5e5 ), color = "red" )
)

fig2b.row <- build.figure2.row(
  data.1 = fig2b.gated[[ 1 ]], data.2 = fig2b.gated[[ 2 ]], data.3 = fig2b.gated[[ 3 ]],
  x.dim.1 = "BUV805", y.dim.1 = "BV605",
  x.dim.2 = "BV650", y.dim.2 = "BV605",
  asp = asp,
  metric.mode = "gated",
  metric.method = "mfi",
  gates = fig2b.gates,
  metric.channel.gated = "BV650",
  x.min = -5e4
)

fig2b <- assemble.figure2.row( fig2b.row, "Skewing" )


# ---- Figure 2C: Spillover spread -------------------------------------------

# PLACEHOLDER: CD19 BUV661 single-stained control file names.
fig2c.classical.file <- "PLACEHOLDER_BUV661_classical.fcs"
fig2c.autospectral.file <- "PLACEHOLDER_BUV661_autospectral.fcs"

fig2c.classical.data <- readFCS( file.path( classical.dir, fig2c.classical.file ) )
fig2c.autospectral.data <- readFCS( file.path( autospectral.dir, fig2c.autospectral.file ) )

fig2c.gate.source <- NULL

fig2c.gated <- apply.row.scatter.gate(
  data.list = list( fig2c.classical.data, fig2c.classical.data, fig2c.autospectral.data ),
  gate.source = fig2c.gate.source,
  asp = asp,
  samp = "Figure_2C_BUV661"
)

# Near-origin (blue) and positive-population (red) gate boxes. Coordinates
# are modifiable; x-axis is deliberately wide (x.min below) so both gates'
# edges remain visible.
fig2c.gates <- list(
  list( x.range = c( -3e3, 3e3 ), y.range = c( -3e3, 3e3 ), color = "blue" ),
  list( x.range = c( -3e4, 3e4 ), y.range = c( 2e4, 3e5 ), color = "red" )
)

fig2c.row <- build.figure2.row(
  data.1 = fig2c.gated[[ 1 ]], data.2 = fig2c.gated[[ 2 ]], data.3 = fig2c.gated[[ 3 ]],
  x.dim.1 = "BUV395", y.dim.1 = "BUV661",
  x.dim.2 = "APC", y.dim.2 = "BUV661",
  asp = asp,
  metric.mode = "gated",
  metric.method = "rsd",
  gates = fig2c.gates,
  metric.channel.gated = "APC",
  x.min = -5e4
)

fig2c <- assemble.figure2.row( fig2c.row, "Spillover spread" )


# ---- Figure 2D: Autofluorescence -------------------------------------------
#
# Two sub-rows: spleen/lung (row 1) and brain/brain (row 2). Each sub-row
# follows the same unaffected/affected/AutoSpectral pattern as A-C, but with
# a whole-plot SD (not rSD) annotation and a tissue label on every panel.

# PLACEHOLDER: unstained tissue file names.
fig2d.spleen.classical.file <- "PLACEHOLDER_unstained_spleen_classical.fcs"
fig2d.lung.classical.file <- "PLACEHOLDER_unstained_lung_classical.fcs"
fig2d.lung.autospectral.file <- "PLACEHOLDER_unstained_lung_autospectral.fcs"
fig2d.brain.classical.file <- "PLACEHOLDER_unstained_brain_classical.fcs"
fig2d.brain.autospectral.file <- "PLACEHOLDER_unstained_brain_autospectral.fcs"

fig2d.spleen.data <- readFCS( file.path( classical.dir, fig2d.spleen.classical.file ) )
fig2d.lung.classical.data <- readFCS( file.path( classical.dir, fig2d.lung.classical.file ) )
fig2d.lung.autospectral.data <- readFCS( file.path( autospectral.dir, fig2d.lung.autospectral.file ) )
fig2d.brain.classical.data <- readFCS( file.path( classical.dir, fig2d.brain.classical.file ) )
fig2d.brain.autospectral.data <- readFCS( file.path( autospectral.dir, fig2d.brain.autospectral.file ) )

fig2d.row1.gate.source <- NULL
fig2d.row2.gate.source <- NULL

fig2d.row1.gated <- apply.row.scatter.gate(
  data.list = list( fig2d.spleen.data, fig2d.lung.classical.data, fig2d.lung.autospectral.data ),
  gate.source = fig2d.row1.gate.source,
  asp = asp,
  samp = "Figure_2D_row1_spleen_lung"
)

fig2d.row2.gated <- apply.row.scatter.gate(
  data.list = list(
    fig2d.brain.classical.data, fig2d.brain.classical.data, fig2d.brain.autospectral.data
  ),
  gate.source = fig2d.row2.gate.source,
  asp = asp,
  samp = "Figure_2D_row2_brain"
)

fig2d.row1 <- build.figure2.row(
  data.1 = fig2d.row1.gated[[ 1 ]], data.2 = fig2d.row1.gated[[ 2 ]], data.3 = fig2d.row1.gated[[ 3 ]],
  x.dim.1 = "Real Blue 744", y.dim.1 = "APC-Fire 810",
  x.dim.2 = "BUV496", y.dim.2 = "BUV563",
  asp = asp,
  metric.mode = "whole",
  metric.method = "sd"
)

# panel.label differs between the unaffected panel (Spleen) and the
# affected/AutoSpectral panels (Lung), so it is added per-panel rather than
# through build.figure2.row()'s single panel.label argument.
subtitle.theme <- theme( plot.subtitle = element_text( size = asp$figure.axis.title.size ) )
fig2d.row1$unaffected <- fig2d.row1$unaffected + labs( subtitle = "Spleen" ) + subtitle.theme
fig2d.row1$affected <- fig2d.row1$affected + labs( subtitle = "Lung" ) + subtitle.theme
fig2d.row1$autospectral <- fig2d.row1$autospectral + labs( subtitle = "Lung" ) + subtitle.theme

fig2d.row2 <- build.figure2.row(
  data.1 = fig2d.row2.gated[[ 1 ]], data.2 = fig2d.row2.gated[[ 2 ]], data.3 = fig2d.row2.gated[[ 3 ]],
  x.dim.1 = "PE-Vio 770", y.dim.1 = "Real Blue 780",
  x.dim.2 = "BUV615", y.dim.2 = "BV711",
  asp = asp,
  metric.mode = "whole",
  metric.method = "sd",
  panel.label = "Brain"
)

fig2d <- plot_grid(
  assemble.figure2.row( fig2d.row1, "Autofluorescence" ),
  assemble.figure2.row( fig2d.row2, "" ),
  ncol = 1
)


# ---- composite assembly -----------------------------------------------------

column.headers <- ggdraw() +
  draw_label( "Classical Unmixing", x = 0.35, size = 14, fontface = "bold" ) +
  draw_label( "AutoSpectral", x = 0.85, size = 14, fontface = "bold" )

figure.2 <- plot_grid(
  column.headers,
  fig2a, fig2b, fig2c, fig2d,
  ncol = 1,
  rel_heights = c( 0.05, 1, 1, 1, 2 )
)

ggsave(
  file.path( output.dir, "Figure_2.jpg" ),
  plot = figure.2,
  device = ragg::agg_jpeg,
  width = 12,
  height = 20,
  limitsize = FALSE
)
