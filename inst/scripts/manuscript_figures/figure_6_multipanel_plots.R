# figure_6_multipanel_plots.R
#
# Figure 6: representative biplots (A-G) and summary metrics (H, I) across
# five spectral flow datasets (Lung, Spleen, ID7000 PBMC, OMIP-102 PBMC), all
# already-unmixed except panel I, which unmixes raw unstained data five ways
# to compute reconstruction RMSE directly.
#
#   A/B - Lung row:  unstained (A) and fully stained (B), OLS no AF vs
#                     AutoSpectral per-cell AF.
#   C/D - Spleen row: unstained (C) and fully stained (D), same comparison.
#   E   - Lung row:   fully stained, four-way comparison (no AF, one AF,
#                      multiple AF, per-cell AF), files supplied manually.
#   F/G - PBMC row:   ID7000 (WLS no AF vs WLS per-cell AF) and OMIP-102
#                      (OLS no AF vs OLS per-cell AF).
#   H   - Three scatter plots of per-fluorophore rSD (MAD), no-AF extraction
#         vs per-cell AF extraction, one per dataset (Lung, ID7000 PBMC,
#         OMIP-102 PBMC), with a log-linear fit overlay.
#   I   - Grouped bar chart of percent RMSE reduction (relative to no AF
#         correction) across four AF-extraction methods, for five unstained
#         datasets, computed directly via compare.af.methods.rmse().
#
# Every biplot panel (A-G) accepts an optional large.gate pre-gate, computed
# once per panel/block from its own no-AF file and applied to every dataset
# in that block, since they are unmixings of the same underlying events.
#
# Panels A-D, F and G each pair a no-AF export with an independently,
# externally-unmixed per-cell-AF export - like Figure 2's classical vs
# AutoSpectral comparison, these are not guaranteed to share a
# channel-naming convention, so each dataset's own column name for a
# requested fluorophore is resolved independently via
# resolve.fluorophore.channels() (see build.figure6.panel.pair() below).
# Panel H compares a per-fluorophore metric across the same kind of pair, so
# it resolves each dataset's channel names separately too, and passes them
# to compare.channel.metric()'s channels/channels.2 arguments rather than
# assuming both files share names. Panel I computes reconstruction RMSE
# directly from raw data and a single spectra object per sample
# (compare.af.methods.rmse()), with no cross-file column matching involved,
# so no channel resolution step is needed there.
#
# Note on RMSE: current AutoSpectral no longer writes an "$RMSE" FCS header
# keyword (unmix.fcs()'s prior "calculate.error" argument was removed), so
# panel I computes RMSE directly from the raw data and each method's
# reconstruction, via compare.af.methods.rmse() / compute.unmix.rmse(),
# rather than reading it back out of a saved file.

library( AutoSpectral )
library( cowplot )
library( ggplot2 )
library( tidyr )
source("~/Bioinformatics/AutoSpectral/inst/scripts/plot_biplot_annotations.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/resolve_fluorophore_channels.R")


# ---- setup ------------------------------------------------------------------

asp <- get.autospectral.param()

output.dir <- "./figure_6"
if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

# Fluorophore database used to resolve fluorophore identity to each
# dataset's own column names. NULL loads the bundled fluorophore_database.csv.
fig6.fluorophore.database <- NULL

# Point size for the bold panel letters (A-I).
fig6.panel.label.size <- 18

# Point size for each panel's method-descriptive subtitle (e.g.
# "OLS, no AF"). NULL uses asp$figure.axis.title.size.
fig6.subtitle.text.size <- NULL

# Point size for panel H's R^2 annotation. NULL uses that panel's own
# ggplot default.
fig6.h.annotation.text.size <- NULL


# ---- shared helpers -----------------------------------------------------

# Computes a large.gate boundary from the first dataset in data.list (if
# use.large.gate is TRUE) and applies it identically to every dataset in the
# list. Pass use.large.gate = FALSE to skip gating.
apply.optional.large.gate <- function(
    data.list,
    asp,
    use.large.gate = FALSE,
    scatter.param = asp$default.scatter.parameter,
    samp = "panel"
) {

  if ( !use.large.gate ) return( data.list )

  gate.boundary <- compute.scatter.gate(
    flow.data = data.list[[ 1 ]],
    asp = asp,
    scatter.param = scatter.param,
    large.gate = TRUE,
    samp = samp,
    output.dir = file.path( output.dir, "gate_definitions" )
  )

  lapply( data.list, function( d ) {
    apply.gate( d, gate.boundary, scatter.param = scatter.param, asp = asp )
  } )
}

# Square aspect ratio + clip off, matching the Figure 2 panel treatment: keeps
# every biplot square in whatever grid cell cowplot allocates it, and stops
# subtitle text from being clipped at the panel edge when a panel renders
# small.
fig6.square.theme <- theme( aspect.ratio = 1, plot.margin = margin( 10, 10, 10, 10 ) )

apply.fig6.panel.style <- function( plot ) {
  plot <- plot + fig6.square.theme
  plot$coordinates$clip <- "off"
  plot
}

fig6.subtitle.theme <- theme( plot.subtitle = element_text(
  size = if ( is.null( fig6.subtitle.text.size ) ) asp$figure.axis.title.size else fig6.subtitle.text.size
) )

# Builds a no-AF/per-cell-AF panel pair for one block: resolves fluor.x/
# fluor.y to each dataset's own column name independently (see the file
# header), defaults axis labels to the clean fluorophore names (never a raw,
# possibly "-A"-suffixed, channel name), and applies the shared panel style
# and subtitle theme. subtitle.1/subtitle.2 are required so each block's
# wording (e.g. "OLS, no AF" vs "WLS, no AF") stays explicit rather than
# defaulting to something that might not match.
build.figure6.panel.pair <- function(
    data.1,
    data.2,
    fluor.x,
    fluor.y,
    asp,
    subtitle.1,
    subtitle.2,
    x.lab = NULL,
    y.lab = NULL,
    fluorophore.database = fig6.fluorophore.database
) {

  if ( is.null( x.lab ) ) x.lab <- fluor.x
  if ( is.null( y.lab ) ) y.lab <- fluor.y

  channels.1 <- resolve.fluorophore.channels( data.1, c( fluor.x, fluor.y ), fluorophore.database )
  channels.2 <- resolve.fluorophore.channels( data.2, c( fluor.x, fluor.y ), fluorophore.database )

  plot.1 <- create.biplot(
    data.1, x.dim = channels.1[[ fluor.x ]], y.dim = channels.1[[ fluor.y ]], asp = asp,
    x.lab = x.lab, y.lab = y.lab, save = FALSE
  )
  plot.2 <- create.biplot(
    data.2, x.dim = channels.2[[ fluor.x ]], y.dim = channels.2[[ fluor.y ]], asp = asp,
    x.lab = x.lab, y.lab = y.lab, save = FALSE
  )

  plot.1 <- apply.fig6.panel.style( plot.1 ) + labs( subtitle = subtitle.1 ) + fig6.subtitle.theme
  plot.2 <- apply.fig6.panel.style( plot.2 ) + labs( subtitle = subtitle.2 ) + fig6.subtitle.theme

  list( left = plot.1, right = plot.2 )
}

# Assembles a two-plot block (no-AF / per-cell-AF, or similar pairing) with
# an optional title above both plots, and an optional bold panel letter at
# the top-left corner of the block (spanning both of its panels as one
# lettered unit, matching this figure's own A/B, C/D, F/G block scheme).
build.figure6.block <- function( plot.left, plot.right, block.title = NULL, block.letter = NULL ) {

  block.grid <- plot_grid(
    plot.left, plot.right, nrow = 1,
    labels = if ( is.null( block.letter ) ) NULL else c( block.letter, "" ),
    label_size = fig6.panel.label.size,
    label_fontface = "bold"
  )

  if ( is.null( block.title ) ) return( block.grid )

  title.grob <- ggdraw() + draw_label( block.title, size = 14, fontface = "bold" )
  plot_grid( title.grob, block.grid, ncol = 1, rel_heights = c( 0.1, 1 ) )
}

# Assembles a full row of blocks/plots with a row title on the right, and an
# optional bold panel letter on the row's first entry only.
build.figure6.row <- function( panel.plots, row.title, row.letter = NULL ) {

  panel.labels <- if ( is.null( row.letter ) ) NULL else
    c( row.letter, rep( "", length( panel.plots ) - 1 ) )

  row.grid <- plot_grid(
    plotlist = panel.plots, nrow = 1,
    labels = panel.labels, label_size = fig6.panel.label.size, label_fontface = "bold"
  )

  title.grob <- ggdraw() +
    draw_label( row.title, angle = -90, size = 12, fontface = "bold" )

  plot_grid( row.grid, title.grob, ncol = 2, rel_widths = c( 1, 0.06 ) )
}


# ---- Figure 6A/B: Lung -------------------------------------------------------

fig6a.large.gate <- FALSE
fig6b.large.gate <- FALSE

# PLACEHOLDER: unstained lung, no-AF and per-cell-AF unmixed exports.
fig6a.noaf.file <- "PLACEHOLDER_unstained_lung_noAF.fcs"
fig6a.percellaf.file <- "PLACEHOLDER_unstained_lung_percellAF.fcs"

fig6a.data <- apply.optional.large.gate(
  list(
    readFCS( fig6a.noaf.file ),
    readFCS( fig6a.percellaf.file )
  ),
  asp = asp, use.large.gate = fig6a.large.gate, samp = "Figure_6A_unstained_lung"
)

fig6a.pair <- build.figure6.panel.pair(
  fig6a.data[[ 1 ]], fig6a.data[[ 2 ]],
  fluor.x = "BUV496", fluor.y = "KIRAVIA Blue 520",
  asp = asp, subtitle.1 = "OLS, no AF", subtitle.2 = "OLS, per-cell AF"
)

block.6a <- build.figure6.block( fig6a.pair$left, fig6a.pair$right, "Unstained", block.letter = "A" )

# PLACEHOLDER: fully stained lung, no-AF and per-cell-AF unmixed exports.
fig6b.noaf.file <- "PLACEHOLDER_stained_lung_noAF.fcs"
fig6b.percellaf.file <- "PLACEHOLDER_stained_lung_percellAF.fcs"

# PLACEHOLDER: marker names for the axis labels.
fig6b.x.lab <- "PLACEHOLDER_marker BUV496"
fig6b.y.lab <- "PLACEHOLDER_marker KIRAVIA Blue 520"

fig6b.data <- apply.optional.large.gate(
  list(
    readFCS( fig6b.noaf.file ),
    readFCS( fig6b.percellaf.file )
  ),
  asp = asp, use.large.gate = fig6b.large.gate, samp = "Figure_6B_stained_lung"
)

fig6b.pair <- build.figure6.panel.pair(
  fig6b.data[[ 1 ]], fig6b.data[[ 2 ]],
  fluor.x = "BUV496", fluor.y = "KIRAVIA Blue 520",
  asp = asp, subtitle.1 = "OLS, no AF", subtitle.2 = "OLS, per-cell AF",
  x.lab = fig6b.x.lab, y.lab = fig6b.y.lab
)

block.6b <- build.figure6.block( fig6b.pair$left, fig6b.pair$right, "Fully stained", block.letter = "B" )

row.6ab <- build.figure6.row( list( block.6a, block.6b ), "Lung" )


# ---- Figure 6C/D: Spleen -----------------------------------------------------

fig6c.large.gate <- FALSE
fig6d.large.gate <- FALSE

# PLACEHOLDER: unstained spleen, no-AF and per-cell-AF unmixed exports.
fig6c.noaf.file <- "PLACEHOLDER_unstained_spleen_noAF.fcs"
fig6c.percellaf.file <- "PLACEHOLDER_unstained_spleen_percellAF.fcs"

fig6c.data <- apply.optional.large.gate(
  list(
    readFCS( fig6c.noaf.file ),
    readFCS( fig6c.percellaf.file )
  ),
  asp = asp, use.large.gate = fig6c.large.gate, samp = "Figure_6C_unstained_spleen"
)

fig6c.pair <- build.figure6.panel.pair(
  fig6c.data[[ 1 ]], fig6c.data[[ 2 ]],
  fluor.x = "PE-Fire 700", fluor.y = "BV510",
  asp = asp, subtitle.1 = "OLS, no AF", subtitle.2 = "OLS, per-cell AF"
)

block.6c <- build.figure6.block( fig6c.pair$left, fig6c.pair$right, "Unstained", block.letter = "C" )

# PLACEHOLDER: fully stained spleen, no-AF and per-cell-AF unmixed exports.
fig6d.noaf.file <- "PLACEHOLDER_stained_spleen_noAF.fcs"
fig6d.percellaf.file <- "PLACEHOLDER_stained_spleen_percellAF.fcs"

# PLACEHOLDER: marker names for the axis labels.
fig6d.x.lab <- "PLACEHOLDER_marker PE-Fire 700"
fig6d.y.lab <- "PLACEHOLDER_marker BV510"

fig6d.data <- apply.optional.large.gate(
  list(
    readFCS( fig6d.noaf.file ),
    readFCS( fig6d.percellaf.file )
  ),
  asp = asp, use.large.gate = fig6d.large.gate, samp = "Figure_6D_stained_spleen"
)

fig6d.pair <- build.figure6.panel.pair(
  fig6d.data[[ 1 ]], fig6d.data[[ 2 ]],
  fluor.x = "PE-Fire 700", fluor.y = "BV510",
  asp = asp, subtitle.1 = "OLS, no AF", subtitle.2 = "OLS, per-cell AF",
  x.lab = fig6d.x.lab, y.lab = fig6d.y.lab
)

block.6d <- build.figure6.block( fig6d.pair$left, fig6d.pair$right, "Fully stained", block.letter = "D" )

row.6cd <- build.figure6.row( list( block.6c, block.6d ), "Spleen" )


# ---- Figure 6E: Lung, four-way comparison ------------------------------------

fig6e.large.gate <- FALSE

# PLACEHOLDER: fully stained lung, four separately-unmixed exports.
fig6e.noaf.file <- "PLACEHOLDER_stained_lung_noAF.fcs"
fig6e.oneaf.file <- "PLACEHOLDER_stained_lung_oneAF.fcs"
fig6e.multiaf.file <- "PLACEHOLDER_stained_lung_multiAF.fcs"
fig6e.percellaf.file <- "PLACEHOLDER_stained_lung_percellAF.fcs"

fig6e.data <- apply.optional.large.gate(
  list(
    readFCS( fig6e.noaf.file ),
    readFCS( fig6e.oneaf.file ),
    readFCS( fig6e.multiaf.file ),
    readFCS( fig6e.percellaf.file )
  ),
  asp = asp, use.large.gate = fig6e.large.gate, samp = "Figure_6E_stained_lung"
)

fig6e.titles <- c( "OLS, no AF", "OLS, one AF", "OLS, multiple AF", "OLS, per-cell AF" )

fig6e.plots <- Map(
  function( d, plot.title ) {
    channels <- resolve.fluorophore.channels(
      d, c( "BUV563", "BUV615" ), fig6.fluorophore.database
    )
    biplot <- create.biplot(
      d, x.dim = channels[[ "BUV563" ]], y.dim = channels[[ "BUV615" ]], asp = asp,
      x.lab = "CD103 BUV563", y.lab = "Siglec F BUV615", save = FALSE
    )
    apply.fig6.panel.style( biplot ) + labs( subtitle = plot.title ) + fig6.subtitle.theme
  },
  fig6e.data, fig6e.titles
)

row.6e <- build.figure6.row( fig6e.plots, "Lung", row.letter = "E" )


# ---- Figure 6F: ID7000 -------------------------------------------------------

fig6f.large.gate <- FALSE

# PLACEHOLDER: fully stained 40C ID7000 file, WLS no-AF and WLS per-cell-AF
# unmixed exports.
fig6f.noaf.file <- "PLACEHOLDER_ID7000_40C_noAF.fcs"
fig6f.percellaf.file <- "PLACEHOLDER_ID7000_40C_percellAF.fcs"

fig6f.data <- apply.optional.large.gate(
  list(
    readFCS( fig6f.noaf.file ),
    readFCS( fig6f.percellaf.file )
  ),
  asp = asp, use.large.gate = fig6f.large.gate, samp = "Figure_6F_ID7000_40C"
)

fig6f.pair <- build.figure6.panel.pair(
  fig6f.data[[ 1 ]], fig6f.data[[ 2 ]],
  fluor.x = "FITC", fluor.y = "BV510",
  asp = asp, subtitle.1 = "WLS, no AF", subtitle.2 = "WLS, per-cell AF",
  x.lab = "CD45RA FITC", y.lab = "CD16 BV510"
)

block.6f <- build.figure6.block( fig6f.pair$left, fig6f.pair$right, "ID7000", block.letter = "F" )


# ---- Figure 6G: OMIP-102 -----------------------------------------------------

fig6g.large.gate <- FALSE

# PLACEHOLDER: OMIP-102 fully stained file, original (no-AF) and AutoSpectral
# (per-cell-AF) unmixed exports.
fig6g.noaf.file <- "PLACEHOLDER_OMIP-102_stained_noAF.fcs"
fig6g.percellaf.file <- "PLACEHOLDER_OMIP-102_stained_percellAF.fcs"

fig6g.data <- apply.optional.large.gate(
  list(
    readFCS( fig6g.noaf.file ),
    readFCS( fig6g.percellaf.file )
  ),
  asp = asp, use.large.gate = fig6g.large.gate, samp = "Figure_6G_OMIP-102"
)

fig6g.pair <- build.figure6.panel.pair(
  fig6g.data[[ 1 ]], fig6g.data[[ 2 ]],
  fluor.x = "Alexa Fluor 647", fluor.y = "Spark NIR 685",
  asp = asp, subtitle.1 = "OLS, no AF", subtitle.2 = "OLS, per-cell AF",
  x.lab = "Va7.2 Alexa Fluor 647", y.lab = "KLRG1 Spark NIR 685"
)

block.6g <- build.figure6.block( fig6g.pair$left, fig6g.pair$right, "OMIP-102", block.letter = "G" )

row.6fg <- build.figure6.row( list( block.6f, block.6g ), "PBMC" )


# ---- Figure 6H: rSD comparison, no AF vs per-cell AF -------------------------
#
# One scatter plot per dataset, comparing the robust SD (MAD) of every
# non-AF fluorophore channel between the no-AF and per-cell-AF unmixed
# exports, with a log-linear (y ~ a * log(x) + b) fit overlaid. Each
# dataset's own column name for the shared fluorophore list is resolved
# separately (fluorophores may not share a channel-naming convention
# between the two exports - see the file header), and passed to
# compare.channel.metric()'s channels/channels.2 arguments.

# Fits y ~ a * log(x) + b by ordinary least squares on log(x), and returns
# the input data augmented with the fitted curve and its R-squared.
fit.log.curve <- function( rsd.compare ) {

  fit <- lm( PerCellAF ~ log( NoAF ), data = rsd.compare )
  r.squared <- summary( fit )$r.squared

  curve.x <- seq( min( rsd.compare$NoAF ), max( rsd.compare$NoAF ), length.out = 200 )
  curve.data <- data.frame( NoAF = curve.x, PerCellAF = predict( fit, data.frame( NoAF = curve.x ) ) )

  list( curve.data = curve.data, r.squared = r.squared )
}

plot.rsd.comparison <- function( rsd.compare, plot.title, annotation.text.size = NULL ) {

  fitted.curve <- fit.log.curve( rsd.compare )

  text.size <- if ( is.null( annotation.text.size ) ) 3.88 else annotation.text.size

  ggplot( rsd.compare, aes( NoAF, PerCellAF ) ) +
    geom_point( size = 0.5 ) +
    geom_line( data = fitted.curve$curve.data, aes( NoAF, PerCellAF ), color = "red" ) +
    annotate(
      "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
      label = sprintf( "R^2 == %.3f", fitted.curve$r.squared ), parse = TRUE,
      size = text.size
    ) +
    xlab( "No AF" ) +
    ylab( "Per-cell AF" ) +
    ggtitle( plot.title ) +
    theme_classic() +
    theme( axis.text.x = element_text( angle = 45, hjust = 1 ) )
}

fig6h.large.gate <- FALSE

# PLACEHOLDER: spectra CSVs (to identify the non-AF fluorophore channels) and
# no-AF / per-cell-AF unmixed exports, for each of the three datasets.
fig6h.lung.spectra.file <- "PLACEHOLDER_lung_spectra.csv"
fig6h.lung.noaf.file <- "PLACEHOLDER_unstained_lung_noAF.fcs"
fig6h.lung.percellaf.file <- "PLACEHOLDER_unstained_lung_percellAF.fcs"

fig6h.id7000.spectra.file <- "PLACEHOLDER_ID7000_spectra.csv"
fig6h.id7000.noaf.file <- "PLACEHOLDER_ID7000_40C_PBMC_unstained_noAF.fcs"
fig6h.id7000.percellaf.file <- "PLACEHOLDER_ID7000_40C_PBMC_unstained_percellAF.fcs"

fig6h.omip.spectra.file <- "PLACEHOLDER_OMIP-102_spectra.csv"
fig6h.omip.noaf.file <- "PLACEHOLDER_OMIP-102_50C_PBMC_unstained_noAF.fcs"
fig6h.omip.percellaf.file <- "PLACEHOLDER_OMIP-102_50C_PBMC_unstained_percellAF.fcs"

build.fig6h.panel <- function(
    spectra.file, noaf.file, percellaf.file, plot.title, use.large.gate, samp
) {

  fluorophores <- rownames( read.spectra( spectra.file, remove.af = TRUE ) )

  gated <- apply.optional.large.gate(
    list( readFCS( noaf.file ), readFCS( percellaf.file ) ),
    asp = asp, use.large.gate = use.large.gate, samp = samp
  )

  channels.noaf <- resolve.fluorophore.channels( gated[[ 1 ]], fluorophores, fig6.fluorophore.database )
  channels.percellaf <- resolve.fluorophore.channels( gated[[ 2 ]], fluorophores, fig6.fluorophore.database )

  rsd.compare <- compare.channel.metric(
    gated[[ 1 ]], gated[[ 2 ]],
    channels = channels.noaf[ fluorophores ],
    channels.2 = channels.percellaf[ fluorophores ],
    channel.labels = fluorophores,
    method = "rsd", label.1 = "NoAF", label.2 = "PerCellAF"
  )

  plot.rsd.comparison( rsd.compare, plot.title, annotation.text.size = fig6.h.annotation.text.size )
}

fig6h.lung <- build.fig6h.panel(
  fig6h.lung.spectra.file, fig6h.lung.noaf.file, fig6h.lung.percellaf.file,
  "Lung: data rSD", fig6h.large.gate, "Figure_6H_lung"
)
fig6h.id7000 <- build.fig6h.panel(
  fig6h.id7000.spectra.file, fig6h.id7000.noaf.file, fig6h.id7000.percellaf.file,
  "ID7000: data rSD", fig6h.large.gate, "Figure_6H_ID7000"
)
fig6h.omip <- build.fig6h.panel(
  fig6h.omip.spectra.file, fig6h.omip.noaf.file, fig6h.omip.percellaf.file,
  "OMIP-102: data rSD", fig6h.large.gate, "Figure_6H_OMIP-102"
)

row.6h <- plot_grid(
  fig6h.lung, fig6h.id7000, fig6h.omip, nrow = 1,
  labels = c( "H", "", "" ), label_size = fig6.panel.label.size, label_fontface = "bold"
)


# ---- Figure 6I: percent RMSE reduction by AF-extraction method --------------
#
# Unmixes each unstained sample five ways and computes the percent reduction
# in reconstruction RMSE, relative to no AF correction, for the four
# AF-extraction methods, via compare.af.methods.rmse(). Each sample is
# unmixed here against its own single `spectra` object, not read from two
# independently-unmixed exports, so no channel resolution step applies.

# PLACEHOLDER: raw unstained FCS files, fluorophore spectra CSVs, and AF
# spectra CSVs, for each of the five samples.
fig6i.samples <- list(
  ID7000 = list(
    raw.file = "PLACEHOLDER_ID7000_unstained_raw.fcs",
    spectra.file = "PLACEHOLDER_ID7000_spectra.csv",
    af.spectra.file = "PLACEHOLDER_ID7000_AF_spectra.csv"
  ),
  `OMIP-102` = list(
    raw.file = "PLACEHOLDER_OMIP-102_unstained_raw.fcs",
    spectra.file = "PLACEHOLDER_OMIP-102_spectra.csv",
    af.spectra.file = "PLACEHOLDER_OMIP-102_AF_spectra.csv"
  ),
  Brain = list(
    raw.file = "PLACEHOLDER_brain_unstained_raw.fcs",
    spectra.file = "PLACEHOLDER_brain_spectra.csv",
    af.spectra.file = "PLACEHOLDER_brain_AF_spectra.csv"
  ),
  Lung = list(
    raw.file = "PLACEHOLDER_lung_unstained_raw.fcs",
    spectra.file = "PLACEHOLDER_lung_spectra.csv",
    af.spectra.file = "PLACEHOLDER_lung_AF_spectra.csv"
  ),
  Spleen = list(
    raw.file = "PLACEHOLDER_spleen_unstained_raw.fcs",
    spectra.file = "PLACEHOLDER_spleen_spectra.csv",
    af.spectra.file = "PLACEHOLDER_spleen_AF_spectra.csv"
  )
)

fig6i.results <- lapply( names( fig6i.samples ), function( sample.name ) {

  sample.info <- fig6i.samples[[ sample.name ]]

  spectra <- read.spectra( sample.info$spectra.file, remove.af = TRUE )
  af.spectra <- as.matrix( read.spectra( sample.info$af.spectra.file, check.collinearity = FALSE ) )

  raw.data <- readFCS( sample.info$raw.file, columns = colnames( spectra ) )

  rmse.result <- compare.af.methods.rmse( raw.data, spectra, af.spectra )
  rmse.result$Sample <- sample.name

  rmse.result
} )

fig6i.results <- do.call( rbind, fig6i.results )

fig6i.results$Sample <- factor(
  fig6i.results$Sample, levels = c( "ID7000", "OMIP-102", "Brain", "Lung", "Spleen" )
)

# only the three original methods plus the new joint-assignment method are
# shown; the No AF row is the 0% baseline every other method is measured
# against, so it is not itself plotted
fig6i.plot.data <- fig6i.results[ fig6i.results$Method != "NoAF", ]
fig6i.plot.data$Method <- factor(
  fig6i.plot.data$Method,
  levels = c(
    "SingleAF", "PerCellAF_Residuals", "PerCellAF_Fluorophores", "PerCellAF_Joint"
  )
)

ggplot( fig6i.plot.data, aes( x = Sample, y = PercentReduction, fill = Method ) ) +
  geom_col( position = "dodge", color = "black" ) +
  scale_fill_manual(
    values = c( "white", "black", "darkgrey", "grey40" ),
    labels = c(
      "SingleAF", "PerCellAF Residuals", "PerCellAF Fluorophores", "PerCellAF Joint"
    ),
    guide = guide_legend( nrow = 4 )
  ) +
  theme_minimal() +
  scale_y_reverse() +
  labs( x = "Sample", y = "Percent Error Reduction" ) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text( face = "bold", size = 18 ),
    legend.text = element_text( size = 18 ),
    axis.text.x = element_text( size = 14, angle = 45, hjust = 1 ),
    axis.text.y = element_text( size = 14 ),
    axis.title.x = element_text( size = 14 ),
    axis.title.y = element_text( size = 14 )
  ) -> row.6i

row.6i <- plot_grid(
  row.6i, labels = "I", label_size = fig6.panel.label.size, label_fontface = "bold"
)

ggsave(
  file.path( output.dir, "Figure_6I_PercentErrorReduction.jpg" ),
  plot = row.6i, width = 6, height = 6
)


# ---- composite assembly (panels A-H) -----------------------------------------

figure.6.abcdefgh <- plot_grid(
  row.6ab, row.6cd, row.6e, row.6fg, row.6h,
  ncol = 1,
  rel_heights = c( 1, 1, 1, 1, 0.7 )
)

ggsave(
  file.path( output.dir, "Figure_6_A-H.jpg" ),
  plot = figure.6.abcdefgh,
  device = ragg::agg_jpeg,
  width = 16,
  height = 22,
  limitsize = FALSE
)
