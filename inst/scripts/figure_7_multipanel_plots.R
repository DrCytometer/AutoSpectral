# figure_7_multipanel_plots.R
#
# Figure 7: representative spectral variant and per-cell fluorophore
# optimization biplots across several tissues and instruments.
#
# A) BUV661 spectral variant plot, from a cached get.spectral.variants()
#    .rds file for the spectral symposium data.
# B) (right of A) Fully-stained spleen, ungated: OLS vs per-cell
#    fluorophore optimization (AutoSpectral), APC vs BUV661.
# C) (below A) Mouse lung, ungated (optional large.gate), already-unmixed
#    files: OLS vs AutoSpectral, BUV615 vs SBB615.
# D) (below E) OMIP-102, ungated (optional large.gate), already-unmixed
#    files: Chorus original vs AutoSpectral, PE vs PE-Fire 810.
# E) (right of C/D, below B) 2x2 grid, mouse lung WT vs GFP mouse x
#    OLS vs AutoSpectral, GFP vs BV510, with power-of-ten gridlines.
#
# Panels C, D and E each compare two or more independently, externally
# unmixed export files - like Figure 2's classical/AutoSpectral exports,
# these are not guaranteed to share a channel-naming convention, so each
# file's own column name for a requested fluorophore is resolved
# independently via resolve.fluorophore.channels(). Panel B instead unmixes
# both arms (OLS and AutoSpectral) in this script against the same `spectra`
# object, so its two results are guaranteed to share column names already -
# no resolution is needed or performed there.
#
# Every file path below is a PLACEHOLDER for you to fill in.

library( AutoSpectral )
library( cowplot )
library( ggplot2 )
source("~/Bioinformatics/AutoSpectral/inst/scripts/plot_biplot_annotations.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/resolve_fluorophore_channels.R")

asp <- get.autospectral.param()
output.dir <- "./figure_7"
if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

# Fluorophore database used to resolve fluorophore identity to each
# dataset's own column names in panels C, D and E. NULL loads the bundled
# fluorophore_database.csv.
fig7.fluorophore.database <- NULL

# Point size for the bold panel letters (A-E).
fig7.panel.label.size <- 16

# Point size for each panel's method-descriptive subtitle (e.g. "OLS",
# "AutoSpectral"). NULL uses asp$figure.axis.title.size.
fig7.subtitle.text.size <- NULL

# ---------------------------------------------------------------------------
# Shared spectral variants (Panel A source, and Panel B's AutoSpectral unmix)
# ---------------------------------------------------------------------------

spectra.variants.rds.path <- "PLACEHOLDER_SPECTRAL_VARIANTS.rds"
spectra.variants <- readRDS( spectra.variants.rds.path )

spectra.file <- "PLACEHOLDER_SPECTRA_FILE.csv"
af.spectra.file <- "PLACEHOLDER_AF_SPECTRA_FILE.csv"

spectra <- read.spectra( spectra.file )
spectral.channel <- colnames( spectra )
af.spectra <- read.spectra( af.spectra.file, check.collinearity = FALSE )

# ---------------------------------------------------------------------------
# Local layout helpers
# ---------------------------------------------------------------------------

# Square aspect ratio + clip off, matching the Figure 2 panel treatment: keeps
# every biplot square in whatever grid cell cowplot allocates it, and stops
# subtitle/annotation text from being clipped at the panel edge when a panel
# renders small.
fig7.square.theme <- theme( aspect.ratio = 1, plot.margin = margin( 10, 10, 10, 10 ) )

apply.fig7.panel.style <- function( plot ) {
  plot <- plot + fig7.square.theme
  plot$coordinates$clip <- "off"
  plot
}

# Adds a method-descriptive subtitle at a configurable text size, replacing
# a plain ggtitle()/labs(subtitle=) call with no size control.
apply.fig7.subtitle <- function( plot, subtitle.text, asp, text.size = fig7.subtitle.text.size ) {
  subtitle.size <- if ( is.null( text.size ) ) asp$figure.axis.title.size else text.size
  plot +
    labs( subtitle = subtitle.text ) +
    theme( plot.subtitle = element_text( size = subtitle.size ) )
}

build.figure7.titled.pair <- function( plot.left, plot.right, block.title = NULL, block.letter = NULL ) {

  pair.row <- plot_grid(
    plot.left, plot.right, nrow = 1,
    labels = if ( is.null( block.letter ) ) NULL else c( block.letter, "" ),
    label_size = fig7.panel.label.size,
    label_fontface = "bold"
  )

  if ( is.null( block.title ) ) return( pair.row )

  title.row <- ggdraw() + draw_label( block.title, fontface = "bold", size = 12 )

  plot_grid( title.row, pair.row, ncol = 1, rel_heights = c( 0.12, 1 ) )
}

# ---------------------------------------------------------------------------
# Panel A: BUV661 spectral variant plot
# ---------------------------------------------------------------------------

fig7a.target.fluor <- "BUV661"
fig7a.variant.matrix <- spectra.variants$variants[[ fig7a.target.fluor ]]

panel.A <- spectral.variant.plot.dens(
  spectra.variants = fig7a.variant.matrix,
  median.spectrum = fig7a.variant.matrix[ 1, ],
  title = paste0( fig7a.target.fluor, " Spectral Variants" ),
  save = FALSE
)

# ---------------------------------------------------------------------------
# Panel B: fully-stained spleen, ungated, OLS vs AutoSpectral
#
# Both unmixed results are produced here, in this script, against the same
# `spectra` object - so unmixed.ols and unmixed.autospectral are guaranteed
# to share column names, and no channel resolution is needed.
# ---------------------------------------------------------------------------

fig7b.spleen.file <- "PLACEHOLDER_SPLEEN_FULLY_STAINED.fcs"

fig7b.raw.data <- readFCS( fig7b.spleen.file, columns = spectral.channel )

fig7b.unmixed.ols <- unmix.ols( fig7b.raw.data, spectra )

fig7b.unmixed.autospectral <- unmix.autospectral.joint(
  raw.data = fig7b.raw.data,
  spectra = spectra,
  af.spectra = af.spectra,
  asp = asp,
  spectra.variants = spectra.variants,
  verbose = FALSE
)

fig7b.left <- create.biplot(
  fig7b.unmixed.ols, x.dim = "APC", y.dim = "BUV661", asp,
  x.lab = "Foxp3 APC", y.lab = "CD19 BUV661", save = FALSE
)
fig7b.left <- apply.fig7.panel.style( fig7b.left )
fig7b.left <- apply.fig7.subtitle( fig7b.left, "Standard unmixing (OLS)", asp )

fig7b.right <- create.biplot(
  fig7b.unmixed.autospectral, x.dim = "APC", y.dim = "BUV661", asp,
  x.lab = "Foxp3 APC", y.lab = "CD19 BUV661", save = FALSE
)
fig7b.right <- apply.fig7.panel.style( fig7b.right )
fig7b.right <- apply.fig7.subtitle( fig7b.right, "Per-cell fluorophore (AutoSpectral)", asp )

panel.B <- plot_grid( fig7b.left, fig7b.right, nrow = 1 )

# ---------------------------------------------------------------------------
# Panel C: mouse lung, ungated, already-unmixed OLS vs AutoSpectral files
# ---------------------------------------------------------------------------

fig7c.large.gate <- FALSE

fig7c.ols.file <- "PLACEHOLDER_LUNG_OLS_UNMIXED.fcs"
fig7c.autospectral.file <- "PLACEHOLDER_LUNG_AUTOSPECTRAL_UNMIXED.fcs"

fig7c.ols.data <- readFCS( fig7c.ols.file )
fig7c.autospectral.data <- readFCS( fig7c.autospectral.file )

if ( fig7c.large.gate ) {
  fig7c.ols.data <- gate.large.sample(
    fig7c.ols.data, asp, samp = "Lung_OLS",
    output.dir = file.path( output.dir, "gate_definitions" )
  )
  fig7c.autospectral.data <- gate.large.sample(
    fig7c.autospectral.data, asp, samp = "Lung_AutoSpectral",
    output.dir = file.path( output.dir, "gate_definitions" )
  )
}

fig7c.channels.ols <- resolve.fluorophore.channels(
  fig7c.ols.data, c( "BUV615", "SBB615" ), fig7.fluorophore.database
)
fig7c.channels.autospectral <- resolve.fluorophore.channels(
  fig7c.autospectral.data, c( "BUV615", "SBB615" ), fig7.fluorophore.database
)

fig7c.left <- create.biplot(
  fig7c.ols.data,
  x.dim = fig7c.channels.ols[[ "BUV615" ]], y.dim = fig7c.channels.ols[[ "SBB615" ]], asp,
  x.lab = "Siglec F BUV615", y.lab = "Ly-6C SBB615", save = FALSE
)
fig7c.left <- apply.fig7.panel.style( fig7c.left )
fig7c.left <- apply.fig7.subtitle( fig7c.left, "OLS", asp )

fig7c.right <- create.biplot(
  fig7c.autospectral.data,
  x.dim = fig7c.channels.autospectral[[ "BUV615" ]], y.dim = fig7c.channels.autospectral[[ "SBB615" ]], asp,
  x.lab = "Siglec F BUV615", y.lab = "Ly-6C SBB615", save = FALSE
)
fig7c.right <- apply.fig7.panel.style( fig7c.right )
fig7c.right <- apply.fig7.subtitle( fig7c.right, "AutoSpectral", asp )

panel.C <- build.figure7.titled.pair( fig7c.left, fig7c.right, "Mouse Lung", block.letter = "C" )

# ---------------------------------------------------------------------------
# Panel D: OMIP-102, ungated, already-unmixed Chorus vs AutoSpectral files
# ---------------------------------------------------------------------------

fig7d.large.gate <- FALSE

fig7d.chorus.file <- "PLACEHOLDER_OMIP102_CHORUS_ORIGINAL_UNMIXED.fcs"
fig7d.autospectral.file <- "PLACEHOLDER_OMIP102_AUTOSPECTRAL_UNMIXED.fcs"

fig7d.chorus.data <- readFCS( fig7d.chorus.file )
fig7d.autospectral.data <- readFCS( fig7d.autospectral.file )

if ( fig7d.large.gate ) {
  fig7d.chorus.data <- gate.large.sample(
    fig7d.chorus.data, asp, samp = "OMIP102_Chorus",
    output.dir = file.path( output.dir, "gate_definitions" )
  )
  fig7d.autospectral.data <- gate.large.sample(
    fig7d.autospectral.data, asp, samp = "OMIP102_AutoSpectral",
    output.dir = file.path( output.dir, "gate_definitions" )
  )
}

fig7d.channels.chorus <- resolve.fluorophore.channels(
  fig7d.chorus.data, c( "PE", "PE-Fire 810" ), fig7.fluorophore.database
)
fig7d.channels.autospectral <- resolve.fluorophore.channels(
  fig7d.autospectral.data, c( "PE", "PE-Fire 810" ), fig7.fluorophore.database
)

fig7d.left <- create.biplot(
  fig7d.chorus.data,
  x.dim = fig7d.channels.chorus[[ "PE" ]], y.dim = fig7d.channels.chorus[[ "PE-Fire 810" ]], asp,
  x.lab = "CD57 PE", y.lab = "HLA-DR PE-Fire 810", save = FALSE
)
fig7d.left <- apply.fig7.panel.style( fig7d.left )
fig7d.left <- apply.fig7.subtitle( fig7d.left, "Chorus Original", asp )

fig7d.right <- create.biplot(
  fig7d.autospectral.data,
  x.dim = fig7d.channels.autospectral[[ "PE" ]], y.dim = fig7d.channels.autospectral[[ "PE-Fire 810" ]], asp,
  x.lab = "CD57 PE", y.lab = "HLA-DR PE-Fire 810", save = FALSE
)
fig7d.right <- apply.fig7.panel.style( fig7d.right )
fig7d.right <- apply.fig7.subtitle( fig7d.right, "AutoSpectral", asp )

panel.D <- build.figure7.titled.pair( fig7d.left, fig7d.right, "OMIP-102", block.letter = "D" )

# ---------------------------------------------------------------------------
# Panel E: 2x2 grid, mouse lung WT vs GFP mouse x OLS vs AutoSpectral,
# with power-of-ten gridlines
# ---------------------------------------------------------------------------

fig7e.large.gate <- FALSE

fig7e.files <- list(
  wt.ols            = "PLACEHOLDER_WT_LUNG_OLS_UNMIXED.fcs",
  wt.autospectral    = "PLACEHOLDER_WT_LUNG_AUTOSPECTRAL_UNMIXED.fcs",
  gfp.ols            = "PLACEHOLDER_GFPMOUSE_LUNG_OLS_UNMIXED.fcs",
  gfp.autospectral   = "PLACEHOLDER_GFPMOUSE_LUNG_AUTOSPECTRAL_UNMIXED.fcs"
)

fig7e.data <- lapply( fig7e.files, readFCS )

if ( fig7e.large.gate ) {
  fig7e.data <- Map(
    function( data, samp.name )
      gate.large.sample(
        data, asp, samp = samp.name,
        output.dir = file.path( output.dir, "gate_definitions" )
      ),
    fig7e.data,
    names( fig7e.data )
  )
}

build.fig7e.panel <- function( data ) {
  channels <- resolve.fluorophore.channels( data, c( "GFP", "BV510" ), fig7.fluorophore.database )
  biplot <- create.biplot(
    data, x.dim = channels[[ "GFP" ]], y.dim = channels[[ "BV510" ]], asp,
    x.lab = "iCas9 GFP", y.lab = "CD4 BV510", save = FALSE
  )
  biplot <- apply.fig7.panel.style( biplot )
  add.biplot.gridlines( biplot, asp )
}

fig7e.panels <- lapply( fig7e.data, build.fig7e.panel )

fig7e.col.titles <- plot_grid(
  ggdraw() + draw_label( "OLS", fontface = "bold" ),
  ggdraw() + draw_label( "AutoSpectral", fontface = "bold" ),
  nrow = 1
)

fig7e.row.wt <- plot_grid(
  fig7e.panels$wt.ols, fig7e.panels$wt.autospectral, nrow = 1
)
fig7e.row.gfp <- plot_grid(
  fig7e.panels$gfp.ols, fig7e.panels$gfp.autospectral, nrow = 1
)

fig7e.row.wt.labeled <- plot_grid(
  fig7e.row.wt, ggdraw() + draw_label( "WT mouse", angle = -90 ),
  ncol = 2, rel_widths = c( 1, 0.08 )
)
fig7e.row.gfp.labeled <- plot_grid(
  fig7e.row.gfp, ggdraw() + draw_label( "GFP mouse", angle = -90 ),
  ncol = 2, rel_widths = c( 1, 0.08 )
)

fig7e.title.row <- ggdraw() + draw_label( "Mouse Lung", fontface = "bold", size = 12 )

panel.E <- plot_grid(
  fig7e.title.row, fig7e.col.titles, fig7e.row.wt.labeled, fig7e.row.gfp.labeled,
  ncol = 1, rel_heights = c( 0.1, 0.1, 1, 1 )
)

# ---------------------------------------------------------------------------
# Composite
# ---------------------------------------------------------------------------
#
# Column 1: A, C, D, stacked (three equal rows).
# Column 2: B (top), E spanning the height of C and D (right of C/D, below
# B), since E is a single 2x2 grid rather than a separate panel per row.
# Panel letters for C and D are already set via build.figure7.titled.pair()
# above; A, B and E are labeled here, at composite assembly, since they
# don't go through that helper.

panel.A <- plot_grid( panel.A, labels = "A", label_size = fig7.panel.label.size, label_fontface = "bold" )
panel.B <- plot_grid( panel.B, labels = "B", label_size = fig7.panel.label.size, label_fontface = "bold" )
panel.E <- plot_grid( panel.E, labels = "E", label_size = fig7.panel.label.size, label_fontface = "bold" )

column.1 <- plot_grid( panel.A, panel.C, panel.D, ncol = 1, rel_heights = c( 1, 1, 1 ) )
column.2 <- plot_grid( panel.B, panel.E, ncol = 1, rel_heights = c( 1, 2 ) )

figure.7 <- plot_grid( column.1, column.2, ncol = 2, rel_widths = c( 1, 1.3 ) )

ggsave(
  file.path( output.dir, "Figure_7.jpg" ),
  plot = figure.7, width = 16, height = 16
)
