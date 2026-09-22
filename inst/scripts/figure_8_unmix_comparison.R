# figure_8_unmix_comparison.R
#
# Figure 8: cross-operator/method unmixing comparison, across three panels
# (7-colour, 42-colour, OMIP-95), each contributing four rows:
#
# A/E/I) Unstained noise (robust SD / MAD), one box per operator/method.
# B/F/J) Total spillover error (sum of |delta MFI| across every off-target
#         channel), one box per operator/method.
# C/G/K) Spillover (delta MFI) for one representative off-target channel of
#         one fluorophore.
# D/H/L) Representative single-stained-control biplots, one per
#         operator/method, with a larger shared axis label (arrow) on each
#         side.
#
# A, B, C, G, K are captured directly from the JPEG files already written by
# unmix.comparison.plot()/compare.unmix.folders() when the attached
# run_unmix_comparison_*.R scripts were run. E, F, I, J need the same
# unstained-MAD / summed-spillover-error plots, but with significance
# brackets (adjusted p < 0.05 vs AutoSpectral) added - this needs the
# ggpubr-based diff to .plot.metric.boxplot()/unmix.comparison.plot()
# (presented separately, not applied here); this script assumes that diff
# has been applied, regenerates those two plots per panel with it, and then
# captures the resulting JPEGs exactly like the other image panels. D, H, L
# are built natively as ggplot/cowplot objects via build.operator.biplot.row()
# (see plot_unmix_comparison_biplots.R), reusing the fluorophore/channel
# mapping already resolved by setup.unmix.comparison() for each folder - the
# per-folder setup CSVs already store one `match.fluorophores()` resolution
# per folder, so no further fluorophore matching is needed here (see the
# note at the top of plot_unmix_comparison_biplots.R).
#
# Panel letters (A-L) are drawn once per row/block via cowplot's own
# plot_grid(labels = ...): A-C/E-G/I-K on the three image panels of each
# row, D/H/L on the biplot row below. The A-C/E-G/I-K image panels are
# captured JPEGs from an earlier run, so their own internal axis/annotation
# text size is fixed at capture time and cannot be resized here; only the
# D/H/L biplot rows' subtitle and shared-axis-label text sizes are
# controllable from this script, via fig8.subtitle.text.size and
# fig8.axis.label.text.size below.

library( AutoSpectral )
library( cowplot )
library( ggplot2 )

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral/R"
devtools::load_all(asp.dir)
source("~/Bioinformatics/AutoSpectral/inst/scripts/plot_unmix_comparison_biplots.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/plot_biplot_annotations.R")

asp <- get.autospectral.param()
output.dir <- "./Figure 8"
if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

# Point size for the bold panel letters (A-L).
fig8.panel.label.size <- 18

# Point size for each D/H/L panel's operator-name subtitle. NULL uses
# asp$figure.axis.title.size.
fig8.subtitle.text.size <- NULL

# Point size for the larger shared axis-label strip on each D/H/L row.
fig8.axis.label.text.size <- 16

# Size (inches) of each square biplot panel in the D/H/L rows, in the final
# composite - the same role fig2.panel.size plays in
# figure_2_multipanel_plots.R. The composite's overall width/height are
# derived from it below, together with this figure's fixed structure: three
# blocks (7-Colour, 42-Colour, OMIP-95), each with a title strip and two
# content rows (a 3-wide row of captured JPEG panels, then a D/H/L row of
# fig8.n.biplots.per.row square biplots wrapped in a left-hand label strip
# by add.biplot.shared.axis.labels()). Unlike Figure 2's uniform panel
# grid, the captured JPEG panels (A/B/C etc.) are not forced square and so
# do not scale with this value the way the D/H/L biplots do - see the note
# at the composite assembly below.
fig8.panel.size <- 3.5

# Structural constants used only to derive the composite's canvas size
# below; change these only if the block/row layout itself changes.
fig8.n.blocks <- 3
fig8.rows.per.block <- 2
fig8.block.title.rel.height <- 0.08
fig8.n.biplots.per.row <- 6
fig8.biplot.label.strip.width <- 0.08

# ---------------------------------------------------------------------------
# Dataset root directories - each is the working directory one of the three
# attached run_unmix_comparison_*.R scripts was (or will be) run from, so
# that script's relative paths ("./unmix_comparison_setup", "comparison_log",
# "unmix_comparison_results.csv", etc.) resolve underneath it.
# ---------------------------------------------------------------------------

fig8.7c.dir     <- "../Unmixing comparison/7C analysis"
fig8.42c.dir    <- "../Unmixing comparison/SpecSymp analysis"
fig8.omip95.dir <- "../Unmixing comparison/OMIP95 analysis"

# ---------------------------------------------------------------------------
# Local layout helpers
# ---------------------------------------------------------------------------

embed.image.panel <- function( image.path ) {
  ggdraw() + draw_image( image.path )
}

build.figure8.block <- function( row.1, row.2, block.title ) {

  title.row <- ggdraw() + draw_label( block.title, fontface = "bold", size = 14 )

  plot_grid( title.row, row.1, row.2, ncol = 1, rel_heights = c( 0.08, 1, 1 ) )
}

# ---------------------------------------------------------------------------
# 7-colour panel (A-D)
# ---------------------------------------------------------------------------

fig8.7c.setup.dir <- file.path( fig8.7c.dir, "unmix_comparison_setup" )

fig8.7c.folders <- c(
  `FlowJo Unmixing Wizard` = file.path(
    fig8.7c.dir, "../7C comparison/FlowJo/Controls_renamed"
  ),
  Operator1    = file.path(
    fig8.7c.dir, "../7C comparison/Emily/7C_Panel_woAF_ES/Unmixed/Plate_001/Reference Group"
  ),
  Operator2    = file.path(
    fig8.7c.dir, "../7C comparison/Jay/7C_NO AF Extraction/Unmixed/Plate_001/Reference Group"
  ),
  Operator3    = file.path(
    fig8.7c.dir, "../7C comparison/Oliver/7C panel Oliver unmix no AF/Unmixed/Reference Group"
  ),
  Operator4    = file.path(
    fig8.7c.dir, "../7C comparison/Sameen/7C_Analysis_Unmixing wo AF/Unmixed/Reference Group"
  ),
  AutoSpectral = file.path(
    fig8.7c.dir, "../7C comparison/AutoSpectral/AutoSpectral_unmixed/SSC"
  )
)

panel.A <- embed.image.panel( file.path( fig8.7c.dir, "comparison_linear/Unstained_rSD.jpg" ) )
panel.B <- embed.image.panel( file.path( fig8.7c.dir, "comparison_linear/Summary_Delta.MFI.jpg" ) )
panel.C <- embed.image.panel(
  file.path( fig8.7c.dir, "comparison_linear/Delta.MFI_eFluor_780.jpg" )
)

panel.D <- build.operator.biplot.row(
  setup.dir        = fig8.7c.setup.dir,
  folders          = fig8.7c.folders,
  panel.order      = c(
    "FlowJo Unmixing Wizard", "Operator1", "Operator2", "Operator3", "Operator4", "AutoSpectral"
  ),
  on.target.fluor  = "eFluor 780",
  off.target.fluor = "PE-Cy7",
  asp              = asp,
  x.lab            = "PE-Cy7 (empty)",
  y.lab            = "Fixable viability dye eFluor 780",
  subtitle.text.size = fig8.subtitle.text.size,
  axis.label.text.size = fig8.axis.label.text.size,
  x.min = -20000, y.min = -20000
)

row.ABC <- plot_grid(
  panel.A, panel.B, panel.C, nrow = 1,
  labels = c( "A", "B", "C" ), label_size = fig8.panel.label.size, label_fontface = "bold"
)

panel.D.labeled <- plot_grid(
  panel.D, labels = "D", label_size = fig8.panel.label.size, label_fontface = "bold"
)

block.7C <- build.figure8.block( row.ABC, panel.D.labeled, "7-Colour Panel" )

# ---------------------------------------------------------------------------
# 42-colour panel (E-H)
# ---------------------------------------------------------------------------

fig8.42c.setup.dir <- file.path( fig8.42c.dir, "unmix_comparison_setup" )

fig8.42c.folders <- c(
  `FlowJo Unmixing Wizard` = file.path(
    fig8.42c.dir, "FlowJo Low Effort Unmix/Reference Group_renamed"
  ),
  Operator1           = file.path( fig8.42c.dir, "SpectroFlo_cells_poor/Unmixed/Reference Group" ),
  Operator2           = file.path( fig8.42c.dir, "SpectroFlo_cells_better/Unmixed/Reference Group" ),
  `Operator2 beads`     = file.path( fig8.42c.dir, "SpectroFlo_beads/Unmixed/SSC_cells" ),
  AutoSpectral        = file.path( fig8.42c.dir, "AutoSpectral/AutoSpectral_unmixed" ),
  `AutoSpectral beads`  = file.path( fig8.42c.dir, "AutoSpectral_beads/unmixed_cell_controls" )
)

# fluorophore-level significance vs AutoSpectral, log (ratio) scale - already
# written by the attached run_unmix_comparison_42C.R script
fig8.42c.stats.detail <- utils::read.csv(
  file.path( fig8.42c.dir, "unmix_comparison_stats_channels_detail_log.csv" ),
  stringsAsFactors = FALSE
)

fig8.42c.results <- utils::read.csv(
  file.path( fig8.42c.dir, "unmix_comparison_results.csv" ), stringsAsFactors = FALSE
)
fig8.42c.summary <- utils::read.csv(
  file.path( fig8.42c.dir, "unmix_comparison_summary.csv" ), stringsAsFactors = FALSE
)

# regenerate only the two plots panels E/F need, with significance brackets;
# requires the ggpubr diff to unmix.comparison.plot()/.plot.metric.boxplot()
unmix.comparison.plot(
  results          = fig8.42c.results,
  summary          = fig8.42c.summary,
  setup.files      = fig8.42c.folders,
  log.scale        = TRUE,
  point.size       = 2,
  plot.dir         = file.path( fig8.42c.dir, "comparison_log_stats" ),
  stats.detail     = fig8.42c.stats.detail,
  reference.folder = "AutoSpectral",
  stats.text.size = 4.5,
  stats.step.increase = 0.25,
  stats.label.style = "stars"
)

panel.E <- embed.image.panel(
  file.path( fig8.42c.dir, "comparison_log_stats/Unstained_rSD.jpg" )
)
panel.F <- embed.image.panel(
  file.path( fig8.42c.dir, "comparison_log_stats/Summary_Delta.MFI.jpg" )
)
panel.G <- embed.image.panel(
  file.path( fig8.42c.dir, "comparison_log/Delta.MFI_BV650.jpg" )
)

panel.H <- build.operator.biplot.row(
  setup.dir        = fig8.42c.setup.dir,
  folders          = fig8.42c.folders,
  panel.order      = c(
    "FlowJo Unmixing Wizard", "Operator1", "Operator2",
    "Operator2 beads", "AutoSpectral", "AutoSpectral beads"
  ),
  on.target.fluor  = "BV650",
  off.target.fluor = "BB660",
  asp              = asp,
  x.lab            = "BB660 (empty)",
  y.lab            = "XCR1 BV650",
  subtitle.text.size = fig8.subtitle.text.size,
  axis.label.text.size = fig8.axis.label.text.size,
  x.min = -10000,
  y.min = -10000,
  x.width.basis = -500,
  y.width.basis = -500
)

row.EFG <- plot_grid(
  panel.E, panel.F, panel.G, nrow = 1,
  labels = c( "E", "F", "G" ), label_size = fig8.panel.label.size, label_fontface = "bold"
)

panel.H.labeled <- plot_grid(
  panel.H, labels = "H", label_size = fig8.panel.label.size, label_fontface = "bold"
)

block.42C <- build.figure8.block( row.EFG, panel.H.labeled, "42-Colour Panel" )

# ---------------------------------------------------------------------------
# OMIP-95 (I-L)
# ---------------------------------------------------------------------------

fig8.omip95.setup.dir <- file.path( fig8.omip95.dir, "unmix_comparison_setup2" )

fig8.omip95.folders <- c(
  `FlowJo Unmixing Wizard` = file.path(
    fig8.omip95.dir, "../OMIP95 comparison/low effort FlowJo/SSC_renamed"
  ),
  Operator1 = file.path( fig8.omip95.dir, "../OMIP95 comparison/Jay/Reference Group" ),
  `Operator1 AF`  = file.path( fig8.omip95.dir, "../OMIP95 comparison/Jay AF/Reference Group" ),
  `Operator2 AF`  = file.path( fig8.omip95.dir, "../OMIP95 comparison/Emily AF" ),
  `Operator3`  = file.path(
    fig8.omip95.dir,
    "../OMIP95 comparison/Oliver/OMIP95 Oliver pregated controls/Unmixed/Reference Group"
  ),
  `Operator3 AF`  = file.path(
    fig8.omip95.dir,
    "../OMIP95 comparison/Oliver/OMIP95 Oliver pregated controls with AF/Unmixed/Reference Group"
  ),
  `Operator4`      = file.path( fig8.omip95.dir, "../OMIP95 comparison/Sameen/Reference Group" ),
  `Operator4 AF`      = file.path( fig8.omip95.dir, "../OMIP95 comparison/Sameen AF/Reference Group" ),
  AutoSpectral   = file.path(
    fig8.omip95.dir, "../OMIP95 comparison/AutoSpectral/AutoSpectral_unmixed"
  )
)

fig8.omip95.stats.detail <- utils::read.csv(
  file.path( fig8.omip95.dir, "unmix_comparison_stats_channels_detail_log.csv" ),
  stringsAsFactors = FALSE
)

fig8.omip95.results <- utils::read.csv(
  file.path( fig8.omip95.dir, "unmix_comparison_results.csv" ), stringsAsFactors = FALSE
)
fig8.omip95.summary <- utils::read.csv(
  file.path( fig8.omip95.dir, "unmix_comparison_summary.csv" ), stringsAsFactors = FALSE
)

# unmix_comparison_results.csv/summary.csv carry every folder unmixed
# across every OMIP-95 comparison run, including a plain "OperatorN" entry
# alongside "OperatorN AF" for operators 1, 3 and 4 (only the AF version
# exists for Operator2). unmix.comparison.plot() requires setup.files to
# name every folder present in results/summary, so results/summary are
# filtered down here to just the folders this panel plots
# (fig8.omip95.folders' own names), rather than widening fig8.omip95.folders
# to cover every extra folder those two files happen to carry.
fig8.omip95.plotted.folders <- names( fig8.omip95.folders )

fig8.omip95.results <- fig8.omip95.results[
  fig8.omip95.results$folder %in% fig8.omip95.plotted.folders, ]
fig8.omip95.summary <- fig8.omip95.summary[
  fig8.omip95.summary$folder %in% fig8.omip95.plotted.folders, ]

unmix.comparison.plot(
  results          = fig8.omip95.results,
  summary          = fig8.omip95.summary,
  setup.files      = fig8.omip95.folders,
  log.scale        = TRUE,
  point.size       = 2,
  plot.dir         = file.path( fig8.omip95.dir, "comparison_log_stats" ),
  stats.detail     = fig8.omip95.stats.detail,
  reference.folder = "AutoSpectral",
  stats.text.size = 4.5,
  stats.step.increase = 0.25,
  stats.label.style = "stars"
)

panel.I <- embed.image.panel(
  file.path( fig8.omip95.dir, "comparison_log_stats/Unstained_rSD.jpg" )
)
panel.J <- embed.image.panel(
  file.path( fig8.omip95.dir, "comparison_log_stats/Summary_Delta.MFI.jpg" )
)
panel.K <- embed.image.panel(
  file.path( fig8.omip95.dir, "comparison_log/Delta.MFI_APC-eFluor_780.jpg" )
)

panel.L <- build.operator.biplot.row(
  setup.dir        = fig8.omip95.setup.dir,
  folders          = fig8.omip95.folders,
  panel.order      = c(
    "FlowJo Unmixing Wizard",
    "Operator1",
    "Operator2 AF",
    "Operator3",
    "Operator4",
    "AutoSpectral"
  ),
  panel.labels     = c(
    "Operator3" = "Operator 3 (beads)"
  ),
  on.target.fluor  = "APC-eFluor 780",
  off.target.fluor = "APC",
  asp              = asp,
  x.lab            = "APC (empty)",
  y.lab            = "FcεRIα APC-eFluor 780",
  subtitle.text.size = fig8.subtitle.text.size,
  axis.label.text.size = fig8.axis.label.text.size,
  x.min = -40000,
  y.min = -40000,
  x.width.basis = -3000,
  y.width.basis = -3000
)

row.IJK <- plot_grid(
  panel.I, panel.J, panel.K, nrow = 1,
  labels = c( "I", "J", "K" ), label_size = fig8.panel.label.size, label_fontface = "bold"
)

panel.L.labeled <- plot_grid(
  panel.L, labels = "L", label_size = fig8.panel.label.size, label_fontface = "bold"
)

block.OMIP95 <- build.figure8.block( row.IJK, panel.L.labeled, "OMIP-95" )

# ---------------------------------------------------------------------------
# Composite
# ---------------------------------------------------------------------------

figure.8 <- plot_grid( block.7C, block.42C, block.OMIP95, ncol = 1 )

# Canvas derived from fig8.panel.size: width is set so each of the 6 D/H/L
# biplots renders at fig8.panel.size inches square (the label strip on the
# left of each of those rows takes fig8.biplot.label.strip.width of the
# row's width, so the full row - and so the figure - must be wider than
# fig8.panel.size * fig8.n.biplots.per.row by that fraction). Height is set
# so each block's two content rows (image row + biplot row) are each
# fig8.panel.size tall, plus its title strip at
# fig8.block.title.rel.height * fig8.panel.size, summed across
# fig8.n.blocks stacked blocks. The row 1 JPEG panels are not forced square
# and will fill whatever width/height this allocates them, same as before
# this derivation was added.
figure.8.width <- fig8.panel.size * fig8.n.biplots.per.row / ( 1 - fig8.biplot.label.strip.width )
figure.8.height <- fig8.n.blocks * fig8.panel.size * ( fig8.rows.per.block + fig8.block.title.rel.height )

ggsave(
  file.path( output.dir, "Figure_8.jpg" ),
  plot = figure.8,
  width = figure.8.width,
  height = figure.8.height,
  limitsize = FALSE
)
