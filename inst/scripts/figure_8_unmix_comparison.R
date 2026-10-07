# figure_8_unmix_comparison.R
#
# Figure 8: cross-operator/method unmixing comparison. One figure is written
# per flow cytometry panel (7-colour, 42-colour, OMIP-95), each with four
# panels:
#
# A) Unstained noise (robust SD / MAD), one box per operator/method.
# B) Total spillover error (sum of |delta MFI| across every off-target
#    channel), one box per operator/method.
# C) Spillover (delta MFI) for one representative off-target channel of
#    one fluorophore.
# D) Representative single-stained-control biplots, one per
#    operator/method, with a larger shared axis label (arrow) on each
#    side.
#
# A and B are captured directly from the JPEG files already written by
# unmix.comparison.plot()/compare.unmix.folders() when the attached
# run_unmix_comparison_*.R scripts were run; C likewise. For the 42-colour
# and OMIP-95 panels, A and B need significance brackets (adjusted p < 0.05
# vs AutoSpectral) - this needs the ggpubr-based diff to
# .plot.metric.boxplot()/unmix.comparison.plot() (presented separately, not
# applied here); this script assumes that diff has been applied, regenerates
# those two plots per panel with it, and then captures the resulting JPEGs
# exactly like the other image panels. D is built natively as ggplot/cowplot
# objects via build.operator.biplot.row() (see
# plot_unmix_comparison_biplots.R), reusing the fluorophore/channel mapping
# already resolved by setup.unmix.comparison() for each folder - the
# per-folder setup CSVs already store one `match.fluorophores()` resolution
# per folder, so no further fluorophore matching is needed here (see the
# note at the top of plot_unmix_comparison_biplots.R).
#
# Outputs, all under output.dir:
#   Figure_8_7C.jpg, Figure_8_42C.jpg, Figure_8_OMIP95.jpg
#       the three assembled figures.
#   7C_panels/, 42C_panels/, OMIP95_panels/
#       one subfolder per figure holding each individual image: every
#       operator/method biplot (for example 7C_FlowJo_biplot.jpg,
#       42C_Operator1_biplot.jpg) and a copy of each captured JPEG panel,
#       prefixed with the panel name (for example 7C_Unstained_rSD.jpg).
#
# Panel letters (A-D) are drawn once per row via cowplot's own
# plot_grid(labels = ...): A-C on the three image panels of the top row, D
# on the biplot row below. The A-C image panels are captured JPEGs from an
# earlier run, so their own internal axis/annotation text size is fixed at
# capture time and cannot be resized here; only the D biplot row's subtitle
# and shared-axis-label text sizes are controllable from this script, via
# fig8.subtitle.text.size and fig8.axis.label.text.size below.

library( AutoSpectral )
library( cowplot )
library( ggplot2 )

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
scripts.dir <- "inst/scripts"
devtools::load_all(asp.dir)
source(file.path(asp.dir, scripts.dir, "plot_unmix_comparison_biplots.R"))
source(file.path(asp.dir, scripts.dir, "plot_biplot_annotations.R"))
source(file.path(asp.dir, scripts.dir, "plot_unmix_comparison.R"))
source(file.path(asp.dir, scripts.dir, "test_unmix_comparison.R"))
source(file.path(asp.dir, scripts.dir, "compare_unmix_folders.R"))

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

# Size (inches) of each square biplot panel in the D rows, in the final
# figures - the same role fig2.panel.size plays in
# figure_2_multipanel_plots.R. Each figure's overall width/height are
# derived from it below, together with the fixed structure shared by all
# three figures: two content rows (a 3-wide row of captured JPEG panels,
# then a D row of fig8.n.biplots.per.row square biplots wrapped in a
# left-hand label strip by add.biplot.shared.axis.labels()). Unlike Figure
# 2's uniform panel grid, the captured JPEG panels (A/B/C) are not forced
# square and so do not scale with this value the way the D biplots do.
fig8.panel.size <- 3.5

# Size (inches) of each individually saved biplot image.
fig8.individual.biplot.size <- 4

# Structural constants used only to derive each figure's canvas size below;
# change these only if the row layout itself changes.
fig8.rows.per.figure <- 2
fig8.n.biplots.per.row <- 6
fig8.biplot.label.strip.width <- 0.08

# Canvas derived from fig8.panel.size: width is set so each of the 6 D
# biplots renders at fig8.panel.size inches square (the label strip on the
# left of that row takes fig8.biplot.label.strip.width of the row's width,
# so the full row - and so the figure - must be wider than
# fig8.panel.size * fig8.n.biplots.per.row by that fraction). Height is set
# so each of the figure's two content rows is fig8.panel.size tall. The top
# row JPEG panels are not forced square and will fill whatever width/height
# this allocates them.
figure.8.width <- fig8.panel.size * fig8.n.biplots.per.row / ( 1 - fig8.biplot.label.strip.width )
figure.8.height <- fig8.rows.per.figure * fig8.panel.size

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

# Copies each captured JPEG panel into the figure's panel subfolder, with the
# panel name prepended to the original filename (e.g. 7C_Unstained_rSD.jpg).
copy.figure8.images <- function( image.paths, panel.dir, file.prefix ) {

  if ( !dir.exists( panel.dir ) ) dir.create( panel.dir, recursive = TRUE )

  file.copy(
    image.paths,
    file.path( panel.dir, paste0( file.prefix, "_", basename( image.paths ) ) ),
    overwrite = TRUE
  )
}

# ---------------------------------------------------------------------------
# 7-colour panel
# ---------------------------------------------------------------------------

fig8.7c.setup.dir <- file.path( fig8.7c.dir, "unmix_comparison_setup" )
fig8.7c.panel.dir <- file.path( output.dir, "7C_panels" )

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

panel.7C.A <- embed.image.panel( file.path( fig8.7c.dir, "comparison_linear/Unstained_rSD.jpg" ) )
panel.7C.B <- embed.image.panel( file.path( fig8.7c.dir, "comparison_linear/Summary_Delta.MFI.jpg" ) )
panel.7C.C <- embed.image.panel(
  file.path( fig8.7c.dir, "comparison_linear/Delta.MFI_eFluor_780.jpg" )
)

panel.7C.D <- build.operator.biplot.row(
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
  x.min = -20000, y.min = -20000,
  save.dir = fig8.7c.panel.dir,
  file.prefix = "7C",
  file.labels = c( "FlowJo Unmixing Wizard" = "FlowJo" ),
  save.size = fig8.individual.biplot.size
)

row.7C <- plot_grid(
  panel.7C.A, panel.7C.B, panel.7C.C, nrow = 1,
  labels = c( "A", "B", "C" ), label_size = fig8.panel.label.size, label_fontface = "bold"
)

panel.7C.D.labeled <- plot_grid(
  panel.7C.D, labels = "D", label_size = fig8.panel.label.size, label_fontface = "bold"
)

figure.8.7C <- plot_grid( row.7C, panel.7C.D.labeled, ncol = 1 )

copy.figure8.images(
  c(
    file.path( fig8.7c.dir, "comparison_linear/Unstained_rSD.jpg" ),
    file.path( fig8.7c.dir, "comparison_linear/Summary_Delta.MFI.jpg" ),
    file.path( fig8.7c.dir, "comparison_linear/Delta.MFI_eFluor_780.jpg" )
  ),
  fig8.7c.panel.dir, "7C"
)

ggsave(
  file.path( output.dir, "Figure_8_7C.jpg" ),
  plot = figure.8.7C,
  width = figure.8.width,
  height = figure.8.height,
  limitsize = FALSE
)

# ---------------------------------------------------------------------------
# 42-colour panel
# ---------------------------------------------------------------------------

fig8.42c.setup.dir <- file.path( fig8.42c.dir, "unmix_comparison_setup" )
fig8.42c.panel.dir <- file.path( output.dir, "42C_panels" )

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

# regenerate only the two plots panels A and B need, with significance brackets;
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

panel.42C.A <- embed.image.panel(
  file.path( fig8.42c.dir, "comparison_log_stats/Unstained_rSD.jpg" )
)
panel.42C.B <- embed.image.panel(
  file.path( fig8.42c.dir, "comparison_log_stats/Summary_Delta.MFI.jpg" )
)
panel.42C.C <- embed.image.panel(
  file.path( fig8.42c.dir, "comparison_log/Delta.MFI_BV650.jpg" )
)

panel.42C.D <- build.operator.biplot.row(
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
  y.width.basis = -500,
  save.dir = fig8.42c.panel.dir,
  file.prefix = "42C",
  file.labels = c( "FlowJo Unmixing Wizard" = "FlowJo" ),
  save.size = fig8.individual.biplot.size
)

row.42C <- plot_grid(
  panel.42C.A, panel.42C.B, panel.42C.C, nrow = 1,
  labels = c( "A", "B", "C" ), label_size = fig8.panel.label.size, label_fontface = "bold"
)

panel.42C.D.labeled <- plot_grid(
  panel.42C.D, labels = "D", label_size = fig8.panel.label.size, label_fontface = "bold"
)

figure.8.42C <- plot_grid( row.42C, panel.42C.D.labeled, ncol = 1 )

copy.figure8.images(
  c(
    file.path( fig8.42c.dir, "comparison_log_stats/Unstained_rSD.jpg" ),
    file.path( fig8.42c.dir, "comparison_log_stats/Summary_Delta.MFI.jpg" ),
    file.path( fig8.42c.dir, "comparison_log/Delta.MFI_BV650.jpg" )
  ),
  fig8.42c.panel.dir, "42C"
)

ggsave(
  file.path( output.dir, "Figure_8_42C.jpg" ),
  plot = figure.8.42C,
  width = figure.8.width,
  height = figure.8.height,
  limitsize = FALSE
)

# ---------------------------------------------------------------------------
# OMIP-95
# ---------------------------------------------------------------------------

fig8.omip95.setup.dir <- file.path( fig8.omip95.dir, "unmix_comparison_setup2" )
fig8.omip95.panel.dir <- file.path( output.dir, "OMIP95_panels" )

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
  file.path( fig8.omip95.dir, "unmix_comparison_results",
             "unmix_comparison_stats_channels_detail_log.csv" ),
  stringsAsFactors = FALSE
)

fig8.omip95.results <- utils::read.csv(
  file.path( fig8.omip95.dir,  "unmix_comparison_results",
             "unmix_comparison_results.csv" ), stringsAsFactors = FALSE
)
fig8.omip95.summary <- utils::read.csv(
  file.path( fig8.omip95.dir,  "unmix_comparison_results",
             "unmix_comparison_summary.csv" ), stringsAsFactors = FALSE
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

panel.OMIP95.A <- embed.image.panel(
  file.path( fig8.omip95.dir, "comparison_log_stats/Unstained_rSD.jpg" )
)
panel.OMIP95.B <- embed.image.panel(
  file.path( fig8.omip95.dir, "comparison_log_stats/Summary_Delta.MFI.jpg" )
)
panel.OMIP95.C <- embed.image.panel(
  file.path( fig8.omip95.dir, "comparison_log/Delta.MFI_APC-eFluor_780.jpg" )
)

panel.OMIP95.D <- build.operator.biplot.row(
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
  y.width.basis = -3000,
  save.dir = fig8.omip95.panel.dir,
  file.prefix = "OMIP95",
  file.labels = c( "FlowJo Unmixing Wizard" = "FlowJo" ),
  save.size = fig8.individual.biplot.size
)

row.OMIP95 <- plot_grid(
  panel.OMIP95.A, panel.OMIP95.B, panel.OMIP95.C, nrow = 1,
  labels = c( "A", "B", "C" ), label_size = fig8.panel.label.size, label_fontface = "bold"
)

panel.OMIP95.D.labeled <- plot_grid(
  panel.OMIP95.D, labels = "D", label_size = fig8.panel.label.size, label_fontface = "bold"
)

figure.8.OMIP95 <- plot_grid( row.OMIP95, panel.OMIP95.D.labeled, ncol = 1 )

copy.figure8.images(
  c(
    file.path( fig8.omip95.dir, "comparison_log_stats/Unstained_rSD.jpg" ),
    file.path( fig8.omip95.dir, "comparison_log_stats/Summary_Delta.MFI.jpg" ),
    file.path( fig8.omip95.dir, "comparison_log/Delta.MFI_APC-eFluor_780.jpg" )
  ),
  fig8.omip95.panel.dir, "OMIP95"
)

ggsave(
  file.path( output.dir, "Figure_8_OMIP95.jpg" ),
  plot = figure.8.OMIP95,
  width = figure.8.width,
  height = figure.8.height,
  limitsize = FALSE
)
