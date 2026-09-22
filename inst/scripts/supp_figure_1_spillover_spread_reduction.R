# supp_figure_1_spillover_spread_reduction.R
#
# Supplementary Figure 1: spillover spread reduction from per-cell
# fluorophore spectrum optimization.
#
# A) Synthetic ground-truth mixtures of BUV661+ and APC+ single-stained
#    control events (and their sum, as a synthetic double-positive
#    population) are unmixed two ways - OLS and per-cell fluorophore
#    optimization - and shown as a 2 (method) x 3 (ground-truth population)
#    grid of biplots, with positivity threshold lines and the percentage of
#    events in each quadrant annotated.
# B) The per-cell cosine similarity of the BUV661 single-stained control's
#    raw signal to its assigned spectrum is compared between OLS (a single
#    fixed reference spectrum) and per-cell fluorophore optimization (the
#    best-fitting spectral variant per event), as overlaid density
#    distributions with a Kolmogorov-Smirnov test p-value.
#
# synthetic.ols and synthetic.percell (panel A) are both unmixed here, in
# this script, against the same `spectra` object, so they are guaranteed to
# share column names and no resolve.fluorophore.channels() step is needed -
# unlike the Figure 2/6/7/8 comparisons, which read two independently,
# externally unmixed exports.
#
# Every file path below is a PLACEHOLDER for you to fill in.

library( AutoSpectral )
library( cowplot )
library( ggplot2 )

asp <- get.autospectral.param()
output.dir <- "./supp_figure_1"
if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

# Point size for the bold panel letters (A, B).
supp1.panel.label.size <- 16

# Point size for panel B's KS-test p-value annotation, passed to
# plot.similarity.comparison(). NULL uses that function's own default (3.5).
supp1.annotation.text.size <- NULL

# ---------------------------------------------------------------------------
# File paths (PLACEHOLDERS)
# ---------------------------------------------------------------------------

spectra.file        <- "PLACEHOLDER_SPECTRA_FILE.csv"
af.spectra.file      <- "PLACEHOLDER_AF_SPECTRA_FILE.csv"
control.dir          <- "PLACEHOLDER_CONTROL_DIR"
control.def.file     <- "PLACEHOLDER_CONTROL_DEF_FILE.csv"
buv661.fcs.file       <- "PLACEHOLDER_BUV661_CONTROL.fcs"
apc.fcs.file           <- "PLACEHOLDER_APC_CONTROL.fcs"
unstained.fcs.file     <- "PLACEHOLDER_UNSTAINED_CONTROL.fcs"

# optional cache: set to a previously-saved get.spectral.variants() .rds
# file to skip recomputing spectral variants; leave NULL to compute fresh
variants.rds.path <- NULL

# ---------------------------------------------------------------------------
# Spectra, autofluorescence, and spectral variants
# ---------------------------------------------------------------------------

spectra <- read.spectra( spectra.file )
spectral.channel <- colnames( spectra )
scatter.channel <- read.scatter.parameter( asp )

af.spectra <- read.spectra( af.spectra.file, check.collinearity = FALSE )

spectra.variants <- if ( !is.null( variants.rds.path ) && file.exists( variants.rds.path ) ) {
  readRDS( variants.rds.path )
} else {
  get.spectral.variants(
    control.dir = control.dir,
    control.def.file = control.def.file,
    asp = asp,
    spectra = spectra,
    figures = FALSE
  )
}

# ---------------------------------------------------------------------------
# Read and gate single-stained controls
# ---------------------------------------------------------------------------

buv661.data <- readFCS( buv661.fcs.file, columns = union( spectral.channel, scatter.channel ) )
apc.data <- readFCS( apc.fcs.file, columns = union( spectral.channel, scatter.channel ) )
unstained.data <- readFCS( unstained.fcs.file, columns = union( spectral.channel, scatter.channel ) )

gate.boundary <- compute.scatter.gate(
  unstained.data, asp, scatter.param = scatter.channel,
  large.gate = FALSE, samp = "lymphocytes",
  output.dir = file.path( output.dir, "gates" )
)

gated.unstained <- apply.gate( unstained.data, gate.boundary, scatter.param = scatter.channel, asp = asp )
gated.buv661 <- apply.gate( buv661.data, gate.boundary, scatter.param = scatter.channel, asp = asp )
gated.apc <- apply.gate( apc.data, gate.boundary, scatter.param = scatter.channel, asp = asp )

spect.unstained <- gated.unstained[ , spectral.channel ]
spect.buv661 <- gated.buv661[ , spectral.channel ]
spect.apc <- gated.apc[ , spectral.channel ]

# ---------------------------------------------------------------------------
# Positivity thresholds and single-color positive event selection
# ---------------------------------------------------------------------------

unmixed.unstained.ols <- unmix.ols( spect.unstained, spectra )
unmixed.buv661.ols <- unmix.ols( spect.buv661, spectra )
unmixed.apc.ols <- unmix.ols( spect.apc, spectra )

apc.pos.thresh <- stats::quantile( unmixed.unstained.ols[ , "APC" ], 0.995 ) * 2
buv661.pos.thresh <- stats::quantile( unmixed.unstained.ols[ , "BUV661" ], 0.995 ) * 2

apc.pos.idx <- which( unmixed.apc.ols[ , "APC" ] > apc.pos.thresh )
buv661.pos.idx <- which( unmixed.buv661.ols[ , "BUV661" ] > buv661.pos.thresh )

n.synthetic.cells <- min( 10000, length( apc.pos.idx ), length( buv661.pos.idx ) )

buv661.part.raw <- spect.buv661[ buv661.pos.idx[ seq_len( n.synthetic.cells ) ], ]
apc.part.raw <- spect.apc[ apc.pos.idx[ seq_len( n.synthetic.cells ) ], ]

ground.truth.buv661.idx <- seq_len( n.synthetic.cells )
ground.truth.apc.idx <- n.synthetic.cells + seq_len( n.synthetic.cells )
ground.truth.doublepos.idx <- 2 * n.synthetic.cells + seq_len( n.synthetic.cells )

synthetic.data <- rbind(
  buv661.part.raw,
  apc.part.raw,
  buv661.part.raw + apc.part.raw
)

# ---------------------------------------------------------------------------
# Unmix synthetic mixtures: OLS vs per-cell fluorophore optimization
# ---------------------------------------------------------------------------

synthetic.ols <- unmix.ols( synthetic.data, spectra )

synthetic.percell <- unmix.autospectral.joint(
  raw.data = synthetic.data,
  spectra = spectra,
  af.spectra = af.spectra,
  asp = asp,
  spectra.variants = spectra.variants,
  verbose = FALSE
)

# ---------------------------------------------------------------------------
# Panel A: 2 (method) x 3 (ground-truth population) biplot grid
# ---------------------------------------------------------------------------

populations <- list(
  list( key = "buv661", idx = ground.truth.buv661.idx, title = "Ground truth:\nBUV661+APC-" ),
  list( key = "apc", idx = ground.truth.apc.idx, title = "Ground truth:\nBUV661-APC+" ),
  list( key = "double", idx = ground.truth.doublepos.idx, title = "Ground truth:\nBUV661+APC+" )
)

methods <- list(
  list( key = "ols", label = "OLS", data = synthetic.ols ),
  list( key = "percell", label = "Per-cell fluorophore\noptimization", data = synthetic.percell )
)

build.supp1a.panel <- function( plot.data, panel.title ) {

  biplot <- create.biplot(
    plot.data, x.dim = "APC", y.dim = "BUV661", asp,
    x.min = -2e4, save = FALSE, title = panel.title
  )
  biplot <- biplot + theme( aspect.ratio = 1, plot.margin = margin( 8, 8, 8, 8 ) )
  biplot$coordinates$clip <- "off"

  annotate.biplot.quadrants(
    plot = biplot,
    data = plot.data,
    x.dim = "APC", y.dim = "BUV661",
    x.thresh = apc.pos.thresh, y.thresh = buv661.pos.thresh,
    asp = asp
  )
}

panel.grid <- list()
for ( method in methods ) {
  for ( population in populations ) {
    panel.grid[[ paste( method$key, population$key, sep = "_" ) ]] <- build.supp1a.panel(
      method$data[ population$idx, ], NULL
    )
  }
}

col.title.row <- plot_grid(
  plotlist = lapply( populations, function( population )
    ggdraw() + draw_label( population$title, fontface = "bold", size = 11 ) ),
  nrow = 1
)

build.supp1a.row <- function( method ) {

  row.plots <- plot_grid(
    plotlist = lapply( populations, function( population )
      panel.grid[[ paste( method$key, population$key, sep = "_" ) ]] ),
    nrow = 1
  )

  row.label <- ggdraw() + draw_label( method$label, angle = 90, size = 11 )

  plot_grid( row.label, row.plots, ncol = 2, rel_widths = c( 0.08, 1 ) )
}

panel.A <- plot_grid(
  col.title.row,
  build.supp1a.row( methods[[ 1 ]] ),
  build.supp1a.row( methods[[ 2 ]] ),
  ncol = 1,
  rel_heights = c( 0.15, 1, 1 )
)

panel.A <- plot_grid(
  panel.A, labels = "A", label_size = supp1.panel.label.size, label_fontface = "bold"
)

ggsave(
  file.path( output.dir, "Supp_Figure_1A_SpilloverSpread.jpg" ),
  plot = panel.A, width = 9, height = 6.5
)

# ---------------------------------------------------------------------------
# Panel B: per-cell cosine similarity, OLS vs per-cell fluorophore
# optimization
# ---------------------------------------------------------------------------

ols.similarity <- per.cell.cosine.similarity( buv661.part.raw, spectra[ "BUV661", ] )

percell.similarity <- compute.percell.variant.cosine.similarity(
  raw.data = buv661.part.raw,
  spectra = spectra,
  af.spectra = af.spectra,
  spectra.variants = spectra.variants,
  target.fluor = "BUV661",
  asp = asp
)

panel.B <- plot.similarity.comparison(
  values.1 = ols.similarity,
  values.2 = percell.similarity,
  label.1 = "OLS",
  label.2 = "PerCell",
  color.1 = "lightblue",
  color.2 = "gold",
  title = "BUV661 Cosine Similarity",
  text.size = supp1.annotation.text.size %||% 3.5
)

panel.B <- plot_grid(
  panel.B, labels = "B", label_size = supp1.panel.label.size, label_fontface = "bold"
)

ggsave(
  file.path( output.dir, "Supp_Figure_1B_CosineSimilarity.jpg" ),
  plot = panel.B, width = 6, height = 3
)

# ---------------------------------------------------------------------------
# Composite
# ---------------------------------------------------------------------------

supp.figure.1 <- plot_grid( panel.A, panel.B, ncol = 1, rel_heights = c( 2, 1 ) )

ggsave(
  file.path( output.dir, "Supp_Figure_1_SpilloverSpreadReduction.jpg" ),
  plot = supp.figure.1, width = 9, height = 9.5
)
