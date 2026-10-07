# =============================================================================
# Figure_1.R  --  Figure 1: Describe the problem
# =============================================================================
#
# AIM
#   Show that compensation/unmixing beads and cells give different spectra for
#   the same fluorophore, that this happens on several cytometers, and that
#   the resulting unmixing error grows with the brightness of the stain.
#
# PANELS (one image file per letter, plus the stitched figure)
#   a  Spectral overlay (Cells vs BD Beads) and difference, PerCP-Cy5.5,
#      FACSDiscover A8, with cosine similarity; detector labels coloured by laser.
#      Spectra L2-normalised (all detectors) before the difference is taken.
#   b  The same for Aurora, A8, ID7000 and Xenith (2 x 2 grid, one shared laser
#      key). ID7000 deep UV detectors are omitted from the display (not from the
#      normalisation or the cosine similarity).
#   c  Spider chart: cell -> bead shift in "laser space" (each laser's detector
#      bank = one axis; centred log2 ratio of bank L2 norms), one unfilled trace
#      per cytometer; directional concordance (mean pairwise cosine of the shift
#      vectors) with a laser-label permutation test.
#   d  Biplot of single-stained cells (titration series, one colour per
#      antibody amount, negative = black) unmixed in a 35-colour space with the
#      CORRECT cell-derived spectrum (from the brightest, 20 ul, cell file).
#   e  Same events unmixed with the INCORRECT bead-derived spectrum (from the
#      20 ul bead file). x = stained fluorophore, y = channel with the largest
#      |SSI| (secondary stain index) under the bead-derived unmixing. Red =
#      spillover-spread boundaries (Residual Model, brightest cell file);
#      secondary SSI for every population is given in the legend.
#
# STRUCTURE OF THE ANALYSIS
#   a, b: read the cell and BD-bead PerCP-Cy5.5 spectra of each cytometer from
#         the 20260901 run of run.bead.cell.comparison() (legacy pipeline,
#         density gating), L2-normalise, compare them (difference, cosine).
#   c   : same spectra, restricted to the lasers present on every cytometer;
#         per-laser bank norms -> centred log2 ratio (beads vs cells) per
#         cytometer (laser_space_shift.R); angles: total, between-laser,
#         within-laser; concordance permutation test across cytometers.
#   d, e: 1. find the titration files (exact "Cells <fluor>" / "Beads <fluor>"
#            folders; "unlabelled" folders are NOT used);
#         2. density FSC/SSC gate from the pooled cell series (+ unstained);
#            separate density gate for the bead series;
#         3. cell and bead reference spectra from the 20 ul files with the
#            legacy AutoSpectral pipeline (define.flow.control ->
#            clean.controls -> get.fluorophore.spectra);
#         4. 35-colour panel: library spectra, with the stained fluorophore
#            replaced by the cell (correct) or bead (incorrect) spectrum;
#         5. unmix every gated file (OLS); keep events above
#            positive.multiplier x the 99.5th percentile of the unmixed
#            unstained control in the target channel (internal negatives are
#            present in the stained files);
#         6. y channel = largest |SSI| in the brightest file under the
#            bead-derived unmixing; SSI per population from calculate.ssi();
#         7. spillover-spread boundaries from the Residual Model replicate
#            (calculate.residual.spread), brightest cell file, correct matrix.
#
# DATA SOURCES (edit below)
#   a-c  : <proj.dir>/20260130 beads cells/results_density_20260901/
#            <aurora|a8|id7000|xenith>/<Cells|BD>/table_spectra/*_autospectral_spectra.csv
#   d, e : <proj.dir>/20250709 Bead Cells with Laura/
#            Aurora Beads and cells072025/Raw/Plate_001 (PerCP-Cy5.5)
#                                         /Plate_002 (BV711)
#          panel spectra: Aurora cell spectra of the 20260901 run (35-colour list)
#   laser assignment: AutoSpectral inst/extdata/cytometer_database.csv
#
# OUTPUT
#   <output.dir>/Fig1a_*.png ... Fig1e_*.png (300 dpi), Figure_1.png,
#   <output.dir>/tables/*.csv, <output.dir>/cache, <output.dir>/work
# =============================================================================

# ---- user settings ----------------------------------------------------------

# Package: devtools::load_all() on the dev copy (TRUE) or library(AutoSpectral)
use.load.all <- TRUE
asp.dir      <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"

proj.dir    <- "/Volumes/T7 Shield/Beads vs Cells"
script.dir  <- file.path( proj.dir, "Figures" )          # holds figure_functions/
output.dir  <- file.path( proj.dir, "Figures", "Figure_1" )

results.dir <- file.path( proj.dir, "20260130 beads cells", "results_density_20260901" )
raw.dir.1c  <- file.path( proj.dir, "20250709 Bead Cells with Laura",
                          "Aurora Beads and cells072025", "Raw" )

# Fluorophore for panels d and e: "PerCP-Cy5.5" or "BV711" (the whole of d/e
# is built for the one chosen). Panels a-c are always PerCP-Cy5.5.
fig1de.fluor <- "PerCP-Cy5.5"

# Panel sizes (width, height) in `size.units`
size.units  <- "in"
panel.size  <- list(
  a = c( width = 6.0,  height = 5.0 ),
  b = c( width = 12.0, height = 10.0 ),
  c = c( width = 4.5,  height = 4.5 ),
  d = c( width = 4.5,  height = 4.5 ),
  e = c( width = 4.5,  height = 4.5 )
)
figure.size <- c( width = 14, height = 15 )               # stitched figure
dpi         <- 300
base.size   <- 9                                          # base font size

# Panel d/e options
positive.quantile   <- 0.995     # plot events above this quantile of unmixed unstained ...
positive.multiplier <- 2         # ... multiplied by this factor (e.g. 2 = twice the quantile)
positive.space    <- "cells"     # "cells": same events in c and d; "each": own threshold
viridis.option    <- "viridis"
viridis.range     <- c( 0.15, 0.95 )
reuse.cache       <- TRUE        # reuse extracted reference spectra if present

# Panel a/b options
fig1ab.fluor.regex <- "^PerCP-Cy5\\.5$"   # row name in the *_autospectral_spectra.csv
panel.a.normalization <- "L2"    # "L2", "max" (peak = 1) or "none"
panel.b.normalization <- "L2"
laser.key.height      <- 0.04    # relative height of the shared laser key under panel b
cytometers <- list(
  Aurora = list( results = "aurora", db.col = "Aurora",  exclude.lasers = NULL ),
  A8     = list( results = "a8",     db.col = "Discover", exclude.lasers = NULL ),
  ID7000 = list( results = "id7000", db.col = "ID7000",  exclude.lasers = "DeepUV" ),
  Xenith = list( results = "xenith", db.col = "Xenith",  exclude.lasers = NULL )
)
panel.a.cytometer <- "A8"
panel.b.order     <- c( "Aurora", "A8", "ID7000", "Xenith" )
panel.b.ncol      <- 2

# Panel c options
radar.instruments <- c( "Aurora", "A8", "ID7000", "Xenith" )
radar.lasers      <- NULL        # NULL = lasers present on every instrument in radar.instruments
radar.metric      <- "clr"       # "clr" (centred log2 ratio of bank norms) or "delta"
radar.n.perm      <- 10000       # permutations for the concordance test
radar.seed        <- 42
instrument.colors <- c( Aurora = "#E69F00", A8 = "#56B4E9", ID7000 = "#009E73",
                        Xenith = "#CC79A7" )   # Okabe-Ito

# 35-colour panel (Aurora names)
fluors.35C <- c(
  "BUV395", "BUV496", "BUV563", "BUV615", "BUV661", "BUV737", "BUV805",
  "BV421", "BV480", "BV510", "BV570", "BV605", "BV650", "BV711", "BV750", "BV786",
  "BB515", "RB545", "BB630", "RB670", "BB700", "PerCP-Cy5.5", "RB705", "RB744", "RB780",
  "PE", "RY586", "PE-CF594", "PE-Cy5", "RY743", "PE-Cy7",
  "APC", "RR688", "R718", "APC-H7"
)
panel.spectra.file <- file.path( results.dir, "aurora", "Cells", "table_spectra",
                                 "Cells_autospectral_spectra.csv" )

# Files and folder names of the 2025-07-09 Aurora titrations
fluor.info.1c <- list(
  "PerCP-Cy5.5" = list( plate = "Plate_001", cell.dir = "Cells PerCP-Cy5.5",
                        bead.dir = "Beads PerCPCy5.5" ),
  "BV711"       = list( plate = "Plate_002", cell.dir = "Cells BV711",
                        bead.dir = "Beads BV711" )
)

# ---- packages and functions -------------------------------------------------

required.packages <- c( "ggplot2", "cowplot", "ragg", "viridis", "scattermore", "FNN" )
invisible( lapply( required.packages, library, character.only = TRUE ) )

if ( use.load.all ) {
  devtools::load_all( asp.dir )
  if ( dir.exists( asp.rcpp.dir ) ) devtools::load_all( asp.rcpp.dir )
  cytometer.db <- file.path( asp.dir, "inst", "extdata", "cytometer_database.csv" )
} else {
  library( AutoSpectral )
  if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) ) library( AutoSpectralRcpp )
  cytometer.db <- system.file( "extdata", "cytometer_database.csv", package = "AutoSpectral" )
}

func.dir <- file.path( script.dir, "figure_functions" )
if ( !dir.exists( func.dir ) ) stop( "figure_functions folder not found: ", func.dir )
invisible( lapply( list.files( func.dir, pattern = "\\.R$", full.names = TRUE ), source ) )

for ( d in c( output.dir, file.path( output.dir, c( "tables", "cache", "work" ) ) ) )
  if ( !dir.exists( d ) ) dir.create( d, recursive = TRUE )

save.fig.panel <- function( plot, name, size )
  save.panel( plot, name, output.dir, size[[ "width" ]], size[[ "height" ]],
              units = size.units, dpi = dpi )

# =============================================================================
# Panels a and b: spectral mismatch, PerCP-Cy5.5, BD beads vs cells
# =============================================================================

cyto.data <- lapply( names( cytometers ), function( cy ) {
  cfg <- cytometers[[ cy ]]
  list(
    spectra = load.cell.bead.spectra( results.dir, cfg$results, fig1ab.fluor.regex ),
    lasers  = get.detector.lasers( cytometer.db, cfg$db.col ),
    exclude = cfg$exclude.lasers
  )
} )
names( cyto.data ) <- names( cytometers )

make.mismatch.plot <- function( cy, normalization, show.laser.legend ) {
  d <- cyto.data[[ cy ]]
  plot.spectral.mismatch(
    ref.spectrum      = d$spectra$cells,
    test.spectrum     = d$spectra$beads,
    detector.lasers   = d$lasers,
    ref.label         = "Cells",
    test.label        = "BD Beads",
    normalization     = normalization,
    exclude.lasers    = d$exclude,
    title             = paste( cy, "PerCP-Cy5.5" ),
    base.size         = base.size,
    show.laser.legend = show.laser.legend
  )
}

# panel a: single instrument, own laser key
p1a <- make.mismatch.plot( panel.a.cytometer, panel.a.normalization, TRUE )
save.fig.panel( p1a, paste0( "Fig1a_", panel.a.cytometer, "_PerCP-Cy5.5" ), panel.size$a )

# panel b: all instruments, one shared single-row laser key
mismatch.plots <- lapply( panel.b.order, make.mismatch.plot,
                          normalization = panel.b.normalization, show.laser.legend = FALSE )
names( mismatch.plots ) <- panel.b.order

cosines <- vapply( mismatch.plots, function( p ) attr( p, "cosine.similarity" ), numeric( 1 ) )
utils::write.csv(
  data.frame( Cytometer = names( cosines ), Cosine.similarity = cosines ),
  file.path( output.dir, "tables", "Fig1ab_cosine_similarity.csv" ), row.names = FALSE
)

shown.lasers <- unique( unlist( lapply( mismatch.plots, attr, "lasers" ) ) )
shown.lasers <- c( intersect( laser.order, shown.lasers ), setdiff( shown.lasers, laser.order ) )

p1b <- cowplot::plot_grid(
  cowplot::plot_grid( plotlist = mismatch.plots, ncol = panel.b.ncol, align = "none" ),
  make.laser.key( shown.lasers, base.size = base.size ),
  ncol = 1, rel_heights = c( 1, laser.key.height )
)
save.fig.panel( p1b, "Fig1b_all_cytometers_PerCP-Cy5.5", panel.size$b )

# =============================================================================
# Panel c: cell -> bead shift in laser space, all instruments
# =============================================================================

if ( is.null( radar.lasers ) ) {
  radar.lasers <- Reduce( intersect, lapply( radar.instruments, function( cy )
    lasers.present( cyto.data[[ cy ]]$spectra$cells, cyto.data[[ cy ]]$lasers ) ) )
  radar.lasers <- c( intersect( laser.order, radar.lasers ), setdiff( radar.lasers, laser.order ) )
}
message( "Laser-space axes: ", paste( radar.lasers, collapse = ", " ) )

laser.shifts <- lapply( radar.instruments, function( cy ) {
  d <- cyto.data[[ cy ]]
  laser.space.shift( d$spectra$cells, d$spectra$beads, d$lasers, radar.lasers,
                     metric = radar.metric )
} )
names( laser.shifts ) <- radar.instruments

shift.matrix <- do.call( rbind, lapply( laser.shifts, `[[`, "shift" ) )
rownames( shift.matrix ) <- radar.instruments

concordance <- laser.shift.concordance( shift.matrix, n.perm = radar.n.perm, seed = radar.seed )
message( sprintf( "Laser-space concordance: mean pairwise cosine = %.3f, permutation p = %s (%d permutations)",
                  concordance$mean.cosine,
                  format.pval( concordance$p.value, digits = 2, eps = 1 / ( radar.n.perm + 1 ) ),
                  radar.n.perm ) )

r.label <- if ( radar.metric == "clr" ) "Centred log2 ratio of laser-bank norm\n(beads / cells)" else
  "Change in laser-bank norm (beads \u2212 cells)"
p1c <- plot.laser.radar(
  shift.matrix,
  instrument.colors = instrument.colors,
  r.label = r.label,
  title = "Cell \u2192 bead shift by laser, PerCP-Cy5.5",
  subtitle = sprintf( "Mean pairwise cosine %.2f; permutation p %s",
                      concordance$mean.cosine,
                      format.pval( concordance$p.value, digits = 2, eps = 1 / ( radar.n.perm + 1 ) ) ),
  base.size = base.size
)
save.fig.panel( p1c, "Fig1c_laser_space_shift_PerCP-Cy5.5", panel.size$c )

# tables for panel c
utils::write.csv( shift.matrix, file.path( output.dir, "tables", "Fig1c_laser_shift.csv" ) )
utils::write.csv(
  data.frame(
    Cytometer = radar.instruments,
    Angle.total.deg   = vapply( laser.shifts, `[[`, numeric( 1 ), "angle.total" ),
    Angle.between.deg = vapply( laser.shifts, `[[`, numeric( 1 ), "angle.between" ),
    do.call( rbind, lapply( laser.shifts, function( x )
      stats::setNames( x$angle.within, paste0( "Angle.within.", names( x$angle.within ), ".deg" ) ) ) ),
    Cosine.leave.one.out = concordance$leave.one.out[ radar.instruments ],
    check.names = FALSE
  ),
  file.path( output.dir, "tables", "Fig1c_laser_space_angles.csv" ), row.names = FALSE
)
utils::write.csv( concordance$pairwise, file.path( output.dir, "tables", "Fig1c_pairwise_cosine.csv" ) )
utils::write.csv(
  data.frame( Mean.pairwise.cosine = concordance$mean.cosine, P.value = concordance$p.value,
              N.permutations = radar.n.perm, Metric = radar.metric,
              Lasers = paste( radar.lasers, collapse = ";" ) ),
  file.path( output.dir, "tables", "Fig1c_concordance_test.csv" ), row.names = FALSE
)
utils::write.csv(
  do.call( rbind, lapply( radar.instruments, function( cy ) data.frame(
    Cytometer = cy, Laser = radar.lasers,
    Cell.norm = laser.shifts[[ cy ]]$cell.norms, Bead.norm = laser.shifts[[ cy ]]$bead.norms ) ) ),
  file.path( output.dir, "tables", "Fig1c_laser_bank_norms.csv" ), row.names = FALSE
)

# =============================================================================
# Panels d and e: unmixing error versus brightness (Aurora, 2025-07-09 data)
# =============================================================================

if ( !fig1de.fluor %in% names( fluor.info.1c ) )
  stop( "`fig1de.fluor` must be one of: ", paste( names( fluor.info.1c ), collapse = ", " ) )
fi <- fluor.info.1c[[ fig1de.fluor ]]
fluor.tag <- gsub( "[^A-Za-z0-9.-]", "", fig1de.fluor )

asp.1c <- get.autospectral.param( cytometer = "Aurora" )

# 1. files -------------------------------------------------------------------
cell.files <- find.titration.files( raw.dir.1c, fi$plate, fi$cell.dir, "Unstained cells" )
bead.files <- find.titration.files( raw.dir.1c, fi$plate, fi$bead.dir, "Unstained beads" )
amounts <- names( cell.files$stained )
brightest <- amounts[ length( amounts ) ]
message( "Cell titration files (ul): ", paste( amounts, collapse = ", " ) )
message( "Bead titration files (ul): ", paste( names( bead.files$stained ), collapse = ", " ) )

# 2. 35-colour panel ----------------------------------------------------------
library.spectra <- as.matrix( utils::read.csv( panel.spectra.file, check.names = FALSE, row.names = 1 ) )
missing.fluors <- setdiff( fluors.35C, rownames( library.spectra ) )
if ( length( missing.fluors ) > 0 )
  stop( "Fluorophores missing from ", panel.spectra.file, ": ",
        paste( missing.fluors, collapse = ", " ) )
panel <- library.spectra[ fluors.35C, , drop = FALSE ]
peak.channel <- colnames( panel )[ which.max( panel[ fig1de.fluor, ] ) ]
message( "Peak channel for ", fig1de.fluor, ": ", peak.channel )

# 3. reference spectra (cache) ------------------------------------------------
cache.file <- file.path( output.dir, "cache", sprintf( "fig1de_%s_reference_spectra.rds", fluor.tag ) )
if ( reuse.cache && file.exists( cache.file ) ) {
  message( "Reading cached reference spectra: ", cache.file )
  ref <- readRDS( cache.file )
} else {
  work.dir <- file.path( output.dir, "work", fluor.tag )
  dir.create( work.dir, recursive = TRUE, showWarnings = FALSE )
  ref <- list(
    cells = extract.reference.spectrum(
      stained.files = cell.files$stained, unstained.file = cell.files$unstained,
      spectrum.from = brightest, fluorophore = fig1de.fluor, peak.channel = peak.channel,
      particle.type = "cells", asp = asp.1c, work.dir = work.dir, label = "Cells"
    ),
    beads = extract.reference.spectrum(
      stained.files = bead.files$stained, unstained.file = bead.files$unstained,
      spectrum.from = names( bead.files$stained )[ length( bead.files$stained ) ],
      fluorophore = fig1de.fluor, peak.channel = peak.channel,
      particle.type = "beads", asp = asp.1c, work.dir = work.dir, label = "Beads"
    )
  )
  saveRDS( ref, cache.file )
}

cell.spectrum <- ref$cells$spectrum[ , colnames( panel ), drop = FALSE ]
bead.spectrum <- ref$beads$spectrum[ , colnames( panel ), drop = FALSE ]
spectrum.cosine <- cosine.similarity( rbind( cell.spectrum, bead.spectrum ) )[ 1, 2 ]
message( sprintf( "Cosine similarity, cell vs bead %s spectrum: %.4f", fig1de.fluor, spectrum.cosine ) )

panel.cells <- panel
panel.cells[ fig1de.fluor, ] <- cell.spectrum
panel.beads <- panel
panel.beads[ fig1de.fluor, ] <- bead.spectrum

utils::write.csv( rbind( Cells = cell.spectrum, Beads = bead.spectrum ),
                  file.path( output.dir, "tables", sprintf( "Fig1de_%s_reference_spectra.csv", fluor.tag ) ) )

# 4. unmix the cell titration (gated), select positive events, SSI ------------
tit <- unmix.titration(
  stained.files = cell.files$stained,
  unstained.file = cell.files$unstained,
  spectra.list = list( cells = panel.cells, beads = panel.beads ),
  fluorophore = fig1de.fluor,
  asp = asp.1c,
  gate = ref$cells$gate,
  positive.quantile = positive.quantile,
  positive.multiplier = positive.multiplier,
  positive.space = positive.space,
  y.from = "beads"
)
y.channel <- tit$y.channel
message( "Condition numbers: cells = ", round( tit$condition.number[[ "cells" ]], 2 ),
         ", beads = ", round( tit$condition.number[[ "beads" ]], 2 ) )
utils::write.csv( tit$ssi, file.path( output.dir, "tables", sprintf( "Fig1de_%s_SSI.csv", fluor.tag ) ),
                  row.names = FALSE )

# 5. spillover spread: Residual Model, brightest cell file, correct matrix ----
spread.fit <- calculate.residual.spread(
  stained.data = tit$raw[[ brightest ]],
  unstained.data = tit$raw[[ "Unstained" ]],
  spectra = panel.cells,
  fluorophore = fig1de.fluor,
  peak.channel = peak.channel,
  asp = asp.1c
)
message( "Residual-model spread ", fig1de.fluor, " -> ", y.channel, ": ",
         signif( spread.fit$spread[[ y.channel ]], 4 ),
         " (", spread.fit$n.events, " events; AF collinear: ", spread.fit$af.collinear, ")" )
utils::write.csv(
  data.frame( Channel = names( spread.fit$spread ), Spread = spread.fit$spread,
              Intercept = spread.fit$intercept ),
  file.path( output.dir, "tables", sprintf( "Fig1de_%s_residual_spread.csv", fluor.tag ) ),
  row.names = FALSE
)

unstained.cells.unmixed <- tit$unmixed$cells[[ "Unstained" ]]
spread.params <- list(
  ss       = spread.fit$spread[[ y.channel ]],
  flat.pos = unname( stats::quantile( unstained.cells.unmixed[ , y.channel ], 0.995 ) ),
  flat.neg = unname( stats::quantile( unstained.cells.unmixed[ , y.channel ], 0.005 ) )
)

# 6. axis limits (as in compare.unmix()) ---------------------------------------
all.pos.y <- unlist( lapply( names( tit$unmixed ), function( e )
  unlist( lapply( amounts, function( a ) tit$unmixed[[ e ]][[ a ]][ tit$keep[[ e ]][[ a ]], y.channel ] ) ) ) )
all.neg.x <- unlist( lapply( names( tit$unmixed ), function( e ) tit$unmixed[[ e ]][[ "Unstained" ]][ , fig1de.fluor ] ) )
x.plot.min <- min( -1000, 2 * stats::quantile( all.neg.x, 0.01 ) )
y.plot.min <- min( -1000, 2 * stats::quantile( all.pos.y, 0.01 ),
                   2 * stats::quantile( unstained.cells.unmixed[ , y.channel ], 0.01 ) )

# 7. colours and legend labels -------------------------------------------------
pop.colors <- c(
  Negative = "black",
  stats::setNames(
    viridis::viridis( length( amounts ), option = viridis.option,
                      begin = viridis.range[ 1 ], end = viridis.range[ 2 ] ),
    amounts )
)

make.populations <- function( entry )
  c( list( Negative = tit$unmixed[[ entry ]][[ "Unstained" ]] ),
     stats::setNames(
       lapply( amounts, function( a )
         tit$unmixed[[ entry ]][[ a ]][ tit$keep[[ entry ]][[ a ]], , drop = FALSE ] ),
       amounts ) )

make.legend.labels <- function( entry ) {
  s <- tit$ssi[ tit$ssi$entry == entry, ]
  c( Negative = "Unstained",
     stats::setNames(
       sprintf( "%s µl  (SSI %s)", s$amount, formatC( s$SSI, digits = 1, format = "f" ) ),
       s$amount ) )
}

make.biplot <- function( entry, title ) {
  create.brightness.biplot(
    populations = make.populations( entry ),
    x.dim = fig1de.fluor,
    y.dim = y.channel,
    asp = asp.1c,
    population.colors = pop.colors,
    legend.labels = make.legend.labels( entry ),
    spread = spread.params,
    x.min = x.plot.min,
    y.min = y.plot.min,
    title = title
  )
}

p1d <- make.biplot( "cells", paste0( "Cell-derived spectrum (correct): ", fig1de.fluor ) )
p1e <- make.biplot( "beads", paste0( "Bead-derived spectrum (incorrect): ", fig1de.fluor ) )

save.fig.panel( p1d, sprintf( "Fig1d_%s_cell_spectrum_unmixing", fluor.tag ), panel.size$d )
save.fig.panel( p1e, sprintf( "Fig1e_%s_bead_spectrum_unmixing", fluor.tag ), panel.size$e )

# =============================================================================
# Stitched figure
# =============================================================================

top <- cowplot::plot_grid(
  p1a, p1b, labels = c( "a", "b" ), ncol = 2, rel_widths = c( 1, 2 ), label_size = 14
)
bottom <- cowplot::plot_grid(
  p1c,
  cowplot::plot_grid( p1d, p1e, labels = c( "d", "e" ), ncol = 2, align = "h", axis = "tb",
                      label_size = 14 ),
  labels = c( "c", "" ), ncol = 2, rel_widths = c( 1, 2 ), label_size = 14
)
figure.1 <- cowplot::plot_grid( top, bottom, ncol = 1, rel_heights = c( 1.4, 1 ) )

save.panel( figure.1, "Figure_1", output.dir, figure.size[[ "width" ]], figure.size[[ "height" ]],
            units = size.units, dpi = dpi )
