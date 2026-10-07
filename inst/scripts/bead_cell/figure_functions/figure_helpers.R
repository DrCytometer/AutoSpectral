# figure_helpers.R
#
# Small shared helpers for the bead-vs-cell manuscript figures:
#  - spectrum normalisation (L2 via AutoSpectral::l2.normalize.spectra())
#  - detector -> laser lookup from AutoSpectral's cytometer_database.csv
#  - loading cell/bead spectra from a run.bead.cell.comparison() results folder
#  - saving individual panels as 300 dpi PNG with controllable size

# Default laser colours (used to colour detector labels) ---------------------
default.laser.colors <- c(
  DeepUV      = "#4B2A7B",
  UV          = "#8E44AD",
  Violet      = "#3F51B5",
  Blue        = "#1E88E5",
  YellowGreen = "#7CB342",
  Red         = "#D32F2F",
  IR          = "#8D6E63",
  Unknown     = "grey50"
)

#' Normalise a spectrum
#'
#' @param x Named numeric vector.
#' @param method "max" (peak = 1), "L2" (unit Euclidean norm) or "none".
normalize.spectrum <- function( x, method = c( "max", "L2", "none" ) ) {
  method <- match.arg( method )
  switch( method,
    max  = x / max( x ),
    L2   = stats::setNames( as.numeric( l2.normalize.spectra( matrix( x, nrow = 1 ) ) ), names( x ) ),
    none = x )
}

#' Detector -> laser lookup from cytometer_database.csv
#'
#' @param cytometer.db Path to cytometer_database.csv
#'   (e.g. system.file("extdata", "cytometer_database.csv", package = "AutoSpectral")).
#' @param detector.col Column holding the detector names for the cytometer
#'   ("Aurora", "Discover", "ID7000", "Xenith", ...). The laser is read from
#'   the matching `<detector.col>_laser` column.
#' @return data.frame with columns Detector, Laser, in database (wavelength) order.
get.detector.lasers <- function( cytometer.db, detector.col ) {

  db <- utils::read.csv( cytometer.db, check.names = FALSE, stringsAsFactors = FALSE )
  laser.col <- paste0( detector.col, "_laser" )

  if ( !all( c( detector.col, laser.col ) %in% colnames( db ) ) )
    stop( "Columns '", detector.col, "' and '", laser.col, "' not found in ",
          cytometer.db, call. = FALSE )

  keep <- !is.na( db[[ detector.col ]] ) & trimws( db[[ detector.col ]] ) != ""

  data.frame(
    Detector = trimws( db[[ detector.col ]][ keep ] ),
    Laser    = trimws( db[[ laser.col ]][ keep ] ),
    stringsAsFactors = FALSE
  )
}

#' Load the cell and bead spectrum of one fluorophore from a results folder
#'
#' Reads `<results.dir>/<cytometer.dir>/<particle>/table_spectra/<particle>_autospectral_spectra.csv`
#' as written by run.bead.cell.comparison().
#'
#' @param results.dir Root results folder (e.g. ".../results_density_20260901").
#' @param cytometer.dir Sub-folder for the cytometer ("aurora", "a8", ...).
#' @param fluor.regex Regular expression matching exactly one row name.
#' @param cell.particle,bead.particle Folder names of the particle types.
#' @return list(cells = 1-row matrix, beads = 1-row matrix), both with the
#'   same detector columns.
load.cell.bead.spectra <- function(
    results.dir,
    cytometer.dir,
    fluor.regex,
    cell.particle = "Cells",
    bead.particle = "BD"
) {

  read.one <- function( particle ) {
    f <- file.path(
      results.dir, cytometer.dir, particle, "table_spectra",
      paste0( particle, "_autospectral_spectra.csv" )
    )
    if ( !file.exists( f ) ) stop( "File not found: ", f, call. = FALSE )
    spec <- as.matrix( utils::read.csv( f, check.names = FALSE, row.names = 1 ) )
    idx <- grep( fluor.regex, rownames( spec ) )
    if ( length( idx ) != 1 )
      stop( "Pattern '", fluor.regex, "' matched ", length( idx ),
            " rows in ", f, ": ", paste( rownames( spec )[ idx ], collapse = ", " ),
            call. = FALSE )
    spec[ idx, , drop = FALSE ]
  }

  cells <- read.one( cell.particle )
  beads <- read.one( bead.particle )

  if ( !identical( colnames( cells ), colnames( beads ) ) )
    stop( "Cell and bead spectra have different detector columns for ",
          cytometer.dir, call. = FALSE )

  list( cells = cells, beads = beads )
}

#' Save one panel as a PNG (300 dpi by default)
#'
#' @param plot ggplot / cowplot object.
#' @param file.name File stem (no extension).
#' @param output.dir Folder (created if needed).
#' @param width,height Size in `units`.
#' @param units "in", "cm" or "mm".
#' @param dpi Resolution. Default 300.
#' @return The saved file path, invisibly.
save.panel <- function(
    plot, file.name, output.dir, width, height, units = "in", dpi = 300
) {
  if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
  path <- file.path( output.dir, paste0( file.name, ".png" ) )
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  message( "Saved ", path )
  invisible( path )
}
