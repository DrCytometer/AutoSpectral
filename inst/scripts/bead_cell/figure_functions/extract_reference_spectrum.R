# extract_reference_spectrum.R
#
# Extracts a single-colour reference spectrum from a titration series using
# the same legacy AutoSpectral pipeline as 20260901_bead_cell_analysis.R /
# run.bead.cell.comparison():
#   define.gate.density() -> define.flow.control() -> clean.controls()
#   -> get.fluorophore.spectra()
# i.e. density FSC/SSC gate, scatter-matched universal-negative (unstained)
# subtraction, spectrum from the brightest positive events.
#
# The raw FCS files are linked (or copied) into a staging folder so that a
# control file can be written next to them without touching the raw data.

#' Find the FCS files of one titration series
#'
#' @param raw.dir Folder holding the plate folders (".../Raw").
#' @param plate Plate folder ("Plate_001").
#' @param particle.dir Exact name of the folder with the stained series
#'   ("Cells PerCP-Cy5.5"; NOT the "... unlabelled" folders).
#' @param unstained.pattern Regular expression identifying the unstained file
#'   in `<plate>/Unstained` ("Unstained cells" / "Unstained beads").
#' @return list(stained = named character vector of paths, names = antibody
#'   amount in ul, ordered ascending; unstained = path; amount = numeric vector).
find.titration.files <- function( raw.dir, plate, particle.dir, unstained.pattern ) {

  stained.dir <- file.path( raw.dir, plate, particle.dir )
  if ( !dir.exists( stained.dir ) ) stop( "Folder not found: ", stained.dir, call. = FALSE )

  stained <- list.files( stained.dir, pattern = "\\.fcs$", full.names = TRUE )
  if ( length( stained ) == 0 ) stop( "No FCS files in ", stained.dir, call. = FALSE )

  # "G2 6.7_Plate_001.fcs": plate location, then antibody amount in ul
  amount <- suppressWarnings(
    as.numeric( sub( "^[A-H][0-9]+ ([0-9.]+)_Plate.*$", "\\1", basename( stained ) ) )
  )
  if ( anyNA( amount ) )
    stop( "Could not read the antibody amount from: ",
          paste( basename( stained )[ is.na( amount ) ], collapse = ", " ), call. = FALSE )

  o <- order( amount )
  stained <- stained[ o ]
  amount <- amount[ o ]
  names( stained ) <- as.character( amount )

  unstained.dir <- file.path( raw.dir, plate, "Unstained" )
  unstained <- list.files( unstained.dir, pattern = "\\.fcs$", full.names = TRUE )
  unstained <- unstained[ grepl( unstained.pattern, basename( unstained ), ignore.case = TRUE ) ]
  if ( length( unstained ) != 1 )
    stop( "Expected exactly one unstained file matching '", unstained.pattern,
          "' in ", unstained.dir, call. = FALSE )

  list( stained = stained, unstained = unstained, amount = amount )
}

# link (or copy) files into a staging folder under new names
.stage.fcs <- function( files, new.names, dir ) {
  if ( !dir.exists( dir ) ) dir.create( dir, recursive = TRUE )
  dest <- file.path( dir, new.names )
  for ( i in seq_along( files ) ) {
    if ( file.exists( dest[ i ] ) ) next
    ok <- suppressWarnings( file.symlink( normalizePath( files[ i ] ), dest[ i ] ) )
    if ( !isTRUE( ok ) ) file.copy( files[ i ], dest[ i ] )
  }
  dest
}

#' Extract a reference spectrum from a titration series (legacy pipeline)
#'
#' @param stained.files Named character vector, all stained files of the series
#'   (names = antibody amount). All are used to define the density gate.
#' @param unstained.file Path of the matched unstained file (universal negative).
#' @param spectrum.from Name (in `stained.files`) of the file used for the
#'   spectrum, normally the brightest.
#' @param fluorophore Fluorophore name as used in the spectral library.
#' @param peak.channel Peak detector of the fluorophore.
#' @param particle.type "cells" or "beads".
#' @param asp AutoSpectral parameter list for the cytometer.
#' @param work.dir Folder for staging, gates and intermediate files.
#' @param label Label used in file names ("Cells", "Beads").
#' @param gate Optional pre-defined gate (list with x, y). If `NULL` a density
#'   gate is defined from all files in `stained.files` plus the unstained file.
#' @param verbose Logical.
#' @return list(spectrum = 1-row matrix (peak = 1), gate = gate polygon).
extract.reference.spectrum <- function(
    stained.files,
    unstained.file,
    spectrum.from,
    fluorophore,
    peak.channel,
    particle.type = c( "cells", "beads" ),
    asp,
    work.dir,
    label = particle.type,
    gate = NULL,
    include.unstained.in.gate = TRUE,
    verbose = TRUE
) {

  particle.type <- match.arg( particle.type )
  if ( !spectrum.from %in% names( stained.files ) )
    stop( "`spectrum.from` must be one of: ", paste( names( stained.files ), collapse = ", " ),
          call. = FALSE )

  gate.name <- paste0( label, "_gate" )
  stage.dir <- file.path( work.dir, paste0( "staging_", label ) )
  dir.create( stage.dir, recursive = TRUE, showWarnings = FALSE )

  # stage files under simple names
  fluor.tag <- gsub( "[^A-Za-z0-9.-]", "", fluorophore )
  stained.names <- sprintf( "%s_%sul_%s.fcs", fluor.tag, names( stained.files ), label )
  names( stained.names ) <- names( stained.files )
  unstained.name <- sprintf( "Unstained_%s.fcs", label )
  .stage.fcs( stained.files, stained.names, stage.dir )
  .stage.fcs( unstained.file, unstained.name, stage.dir )

  # per-run output folders (the legacy pipeline writes into these)
  asp$figure.gate.dir               <- file.path( work.dir, label, "figure_gate" )
  asp$figure.clean.control.dir      <- file.path( work.dir, label, "figure_clean_controls" )
  asp$figure.spectral.ribbon.dir    <- file.path( work.dir, label, "figure_spectral_ribbon" )
  asp$figure.spectra.dir            <- file.path( work.dir, label, "figure_spectra" )
  asp$figure.similarity.heatmap.dir <- file.path( work.dir, label, "figure_similarity_heatmap" )
  asp$figure.scatter.dir.base       <- file.path( work.dir, label, "figure_scatter" )
  asp$table.spectra.dir             <- file.path( work.dir, label, "table_spectra" )
  asp$variant.dir                   <- file.path( work.dir, label, "figure_spectral_variants" )
  for ( d in c( asp$figure.gate.dir, asp$figure.clean.control.dir,
                asp$figure.spectral.ribbon.dir, asp$figure.spectra.dir,
                asp$figure.similarity.heatmap.dir, asp$figure.scatter.dir.base,
                asp$table.spectra.dir, asp$variant.dir ) )
    dir.create( d, recursive = TRUE, showWarnings = FALSE )

  # density gate from the whole series (pooled), unless supplied
  if ( is.null( gate ) ) {
    gate.files <- c( stained.names, if ( include.unstained.in.gate ) unstained.name )
    gate.table <- data.frame(
      filename           = unname( gate.files ),
      fluorophore        = NA_character_,
      marker             = NA_character_,
      channel            = NA_character_,
      control.type       = particle.type,
      universal.negative = NA_character_,
      large.gate         = NA,
      gate.name          = gate.name,
      gate.define        = TRUE,
      is.viability       = NA,
      stringsAsFactors   = FALSE
    )
    gate <- define.gate.density(
      control.file = NULL,
      control.dir  = stage.dir,
      asp          = asp,
      gate.name    = gate.name,
      output.dir   = asp$figure.gate.dir,
      verbose      = verbose,
      control.table = gate.table
    )
  }

  # control file: the one stained file used for the spectrum + the unstained.
  # Cells: unstained row = "AF"; beads: unstained row = "Negative"
  # (the convention of finalize.control.file() in run_bead_cell_comparison.R)
  unstained.label <- if ( particle.type == "cells" ) "AF" else "Negative"
  ctrl <- data.frame(
    filename           = c( stained.names[[ spectrum.from ]], unstained.name ),
    fluorophore        = c( fluorophore, unstained.label ),
    marker             = NA_character_,
    channel            = c( peak.channel, NA ),
    control.type       = particle.type,
    universal.negative = c( unstained.name, unstained.name ),
    large.gate         = NA,
    gate.name          = gate.name,
    gate.define        = c( TRUE, NA ),
    is.viability       = NA,
    stringsAsFactors   = FALSE
  )
  ctrl.path <- file.path( stage.dir, "fcs_control_file.csv" )
  utils::write.csv( ctrl, ctrl.path, row.names = FALSE, na = "" )

  flow.control <- define.flow.control(
    control.dir      = stage.dir,
    control.def.file = ctrl.path,
    asp              = asp,
    gate             = TRUE,
    gating.system    = "density",
    gate.list        = stats::setNames( list( gate ), gate.name ),
    verbose          = verbose,
    allow.duplicate.controls = TRUE
  )

  flow.control <- clean.controls(
    flow.control = flow.control,
    asp          = asp,
    af.remove    = FALSE,
    main.figures = FALSE,
    verbose      = verbose
  )

  spectra <- get.fluorophore.spectra(
    flow.control   = flow.control,
    asp            = asp,
    use.clean.expr = TRUE,
    title          = label,
    figures        = FALSE
  )

  row.idx <- which( rownames( spectra ) == fluorophore )
  if ( length( row.idx ) != 1 )
    row.idx <- grep( fluorophore, rownames( spectra ), fixed = TRUE )
  if ( length( row.idx ) != 1 )
    stop( "Could not identify the '", fluorophore, "' row in the extracted spectra: ",
          paste( rownames( spectra ), collapse = ", " ), call. = FALSE )

  spectrum <- spectra[ row.idx, , drop = FALSE ]
  rownames( spectrum ) <- fluorophore

  list( spectrum = spectrum, gate = gate )
}
