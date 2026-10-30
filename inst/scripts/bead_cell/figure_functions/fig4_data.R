# fig4_data.R
#
# Data assembly for Figure 4 ("Is it the fluorophores?"):
#   - fluorophore class (the `class` column of AutoSpectral's
#     fluorophore_database.csv, names resolved through the database synonyms)
#   - spectral "messiness" of each fluorophore's cell spectrum (two scores)
#   - excitation-emission separation (emission peak - excitation laser)
#   - within-control spectral variability of cells and of every bead type
#
# Everything is read from the bead_cell_comparison_results.rds files written by
# AutoSpectral's run.bead.cell.comparison() (same files as Figure 2).

# ---- fluorophore names and classes ------------------------------------------

#' Map every fluorophore name / synonym in fluorophore_database.csv to its
#' canonical name (the `fluorophore` column).
#' @param fluor.db data.frame read from fluorophore_database.csv.
#' @return Named character vector: names = any known spelling, value = canonical.
fluor.synonym.map <- function( fluor.db ) {
  syn.cols <- intersect( c( "fluorophore", paste0( "synonym", 1:9 ) ), names( fluor.db ) )
  out <- character( 0 )
  for ( col in syn.cols ) {
    s <- trimws( as.character( fluor.db[[ col ]] ) )
    ok <- !is.na( s ) & s != "" & !( s %in% names( out ) )
    out <- c( out, stats::setNames( fluor.db$fluorophore[ ok ], s[ ok ] ) )
  }
  out
}

#' Canonical fluorophore name for each input name (unknown names unchanged)
canonical.fluor <- function( x, syn.map ) {
  y <- unname( syn.map[ x ] )
  ifelse( is.na( y ), x, y )
}

#' Fluorophore class for a set of fluorophores
#'
#' Read from the `class` column of AutoSpectral's fluorophore_database.csv
#' (names matched through the database synonyms) and translated to the class
#' codes used in the figures.
#'
#' @param fluorophores Fluorophore names as used in the results files.
#' @param fluor.db data.frame from fluorophore_database.csv.
#' @param db.class.map Named vector translating database classes to the class
#'   codes used in the figure.
#' @param overrides Named vector fluorophore -> class code, applied last
#'   (default NULL: none).
#' @return data.frame(Fluorophore, Class, Class.source)
get.fluor.classes <- function(
    fluorophores,
    fluor.db,
    db.class.map = c( "Small molecule" = "Small", "SIRIGEN" = "Brilliant",
                      "SIRIGEN tandem" = "Brilliant_tandem", "Protein" = "Protein",
                      "Protein tandem" = "Protein_tandem", "Real" = "Real" ),
    overrides = NULL
) {
  syn <- fluor.synonym.map( fluor.db )
  canon <- canonical.fluor( fluorophores, syn )

  cls <- unname( db.class.map[ fluor.db$class[ match( canon, fluor.db$fluorophore ) ] ] )
  src <- rep( "fluorophore_database.csv", length( cls ) )

  if ( !is.null( overrides ) ) {
    ov.canon <- canonical.fluor( names( overrides ), syn )
    hit <- match( canon, ov.canon )
    cls[ !is.na( hit ) ] <- overrides[ hit[ !is.na( hit ) ] ]
    src[ !is.na( hit ) ] <- "override"
  }

  if ( any( is.na( cls ) ) )
    warning( "No class for: ", paste( fluorophores[ is.na( cls ) ], collapse = ", " ),
             call. = FALSE )

  data.frame( Fluorophore = fluorophores, Class = cls, Class.source = src,
              stringsAsFactors = FALSE )
}

# ---- spectral messiness -----------------------------------------------------

# spectral.messiness() (Linf and PR scores; Hill numbers of order infinity
# and 2 divided by the number of detectors) is an exported AutoSpectral
# function (R/spectral_messiness.R).

# ---- excitation-emission separation -----------------------------------------

#' Emission peak, excitation laser and their separation for each fluorophore
#'
#' The database holds the nominal emission peak (`nominal.wavelength`) and the
#' primary excitation laser, not the excitation maximum, so the "Stokes shift"
#' is the separation between the laser line and the emission peak. For a
#' tandem this is the donor-excitation to acceptor-emission separation seen by
#' the cytometer.
#'
#' @param fluorophores Names as used in the results.
#' @param fluor.db data.frame from fluorophore_database.csv.
#' @param laser.nm Named vector: laser label (database `excitation.laser`) -> nm.
#' @return data.frame(Fluorophore, Laser, Laser.nm, Emission.nm, Stokes.nm)
get.fluor.stokes <- function(
    fluorophores,
    fluor.db,
    laser.nm = c( UV = 355, Violet = 405, Blue = 488, YellowGreen = 561, YG = 561, Red = 640 )
) {
  syn <- fluor.synonym.map( fluor.db )
  canon <- canonical.fluor( fluorophores, syn )
  i <- match( canon, fluor.db$fluorophore )
  laser <- trimws( as.character( fluor.db$excitation.laser[ i ] ) )
  em <- suppressWarnings( as.numeric( fluor.db$nominal.wavelength[ i ] ) )
  lnm <- unname( laser.nm[ laser ] )
  out <- data.frame( Fluorophore = fluorophores, Laser = laser, Laser.nm = lnm,
                     Emission.nm = em, Stokes.nm = em - lnm, stringsAsFactors = FALSE )
  miss <- is.na( out$Stokes.nm )
  if ( any( miss ) )
    warning( "No emission peak or laser wavelength for: ",
             paste( fluorophores[ miss ], collapse = ", " ), call. = FALSE )
  out
}

# ---- per-cytometer properties and variability -------------------------------

#' Read one cytometer's results file (cached in an environment)
.fig4.cache <- new.env( parent = emptyenv() )
read.bead.cell.results <- function( results.dir, cyto.dir ) {
  f <- file.path( results.dir, cyto.dir, "bead_cell_comparison_results.rds" )
  if ( !file.exists( f ) ) stop( "File not found: ", f, call. = FALSE )
  if ( is.null( .fig4.cache[[ f ]] ) ) .fig4.cache[[ f ]] <- readRDS( f )
  .fig4.cache[[ f ]]
}

#' Messiness of each fluorophore's spectrum, per cytometer
#'
#' @param results.dir,cytometers As for load.bead.angles().
#' @param particle Particle type whose spectrum is scored (default "Cells").
#' @param exclude.fluorophores Fluorophores to drop ("AF" is always dropped).
#' @return data.frame(Cytometer, Fluorophore, N.detectors, Messiness.Linf, Messiness.PR)
load.fluor.messiness <- function( results.dir, cytometers, particle = "Cells",
                                  exclude.fluorophores = NULL ) {
  do.call( rbind, lapply( names( cytometers ), function( cy ) {
    res <- read.bead.cell.results( results.dir, cytometers[[ cy ]]$results )
    sp <- res$Extraction[[ particle ]]$spectra
    if ( is.null( sp ) ) stop( "No '", particle, "' spectra for ", cy, call. = FALSE )
    sp <- sp[ setdiff( rownames( sp ), c( "AF", exclude.fluorophores ) ), , drop = FALSE ]
    m <- t( apply( sp, 1, spectral.messiness ) )
    data.frame( Cytometer = cy, Fluorophore = rownames( sp ), N.detectors = ncol( sp ),
                Messiness.Linf = m[ , "Linf" ], Messiness.PR = m[ , "PR" ],
                stringsAsFactors = FALSE, row.names = NULL )
  } ) )
}

#' Within-control spectral variability of every particle type, per cytometer
#'
#' Uses the stored spectral variants (Extraction[[particle]]$variants$variants),
#' renormalised to unit L2 norm exactly as run.bead.cell.comparison() does
#' before its comparisons, and AutoSpectral::assess.variability(). The bead
#' values are checked against the stored Comparison[[bead]]$Variability.
#'
#' @param particles Particle types to include (NULL = all in the results).
#' @param min.variant.rows Fluorophores with fewer variant rows are dropped
#'   (same rule as run.bead.cell.comparison()).
#' @param recompute.check Stop if a recomputed bead value differs from the
#'   stored one by more than `check.tolerance`.
#' @return Long data.frame(Cytometer, Particle, Fluorophore, rSD, Scaled,
#'   Sim_rSD, Denoised_rSD, N.variants)
load.particle.variability <- function(
    results.dir, cytometers, particles = NULL, exclude.fluorophores = NULL,
    min.variant.rows = 3L, recompute.check = TRUE, check.tolerance = 1e-8
) {
  do.call( rbind, lapply( names( cytometers ), function( cy ) {
    res <- read.bead.cell.results( results.dir, cytometers[[ cy ]]$results )
    pts <- names( res$Extraction )
    if ( !is.null( particles ) ) pts <- intersect( particles, pts )

    do.call( rbind, lapply( pts, function( p ) {
      vl <- res$Extraction[[ p ]]$variants$variants
      n.rows <- vapply( vl, function( m ) if ( is.null( m ) ) 0L else nrow( m ), integer( 1 ) )
      keep <- setdiff( names( vl )[ n.rows >= min.variant.rows ], c( "AF", exclude.fluorophores ) )
      if ( length( keep ) == 0 ) return( NULL )
      vl <- lapply( vl[ keep ], l2.normalize.spectra )
      v <- assess.variability( vl )

      if ( recompute.check && p %in% names( res$Comparison ) ) {
        st <- res$Comparison[[ p ]]$Variability
        common <- intersect( rownames( st ), rownames( v ) )
        dev <- max( abs( st[ common, "Denoised_rSD" ] - v[ common, "Denoised_rSD" ] ) )
        if ( !is.finite( dev ) || dev > check.tolerance )
          stop( sprintf( "%s / %s: recomputed variability differs from stored values by up to %.3g",
                         cy, p, dev ), call. = FALSE )
      }

      data.frame( Cytometer = cy, Particle = p, Fluorophore = rownames( v ),
                  v, N.variants = n.rows[ rownames( v ) ],
                  stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE )
    } ) )
  } ) )
}

#' Geometric mean of the angle over bead types, per cytometer and fluorophore
#' (one point per fluorophore per cytometer for display).
#' @param d data.frame with Cytometer, Fluorophore, Angle (+ extra columns
#'   constant within Cytometer x Fluorophore, carried over from the first row).
fluor.mean.angle <- function( d, carry = character( 0 ) ) {
  key <- paste( d$Cytometer, d$Fluorophore, sep = "\r" )
  out <- do.call( rbind, lapply( split( seq_len( nrow( d ) ), key ), function( i ) {
    data.frame( Cytometer = as.character( d$Cytometer[ i[ 1 ] ] ),
                Fluorophore = d$Fluorophore[ i[ 1 ] ],
                Angle = 2^mean( log2( d$Angle[ i ] ) ),
                N.beads = length( i ),
                d[ i[ 1 ], carry, drop = FALSE ],
                stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE )
  } ) )
  rownames( out ) <- NULL
  out
}
