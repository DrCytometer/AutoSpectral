# bead_angle_data.R
#
# Compile the bead-vs-cell spectral angle (one value per fluorophore, per bead
# type, per cytometer) from the bead_cell_comparison_results.rds files written
# by AutoSpectral's run.bead.cell.comparison().
#
#   Comparison[[ bead ]]$Angle : one-column matrix (rownames = fluorophores)
#                                from assess.mismatch.angle(): the angle, in
#                                degrees, between the cell reference spectrum
#                                and the bead reference spectrum (row 1 of
#                                each fluorophore's spectral-variant matrix).

#' Read the spectral angles of one or more cytometers into one long table
#'
#' @param results.dir Root results folder (e.g. ".../results_density_20260901").
#' @param cytometers Named list; each element has `results` (sub-folder name,
#'   e.g. "aurora"). Names are the display names ("Aurora", "A8", ...).
#' @param beads Optional character vector: bead types to keep (NULL = all
#'   found in the results).
#' @param exclude.fluorophores Optional character vector of fluorophores to drop
#'   (on top of the exclusions already applied by run.bead.cell.comparison()).
#' @param recompute.check Logical. Recompute each angle from the stored
#'   Extraction variants with AutoSpectral::assess.mismatch.angle() and stop if
#'   any differs from the stored value by more than `check.tolerance` degrees.
#'   Requires AutoSpectral to be loaded.
#' @param check.tolerance Degrees.
#' @return data.frame: Cytometer, Bead, Fluorophore, Angle (degrees), Cosine.
load.bead.angles <- function(
    results.dir,
    cytometers,
    beads = NULL,
    exclude.fluorophores = NULL,
    recompute.check = TRUE,
    check.tolerance = 1e-3,
    reference.type = "Cells"
) {

  out <- lapply( names( cytometers ), function( cy ) {

    f <- file.path( results.dir, cytometers[[ cy ]]$results, "bead_cell_comparison_results.rds" )
    if ( !file.exists( f ) ) stop( "File not found: ", f, call. = FALSE )
    res <- readRDS( f )

    bead.types <- names( res$Comparison )
    if ( !is.null( beads ) ) bead.types <- intersect( bead.types, beads )

    do.call( rbind, lapply( bead.types, function( b ) {

      ang <- res$Comparison[[ b ]]$Angle
      cs  <- res$Comparison[[ b ]]$Cosine
      fl  <- rownames( ang )

      if ( recompute.check ) {
        ref.var  <- res$Extraction[[ reference.type ]]$variants$variants[ fl ]
        bead.var <- res$Extraction[[ b ]]$variants$variants[ fl ]
        re <- assess.mismatch.angle( ref.var, bead.var )
        dev <- max( abs( re[ match( fl, rownames( re ) ), 1 ] - ang[ fl, 1 ] ) )
        if ( !is.finite( dev ) || dev > check.tolerance )
          stop( sprintf( "%s / %s: recomputed angles differ from stored values by up to %.3g deg",
                         cy, b, dev ), call. = FALSE )
      }

      data.frame(
        Cytometer   = cy,
        Bead        = b,
        Fluorophore = fl,
        Angle       = as.numeric( ang[ fl, 1 ] ),
        Cosine      = as.numeric( cs[ fl, 1 ] ),
        stringsAsFactors = FALSE
      )
    } ) )
  } )

  d <- do.call( rbind, out )
  if ( !is.null( exclude.fluorophores ) )
    d <- d[ !d$Fluorophore %in% exclude.fluorophores, , drop = FALSE ]
  d <- d[ is.finite( d$Angle ) & d$Angle > 0, , drop = FALSE ]
  rownames( d ) <- NULL
  d
}

#' Descriptive summary of the angles per cytometer and bead type
#'
#' @param d Output of load.bead.angles().
#' @return data.frame with n, median, IQR and geometric mean (degrees).
summarize.bead.angles <- function( d ) {
  keys <- unique( d[ , c( "Cytometer", "Bead" ) ] )
  do.call( rbind, lapply( seq_len( nrow( keys ) ), function( i ) {
    x <- d$Angle[ d$Cytometer == keys$Cytometer[ i ] & d$Bead == keys$Bead[ i ] ]
    q <- stats::quantile( x, c( 0.25, 0.5, 0.75 ) )
    data.frame(
      Cytometer = keys$Cytometer[ i ], Bead = keys$Bead[ i ],
      N.fluorophores = length( x ),
      Median.deg = q[[ 2 ]], Q1.deg = q[[ 1 ]], Q3.deg = q[[ 3 ]],
      Geometric.mean.deg = 2^mean( log2( x ) ),
      Max.deg = max( x ),
      stringsAsFactors = FALSE
    )
  } ) )
}
