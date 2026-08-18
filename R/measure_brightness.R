# get_brightness_automated.R

#' @title Get Brightness (Automated)
#'
#' @description
#' Fast, non-interactive estimate of fluorophore brightness (background-
#' subtracted MFI) for every single-stained control listed in
#' \code{control.def.file}. Unlike \code{get.brightness()}, no scatter gate is
#' drawn: brightness is measured directly from the whole acquired file using
#' the top \code{n.cells} events in each fluorophore's peak detector, with
#' background subtracted using the matched universal-negative control.
#'
#' The peak detector for each fluorophore is taken from \code{spectra}
#' (the output of \code{get.spectra.automated()}) rather than the nominal
#' \code{channel} column in the control table, so brightness is always
#' measured on the same channel the extracted spectrum is actually built
#' from, including any legacy refinement.
#'
#' @param control.dir Character. Path to the single-stained control FCS files.
#' @param control.def.file Character. Path to the control definition CSV.
#' @param asp The AutoSpectral parameter list from \code{get.autospectral.param()}.
#'   Currently unused directly but accepted for interface consistency with
#'   \code{get.spectra.automated()} / \code{get.spectral.variants()}, and to
#'   allow future use (e.g. saturation thresholds).
#' @param spectra Numeric matrix as returned by \code{get.spectra.automated()};
#'   fluorophores in rows (may include \code{"AF"}, which is ignored),
#'   detectors in columns.
#' @param n.cells Integer, default \code{500L}. Number of brightest events
#'   (in the fluorophore's peak channel) averaged to estimate the positive
#'   MFI.
#' @param verbose Logical, default \code{TRUE}.
#'
#' @return A one-column numeric matrix of background-subtracted MFI values;
#'   rownames are fluorophore names, column name \code{"MFI"}.
#'
#' @export

get.brightness.automated <- function(
    control.dir,
    control.def.file,
    asp,
    spectra,
    n.cells = 500L,
    verbose = TRUE
) {

  if ( !dir.exists( control.dir ) )
    stop( "control.dir does not exist: ", control.dir, call. = FALSE )

  ctrl.path <- if ( file.exists( control.def.file ) ) {
    control.def.file
  } else {
    file.path( control.dir, control.def.file )
  }

  ctrl.tbl <- utils::read.csv(
    ctrl.path, stringsAsFactors = FALSE, strip.white = TRUE
  )

  fluorophores <- rownames( spectra )[ rownames( spectra ) != "AF" ]

  row.idx <- match( fluorophores, ctrl.tbl$fluorophore )

  if ( any( is.na( row.idx ) ) )
    stop(
      "Fluorophore(s) present in `spectra` but not found in control file: ",
      paste( fluorophores[ is.na( row.idx ) ], collapse = ", " ),
      call. = FALSE
    )

  fluor.files <- stats::setNames( ctrl.tbl$filename[ row.idx ], fluorophores )

  unstained.sources <- stats::setNames(
    lapply( ctrl.tbl$universal.negative[ row.idx ], .parse.unstained.source ),
    fluorophores
  )

  # cache raw unstained event matrices by filename, each loaded only once
  unstained.cache <- new.env( parent = emptyenv() )

  read.unstained <- function( uf ) {
    if ( !is.null( unstained.cache[[ uf ]] ) ) return( unstained.cache[[ uf ]] )
    uf.path <- file.path( control.dir, uf )
    mat <- if ( file.exists( uf.path ) ) readFCS( uf.path ) else NULL
    unstained.cache[[ uf ]] <- mat
    mat
  }

  brightness.vals <- vapply( fluorophores, function( fluor ) {

    if ( verbose ) message( "Measuring brightness for ", fluor )

    tryCatch(
      {
        channel <- colnames( spectra )[ which.max( spectra[ fluor, ] ) ]

        stained <- readFCS( file.path( control.dir, fluor.files[[ fluor ]] ) )

        if ( !channel %in% colnames( stained ) )
          stop( "Channel '", channel, "' not found in ", fluor.files[[ fluor ]] )

        expr  <- stained[ , channel ]
        top.n <- min( n.cells, length( expr ) )

        if ( top.n < n.cells )
          warning(
            "Fewer than ", n.cells, " events in '", fluor, "'; using all ",
            top.n, " events.", call. = FALSE
          )

        pos.mfi <- mean( sort( expr, decreasing = TRUE )[ seq_len( top.n ) ] )

        src <- unstained.sources[[ fluor ]]

        neg.mfi <- if ( src$type == "file" ) {
          uf.mat <- read.unstained( src$file )
          if ( is.null( uf.mat ) || !channel %in% colnames( uf.mat ) ) {
            warning(
              "Unstained file unavailable for '", fluor,
              "'; using 25th percentile of the stained file instead.",
              call. = FALSE
            )
            stats::quantile( expr, 0.25 )
          } else {
            mean( uf.mat[ , channel ] )
          }
        } else {
          stats::quantile( expr, 0.25 )
        }

        pos.mfi - neg.mfi
      },
      error = function( e ) {
        message(
          sprintf( "[ERROR] Fluorophore %s failed: %s", fluor, conditionMessage( e ) )
        )
        NA_real_
      }
    )

  }, numeric( 1 ) )

  matrix(
    brightness.vals, ncol = 1,
    dimnames = list( fluorophores, "MFI" )
  )
}
