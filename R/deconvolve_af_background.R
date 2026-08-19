# deconvolve_af_background.R

#' @title Get Autofluorescence Basis
#'
#' @description
#' Builds a low-dimensional autofluorescence basis from an unstained control as
#' the leading right singular vectors of the raw, uncentred unstained matrix.
#' The leading component is therefore the mean background direction and the
#' remainder describe how background varies from cell to cell.
#'
#' A single autofluorescence spectrum cannot represent that variation. In a real
#' sample autofluorescence is a family of profiles that tracks cell type, and
#' cell type also determines which markers are expressed, so whatever a single
#' row leaves behind is correlated with marker expression and is readily
#' mistaken for spillover by any method that fits one channel against another.
#'
#' @param unstained Numeric matrix (events x detectors), raw unstained control,
#'   already gated to the population of interest.
#' @param n.pc Integer or `"auto"`. Number of components to retain. `"auto"`
#'   keeps every component whose singular value exceeds the largest value a
#'   matrix of the same shape containing only read noise would produce.
#'   Default `"auto"`.
#' @param max.pc Integer, hard cap on the number of components. Default `6`.
#' @param read.var Numeric, per-detector read noise variance in raw units, the
#'   floor of the retention threshold when `n.pc = "auto"`. Default `125^2`.
#' @param n.permutations Integer, permutations used to build the retention
#'   threshold when `n.pc = "auto"`. Each permutation costs one decomposition.
#'   Default `3`.
#' @param max.events Integer, maximum events used for the decomposition.
#'   Default `50000`.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#'
#' @return Numeric matrix, components x detectors, with unit L2 rows and the
#'   singular values attached as the `singular.values` attribute.
#'
#' @importFrom stats setNames median
#'
#' @export

get.af.basis <- function(
    unstained,
    n.pc       = "auto",
    max.pc         = 6L,
    read.var       = 125^2,
    n.permutations = 3L,
    max.events     = 50000L,
    verbose        = TRUE
) {

  unstained <- as.matrix( unstained )

  if ( nrow( unstained ) < 3 )
    stop( "`unstained` must contain at least three events.", call. = FALSE )

  fit.data <- unstained
  if ( nrow( fit.data ) > max.events )
    fit.data <- fit.data[ sample.int( nrow( fit.data ), max.events ), ,
                          drop = FALSE ]

  max.pc <- min( as.integer( max.pc ), ncol( fit.data ), nrow( fit.data ) - 1L )

  sv <- svd( fit.data, nu = 0L, nv = max.pc )

  if ( identical( n.pc, "auto" ) ) {

    # The threshold comes from the same matrix with each detector column
    # independently permuted. That destroys every cross-detector correlation
    # while preserving each detector's own distribution, so the largest
    # singular value the permuted matrix produces past the mean direction is
    # the largest this data could produce carrying no shared structure at all.
    # It needs no noise model, which matters because detector noise includes
    # shot noise from the background and is not one number across the array.
    # The electronic read noise is kept only as a floor.
    n.row <- nrow( fit.data )
    n.col <- ncol( fit.data )

    permuted.edge <- vapply( seq_len( as.integer( n.permutations ) ),
                             function( i )
                               svd( apply( fit.data, 2, sample ),
                                    nu = 0L, nv = 0L )$d[ 2 ],
                             numeric( 1 ) )

    noise.edge <- max( permuted.edge,
                       sqrt( read.var ) * ( sqrt( n.row ) + sqrt( n.col ) ) )

    n.keep <- sum( sv$d[ seq_len( max.pc ) ] > noise.edge )
    n.keep <- max( 1L, min( n.keep, max.pc ) )

  } else {

    n.keep <- max( 1L, min( as.integer( n.pc ), max.pc ) )
  }

  basis <- t( sv$v[ , seq_len( n.keep ), drop = FALSE ] )
  basis <- basis / sqrt( rowSums( basis^2 ) )

  # sign convention: the largest-magnitude detector of each component reads
  # positive, so the leading component prints as a background spectrum
  flip <- apply( basis, 1, function( v ) v[ which.max( abs( v ) ) ] < 0 )
  if ( any( flip ) )
    basis[ flip, ] <- -basis[ flip, , drop = FALSE ]

  dimnames( basis ) <- list( paste0( "AFPC", seq_len( n.keep ) ),
                             colnames( unstained ) )

  attr( basis, "singular.values" ) <- sv$d[ seq_len( max.pc ) ]

  if ( verbose )
    message( sprintf(
      "\033[34mAutofluorescence basis: %d of %d component(s) retained.\033[0m",
      n.keep, max.pc ) )

  basis
}


#' @title Deconvolve Autofluorescence Background
#'
#' @description
#' Removes autofluorescence from raw detector data by fitting it jointly with
#' the panel rather than subtracting it beforehand. Each event is unmixed
#' against a design of the autofluorescence basis stacked on the fluorophore
#' spectra and only the fitted autofluorescence part is subtracted, so a
#' fluorophore's abundance and the background amount are estimated together
#' instead of competing for the same signal.
#'
#' Fluorophores whose spectra lie close to the span of the autofluorescence
#' basis cannot be separated from it, and projecting them out removes real
#' signal. The returned `hotspot` matrix reports that coupling. When `target` is
#' supplied, components coupling to that fluorophore above `max.hotspot` are
#' dropped from its fit; the leading component is always retained, because every
#' event carries mean background.
#'
#' @param raw.data Numeric matrix (events x detectors), raw detector data.
#' @param spectra Numeric matrix (fluorophores x detectors), reference spectra.
#'   Any autofluorescence row named by `af.name` is removed before fitting, since
#'   the basis supersedes it.
#' @param af.basis Numeric matrix (components x detectors) from
#'   `get.af.basis()`.
#' @param af.name Character or `NULL`, the name of an autofluorescence row in
#'   `spectra` to drop. Default `"AF"`.
#' @param target Character or `NULL`. When supplied, autofluorescence components
#'   that this fluorophore cannot be separated from are dropped before fitting.
#'   Default `NULL`.
#' @param max.hotspot Numeric, the hotspot scale above which a component is
#'   considered inseparable from `target`. Default `5`.
#' @param nonneg Logical, whether to clamp fitted autofluorescence coefficients
#'   at zero before subtraction, so background can only be removed and never
#'   added. Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`residual`}{`raw.data` with the fitted autofluorescence removed.}
#'   \item{`background`}{The fitted autofluorescence contribution.}
#'   \item{`abundance`}{Fluorophore abundances from the joint fit.}
#'   \item{`af.coef`}{Per-event autofluorescence component coefficients, clamped
#'     at zero when `nonneg` is `TRUE`.}
#'   \item{`af.coef.raw`}{The same coefficients before clamping, for use as
#'     regressors.}
#'   \item{`hotspot`}{Hotspot matrix of the joint design.}
#'   \item{`af.used`, `af.dropped`}{Component names kept and dropped.}
#' }
#'
#' @export

deconvolve.af.background <- function(
    raw.data,
    spectra,
    af.basis,
    af.name     = "AF",
    target      = NULL,
    max.hotspot = 5,
    nonneg      = TRUE
) {

  raw.data <- as.matrix( raw.data )
  spectra  <- as.matrix( spectra )
  af.basis <- as.matrix( af.basis )

  if ( !is.null( af.name ) )
    spectra <- spectra[ setdiff( rownames( spectra ), af.name ), , drop = FALSE ]

  if ( ncol( raw.data ) != ncol( spectra ) ||
       ncol( raw.data ) != ncol( af.basis ) )
    stop( "`raw.data`, `spectra` and `af.basis` must have the same detectors.",
          call. = FALSE )

  if ( nrow( af.basis ) + nrow( spectra ) > ncol( spectra ) )
    stop( paste0( "The autofluorescence basis plus the panel exceeds the ",
                  "detector count; reduce `max.pc` in `get.af.basis()`." ),
          call. = FALSE )

  if ( !is.null( target ) && ! target %in% rownames( spectra ) )
    stop( "`target` is not a row of `spectra`.", call. = FALSE )

  keep.pc <- rownames( af.basis )

  hotspot <- calculate.hotspot.matrix(
    rbind( af.basis[ keep.pc, , drop = FALSE ], spectra ) )

  if ( !is.null( target ) ) {

    while ( length( keep.pc ) > 1L &&
            max( hotspot[ target, keep.pc ] ) > max.hotspot ) {

      minor   <- keep.pc[ -1L ]
      keep.pc <- setdiff( keep.pc,
                          minor[ which.max( hotspot[ target, minor ] ) ] )

      hotspot <- calculate.hotspot.matrix(
        rbind( af.basis[ keep.pc, , drop = FALSE ], spectra ) )
    }
  }

  design <- rbind( af.basis[ keep.pc, , drop = FALSE ], spectra )

  coefs <- unmix.ols.fast( raw.data, design )
  colnames( coefs ) <- rownames( design )

  af.coef.raw <- coefs[ , keep.pc, drop = FALSE ]

  af.coef <- af.coef.raw
  if ( nonneg ) af.coef[ af.coef < 0 ] <- 0

  background <- af.coef %*% af.basis[ keep.pc, , drop = FALSE ]

  list(
    residual    = raw.data - background,
    background  = background,
    abundance   = coefs[ , rownames( spectra ), drop = FALSE ],
    af.coef     = af.coef,
    af.coef.raw = af.coef.raw,
    hotspot     = hotspot,
    af.used     = keep.pc,
    af.dropped  = setdiff( rownames( af.basis ), keep.pc )
  )
}
