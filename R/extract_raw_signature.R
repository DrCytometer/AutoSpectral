# extract_raw_signature.R

#' @title Extract Raw Signature
#'
#' @description
#' Estimates one fluorophore's spectral signature directly from background
#' subtracted raw detector data, in the same way a signature is extracted from a
#' single stained control: bin the population by that fluorophore's abundance,
#' remove the contribution of every other fluorophore active in the population
#' at its current spectrum, and take the no-intercept slope of the remaining
#' detector signal against abundance.
#'
#' The estimate is a positively weighted combination of measured detector means,
#' so it is a spectrum in the sense that matters: it has a peak channel, it is
#' non-negative up to noise, and it can be compared with a control-derived row
#' by cosine similarity. This is the difference from back-solving a signature
#' out of a compensation matrix, which returns whichever linear combination of
#' the existing rows reproduces the required correction and need not resemble a
#' spectrum at all.
#'
#' Nothing is decided here. The function returns the candidate row together with
#' every quantity an acceptance gate needs, and the caller chooses whether to
#' adopt it.
#'
#' @param raw.data Numeric matrix (events x detectors), already restricted to
#'   the population for this fluorophore and already background subtracted.
#' @param spectra Numeric matrix (fluorophores x detectors), the current
#'   reference spectra. Rows other than `target` are held fixed.
#' @param abundance Numeric matrix (events x fluorophores), row-matched to
#'   `raw.data`, the current best abundance estimates for the same population.
#' @param target Character, the fluorophore whose row is being estimated.
#' @param active Character vector, `target` plus every fluorophore whose
#'   contribution is to be removed at its current spectrum. Defaults to the
#'   whole panel, which is the correct block-coordinate step: a fluorophore
#'   whose abundance in this population is noise adds variance but no bias,
#'   whereas one that is left in is projected onto the target's abundance and
#'   absorbed into its row.
#' @param extra.abundance Optional numeric matrix (events x components),
#'   row-matched to `abundance`, additional regressors fitted jointly with the
#'   panel, typically the unclamped autofluorescence component coefficients.
#'   Background variation that tracks the abundances is an omitted variable in
#'   the joint fit, and with collinear abundances its bias is amplified onto
#'   every partial slope; giving it its own columns removes the transfer.
#'   Requires `multivariate`. Default `NULL`.
#' @param intercept Logical, whether to fit an intercept. The intercept absorbs
#'   any constant offset the removal left behind, such as the mean level of
#'   other markers multiplied by their own row errors. Default `TRUE`.
#' @param multivariate Logical, whether every active abundance is fitted jointly
#'   rather than the nuisance fluorophores being subtracted at their current
#'   spectra beforehand. Subtracting first books any error in a nuisance row
#'   against the target in proportion to how closely the two abundances track
#'   each other within this population, which is what `vif.target` measures.
#'   The joint fit removes that transfer at the cost of variance, and requires
#'   more abundance bins than active fluorophores. Default `TRUE`.
#' @param ridge Numeric, ridge penalty applied to the scaled joint design,
#'   keeping the solve stable when two fluorophores are nearly collinear within
#'   the population. The intercept is never penalised. Default `1e-6`.
#' @param n.levels Integer, the maximum number of abundance bins. The out-of-span
#'   signal is carried by the covariance between the binned residual and the
#'   binned abundances, so the bin count is the fit's sample size and the joint
#'   design spends one point per regressor; the count actually used is the
#'   largest that keeps `min.bin.events` events per bin. Default `60`.
#' @param min.bin.events Integer, the fewest events an abundance bin may contain.
#'   A population that cannot supply more bins than regressors returns `NULL`.
#'   Default `50`.
#' @param min.events Integer, minimum events in the population. Default `200`.
#' @param background.raw Optional numeric matrix (events x detectors) of
#'   background events after the same subtraction, used to report how much
#'   background the subtraction failed to remove. Default `NULL`.
#'
#' @return `NULL` when the population cannot support a fit, otherwise a named
#'   list:
#' \describe{
#'   \item{`signature`}{The candidate row, non-negative and L-infinity
#'     normalised.}
#'   \item{`signature.raw`}{The unclamped, unnormalised slope.}
#'   \item{`stats`}{One-row data frame of fit and plausibility quantities:
#'     `x.span`, `explained`, `explained.total`, `resid.rel`, `intercept.rel`,
#'     `bg.align`, `clamp.frac`, `anchor.rel`, `vif.target`, `vif.max`,
#'     `vif.partner`, `deg.change`, `peak.curr`, `peak.new`, `peak.new.rel`.
#'     `bg.align` is the cosine between the fitted intercept and the fitted
#'     signature; a value near minus one means the fit is trading a
#'     background floor against the slope, the signature of an unmodelled
#'     background confound. `peak.new.rel` is the current signature's own
#'     value at the candidate's new peak channel, relative to the current
#'     signature's own peak - 1 when the candidate does not propose a shift,
#'     falling toward 0 the less that channel has to do with the dye today.}
#' }
#'
#' @importFrom MASS ginv
#' @importFrom stats cor quantile setNames
#'
#' @export

extract.raw.signature <- function(
    raw.data,
    spectra,
    abundance,
    target,
    active          = rownames( spectra ),
    extra.abundance = NULL,
    intercept       = TRUE,
    multivariate    = TRUE,
    ridge           = 1e-6,
    n.levels        = 60L,
    min.bin.events  = 50L,
    min.events      = 200L,
    background.raw  = NULL
) {

  raw.data  <- as.matrix( raw.data )
  spectra   <- as.matrix( spectra )
  abundance <- as.matrix( abundance )

  if ( ! target %in% rownames( spectra ) )
    stop( "`target` is not a row of `spectra`.", call. = FALSE )

  active <- c( target,
               setdiff( intersect( active, rownames( spectra ) ), target ) )

  if ( ! all( active %in% colnames( abundance ) ) )
    stop( "`abundance` must have a column for every active fluorophore.",
          call. = FALSE )

  n.extra <- if ( is.null( extra.abundance ) ) 0L else
    ncol( as.matrix( extra.abundance ) )

  if ( n.extra > 0L && !multivariate )
    stop( "`extra.abundance` requires `multivariate = TRUE`.", call. = FALSE )

  n.floor <- length( active ) + n.extra + 3L

  if ( multivariate && as.integer( n.levels ) < n.floor )
    stop( paste0( "`multivariate` needs more abundance bins than regressors; ",
                  "raise `n.levels` to at least ", n.floor, "." ),
          call. = FALSE )

  if ( nrow( raw.data ) < min.events ) return( NULL )
  if ( nrow( raw.data ) != nrow( abundance ) )
    stop( "`raw.data` and `abundance` must have the same events.",
          call. = FALSE )

  x <- abundance[ , active, drop = FALSE ]

  if ( n.extra > 0L ) {

    extra.abundance <- as.matrix( extra.abundance )

    if ( nrow( extra.abundance ) != nrow( abundance ) )
      stop( "`extra.abundance` must have the same events as `abundance`.",
            call. = FALSE )

    x <- cbind( x, extra.abundance )
  }

  x.target <- x[ , 1 ]

  # Everything this function can learn beyond the spectra it was handed lives in
  # the covariance between the binned residual and the binned abundances, so the
  # bin count is the sample size of the fit and the joint design spends one
  # point per regressor. Bins are therefore taken as fine as the population
  # supports, with the requested count as a ceiling; a population too small to
  # exceed the regressor count cannot identify the row at all.
  n.levels.use <- max( 3L, min(
    as.integer( n.levels ),
    floor( nrow( raw.data ) / as.integer( min.bin.events ) ) ) )

  # A population too small to give the joint design more bins than regressors
  # cannot identify the row jointly, but it can still support the univariate
  # fit, which needs three. Dropping to that is better than refusing the row,
  # and `joint` records which fit produced the answer.
  fit.joint <- multivariate && n.levels.use >= n.floor

  if ( n.extra > 0L && !fit.joint ) return( NULL )

  brk <- unique( stats::quantile(
    x.target, probs = seq( 0, 1, length.out = n.levels.use + 1 ),
    names = FALSE ) )
  if ( length( brk ) < 3 ) return( NULL )

  bin  <- as.integer( cut( x.target, breaks = brk, include.lowest = TRUE ) )
  bins <- sort( unique( bin ) )
  if ( length( bins ) < 3 ) return( NULL )

  idx.by.bin <- split( seq_len( nrow( raw.data ) ), bin )
  idx.by.bin <- idx.by.bin[ as.character( bins ) ]

  y.bin <- matrix(
    vapply( idx.by.bin, function( idx )
      colMeans( raw.data[ idx, , drop = FALSE ] ),
      numeric( ncol( raw.data ) ) ),
    nrow = length( bins ), byrow = TRUE,
    dimnames = list( NULL, colnames( raw.data ) ) )

  x.bin <- matrix(
    vapply( idx.by.bin, function( idx )
      colMeans( x[ idx, , drop = FALSE ] ), numeric( ncol( x ) ) ),
    nrow = length( bins ), byrow = TRUE,
    dimnames = list( NULL, colnames( x ) ) )

  nuisance <- active[ -1 ]
  xt       <- x.bin[ , 1 ]

  if ( max( xt ) - min( xt ) <= 0 ) return( NULL )

  if ( fit.joint && length( nuisance ) > 0 &&
       nrow( x.bin ) >= ncol( x.bin ) + 2L ) {

    # Every active abundance is fitted jointly rather than subtracting the
    # nuisance dyes at their current spectra beforehand. Subtracting first books
    # any error in a nuisance row against the target in proportion to how
    # closely the two abundances track each other within this population, which
    # is what `vif.target` measures; where that coupling is strong the target's
    # row absorbs its neighbours' errors as readily as its own signal. The joint
    # fit removes the transfer at the cost of variance, and the ridge term keeps
    # the solve stable when two dyes are nearly collinear here.
    design.bin <- if ( intercept ) cbind( 1, x.bin ) else x.bin

    scale.bin <- sqrt( colMeans( design.bin^2 ) )
    scale.bin[ !is.finite( scale.bin ) | scale.bin <= 0 ] <- 1

    z <- sweep( design.bin, 2, scale.bin, "/" )

    penalty <- diag( ridge, ncol( z ) )
    if ( intercept ) penalty[ 1, 1 ] <- 0

    coefficients <- tryCatch(
      solve( crossprod( z ) + penalty, crossprod( z, y.bin ) ),
      error = function( e ) NULL )

    if ( is.null( coefficients ) ) return( NULL )

    coefficients <- coefficients / scale.bin

    offset        <- if ( intercept ) coefficients[ 1, ] else
      rep( 0, ncol( y.bin ) )
    signature.raw <- coefficients[ if ( intercept ) 2L else 1L, ]
    fitted        <- design.bin %*% coefficients
    y.res         <- y.bin

  } else {

    y.res <- y.bin

    if ( length( nuisance ) > 0 )
      y.res <- y.bin - x.bin[ , nuisance, drop = FALSE ] %*%
        spectra[ nuisance, , drop = FALSE ]

    if ( intercept ) {

      fit <- stats::lm.fit( x = cbind( 1, xt ), y = y.res )

      offset        <- stats::coef( fit )[ 1, ]
      signature.raw <- stats::coef( fit )[ 2, ]
      fitted        <- cbind( 1, xt ) %*% stats::coef( fit )

    } else {

      den <- sum( xt^2 )
      if ( !is.finite( den ) || den <= 0 ) return( NULL )

      offset        <- rep( 0, ncol( y.res ) )
      signature.raw <- as.numeric( crossprod( xt, y.res ) ) / den
      fitted        <- outer( xt, signature.raw )
    }
  }

  offset[ !is.finite( offset ) ] <- 0
  signature.raw[ !is.finite( signature.raw ) ] <- 0
  names( signature.raw ) <- colnames( raw.data )

  resid.rel <- sqrt( sum( ( y.res - fitted )^2 ) ) /
    max( sqrt( sum( y.res^2 ) ), .Machine$double.eps )

  intercept.rel <- sqrt( sum( offset^2 ) ) /
    max( sqrt( sum( ( max( xt ) * signature.raw )^2 ) ), .Machine$double.eps )

  clamped    <- pmax( signature.raw, 0 )
  clamp.frac <- sum( pmax( -signature.raw, 0 ) ) /
    max( sum( abs( signature.raw ) ), .Machine$double.eps )

  if ( max( clamped ) <= 0 ) return( NULL )

  signature <- clamped / max( clamped )

  # `explained` is the share of the brightest bin's signal carried by the
  # fluorophore's own term. Under the joint fit that term is a partial slope, so
  # a collinear neighbour can take a share of it and the value drifts below one
  # without the fit being wrong. `explained.total` separates the two: it asks
  # whether the whole fit reproduces the bin, and stays near one whenever the
  # decomposition is sound whatever the collinearity.
  top.bin   <- which.max( xt )
  top.norm  <- sqrt( sum( y.bin[ top.bin, ]^2 ) )

  explained <- if ( top.norm > 0 )
    sqrt( sum( ( xt[ top.bin ] * signature.raw )^2 ) ) / top.norm else 0

  explained.total <- if ( top.norm > 0 )
    sqrt( sum( fitted[ top.bin, ]^2 ) ) / top.norm else 0

  # How much background the subtraction failed to remove, as a fraction of the
  # brightest signal being fitted. A large value means the slope is being
  # measured on top of a floor that moves, and the row should not be updated.
  anchor.rel <- NA_real_

  if ( !is.null( background.raw ) ) {
    background.raw <- as.matrix( background.raw )
    if ( nrow( background.raw ) > 0 ) {
      bg.mean    <- colMeans( background.raw )
      anchor.rel <- sqrt( sum( bg.mean^2 ) ) /
        max( sqrt( sum( ( max( xt ) * signature.raw )^2 ) ),
             .Machine$double.eps )
    }
  }

  vif <- .signature.vif( x )

  vif.partner <- NA_character_

  if ( ncol( x ) > 1 ) {
    r.target <- suppressWarnings(
      stats::cor( x[ , 1 ], x[ , -1, drop = FALSE ] ) )
    r.target[ !is.finite( r.target ) ] <- 0
    vif.partner <- colnames( x )[ -1 ][ which.max( abs( r.target ) ) ]
  }

  bg.align <- NA_real_

  if ( intercept ) {
    offset.norm <- sqrt( sum( offset^2 ) )
    slope.norm  <- sqrt( sum( signature.raw^2 ) )
    if ( offset.norm > 0 && slope.norm > 0 )
      bg.align <- sum( offset * signature.raw ) / ( offset.norm * slope.norm )
  }

  current  <- spectra[ target, ]
  cos.curr <- sum( current * signature ) /
    max( sqrt( sum( current^2 ) ) * sqrt( sum( signature^2 ) ),
         .Machine$double.eps )

  peak.new.rel <- current[ which.max( signature ) ] /
    max( max( current ), .Machine$double.eps )

  stats.out <- data.frame(
    fluorophore      = target,
    n.events         = nrow( raw.data ),
    n.bins           = length( bins ),
    n.active         = length( active ),
    joint            = fit.joint,
    x.span           = max( xt ) - min( xt ),
    explained        = explained,
    explained.total  = explained.total,
    resid.rel        = resid.rel,
    intercept.rel    = intercept.rel,
    bg.align         = bg.align,
    clamp.frac       = clamp.frac,
    anchor.rel       = anchor.rel,
    vif.target       = unname( vif[ 1 ] ),
    vif.max          = max( vif ),
    vif.partner      = vif.partner,
    deg.change       = 180 / pi * acos( pmin( 1, pmax( -1, cos.curr ) ) ),
    peak.curr        = colnames( raw.data )[ which.max( current ) ],
    peak.new         = colnames( raw.data )[ which.max( signature ) ],
    peak.new.rel     = unname( peak.new.rel ),
    row.names        = NULL,
    stringsAsFactors = FALSE
  )

  list( signature     = signature,
        signature.raw = signature.raw,
        stats         = stats.out )
}


#' Variance inflation factors for the columns of an abundance matrix. A
#' fluorophore whose abundance within a population is close to a linear function
#' of the others cannot have its signature separated from theirs there, whether
#' the cause is spectral collinearity or two markers being carried by the same
#' cells. Values near one are safe; large values mean the population carries no
#' information about that row.
#' @noRd
.signature.vif <- function( x ) {

  if ( ncol( x ) < 2 )
    return( stats::setNames( rep( 1, ncol( x ) ), colnames( x ) ) )

  r <- suppressWarnings( stats::cor( x ) )
  r[ !is.finite( r ) ] <- 0
  diag( r ) <- 1

  r.inv <- tryCatch( solve.default( r ),
                     error = function( e ) MASS::ginv( r ) )

  stats::setNames( pmax( diag( r.inv ), 1 ), colnames( x ) )
}

#' Split a spectral change into components inside and outside the row space of
#' a reference spectra matrix. The in-span part is what a compensation matrix
#' can represent with one coefficient per fluorophore; the out-of-span part is
#' colour present in no reference row, which only a detector-space estimate can
#' see.
#' @noRd
.signature.span.split <- function( delta, basis ) {

  basis <- as.matrix( basis )
  delta <- as.matrix( delta )

  projector <- t( basis ) %*% MASS::ginv( basis %*% t( basis ) ) %*% basis

  parallel <- delta %*% projector

  list( parallel = parallel, perpendicular = delta - parallel )
}
