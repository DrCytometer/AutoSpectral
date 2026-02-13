# fit_af_spline.r

#' @title Fit Spline to Autofluorescence Data
#'
#' @description
#' This function fits a spline to autofluorescence data, removing extreme events
#' and defining bounds equally far from zero. It uses robust linear modeling
#' and identifies events within a specified number of standard deviations from
#' the spline.
#'
#' @importFrom MASS rlm
#' @importFrom tripack tri.mesh convex.hull
#'
#' @param af.cells A matrix containing the autofluorescence data.
#' @param non.af.cells A matrix containing the low autofluorescence data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return The boundary of the events within the specified number of standard
#' deviations from the spline.
#'
#' @export

fit.af.spline <- function(
    af.cells,
    non.af.cells,
    asp
) {

  # set up boundaries structure for two gates
  af.boundaries <- list(
    upper = NULL,
    lower = NULL
  )

  # failsafe boundary (triangle outside the minimum range--excludes no events)
  af.boundaries$upper <- list(
    x = c( asp$expr.data.min - 10, asp$expr.data.min - 5, asp$expr.data.min - 10 ),
    y = c( asp$expr.data.min - 10, asp$expr.data.min - 10, asp$expr.data.min - 5 )
  )

  if ( is.null( af.cells) || is.null( non.af.cells ) ||
      nrow( af.cells ) < 10 || nrow( non.af.cells ) < 10 ) {
    warning( "Insufficient cells for AF spline fitting. Returning default gate." )
    return( af.boundaries )
  }

  # fit a spline using rlm
  model.data <- data.frame( rbind( af.cells, non.af.cells ) )
  colnames( model.data ) <- c( "x", "y" )

  # fit a model around the autofluorescence data
  rlm.fit <- tryCatch( {
    MASS::rlm( y ~ x, data = model.data, maxit = asp$af.spline.maxit )
  }, error = function( e ) return( NULL ) )

  # check that it worked, fallback if not
  if ( is.null( rlm.fit ) ) {
    warning( "RLM fit failed for AF spline." )
    return( af.boundaries )
  }
  if ( !rlm.fit$converged ) {
    warning( "The IRLS algorithm employed in 'rlm' did not converge." )
    return( af.boundaries )
  }

  # define events within n standard deviations of the spline
  predicted <- stats::predict( rlm.fit, newdata = model.data )
  residuals <- model.data$y - predicted
  sd.res <- stats::sd( residuals, na.rm = TRUE )
  # protect against zero variance
  if (is.na( sd.res ) || sd.res == 0 ) sd.res <- 1e-6

  # filter to points within n standard deviations of the fit
  model.fit <- model.data[ abs( residuals ) <= asp$af.spline.sd.n * sd.res, ]

  # expand the pool of points if not well represented
  x.bound.low <- stats::quantile(
    non.af.cells[ , 1 ],
    asp$af.density.threshold,
    na.rm = TRUE )
  iter <- 0
  model.fit.data <- model.fit[ model.fit$x > x.bound.low, ]

  while ( nrow( model.fit.data ) < 10 && iter < 10 ) {
    # shift bound lower to include more points
    x.bound.low <- x.bound.low * 0.9
    model.fit.data <- model.fit[ model.fit$x > x.bound.low, ]
    iter <- iter + 1
  }

  # create a boundary gate around the AF signal if we have enough events
  if ( nrow( model.fit.data ) >= 5 ) {
    # expand to catch extreme AF
    q.vals <- apply(
      model.fit.data[ , 1:2 ],
      2,
      stats::quantile,
      probs = 0.75,
      na.rm = TRUE
    )
    upper.points <- model.fit.data[
      model.fit.data$x >= q.vals[ 1 ] | model.fit.data$y >= q.vals[ 2 ], ]

    upper.points$x <- upper.points$x * asp$af.spline.expand
    upper.points$y <- upper.points$y * asp$af.spline.expand

    expanded.points <- unique(
      rbind( model.fit.data[ , 1:2 ], upper.points[ , 1:2 ] )
    )

    # attempt a convex hull
    af.boundaries$upper <- tryCatch( {
      tripack::convex.hull( tripack::tri.mesh( expanded.points$x, expanded.points$y ) )
    }, error = function( e ) {
      warning( "AF hull geometry failed. Returning default gate." )
      return( af.boundaries )
    } )
  }

  return( af.boundaries )
}
