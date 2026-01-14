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
#' @importFrom stats predict sd
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

fit.af.spline <- function( af.cells, non.af.cells, asp ) {

  # set up boundaries structure for two gates
  af.boundaries <- list(
    upper = NULL,
    lower = NULL
  )

  # set lower bound based on non-AF cells
  x.bound.low <- quantile( non.af.cells[ , 1 ], asp$af.density.threshold )
  y.bound.low <- quantile( non.af.cells[ , 2 ], asp$af.density.threshold )

  # fit a spline using rlm
  model.data <- data.frame( rbind( af.cells, non.af.cells ) )
  colnames( model.data ) <- c( "x", "y" )

  if ( nrow( model.data ) < 10 )
    stop( "Failed to identify autofluorescence" )

  rlm.fit <- rlm( y ~ x, data = model.data, maxit = asp$af.spline.maxit )

  if ( !rlm.fit$converged )
    warning( "The IRLS algorithm employed in 'rlm' did not converge." )

  # define events within n standard deviations of the spline
  predicted <- predict( rlm.fit, newdata = model.data )
  model.data$predicted <- predicted
  model.data$residuals <- model.data$y - predicted

  sd.residuals <- sd( model.data$residuals )

  model.fit <- model.data[ abs(
    model.data$residuals ) <= asp$af.spline.sd.n * sd.residuals, ]

  model.fit.data <- model.fit[ which( model.fit$x > x.bound.low ), ]

  iter <- 0

  while ( nrow( model.fit.data ) < 10  & iter < 10 ) {
    x.bound.low <- x.bound.low * 0.9
    model.fit.data <- model.fit.data[ which( model.fit$x > x.bound.low ), ]
    iter <- iter + 1
  }

  if ( nrow( model.fit.data ) < 5 )
    stop( "Failed to identify autofluorescence" )

  # expand upper points outwards to catch more highly autofluorescent events
  x.q75 <- quantile( model.fit.data$x, 0.75 )
  y.q75 <- quantile( model.fit.data$y, 0.75 )

  upper.points <- model.fit.data[
    model.fit.data$x >= x.q75 | model.fit.data$y >= y.q75,
  ]

  upper.points$x <- upper.points$x * asp$af.spline.expand
  upper.points$y <- upper.points$y * asp$af.spline.expand

  # add to existing model data
  expanded.points <- rbind(
    model.fit.data,
    upper.points
  )

  expanded.points <- unique( expanded.points )

  # get the boundary of those events
  af.remove.boundary <- tripack::convex.hull(
    tripack::tri.mesh(
      expanded.points$x, expanded.points$y
    ) )

  if ( length( af.remove.boundary ) == 1 )
    stop( "Failed to identify autofluorescence" )
  else
    af.boundaries$upper <- af.remove.boundary

  return( af.boundaries )
}
