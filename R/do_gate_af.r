# do_gate_af.r

#' @title Perform Gating on Autofluorescence Parameters
#'
#' @description
#' This function returns a vector with the indexes of events inside the initial
#' gate on autofluorescence parameters. It proceeds through several steps to
#' define the gate boundaries and identify density maxima using numerical
#' search and Voronoi tessellations.
#'
#' @importFrom deldir deldir tile.list
#' @importFrom sp point.in.polygon
#' @importFrom MASS kde2d bandwidth.nrd
#'
#' @param gate.data A data frame containing the gate data.
#' @param samp A sample identifier.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param intermediate.figures Logical, if `TRUE` returns additional figures to
#' show the inner workings of the cleaning, including definition of low-AF cell
#' gates on the PCA-unmixed unstained and spectral ribbon plots of the AF
#' exclusion from the unstained.
#'
#' @return A vector with the indexes of events inside the initial gate.
#'
#' @export

do.gate.af <- function( gate.data, samp, asp, intermediate.figures = FALSE ) {

  # change to restrict to lower 50%
  gate.data.trim.factor.x.min <- asp$gate.data.trim.factor.x.min.cells
  gate.data.trim.factor.x.max <- asp$gate.data.trim.factor.x.max.cells
  gate.data.trim.factor.y.min <- asp$gate.data.trim.factor.y.min.cells
  gate.data.trim.factor.y.max <- asp$gate.data.trim.factor.y.max.cells

  gate.data.x.min <- min( gate.data[ , 1 ] )
  gate.data.x.max <- max( gate.data[ , 1 ] )

  gate.data.y.min <- min( gate.data[ , 2 ] )
  gate.data.y.max <- max( gate.data[ , 2 ] )

  gate.data.x.min <- ( 1 - gate.data.trim.factor.x.min ) * gate.data.x.min +
    gate.data.trim.factor.x.min * gate.data.x.max
  gate.data.x.max <- ( 1 - gate.data.trim.factor.x.max ) * gate.data.x.min +
    gate.data.trim.factor.x.max * gate.data.x.max

  gate.data.y.min <- ( 1 - gate.data.trim.factor.y.min ) * gate.data.y.min +
    gate.data.trim.factor.y.min * gate.data.y.max
  gate.data.y.max <- ( 1 - gate.data.trim.factor.y.max ) * gate.data.y.min +
    gate.data.trim.factor.y.max * gate.data.y.max

  grid.n <- asp$af.gate.bound.density.grid.n
  bw <- apply( gate.data, 2, bandwidth.nrd )
  gate.bound.density <- MASS::kde2d(
    gate.data[ , 1 ], gate.data[ , 2 ],
    h = asp$af.gate.density.bw.factor * bw,
    n = grid.n )

  # Identify density maxima
  neighbor.n <- asp$af.gate.bound.density.neigh.size

  gate.bound.neighbor.idx <- list(
    x = - neighbor.n :neighbor.n,
    y = - neighbor.n :neighbor.n )

  gate.bound.density.max.bool <- matrix( FALSE, nrow = grid.n, ncol = grid.n )

  for ( x.idx in 1 : grid.n )
    for ( y.idx in 1 : grid.n )
      gate.bound.density.max.bool[ x.idx, y.idx ] <-
    gate.bound.density$z[ x.idx, y.idx ] >=
    max( gate.bound.density$z[
      pmax( 0, pmin( grid.n, x.idx + gate.bound.neighbor.idx$x ) ),
      pmax( 0, pmin( grid.n, y.idx + gate.bound.neighbor.idx$y ) ) ] )

  gate.bound.density.max.idx <- which( gate.bound.density.max.bool,
                                       arr.ind = TRUE )

  gate.bound.density.max.n <- nrow( gate.bound.density.max.idx )

  if ( gate.bound.density.max.n < 1 ) {
    stop(
      paste0( "gate error: no population found in sample bound ", samp ),
      call. = FALSE
    )
  }

  gate.bound.density.max <- data.frame(
    x = gate.bound.density$x[ gate.bound.density.max.idx[ , 1 ] ],
    y = gate.bound.density$y[ gate.bound.density.max.idx[ , 2 ] ],
    z = gate.bound.density$z[ gate.bound.density.max.idx ] )

  gate.bound.density.max <- gate.bound.density.max[
    order( gate.bound.density.max$z, decreasing = TRUE ), ]

  row.names( gate.bound.density.max ) <- NULL
  gate.bound.density.max$num.label <- paste0(
    " ", row.names( gate.bound.density.max ) )

  # Identify the target maximum
  target.max.idx <- asp$af.gate.target.max
  target.max <- gate.bound.density.max[ target.max.idx, ]

  gate.bound.voronoi <- deldir::deldir(
    gate.bound.density.max,
    rw = c( gate.data.x.min, gate.data.x.max,
            gate.data.y.min, gate.data.y.max ),
    suppressMsge = TRUE
    )

  tiles <- tile.list( gate.bound.voronoi )

  # Identify points in the tile of the global maximum
  target.tile <- tiles[[ target.max.idx ]]
  tile.polygon <- cbind( target.tile$x, target.tile$y )

  gate.region.data.idx <- which(
    sp::point.in.polygon(
    gate.data[, 1], gate.data[, 2],
    tile.polygon[ , 1 ], tile.polygon[ , 2 ]
    ) > 0 )

  if ( length( gate.region.data.idx ) == 0 ) {
    stop(
      paste( "Error: No points in target Voronoi tile for sample", samp ),
      call. = FALSE
    )
  }

  # define region boundaries
  gate.region.data <- gate.data[ gate.region.data.idx, ]

  gate.region.x.low <- min( gate.region.data[ , 1 ] )
  gate.region.x.high <- max( gate.region.data[ , 1 ]  )
  gate.region.y.low <- min( gate.region.data[ , 2 ] )
  gate.region.y.high <- max( gate.region.data[ , 2 ] )

  # define gate.region for plotting
  gate.region <- data.frame(
    x = c(
      gate.region.x.low,
      gate.region.x.high,
      gate.region.x.high,
      gate.region.x.low,
      gate.region.x.low
    ),
    y = c(
      gate.region.y.low,
      gate.region.y.low,
      gate.region.y.high,
      gate.region.y.high,
      gate.region.y.low )
  )

  if ( asp$figures & intermediate.figures )
    gate.af.identify.plot(
      gate.data,
      samp,
      gate.region,
      gate.bound.density,
      asp
      )

  gate.population.idx <- which(
    gate.data[ , 1 ] > gate.region.x.low &
      gate.data[ , 1 ] < gate.region.x.high &
      gate.data[ , 2 ] > gate.region.y.low &
      gate.data[ , 2 ] < gate.region.y.high
    )

  return( gate.population.idx )
}
