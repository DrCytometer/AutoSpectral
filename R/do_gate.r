# do_gate.r

#' @title Do Gate
#'
#' @description
#' Performs gating on scatter parameters and returns a vector with the indexes
#' of events inside the initial gate.
#'
#' The gating proceeds in three steps:
#'
#' - Defines bounds by data trimming
#' - Defines a region around the target maximum found within the bounds
#' - Defines a gate around the target maximum, only within that region
#'
#' The method uses numerical search of maxima over estimated densities and
#' Voronoi tessellations to improve density estimation around maxima.
#'
#' @importFrom MASS kde2d bandwidth.nrd
#' @importFrom deldir deldir tile.list which.tile
#' @importFrom sp point.in.polygon
#' @importFrom tripack tri.mesh convex.hull
#' @importFrom fields interp.surface
#'
#' @param gate.data A data frame containing the gate data.
#' @param viability.gate A logical vector indicating the viability gate.
#' @param large.gate A logical vector indicating the large gate.
#' @param samp A sample identifier.
#' @param scatter.and.channel.label A label for scatter and channel.
#' @param control.type The type of control used, either "beads" or "cells".
#' @param asp The AutoSpectral parameter list, prepared using
#' `get.autospectral.param`.
#'
#' @return A vector with the indexes of events inside the initial gate.
#'
#' @export


do.gate <- function( gate.data, viability.gate, large.gate,
                     samp, scatter.and.channel.label, control.type, asp )
{
  # set parameters for beads or cells
  if ( control.type == "beads" ) {
    default.gate.param <- asp$default.gate.param.beads
    gate.data.trim.factor.x.min <- asp$gate.data.trim.factor.x.min.beads
    gate.data.trim.factor.x.max <- asp$gate.data.trim.factor.x.max.beads
    gate.data.trim.factor.y.min <- asp$gate.data.trim.factor.y.min.beads
    gate.data.trim.factor.y.max <- asp$gate.data.trim.factor.y.max.beads

    gate.bound.density.bw.factor <- asp$gate.bound.density.bw.factor.beads
    gate.bound.density.grid.n <- asp$gate.bound.density.grid.n.beads
    gate.bound.density.neigh.size <- asp$gate.bound.density.neigh.size.beads

    gate.bound.density.max.target <- asp$gate.bound.density.max.target.beads
    gate.bound.density.max.exclusion.x <- asp$gate.bound.density.max.exclusion.x.beads
    gate.bound.density.max.exclusion.y <- asp$gate.bound.density.max.exclusion.y.beads
    gate.bound.density.max.mad.factor <- asp$gate.bound.density.max.mad.factor.beads

    gate.region.density.bw.factor <- asp$gate.region.density.bw.factor.beads
    gate.region.density.grid.n <- asp$gate.region.density.grid.n.beads
    gate.region.density.neigh.size <- asp$gate.region.density.neigh.size.beads

    gate.region.max.density.bw.factor <- asp$gate.region.max.density.bw.factor.beads
    gate.region.max.density.grid.n <- asp$gate.region.max.density.grid.n.beads
    gate.downsample.n <- asp$gate.downsample.n.beads

  } else {
    default.gate.param <- asp$default.gate.param.cells
    gate.data.trim.factor.x.min <- asp$gate.data.trim.factor.x.min.cells
    gate.data.trim.factor.x.max <- asp$gate.data.trim.factor.x.max.cells
    gate.data.trim.factor.y.min <- asp$gate.data.trim.factor.y.min.cells
    gate.data.trim.factor.y.max <- asp$gate.data.trim.factor.y.max.cells

    gate.bound.density.bw.factor <- asp$gate.bound.density.bw.factor.cells
    gate.bound.density.grid.n <- asp$gate.bound.density.grid.n.cells
    gate.bound.density.neigh.size <- asp$gate.bound.density.neigh.size.cells

    gate.bound.density.max.target <- asp$gate.bound.density.max.target.cells
    gate.bound.density.max.exclusion.x <- asp$gate.bound.density.max.exclusion.x.cells
    gate.bound.density.max.exclusion.y <- asp$gate.bound.density.max.exclusion.y.cells
    gate.bound.density.max.mad.factor <- asp$gate.bound.density.max.mad.factor.cells

    gate.region.density.bw.factor <- asp$gate.region.density.bw.factor.cells
    gate.region.density.grid.n <- asp$gate.region.density.grid.n.cells
    gate.region.density.neigh.size <- asp$gate.region.density.neigh.size.cells

    gate.region.max.density.bw.factor <- asp$gate.region.max.density.bw.factor.cells
    gate.region.max.density.grid.n <- asp$gate.region.max.density.grid.n.cells
    gate.downsample.n <- asp$gate.downsample.n.cells
  }

  gate.marker <- colnames( gate.data )

  gate.bound <- NULL
  gate.region <- NULL
  gate.boundary <- NULL

  # trim data
  gate.data.x.min <- max( asp$scatter.data.min.x, min( gate.data[ , 1 ] ) )
  gate.data.x.max <- min( asp$scatter.data.max.x, max( gate.data[ , 1 ] ) )

  gate.data.y.min <- max( asp$scatter.data.min.y, min( gate.data[ , 2 ] ) )
  gate.data.y.max <- min( asp$scatter.data.max.y, max( gate.data[ , 2 ] ) )

  gate.data.trim.x.min <-
    ( 1 - gate.data.trim.factor.x.min ) * gate.data.x.min +
    gate.data.trim.factor.x.min * gate.data.x.max
  gate.data.trim.x.max <-
    ( 1 - gate.data.trim.factor.x.max ) * gate.data.x.min +
    gate.data.trim.factor.x.max * gate.data.x.max

  gate.data.trim.y.min <-
    ( 1 - gate.data.trim.factor.y.min ) * gate.data.y.min +
    gate.data.trim.factor.y.min * gate.data.y.max
  gate.data.trim.y.max <-
    ( 1 - gate.data.trim.factor.y.max ) * gate.data.y.min +
    gate.data.trim.factor.y.max * gate.data.y.max

  # set bound from trimmed data
  gate.bound.x.low <- gate.data.trim.x.min
  gate.bound.x.high <- gate.data.trim.x.max

  gate.bound.y.low <- gate.data.trim.y.min
  gate.bound.y.high <- gate.data.trim.y.max

  gate.bound.data.idx <- which(
    gate.data[ , 1 ] > gate.bound.x.low &
      gate.data[ , 1 ] < gate.bound.x.high &
      gate.data[ , 2 ] > gate.bound.y.low &
      gate.data[ , 2 ] < gate.bound.y.high )

  bw <- apply( gate.data[ gate.bound.data.idx, ], 2, bandwidth.nrd )
  gate.bound.density <- MASS::kde2d(
    gate.data[ gate.bound.data.idx, 1 ],
    gate.data[ gate.bound.data.idx, 2 ],
    h = gate.bound.density.bw.factor * bw,
    n = gate.bound.density.grid.n )

  # get density maxima in bound
  gate.bound.neighbor.idx <- list(
    x = - gate.bound.density.neigh.size :
      gate.bound.density.neigh.size,
    y = - gate.bound.density.neigh.size :
      gate.bound.density.neigh.size )

  gate.bound.density.max.bool <- matrix( FALSE,
                                         nrow = gate.bound.density.grid.n,
                                         ncol = gate.bound.density.grid.n )

  for ( x.idx in 1 : gate.bound.density.grid.n )
    for ( y.idx in 1 : gate.bound.density.grid.n )
      gate.bound.density.max.bool[ x.idx, y.idx ] <-
    gate.bound.density$z[ x.idx, y.idx ] >=
    max( gate.bound.density$z[
      pmax( 0, pmin( gate.bound.density.grid.n,
                     x.idx + gate.bound.neighbor.idx$x ) ),
      pmax( 0, pmin( gate.bound.density.grid.n,
                     y.idx + gate.bound.neighbor.idx$y ) ) ] )

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
  gate.bound.density.max$num.label <- paste0( " ",
                                              row.names( gate.bound.density.max ) )

  # locate target maximum in bound, avoiding maxima located near the
  # bottom left corner
  gate.bound.density.max.offset <- 1

  while ( gate.bound.density.max.offset <= gate.bound.density.max.n )
  {
    if ( ( gate.bound.density.max$x[ gate.bound.density.max.offset ] -
           gate.bound.x.low ) /
         ( gate.bound.x.high - gate.bound.x.low ) >
         gate.bound.density.max.exclusion.x ||
         ( gate.bound.density.max$y[ gate.bound.density.max.offset ] -
           gate.bound.y.low ) /
         ( gate.bound.y.high - gate.bound.y.low ) >
         gate.bound.density.max.exclusion.y )
      break

    gate.bound.density.max.offset <- gate.bound.density.max.offset + 1
  }

  if ( gate.bound.density.max.offset > gate.bound.density.max.n ) {
    stop(
      paste( "gate error: no good maximum found in sample bound", samp ),
      call. = FALSE
    )
  }

  gate.bound.density.max.target <- gate.bound.density.max.target +
    gate.bound.density.max.offset - 1

  if ( gate.bound.density.max.target > gate.bound.density.max.n ) {
    stop(
      paste( "gate error: target maximum not found in sample bound", samp ),
      call. = FALSE
    )
  }

  if ( gate.bound.density.max.n > 1 )
  {
    # get voronoi tesselation for density maxima
    gate.bound.voronoi <- deldir( gate.bound.density.max,
                                  rw = c( gate.bound.x.low, gate.bound.x.high, gate.bound.y.low,
                                          gate.bound.y.high ), suppressMsge = TRUE )

    gate.bound.tile <- tile.list( gate.bound.voronoi )

    # get data in the tile of target maximum
    gate.bound.density.max.data.idx <- gate.bound.data.idx[
      sapply( gate.bound.data.idx, function( gbdi )
        which.tile( gate.data[ gbdi, 1 ], gate.data[ gbdi, 2 ],
                    gate.bound.tile ) == gate.bound.density.max.target )
    ]
  }
  else
  {
    gate.bound.voronoi <- NULL
    gate.bound.density.max.data.idx <- gate.bound.data.idx
  }

  if ( default.gate.param$region.auto )
  {
    # set region from target maximum found in bound
    gate.bound.density.max.x.median <- median( gate.data[
      gate.bound.density.max.data.idx, 1 ] )
    gate.bound.density.max.x.mad <- mad(
      gate.data[ gate.bound.density.max.data.idx, 1 ],
      center = gate.bound.density.max.x.median )

    gate.bound.density.max.y.median <- median( gate.data[
      gate.bound.density.max.data.idx, 2 ] )
    gate.bound.density.max.y.mad <- mad(
      gate.data[ gate.bound.density.max.data.idx, 2 ],
      center = gate.bound.density.max.y.median )

    gate.region.x.low <- max(
      gate.data.trim.x.min,
      gate.bound.density.max.x.median - gate.bound.density.max.mad.factor *
        gate.bound.density.max.x.mad
      )
    gate.region.x.high <- min(
      gate.data.trim.x.max,
      gate.bound.density.max.x.median + gate.bound.density.max.mad.factor *
        gate.bound.density.max.x.mad
      )

    gate.region.y.low <- max(
      gate.data.trim.y.min,
      gate.bound.density.max.y.median - gate.bound.density.max.mad.factor *
        gate.bound.density.max.y.mad
      )
    gate.region.y.high <- min(
      gate.data.trim.y.max,
      gate.bound.density.max.y.median + gate.bound.density.max.mad.factor *
        gate.bound.density.max.y.mad
      )

  }
  else
  {
    gate.region.x.low <-
      ( 1 - default.gate.param$region.factor.x.low ) * gate.bound.x.low +
      default.gate.param$region.factor.x.low * gate.bound.x.high
    gate.region.x.high <-
      ( 1 - default.gate.param$region.factor.x.high ) * gate.bound.x.low +
      default.gate.param$region.factor.x.high * gate.bound.x.high

    gate.region.y.low <-
      ( 1 - default.gate.param$region.factor.y.low ) * gate.bound.y.low +
      default.gate.param$region.factor.y.low * gate.bound.y.high
    gate.region.y.high <-
      ( 1 - default.gate.param$region.factor.y.high ) * gate.bound.y.low +
      default.gate.param$region.factor.y.high * gate.bound.y.high
  }

  gate.bound <- list(
    density = gate.bound.density,
    density.max = gate.bound.density.max,
    density.max.n = gate.bound.density.max.n,
    density.max.data.idx = gate.bound.density.max.data.idx,
    density.max.target = gate.bound.density.max.target,
    voronoi = gate.bound.voronoi,
    x.low = gate.bound.x.low,
    x.high = gate.bound.x.high,
    y.low = gate.bound.y.low,
    y.high = gate.bound.y.high
  )

  gate.region.data.idx <- which(
    gate.data[ , 1 ] > gate.region.x.low &
      gate.data[ , 1 ] < gate.region.x.high &
      gate.data[ , 2 ] > gate.region.y.low &
      gate.data[ , 2 ] < gate.region.y.high )

  if ( viability.gate ){

    gate.region.voronoi <- NULL

    # extend gate to the left by scaling factor (to include dead cells)
    gate.region.x.low <- gate.region.x.low * asp$viability.gate.scaling

    gate.region.data.idx <- which(
      gate.data[ , 1 ] > gate.region.x.low &
        gate.data[ , 1 ] < gate.region.x.high &
        gate.data[ , 2 ] > gate.region.y.low &
        gate.data[ , 2 ] < gate.region.y.high )

    gate.region.density.max.data.idx <- gate.region.data.idx

    gate.region <- list(
      data.idx = gate.region.data.idx,
      density = gate.bound.density,
      density.max = gate.bound.density.max,
      density.max.n = gate.bound.density.max.n,
      density.max.data.idx = gate.region.density.max.data.idx,
      voronoi = gate.region.voronoi,
      x.low = gate.region.x.low,
      x.high = gate.region.x.high,
      y.low = gate.region.y.low,
      y.high = gate.region.y.high
    )

    gate.population.boundary <- tripack::convex.hull(
      tripack::tri.mesh(
      gate.data[ gate.region.data.idx, 1 ],
      gate.data[ gate.region.data.idx, 2 ]
      ) )

    gate.population.pip <- sp::point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      gate.population.boundary$x, gate.population.boundary$y )

    gate.population.idx <- which( gate.population.pip != 0 )

    gate.population <- list( boundary = gate.population.boundary )

  } else {

    # get density maxima in region
    bw <- apply( gate.data[ gate.region.data.idx, ], 2, bandwidth.nrd )

    gate.region.density <- MASS::kde2d(
      gate.data[ gate.region.data.idx, 1 ],
      gate.data[ gate.region.data.idx, 2 ],
      h = gate.region.density.bw.factor * bw,
      n = gate.region.density.grid.n )

    gate.region.neighbor.idx <- list(
      x = - gate.region.density.neigh.size :
        gate.region.density.neigh.size,
      y = - gate.region.density.neigh.size :
        gate.region.density.neigh.size )

    gate.region.density.max.bool <- matrix(
      FALSE,
      nrow = gate.region.density.grid.n,
      ncol = gate.region.density.grid.n
      )

    for ( x.idx in 1 : gate.region.density.grid.n )
      for ( y.idx in 1 : gate.region.density.grid.n )
        gate.region.density.max.bool[ x.idx, y.idx ] <-
      gate.region.density$z[ x.idx, y.idx ] >=
      max( gate.region.density$z[
        pmax( 0, pmin( gate.region.density.grid.n,
                       x.idx + gate.region.neighbor.idx$x ) ),
        pmax( 0, pmin( gate.region.density.grid.n,
                       y.idx + gate.region.neighbor.idx$y ) ) ] )

    gate.region.density.max.idx <- which( gate.region.density.max.bool,
                                          arr.ind = TRUE )

    gate.region.density.max.n <- nrow( gate.region.density.max.idx )

    if ( gate.region.density.max.n < 1 ) {
      stop(
        paste( "gate error: no population found in sample region", samp ),
        call. = FALSE
      )
    }

    gate.region.density.max <- data.frame(
      x = gate.region.density$x[ gate.region.density.max.idx[ , 1 ] ],
      y = gate.region.density$y[ gate.region.density.max.idx[ , 2 ] ],
      z = gate.region.density$z[ gate.region.density.max.idx ] )

    gate.region.density.max <- gate.region.density.max[
      order( gate.region.density.max$z, decreasing = TRUE ), ]

    row.names( gate.region.density.max ) <- NULL
    gate.region.density.max$num.label <- paste0( " ",
                                                 row.names( gate.region.density.max ) )

    if ( gate.region.density.max.n > 1 )
    {
      # get voronoi tesselation for density maxima
      gate.region.voronoi <- deldir(
        gate.region.density.max,
        rw = c( gate.region.x.low, gate.region.x.high, gate.region.y.low,
                gate.region.y.high ), suppressMsge = TRUE
        )

      gate.region.tile <- tile.list( gate.region.voronoi )

      # get data in the tile of largest and second-largest maxima
      gate.region.density.max.data.idx <- gate.region.data.idx[
        sapply( gate.region.data.idx, function( grdi )
          which.tile( gate.data[ grdi, 1 ], gate.data[ grdi, 2 ],
                      gate.region.tile ) == 1 | 2 )
      ]
    }
    else
    {
      gate.region.voronoi <- NULL
      gate.region.density.max.data.idx <- gate.region.data.idx
    }

    gate.region <- list(
      data.idx = gate.region.data.idx,
      density = gate.region.density,
      density.max = gate.region.density.max,
      density.max.n = gate.region.density.max.n,
      density.max.data.idx = gate.region.density.max.data.idx,
      voronoi = gate.region.voronoi,
      x.low = gate.region.x.low,
      x.high = gate.region.x.high,
      y.low = gate.region.y.low,
      y.high = gate.region.y.high
    )

    # threshold data in region around target maximum
    bw <- apply( gate.data[ gate.region.density.max.data.idx, ], 2, bandwidth.nrd )

    gate.region.max.density <- MASS::kde2d(
      gate.data[ gate.region.density.max.data.idx, 1 ],
      gate.data[ gate.region.density.max.data.idx, 2 ],
      h = gate.region.max.density.bw.factor * bw,
      n = gate.region.max.density.grid.n )

    gate.region.max.density.interp <- fields::interp.surface(
      gate.region.max.density,
      gate.data[ gate.region.density.max.data.idx, ]
      )

    gate.region.max.density.threshold <-
      ( 1 - default.gate.param$density.threshold ) *
      min( gate.region.max.density.interp ) +
      default.gate.param$density.threshold * max( gate.region.max.density.interp )

    gate.population.strict.idx <- gate.region.density.max.data.idx[
      gate.region.max.density.interp > gate.region.max.density.threshold ]

    gate.population.strict.idx <- gate.population.strict.idx[
      ! duplicated( data.frame( gate.data[ gate.population.strict.idx, ] ) ) ]

    if ( large.gate ) {

      original.hull <- tripack::convex.hull(
        tripack::tri.mesh(
        gate.data[ gate.population.strict.idx, 1 ],
        gate.data[ gate.population.strict.idx, 2 ]
      ) )

      x.hull <- original.hull$x
      y.hull <- original.hull$y

      threshold.x <- quantile( x.hull, asp$large.gate.quantile )
      threshold.y <- quantile( y.hull, asp$large.gate.quantile )

      upper.idx.x <- which( x.hull > threshold.x )
      upper.idx.y <- which( y.hull > threshold.y & x.hull > threshold.x )

      x.hull.expanded <- x.hull
      y.hull.expanded <- y.hull
      x.hull.expanded[ upper.idx.x ] <- x.hull[ upper.idx.x ] +
        ( x.hull[ upper.idx.x ] - threshold.x ) *
        ( asp$large.gate.scaling.x - 1 )
      y.hull.expanded[ upper.idx.y ] <- y.hull[ upper.idx.y ] +
        ( y.hull[ upper.idx.y ] - threshold.y ) *
        ( asp$large.gate.scaling.y - 1 )

      gate.population.boundary <- list(
        x = x.hull.expanded,
        y = y.hull.expanded,
        i = original.hull$i
      )

    } else {

      gate.population.boundary <- tripack::convex.hull(
        tripack::tri.mesh(
        gate.data[ gate.population.strict.idx, 1 ],
        gate.data[ gate.population.strict.idx, 2 ] ) )

    }

    gate.population.pip <- sp::point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      gate.population.boundary$x, gate.population.boundary$y )

    gate.population.idx <- which( gate.population.pip != 0 )

    gate.population <- list( boundary = gate.population.boundary )

  }

  # prevent warnings about some gate.bound$density.max values being offscale
  suppressWarnings(
    gate.define.plot(
      samp,
      gate.data,
      gate.marker,
      gate.bound,
      gate.region,
      gate.population,
      scatter.and.channel.label,
      asp )
  )

  gate.population.boundary
}
