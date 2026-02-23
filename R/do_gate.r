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

do.gate <- function(
    gate.data,
    viability.gate,
    large.gate,
    samp,
    scatter.and.channel.label,
    control.type,
    asp
) {

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

  # trim data to cytometer's limits
  gate.data.x.min <- max( asp$scatter.data.min.x, min( gate.data[ , 1 ] ) )
  gate.data.x.max <- min( asp$scatter.data.max.x, max( gate.data[ , 1 ] ) )
  gate.data.y.min <- max( asp$scatter.data.min.y, min( gate.data[ , 2 ] ) )
  gate.data.y.max <- min( asp$scatter.data.max.y, max( gate.data[ , 2 ] ) )

  # trim the boundary limits (usually by 5%)
  gate.bound.x.low <- (1 - gate.data.trim.factor.x.min) * gate.data.x.min + gate.data.trim.factor.x.min * gate.data.x.max
  gate.bound.x.high <- (1 - gate.data.trim.factor.x.max) * gate.data.x.min + gate.data.trim.factor.x.max * gate.data.x.max
  gate.bound.y.low <- (1 - gate.data.trim.factor.y.min) * gate.data.y.min + gate.data.trim.factor.y.min * gate.data.y.max
  gate.bound.y.high <- (1 - gate.data.trim.factor.y.max) * gate.data.y.min + gate.data.trim.factor.y.max * gate.data.y.max

  lims <- c(gate.data.x.min, gate.data.x.max, gate.data.y.min, gate.data.y.max)
  gate.bound.data.idx <- which(gate.data[,1] > lims[1] & gate.data[,1] < lims[2] &
                                 gate.data[,2] > lims[3] & gate.data[,2] < lims[4])

  bw <- apply(gate.data[gate.bound.data.idx, ], 2, bandwidth.nrd)
  master.density <- MASS::kde2d(
    gate.data[gate.bound.data.idx, 1],
    gate.data[gate.bound.data.idx, 2],
    h = gate.bound.density.bw.factor * bw,
    n = gate.bound.density.grid.n )

  # locate target maxima
  if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
       "find_local_maxima" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
    gate.bound.density.max.bool <- AutoSpectralRcpp::find_local_maxima(
      master.density$z,
      gate.bound.density.neigh.size
    )
  } else {
    # get density maxima in bound
    gate.bound.neighbor.idx <- list(
      x = - gate.bound.density.neigh.size : gate.bound.density.neigh.size,
      y = - gate.bound.density.neigh.size : gate.bound.density.neigh.size
    )

    gate.bound.density.max.bool <- matrix(
      FALSE,
      nrow = gate.bound.density.grid.n,
      ncol = gate.bound.density.grid.n
    )

    for ( x.idx in 1 : gate.bound.density.grid.n )
      for ( y.idx in 1 : gate.bound.density.grid.n )
        gate.bound.density.max.bool[ x.idx, y.idx ] <-
      master.density$z[ x.idx, y.idx ] >= max( master.density$z[
        pmax( 0, pmin( gate.bound.density.grid.n, x.idx + gate.bound.neighbor.idx$x ) ),
        pmax( 0, pmin( gate.bound.density.grid.n, y.idx + gate.bound.neighbor.idx$y ) ) ] )
  }

  gate.bound.density.max.idx <- which( gate.bound.density.max.bool, arr.ind = TRUE )
  gate.bound.density.max.n <- nrow( gate.bound.density.max.idx )

  if ( gate.bound.density.max.n < 1 ) {
    stop(
      paste0( "gate error: no population found in sample bound ", samp ),
      call. = FALSE
    )
  }

  gate.bound.density.max <- data.frame(
    x = master.density$x[gate.bound.density.max.idx[,1]],
    y = master.density$y[gate.bound.density.max.idx[,2]],
    z = master.density$z[gate.bound.density.max.idx]
  )
  gate.bound.density.max <- gate.bound.density.max[order(gate.bound.density.max$z, decreasing = TRUE), ]
  row.names( gate.bound.density.max ) <- NULL
  gate.bound.density.max$num.label <- paste0(
    " ", row.names( gate.bound.density.max )
  )

  # locate target maximum in bound, avoiding maxima located near the bottom left corner
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

  # define the region, using voronoi if multiple density peaks exist
  if ( gate.bound.density.max.n > 1) {
    gate.bound.voronoi <- deldir(gate.bound.density.max, rw = lims, suppressMsge = TRUE)
    gate.bound.tile <- tile.list(gate.bound.voronoi)
    gate.bound.density.max.data.idx <- gate.bound.data.idx[
      sapply(gate.bound.data.idx, function(gbdi)
        which.tile(gate.data[gbdi,1], gate.data[gbdi,2], gate.bound.tile) == gate.bound.density.max.target)
    ]
  } else {
    gate.bound.voronoi <- NULL
    gate.bound.density.max.data.idx <- gate.bound.data.idx
  }

  # calculate region limits (used for plotting the gate box)
  gate.region.x.median <- stats::median(gate.data[gate.bound.density.max.data.idx, 1])
  gate.region.x.mad <- stats::mad(gate.data[gate.bound.density.max.data.idx, 1], center = gate.region.x.median)
  gate.region.y.median <- stats::median(gate.data[gate.bound.density.max.data.idx, 2])
  gate.region.y.mad <- stats::mad(gate.data[gate.bound.density.max.data.idx, 2], center = gate.region.y.median)

  gate.region.x.low <- max(gate.data.x.min, gate.region.x.median - gate.bound.density.max.mad.factor * gate.region.x.mad)
  gate.region.x.high <- min(gate.data.x.max, gate.region.x.median + gate.bound.density.max.mad.factor * gate.region.x.mad)
  gate.region.y.low <- max(gate.data.y.min, gate.region.y.median - gate.bound.density.max.mad.factor * gate.region.y.mad)
  gate.region.y.high <- min(gate.data.y.max, gate.region.y.median + gate.bound.density.max.mad.factor * gate.region.y.mad)

  if ( viability.gate ) {
    # extend gate to the left by scaling factor (to include dead cells)
    gate.region.x.low <- gate.region.x.low * asp$viability.gate.scaling
  }

  spatial.mask <- gate.data[gate.bound.density.max.data.idx, 1] > gate.region.x.low &
    gate.data[gate.bound.density.max.data.idx, 1] < gate.region.x.high &
    gate.data[gate.bound.density.max.data.idx, 2] > gate.region.y.low &
    gate.data[gate.bound.density.max.data.idx, 2] < gate.region.y.high

  gate.region.data.idx <- gate.bound.density.max.data.idx[ spatial.mask ]

  # re-define density in the target region
  bw <- apply( gate.data[ gate.region.data.idx, ], 2, bandwidth.nrd )
  region.density <- MASS::kde2d(
    gate.data[gate.region.data.idx, 1],
    gate.data[gate.region.data.idx, 2],
    h = gate.region.density.bw.factor * bw,
    n = gate.region.density.grid.n )

  gate.region.max.density.interp <- fields::interp.surface(
    region.density,
    gate.data[gate.region.data.idx, ]
  )

  # Calculate threshold based on the interpolated values of the target cluster
  gate.region.max.density.threshold <-
    (1 - default.gate.param$density.threshold) * min(gate.region.max.density.interp) +
    (default.gate.param$density.threshold) * max(gate.region.max.density.interp)

  gate.population.strict.idx <- gate.region.data.idx[
    gate.region.max.density.interp > gate.region.max.density.threshold
  ]

  gate.population.strict.idx <- gate.population.strict.idx[
    ! duplicated( data.frame( gate.data[ gate.population.strict.idx, ] ) ) ]

  # extend the gate upwards and outwards for larger cells
  if ( large.gate ) {
    original.hull <- tripack::convex.hull(
      tripack::tri.mesh(
        gate.data[ gate.population.strict.idx, 1 ],
        gate.data[ gate.population.strict.idx, 2 ]
      ) )

    x.hull <- original.hull$x
    y.hull <- original.hull$y

    threshold.x <- stats::quantile( x.hull, asp$large.gate.quantile )
    threshold.y <- stats::quantile( y.hull, asp$large.gate.quantile )

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

  # format reginos for plotting
  gate.bound <- list(
    density = master.density,
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
  gate.region <- list(
    data.idx = gate.bound.density.max.data.idx,
    density = master.density,
    density.max = gate.bound.density.max,
    density.max.n = gate.bound.density.max.n,
    density.max.data.idx = gate.bound.density.max.data.idx,
    voronoi = NULL,
    x.low = gate.region.x.low,
    x.high = gate.region.x.high,
    y.low = gate.region.y.low,
    y.high = gate.region.y.high
  )

  gate.population <- list( boundary = gate.population.boundary )

  # plotting
  gate.define.plot(
    samp,
    gate.data,
    gate.marker,
    gate.bound,
    gate.region,
    gate.population,
    scatter.and.channel.label,
    asp )

  return( gate.population.boundary )
}
