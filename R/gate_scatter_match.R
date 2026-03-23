# gate_scatter_match.r

#' @title Gate Matching Scatter
#' 
#' @description
#' Calculates the density of the input scatter coordinates, cuts at the defined
#' density percentile and calculates a convex hull boundary around the points
#' within the dense region. Permits selection of events with similar forward
#' and side scatter profiles in a paired sample, e.g., an unstained sample.
#' 
#' @importFrom MASS kde2d bandwidth.nrd
#' @importFrom tripack tri.mesh convex.hull
#' 
#' @param scatter.data Matrix of forward and side scatter data
#' @param percentile Numeric 0-1. Percentile of density to use for defining the 
#' contour for the gate boundary. Default is `0.5`.
#' @param bw.factor Numeric, scalar bandwidth multiplier. Larger numbers smooth
#' the density, expanding the gate outwards. Default is `1`.
#' @param grid.n Numeric, size of the grid for defining the kernel density
#' estimate. Default is `10`.
#' 
#' @return A trimesh object gate boundary

gate.scatter.match <- function(
    scatter.data,
    percentile = 0.5,
    bw.factor = 1,
    grid.n = 10
) {
  # calculate bandwidth for kernel density estimation
  bw <- apply( scatter.data, 2, bandwidth.nrd )
  
  # estimate density of scatter data
  dens <- MASS::kde2d(
    scatter.data[, 1 ],
    scatter.data[, 2 ],
    h = bw * bw.factor,
    n = grid.n
  )
  
  # find the density threshold
  z.sort <- sort( dens$z, decreasing = TRUE )
  cumulative.dens <- cumsum( z.sort ) / sum( z.sort )
  threshold <- z.sort[ which.min( abs( cumulative.dens - percentile ) ) ]
  
  # extract the contour line at that threshold
  contour.lines <- grDevices::contourLines( dens$x, dens$y, dens$z, levels = threshold )
  
  # take the longest contour (in case there are multiple)
  main.contour <- contour.lines[[ which.max( sapply( contour.lines, function(l) length(l$x) ) ) ]]
  gate.coords <- cbind( main.contour$x, main.contour$y )
  
  # define a smooth convex hull around those coordinates
  gate.hull <- unique( gate.coords[ grDevices::chull( gate.coords ), ] )
  scatter.gate <- tripack::convex.hull(
    tripack::tri.mesh( gate.hull[, 1], gate.hull[, 2] )
  )
  
  return( scatter.gate )
}