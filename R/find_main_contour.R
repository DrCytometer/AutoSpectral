# find_main_contour.r

## Extracts the single largest closed density contour at (or near) a requested
## percentile threshold from a 2D kernel density estimate (as returned by
## MASS::kde2d() or AutoSpectralRcpp::fast_kde2d_cpp()). If the requested
## percentile produces a degenerate result -- i.e., grDevices::contourLines()
## returns no contour at that exact density level, which happens when the
## input population is tightly clustered relative to the density grid
## resolution (e.g., highly monodisperse beads) -- the threshold is relaxed
## towards progressively lower density levels until a valid contour is found.
## Used by define.gate.landmarks(), gate.scatter.match() and tune.gate().
##
## @param dens A list with components `x`, `y`, `z` as returned by
## MASS::kde2d() (or a compatible grid-based 2D density estimate).
## @param percentile Numeric 0-1. The target density percentile for the
## contour, e.g. `0.7` keeps roughly the 70% of density mass closest to the
## peak.
## @return A list with components `x` and `y` describing the boundary of the
## largest contour at the (possibly relaxed) threshold.
## @keywords internal
.find.main.contour <- function( dens, percentile ) {

  z.sort <- sort( dens$z, decreasing = TRUE )

  if ( length( z.sort ) == 0 || !is.finite( sum( z.sort ) ) || sum( z.sort ) <= 0 ) {
    stop(
      "Density grid contains no usable (finite, non-zero) values, so no ",
      "contour can be extracted. This typically means the bandwidth used ",
      "for the kernel density estimate collapsed to zero because the ",
      "scatter data are too tightly clustered (e.g., highly monodisperse ",
      "beads, or saturated/quantized channel values). Check `bw`/`h` used ",
      "for the density estimate.",
      call. = FALSE
    )
  }

  cumulative.dens <- cumsum( z.sort ) / sum( z.sort )
  start.idx <- which.min( abs( cumulative.dens - percentile ) )

  # unique candidate threshold levels, starting at the requested percentile
  # and relaxing towards lower density (a looser gate) if no contour is found
  candidate.levels <- unique( z.sort[ start.idx : length( z.sort ) ] )

  for ( i in seq_along( candidate.levels ) ) {
    contour.lines <- grDevices::contourLines(
      dens$x, dens$y, dens$z, levels = candidate.levels[ i ]
    )

    if ( length( contour.lines ) > 0 ) {
      if ( i > 1 ) {
        warning(
          "Requested density percentile produced a degenerate contour ",
          "(commonly caused by a very tightly clustered population, e.g. ",
          "beads). Automatically relaxed the threshold to the nearest ",
          "workable density level.",
          call. = FALSE
        )
      }
      return( contour.lines[[
        which.max( sapply( contour.lines, function( l ) length( l$x ) ) )
      ]] )
    }
  }

  stop(
    "Could not extract a density contour at any threshold from the requested ",
    "percentile down to the minimum density in the grid. This usually means ",
    "the scatter data are too tightly clustered relative to the kernel ",
    "density grid resolution (common with highly monodisperse bead ",
    "populations). Try increasing `grid.n`, increasing the bandwidth ",
    "multiplier, or using `gating.system = \"density\"`.",
    call. = FALSE
  )
}
