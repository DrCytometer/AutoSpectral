# safe_bandwidth.r

## Computes MASS::bandwidth.nrd() per column, with a fallback for degenerate
## columns where the interquartile range is zero (e.g., highly monodisperse
## or saturated/quantized bead scatter data), which would otherwise silently
## return a bandwidth of 0. A bandwidth of 0 passed to MASS::kde2d() throws
## "bandwidths must be strictly positive", but passed to
## AutoSpectralRcpp::fast_kde2d_cpp() it silently produces a density surface
## of all zeros (0/0 comparisons in the Gaussian smoothing kernel fall through
## to the zero branch rather than propagating NaN), which then breaks
## downstream contour/threshold extraction with no informative error message
## at the point of failure. When bandwidth.nrd() returns a non-positive or
## non-finite value, this falls back to a standard-deviation-based bandwidth,
## and if that is also degenerate, to a fixed fraction of the observed data
## range, with a warning either way so the degenerate input doesn't pass
## silently.
## Used by define.gate.landmarks(), gate.scatter.match() and tune.gate().
##
## @param x Numeric matrix or data.frame of scatter data columns.
## @return Named numeric vector of per-column bandwidths, all strictly
## positive and finite.
## @keywords internal
.safe.bandwidth <- function( x ) {
  
  apply( x, 2, function( col ) {
    
    bw <- tryCatch( MASS::bandwidth.nrd( col ), error = function( e ) NA_real_ )
    
    if ( is.na( bw ) || !is.finite( bw ) || bw <= 0 ) {
      # IQR-based estimate is degenerate (many tied values); fall back to a
      # standard-deviation-based bandwidth of the same functional form
      bw <- 4 * 1.06 * stats::sd( col ) * length( col ) ^ ( -1 / 5 )
    }
    
    if ( is.na( bw ) || !is.finite( bw ) || bw <= 0 ) {
      # spread is (near) zero in both IQR and SD; fall back to a small
      # fraction of the observed range so the density estimate stays defined
      col.range <- diff( range( col ) )
      bw <- if ( col.range > 0 ) col.range * 0.01 else 1
      
      warning(
        "Scatter data has (near) zero variance in one gating dimension; ",
        "using a fallback bandwidth. This usually indicates a highly ",
        "monodisperse or saturated population (e.g., beads); the resulting ",
        "gate may be unreliable and should be checked visually.",
        call. = FALSE
      )
    }
    
    bw
  } )
}