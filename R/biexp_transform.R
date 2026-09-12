# biexp_transform.R

#' Fast, uncapped biexponential (logicle) transform
#'
#' A drop-in replacement for \code{flowWorkspace::flowjo_biexp()} built
#' directly on flowCore's compiled logicle transform (Moore & Parks,
#' Cytometry A 2012), rather than flowWorkspace's spline-interpolated
#' re-implementation of FlowJo's own (undocumented) display algorithm.
#'
#' FlowJo's \code{widthBasis} is a legacy re-parameterisation of the
#' logicle transform's width parameter \code{W} (in asymptotic decades):
#' \deqn{W = \log_{10}(-\mathrm{widthBasis}) / 2}
#' The logicle transform requires \eqn{0 \le W \le M/2} and
#' \eqn{-W \le A \le M - 2W}, where \code{M} is \code{pos} and \code{A} is
#' \code{neg}. FlowJo's UI, and \code{flowWorkspace::flowjo_biexp()}, impose
#' a hard floor of \code{widthBasis = -1000} (\code{W = 1.5}) regardless of
#' \code{pos}. This  function checks the real constraint against whatever
#' \code{pos}/\code{neg} you supply and errors with the exact numbers if that's
#' violated, rather than silently clipping or producing a degenerate curve.
#'
#' @param channelRange numeric. Maximum value of the transformed (display)
#'   scale. Default \code{4096}.
#' @param maxValue numeric. Maximum value of the input (raw) scale.
#'   Default \code{262144}.
#' @param pos numeric. Number of positive decades spanned by the transform
#'   (\code{M} in the logicle parameterisation). Default \code{4.5}.
#' @param neg numeric. Additional negative decades brought on scale
#'   (\code{A} in the logicle parameterisation). Default \code{0}.
#' @param widthBasis numeric, negative. FlowJo-style width parameter; see
#'   Details for its relationship to the logicle \code{W}. Default \code{-10}.
#' @param inverse logical. If \code{TRUE}, returns the inverse transform
#'   (display scale back to raw values) instead of the forward transform.
#'
#' @return A function mapping a numeric vector on the input scale to the
#'   transformed display scale (or the reverse, if \code{inverse = TRUE}).
#'
#' @seealso \code{\link[flowCore]{logicleTransform}},
#'   \code{\link[flowCore]{inverseLogicleTransform}}
#' @importFrom flowCore logicleTransform inverseLogicleTransform
#' @export
biexp.transform <- function( channelRange = 4096, maxValue = 262144, pos = 4.5,
                        neg = 0, widthBasis = -10, inverse = FALSE ) {
  
  if ( widthBasis >= 0 ) {
    stop( "widthBasis must be negative (got ", widthBasis, ")." )
  }
  
  w <- log10( -widthBasis ) / 2
  
  if ( w > pos / 2 ) {
    stop(
      "widthBasis = ", widthBasis, " implies a logicle width of W = ",
      round( w, 3 ), " decades, which exceeds the valid maximum of ",
      "pos / 2 = ", round( pos / 2, 3 ), " decades (pos = ", pos, "). ",
      "Either raise pos to at least ", round( 2 * w, 3 ),
      ", or use a widthBasis no more negative than ",
      round( -10 ^ pos, 1 ), "."
    )
  }
  
  if ( neg < -w || neg > pos - 2 * w ) {
    stop(
      "neg = ", neg, " is outside the valid range [", round( -w, 3 ), ", ",
      round( pos - 2 * w, 3 ), "] for widthBasis = ", widthBasis,
      " and pos = ", pos, "."
    )
  }
  
  trans <- flowCore::logicleTransform( w = w, t = maxValue, m = pos, a = neg )
  scale.factor <- channelRange / pos
  
  if ( inverse ) {
    inv <- flowCore::inverseLogicleTransform( trans = trans )
    function( x ) inv( x / scale.factor )
  } else {
    function( x ) trans( x ) * scale.factor
  }
}