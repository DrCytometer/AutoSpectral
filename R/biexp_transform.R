# biexp_transform.R

#' Fast, uncapped biexponential (logicle) transform
#'
#' @description
#' A self-contained implementation of the logicle transform (Parks,
#' Roederer & Moore, Cytometry A 2006), used as a replacement for
#' \code{flowWorkspace::flowjo_biexp()} with the same
#' \code{channelRange}/\code{maxValue}/\code{pos}/\code{neg}/\code{widthBasis}
#' argument names, but with no hard-coded \code{-1000} floor on
#' \code{widthBasis}. Uses only base R (\code{stats::uniroot()} for a
#' one-time scalar solve).
#'
#' @details
#' FlowJo's \code{widthBasis} is a legacy re-parameterisation of the
#' logicle transform's width parameter \code{W} (in asymptotic decades):
#' \deqn{W = \log_{10}(-\mathrm{widthBasis}) / 2}
#' The logicle transform requires \eqn{0 < W \le M/2}, where \code{M} is
#' \code{pos}. FlowJo's UI, and \code{flowWorkspace::flowjo_biexp()}, impose
#' a hard floor of \code{widthBasis = -1000} (\code{W = 1.5}) regardless of
#' \code{pos} -- a historical GUI choice, not a mathematical one. This
#' function checks the real constraint against whatever \code{pos} you
#' supply and errors with the exact numbers if it's violated.
#'
#' Only \code{neg = 0} is supported (the value used by every cytometer
#' profile in this package). Extending to \code{neg > 0} requires the
#' shifted-origin form of the transform (GatingML 2.0's \code{A} parameter)
#' and is deliberately left unimplemented until there's a concrete need for
#' it and data to verify it against.
#'
#' The forward direction (raw intensity -> display decade) has no closed
#' form and is solved by Newton-Raphson, applied to the whole input vector
#' at once for a fixed number of iterations (not a per-point loop), starting
#' from \code{asinh(x / (2 * a.scale))} -- correct in the large-\code{|x|}
#' limit and reasonable near zero. The inverse direction (display decade ->
#' raw intensity) is closed-form.
#'
#' @param channelRange numeric. Maximum value of the transformed (display)
#'   scale. Default \code{4096}.
#' @param maxValue numeric. Maximum value of the input (raw) scale.
#'   Default \code{262144}.
#' @param pos numeric. Number of positive decades spanned by the transform
#'   (\code{M} in the logicle parameterisation). Default \code{4.5}.
#' @param neg numeric. Must be \code{0}; see Details.
#' @param widthBasis numeric, negative. FlowJo-style width parameter; see
#'   Details for its relationship to the logicle \code{W}. Default \code{-10}.
#' @param inverse logical. If \code{TRUE}, returns the inverse transform
#'   (display scale back to raw values) instead of the forward transform.
#' @param newton.iter integer. Fixed number of vectorised Newton iterations
#'   used to solve the forward transform. Default \code{5} -- the smallest
#'   value that reaches double-precision round-trip error across the ranges
#'   swept in testing. Has no effect on the inverse direction, which is
#'   closed-form. If a different \code{pos}/\code{widthBasis} combination
#'   needs more, the built-in residual check will warn rather than silently
#'   under-converge.
#'
#' @return A function mapping a numeric vector on the input scale to the
#'   transformed display scale (or the reverse, if \code{inverse = TRUE}).
#'
#' @references
#' Parks DR, Roederer M, Moore WA (2006). A new "Logicle" display method
#' avoids deceptive effects of logarithmic scaling for low signals and
#' compensated data. Cytometry A, 69(6):541-551.
#'
#' @export
biexp.transform <- function(
    channelRange = 4096,
    maxValue = 262144,
    pos = 4.5,
    neg = 0,
    widthBasis = -10,
    inverse = FALSE,
    newton.iter = 5L
  ) {

  if ( widthBasis >= 0 ) {
    stop( "widthBasis must be negative (got ", widthBasis, ")." )
  }

  if ( neg != 0 ) {
    stop(
      "biexp.fast() currently only supports neg = 0 (got ", neg, "). ",
      "Extending to neg > 0 needs the GatingML 2.0 shifted-origin form -- ",
      "not implemented. Use flowCore::logicleTransform() for that case."
    )
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

  # --- one-time solve for the shape parameter p, from w = 2p*ln(p)/(p+1) ---
  w.nat <- w * log( 10 )
  m.nat <- pos * log( 10 )

  if ( w.nat <= 0 ) {
    p <- 1
  } else {
    p.eqn <- function( p ) 2 * p * log( p ) / ( p + 1 ) - w.nat
    p <- stats::uniroot( p.eqn, interval = c( 1 + 1e-12, 1e12 ),
                         tol = .Machine$double.eps ^ 0.5 )$root
  }

  p2 <- p ^ 2
  a.scale <- maxValue / ( exp( m.nat - w.nat ) - p2 * exp( -( m.nat - w.nat ) / p ) + p2 - 1 )
  scale.factor <- channelRange / pos

  if ( inverse ) {

    function( x ) {
      u <- ( x / scale.factor ) * log( 10 ) - w.nat
      sign( u ) * a.scale * ( exp( abs( u ) ) - p2 * exp( -abs( u ) / p ) + p2 - 1 )
    }

  } else {

    function( x ) {

      u <- asinh( x / ( 2 * a.scale ) )

      for ( i in seq_len( newton.iter ) ) {
        abs.u   <- abs( u )
        e.pos   <- exp( abs.u )
        e.neg   <- exp( -abs.u / p )
        v       <- sign( u ) * a.scale * ( e.pos - p2 * e.neg + p2 - 1 )
        v.prime <- a.scale * ( e.pos + p * e.neg )
        u <- u - ( v - x ) / v.prime
      }

      abs.u <- abs( u )
      v <- sign( u ) * a.scale * ( exp( abs.u ) - p2 * exp( -abs.u / p ) + p2 - 1 )

      residual <- max( abs( v - x ) / pmax( abs( x ), 1 ) )
      if ( residual > 1e-6 ) {
        warning(
          "biexp.fast(): Newton solve did not fully converge (max relative ",
          "residual = ", signif( residual, 3 ), "). Consider raising ",
          "newton.iter."
        )
      }

      ( u + w.nat ) / log( 10 ) * scale.factor
    }
  }
}
