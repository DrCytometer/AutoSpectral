#' @title Rank Replacement Fluorophores by Condition Number
#'
#' @description
#' Removes `remove` from `panel.fluors` and tries every combination of the
#' same number of fluorophores from `candidates`, returning the combinations
#' that give the lowest condition number ([calculate.condition.number()]).
#'
#' @param spectra.pool Matrix of available spectra (fluorophores x detectors).
#' @param panel.fluors Character vector, the current panel (row names of
#'   `spectra.pool`).
#' @param remove Character vector of fluorophores to take out.
#' @param candidates Fluorophores that may be added. Default `NULL`: every
#'   fluorophore in `spectra.pool` that is not in the panel, not in `avoid`
#'   and not `"AF"`.
#' @param avoid Fluorophores never added (e.g. a group of mutually collinear
#'   dyes).
#' @param extra Optional matrix of extra rows included in every matrix (e.g.
#'   an autofluorescence spectrum).
#' @param n.top Number of combinations returned. Default `10`.
#'
#' @return A data frame with columns `Remove`, `Add` and `Condition.number`,
#'   best first. Attribute `"original"`: the condition number of the
#'   unchanged panel (with `extra`).
#'
#' @export
suggest.panel.swaps <- function(
    spectra.pool,
    panel.fluors,
    remove,
    candidates = NULL,
    avoid = NULL,
    extra = NULL,
    n.top = 10
) {
  base <- setdiff( panel.fluors, remove )
  if ( is.null( candidates ) )
    candidates <- setdiff( rownames( spectra.pool ), c( panel.fluors, avoid, "AF" ) )
  candidates <- setdiff( intersect( candidates, rownames( spectra.pool ) ), c( base, avoid ) )
  k <- length( remove )
  if ( k == 0 || length( candidates ) < k ) stop( "Nothing to swap", call. = FALSE )
  combos <- utils::combn( candidates, k, simplify = FALSE )
  res <- do.call( rbind, lapply( combos, function( add ) {
    m <- rbind( spectra.pool[ c( base, add ), , drop = FALSE ], extra )
    data.frame( Remove = paste( remove, collapse = " + " ), Add = paste( add, collapse = " + " ),
                Condition.number = calculate.condition.number( m ), stringsAsFactors = FALSE )
  } ) )
  res <- res[ order( res$Condition.number ), , drop = FALSE ]
  res <- res[ seq_len( min( n.top, nrow( res ) ) ), , drop = FALSE ]
  rownames( res ) <- NULL
  attr( res, "original" ) <- calculate.condition.number(
    rbind( spectra.pool[ panel.fluors, , drop = FALSE ], extra ) )
  res
}
