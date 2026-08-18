# get_spread_thresholds.R

#' @title Align a Spillover Spreading Matrix to a Fluorophore Set
#'
#' @description
#' Coerces the `spillover.spread` element returned by `get.spectral.variants()`
#' into a full, square, finite, non-negative source x target matrix over a given
#' fluorophore set. The stored matrix is not square (rows cover only controls
#' with enough positive events, columns cover the whole panel), carries `NA` on
#' the diagonal, and can contain small negative entries where the positive
#' population's MAD fell below the unstained MAD.
#'
#' @param spillover.spread Matrix (source fluorophore x target channel), or
#'   `NULL`.
#' @param fluor.names Character vector of fluorophore names to align to.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#'
#' @return Numeric matrix (`fluor.names` x `fluor.names`), zero on the diagonal
#'   and zero wherever no estimate was available.
#'
#' @noRd

.align.spillover.spread <- function( spillover.spread, fluor.names, verbose = TRUE ) {

  spread <- matrix(
    0, nrow = length( fluor.names ), ncol = length( fluor.names ),
    dimnames = list( fluor.names, fluor.names )
  )

  if ( is.null( spillover.spread ) )
    return( spread )

  source.shared <- intersect( rownames( spillover.spread ), fluor.names )
  target.shared <- intersect( colnames( spillover.spread ), fluor.names )

  if ( length( source.shared ) == 0 || length( target.shared ) == 0 ) {
    if ( verbose )
      message( "No fluorophore in `spillover.spread` matches the supplied set; ",
               "positivity boundaries will be flat." )
    return( spread )
  }

  block <- spillover.spread[ source.shared, target.shared, drop = FALSE ]
  block[ !is.finite( block ) ] <- 0
  block[ block < 0 ] <- 0

  spread[ source.shared, target.shared ] <- block
  diag( spread ) <- 0

  missing.source <- setdiff( fluor.names, source.shared )
  if ( length( missing.source ) > 0 && verbose )
    message(
      sprintf( "No spread estimate for %d fluorophore(s), treated as contributing none: %s",
               length( missing.source ), paste( missing.source, collapse = ", " ) )
    )

  spread
}


#' @title Spread-Scaled Positivity Thresholds
#'
#' @description
#' Converts a flat, per-fluorophore positivity threshold into a per-event (or
#' per-cluster) threshold that widens with the spillover spread contributed by
#' every fluorophore present in that event.
#'
#' A flat cut taken from an unstained control describes the width of the
#' negative population only. Once a dye is bright, its photon and unmixing noise
#' inflates the estimate in every channel it spills into, so a bright event can
#' clear a flat cut in a channel it carries no signal in at all. The Spillover
#' Spreading Matrix from `get.spectral.variants()` measures exactly this: entry
#' `[a, b]` is the *variance* added to channel `b` per unit of fluorophore `a`'s
#' own on-channel abundance. The boundary therefore grows as the square root of
#' abundance, not in proportion to it:
#'
#' \deqn{t_{c,b} = m \cdot t_b + \kappa \sqrt{ \sum_a S\!S_{a,b} \max(x_{c,a}, 0) }}
#'
#' Contributions are summed over every source, since in a stained sample several
#' dyes spread into the same channel at once.
#'
#' @param unmixed Numeric matrix (events or clusters x fluorophores), unmixed
#'   abundances. Column names must be the fluorophore names.
#' @param thresholds Named numeric vector of flat positivity thresholds in
#'   unmixed space, covering every column of `unmixed` (typically the 99.5th
#'   percentile of an unstained control, or `get.spectral.variants()$thresholds`).
#' @param spillover.spread Matrix (source fluorophore x target channel), the
#'   Spillover Spreading Matrix from `get.spectral.variants()$spillover.spread`.
#'   `NULL` (default) reproduces a flat threshold, broadcast to matrix shape.
#' @param spread.kappa Numeric, how many spread standard deviations to allow
#'   above the flat threshold. Default `2`.
#' @param margin Numeric multiplier applied to the flat component only, matching
#'   the `unstained.margin` convention in `fix.my.unmix()`. Default `1`.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#'
#' @return Numeric matrix with the same dimensions and dimnames as `unmixed`,
#'   giving the positivity threshold for each row and fluorophore.
#'
#' @export

get.spread.thresholds <- function(
    unmixed,
    thresholds,
    spillover.spread = NULL,
    spread.kappa     = 2,
    margin           = 1,
    verbose          = TRUE
  ) {

  if ( !is.matrix( unmixed ) ) unmixed <- as.matrix( unmixed )

  fluor.names <- colnames( unmixed )

  if ( is.null( fluor.names ) )
    stop( "`unmixed` must have fluorophore column names.", call. = FALSE )

  if ( !all( fluor.names %in% names( thresholds ) ) )
    stop( "`thresholds` must be named and cover every column of `unmixed`.",
          call. = FALSE )

  flat <- thresholds[ fluor.names ] * margin

  spread <- .align.spillover.spread( spillover.spread, fluor.names, verbose = verbose )

  if ( all( spread == 0 ) ) {

    threshold.matrix <- matrix(
      flat, nrow = nrow( unmixed ), ncol = length( fluor.names ), byrow = TRUE
    )

  } else {

    spread.variance  <- pmax( unmixed, 0 ) %*% spread
    threshold.matrix <- sweep( spread.kappa * sqrt( spread.variance ), 2, flat, "+" )

  }

  dimnames( threshold.matrix ) <- dimnames( unmixed )

  threshold.matrix
}
