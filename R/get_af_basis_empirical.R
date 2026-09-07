# get_af_basis_empirical.r

#' @title Build A Continuous Autofluorescence Basis From Raw Events
#'
#' @description
#' As `get.af.basis()`, but derives the basis directly from the panel-oblique
#' residual of real unstained single-cell events rather than from a SOM-derived
#' spectral library. A SOM library is itself a quantisation of the unstained
#' sample's cell-to-cell diversity into a fixed set of node shapes, so building
#' the basis one step further upstream, from the events the SOM was trained on,
#' removes that quantisation before it happens rather than fitting around it.
#'
#' Structurally identical to `get.af.basis()`: both take the panel-oblique
#' residual of a set of reference spectra (there, SOM nodes; here, raw events),
#' take its SVD, and keep the leading out-of-span directions. Ranking by
#' out-of-span variance rather than total variance matters here more than in
#' the library case, because raw single-cell data also carries shot noise and
#' any residual spillover contamination, both of which add variance that has
#' nothing to do with autofluorescence shape.
#'
#' @param unstained.exprs Numeric matrix of raw unstained expression data,
#' cells in rows and detectors in columns. Columns must match `spectra`.
#' @param spectra Fluorophore spectral signatures, fluorophores in rows and
#' detectors in columns.
#' @param n.components Integer, number of basis directions to retain. Default
#' `NULL`, in which case the count is chosen from `var.explained`.
#' @param var.explained Numeric in (0, 1], default `0.99`. Fraction of the
#' out-of-span variance of the trimmed, subsampled event set the retained
#' basis must capture. Ignored when `n.components` is supplied.
#' @param trim.quantile Numeric in (0, 1], default `0.99`. Events whose
#' out-of-span residual norm exceeds this quantile are excluded before the
#' basis is built. Real unstained samples can carry debris, doublets, or
#' residual spillover-contaminated events; the SVD has no built-in robustness
#' to them, so they are trimmed explicitly rather than left to dominate the
#' leading components. Set to `1` to disable.
#' @param max.events Integer, default `5e4`. If more events remain after
#' trimming, a random subsample of this size is used to build the basis. The
#' basis is then applied to every event via `unmix.af.basis()`, so this only
#' bounds the cost of basis *construction*, not of applying it.
#' @param seed Integer, default `1`. Seed for the subsampling step.
#'
#' @return A list with `basis` (detectors x components, real combinations of
#' the sampled events), `directions` (detectors x components, the orthonormal
#' out-of-span directions), `sigma` (singular values), `var.explained`
#' (cumulative fraction captured), `n.components`, and `n.trimmed` (events
#' excluded by `trim.quantile`). Compatible with `unmix.af.basis()`.
#'
#' @seealso `get.af.basis`, `unmix.af.basis`
#'
#' @export

get.af.basis.empirical <- function(
    unstained.exprs,
    spectra,
    n.components  = NULL,
    var.explained = 0.99,
    trim.quantile = 0.99,
    max.events    = 5e4,
    seed          = 1
) {

  unstained.exprs <- as.matrix( unstained.exprs )
  spectra         <- as.matrix( spectra )

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )

  # panel residual of every unstained event: the part the fluorophores cannot
  # explain, which is exactly the space an AF basis should be built from
  unmixed.no.af <- unstained.exprs %*% t( unmixing.matrix )
  resid         <- unstained.exprs - unmixed.no.af %*% spectra   # N x D

  # trim events whose residual norm is anomalously large before they can
  # dominate the SVD
  resid.norm <- sqrt( rowSums( resid^2 ) )
  cutoff     <- stats::quantile( resid.norm, trim.quantile, na.rm = TRUE )
  keep       <- resid.norm <= cutoff
  n.trimmed  <- sum( !keep )

  events.kept <- unstained.exprs[ keep, , drop = FALSE ]

  # subsample for the SVD only; the basis is applied to every event
  # afterward via unmix.af.basis(), so this bounds construction cost only
  if ( nrow( events.kept ) > max.events ) {
    set.seed( seed )
    events.kept <- events.kept[ sample( nrow( events.kept ), max.events ), , drop = FALSE ]
  }

  resid.sub <- events.kept - ( events.kept %*% t( unmixing.matrix ) ) %*% spectra

  # D x N SVD, mirroring get.af.basis()'s r.library orientation exactly:
  # r.emp = Pperp %*% t(events.kept) = t(resid.sub). u gives the orthonormal
  # out-of-span directions ("W"); v %*% events.kept gives the basis as real
  # combinations of the sampled events, matching t(af.spectra) %*% v in
  # get.af.basis().
  r.emp <- t( resid.sub )                       # D x N.sub
  sv    <- svd( r.emp )

  power <- sv$d^2
  cumul <- cumsum( power ) / sum( power )

  q <- if ( !is.null( n.components ) )
    min( as.integer( n.components ), length( sv$d ) ) else
      which( cumul >= var.explained )[ 1 ]
  if ( is.na( q ) ) q <- length( sv$d )

  keep.q <- seq_len( q )
  basis      <- t( events.kept ) %*% sv$v[ , keep.q, drop = FALSE ]   # D x q
  directions <- sv$u[ , keep.q, drop = FALSE ]                        # D x q

  rownames( basis )      <- colnames( unstained.exprs )
  colnames( basis )      <- paste0( "AF", keep.q )
  rownames( directions ) <- colnames( unstained.exprs )
  colnames( directions ) <- paste0( "AF", keep.q )

  list(
    basis         = basis,
    directions    = directions,
    sigma         = sv$d[ keep.q ],
    var.explained = cumul[ q ],
    n.components  = q,
    n.trimmed     = n.trimmed
  )
}
