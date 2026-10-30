# fix_my_unmix_smooth.R

#' Fix My Unmix in rounds on fresh subsamples (`fix.my.unmix()`'s
#' `smooth.updates`).
#'
#' A single pass adopts its own subsample's gate decisions, and those are hard
#' cuts on noisy statistics, so which rows change depends on which events were
#' seen. Each round here draws a fresh `fraction` of the stained events,
#' unmixes them and the unstained events again with the spectra the round
#' starts from (the per-event autofluorescence assignment and every
#' abundance follow the spectra), and runs the core on them.
#'
#' With `mode = "damped"` each round starts from the current spectra and moves
#' the rows it changed a step `1 / round` toward its result (Robbins-Monro), so
#' a correction every round finds is adopted and one found in a single round
#' only partly. With `mode = "consensus"` every round starts from
#' `spectra.fluor`, and a row is adopted as the mean of the rounds that accepted
#' it when at least `consensus.min.rounds` of them did and their mean angle to
#' that mean is at most `consensus.max.spread` of the move it makes. Either way
#' the averaged rows are held to the plausibility box of `spectra.fluor`.
#'
#' @param spectra.fluor Numeric matrix (fluorophores x detectors), the start.
#' @param project `project( spectra.fluor, keep )` from `.fix.projector()`.
#' @param n.stained Integer, the number of stained events `keep` indexes.
#' @param run.core `function( start, projection, seed )` running
#'   `.fix.my.unmix.core()`.
#' @param rounds,fraction,mode,consensus.min.rounds,consensus.max.spread See
#'   `fix.my.unmix()`'s `smooth.*` and `consensus.*` arguments.
#' @param plausibility.ratio,plausibility.floor,plausibility.reset See
#'   `fix.my.unmix()`.
#' @param seed Integer; round `r` draws its events with seed `seed + r - 1`.
#' @param verbose Logical.
#'
#' @return A list: `result`, the last round's core result with `spectra` the
#'   smoothed rows, `smooth.log` one row per round, and in `signature.log`
#'   `rounds.accepted` and (consensus) `consensus.spread` and
#'   `consensus.move`, every row the rounds moved marked accepted; and
#'   `projection`, the last round's projection.
#'
#' @noRd
.fix.smooth <- function( spectra.fluor, project, n.stained, run.core,
                         rounds = 5L, fraction = 0.5,
                         mode = c( "damped", "consensus" ),
                         consensus.min.rounds = 0.6, consensus.max.spread = 0.5,
                         plausibility.ratio = 1, plausibility.floor = 0.003,
                         plausibility.reset = TRUE, seed = 0L, verbose = TRUE ) {

  mode <- match.arg( mode )

  normalise <- function( rows ) {
    rows <- pmax( rows, 0 )
    rows / pmax( apply( rows, 1, max ), .Machine$double.eps )
  }

  start   <- spectra.fluor
  current <- start

  reference <- normalise( start )
  allow <- if ( plausibility.ratio > 0 ) .fix.plausible.change(
    reference, plausibility.ratio, plausibility.floor ) else NULL

  hold <- function( rows, which.rows ) {
    if ( is.null( allow ) ) return( rows )
    for ( f in which.rows )
      rows[ f, ] <- .fix.clip.to.plausible( rows[ f, ], reference[ f, ], allow[ f, ],
                                            plausibility.reset )$row
    rows
  }

  n.keep <- max( 1L, as.integer( round( fraction * n.stained ) ) )

  smooth.log     <- vector( "list", rounds )
  round.rows     <- vector( "list", rounds )
  round.accepted <- vector( "list", rounds )
  result         <- NULL
  projection     <- NULL

  for ( r in seq_len( rounds ) ) {

    if ( verbose )
      message( sprintf( "\033[36mRound %d of %d.\033[0m", r, rounds ) )

    set.seed( seed + r - 1L )
    keep <- if ( n.keep >= n.stained ) NULL else
      sort( sample.int( n.stained, n.keep ) )

    from       <- if ( mode == "consensus" ) start else current
    projection <- project( from, keep )
    result     <- run.core( from, projection, seed + r - 1L )

    moved <- .signature.row.angle( current, result$spectra )

    round.accepted[[ r ]] <- result$accepted[ rownames( start ) ]

    if ( mode == "consensus" ) {

      round.rows[[ r ]] <- result$spectra

    } else {

      # Only the rows this round changed take a step, so the others stay
      # exactly as they were.
      changed <- which( rowSums( result$spectra != current ) > 0 )

      if ( length( changed ) > 0L ) {
        current[ changed, ] <- normalise(
          current[ changed, , drop = FALSE ] +
            ( result$spectra[ changed, , drop = FALSE ] -
                current[ changed, , drop = FALSE ] ) / r )
        current <- hold( current, changed )
      }
    }

    smooth.log[[ r ]] <- data.frame(
      round         = r,
      n.events      = if ( is.null( keep ) ) n.stained else length( keep ),
      stop.reason   = result$stop.reason,
      iterations    = nrow( result$convergence.log ),
      accepted      = paste( names( which( result$accepted ) ), collapse = ", " ),
      max.deg.round = if ( length( moved ) > 0L ) max( moved ) else 0,
      row.names = NULL, stringsAsFactors = FALSE )
  }

  smooth.log <- do.call( rbind, smooth.log )

  fluorophores    <- rownames( start )
  accepted.rounds <- do.call( rbind, round.accepted )
  agreement       <- NULL

  if ( mode == "consensus" ) {

    agreement <- data.frame(
      fluorophore      = fluorophores,
      consensus.spread = NA_real_,
      consensus.move   = NA_real_,
      row.names = NULL, stringsAsFactors = FALSE )

    for ( f in seq_along( fluorophores ) ) {

      chosen.rounds <- which( accepted.rounds[ , f ] )
      if ( length( chosen.rounds ) < consensus.min.rounds * rounds ) next

      chosen <- do.call( rbind, lapply( round.rows[ chosen.rounds ], function( m )
        m[ f, , drop = FALSE ] ) )
      mean.row <- normalise( matrix( colMeans( chosen ), nrow = 1 ) )

      spread <- mean( .signature.row.angle(
        chosen, mean.row[ rep( 1L, nrow( chosen ) ), , drop = FALSE ] ) )
      move <- .signature.row.angle( start[ f, , drop = FALSE ], mean.row )

      agreement$consensus.spread[ f ] <- spread
      agreement$consensus.move[ f ]   <- move

      if ( move > 0 && spread <= consensus.max.spread * move ) {
        current[ f, ] <- mean.row[ 1, ]
        current <- hold( current, f )
      }
    }
  }

  # A row counts as accepted when the rounds moved it; its reason says in how
  # many rounds it was.
  moved.rows <- rowSums( current != start ) > 0

  if ( !is.null( result$signature.log ) ) {

    signature.log <- result$signature.log
    row.of <- match( signature.log$fluorophore, fluorophores )

    rounds.accepted <- colSums( accepted.rounds )[ row.of ]

    signature.log$rounds.accepted <- unname( rounds.accepted )

    if ( !is.null( agreement ) ) {
      signature.log$consensus.spread <- agreement$consensus.spread[ row.of ]
      signature.log$consensus.move   <- agreement$consensus.move[ row.of ]
    }

    for ( i in seq_len( nrow( signature.log ) ) ) {

      f <- row.of[ i ]

      if ( moved.rows[ f ] ) {
        signature.log$accepted[ i ] <- TRUE
        signature.log$reason[ i ] <- if ( rounds.accepted[ i ] == rounds ) "accepted" else
          sprintf( "accepted in %d of %d rounds", rounds.accepted[ i ], rounds )
      } else if ( !is.null( agreement ) && is.finite( agreement$consensus.move[ f ] ) ) {
        signature.log$accepted[ i ] <- FALSE
        signature.log$reason[ i ]   <- "rounds.disagree"
      } else if ( signature.log$accepted[ i ] ) {
        signature.log$accepted[ i ] <- FALSE
        signature.log$reason[ i ]   <- sprintf( "accepted in %d of %d rounds, too few",
                                                rounds.accepted[ i ], rounds )
      }
    }

    result$signature.log <- signature.log
  }

  result$spectra    <- current
  result$accepted   <- stats::setNames( moved.rows, fluorophores )
  result$smooth.log <- smooth.log

  list( result = result, projection = projection )
}
