# unmix_autospectral_joint.R

#' @title Joint AutoSpectral Unmixing (Pure R)
#'
#' @description
#' Pure-R implementation of the joint covariance-weighted AutoSpectral pipeline.
#' Mirrors `unmix_autospectral_joint_cpp` (in
#' `unmix_autospectral_joint_pipeline.cpp`) instruction-for-instruction,
#' without requiring `AutoSpectralRcpp`.
#'
#' The pipeline is staged identically to the C++ version:
#'   1. Global weighted pseudoinverse + AF helper pre-computation (Section 1).
#'   2. Per-fluorophore variant pre-computation, including leakage weights and
#'      rank-1 residual-update helpers (Section 2), plus a structural
#'      collinearity table used later for joint-pair conflict retries
#'      (Section 2B).
#'   3. Per-cell AF selection, with optional multi-pass AF refinement against
#'      the residual (`n.af.passes`).
#'   4. Per-cell fluorophore solve against the AF-subtracted residual (with
#'      optional per-cell Poisson-style weighting, `cell.weight`), followed by
#'      joint variant selection: candidates are scored with an alpha-weighted
#'      composite of residual ratio and covariance-propagated leakage ratio,
#'      conflicting swaps are resolved by residual-delta cosine similarity,
#'      and swaps flagged as a structurally collinear pair are queued for a
#'      retry once the winning partner's own swap has been applied. Each
#'      commit is individually verified against the currently accepted RSS
#'      and reverted if it does not improve it.
#'
#' @importFrom parallelly availableCores
#'
#' @param raw.data Numeric matrix (cells x detectors). Columns must match
#'   those of `spectra`.
#' @param spectra Numeric matrix (fluorophores x detectors, normalised 0-1).
#'   Must not contain an `"AF"` row.
#' @param af.spectra Numeric matrix (n_af x detectors, normalised 0-1).
#'   Prepare using `get.af.spectra()`. At least two rows are required.
#' @param asp The AutoSpectral parameter list, prepared using
#'   `get.autospectral.param()`. Required for parallel backend setup.
#' @param spectra.variants Named list as returned by `get.spectral.variants()`,
#'   containing at minimum `$variants`, `$delta.list`, and `$thresholds`.
#'   The optional `$optimize.recommended` slot is respected: fluorophores
#'   flagged `FALSE` are excluded from variant optimisation. Pass `NULL` for
#'   AF-only mode (no variant optimisation).
#' @param n.passes Integer, default `1L`. Number of joint optimisation passes
#'   per cell. Matches `n_passes` in `unmix_autospectral_joint_cpp`.
#' @param parallel Logical, default `TRUE`. Whether to use parallel processing
#'   across cells. Uses `create.parallel.lapply()` to handle platform
#'   differences.
#' @param threads Numeric, default `NULL`. Number of worker threads. When
#'   `NULL`, `asp$worker.process.n` is used. When `threads = 0`,
#'   `parallelly::availableCores()` is used. Ignored when `parallel = FALSE`.
#' @param cell.weight Logical, default `FALSE`. Enables per-cell Poisson-style
#'   detector weighting (`1 / max(|y_hat|, noise.floor)`), matching
#'   `cell_weight` in the C++ pipeline. When `FALSE`, weighting collapses to
#'   the (all-ones, unless a global weight was requested) `w.global` used in
#'   Section 1.
#' @param noise.floor Numeric scalar or length-D vector, default `NULL`.
#'   Per-detector noise floor used by the weighting scheme; falls back to a
#'   constant `125.0` per detector when `NULL`, matching the C++ default.
#'   Only used when `cell.weight = TRUE`.
#' @param alpha Numeric in `[0, 1]`, default `0.5`. Weighting exponent for the
#'   joint candidate score: `resid.ratio^alpha * leakage.ratio^(1 - alpha)`.
#' @param collinear.thresh Numeric, default `0.5`. Cosine-similarity
#'   threshold (in the pseudoinverse row space) above which two optimisable
#'   fluorophores are flagged as structurally collinear for joint-pair
#'   conflict retries.
#' @param joint.pair.resolution Logical, default `TRUE`. When `TRUE`, a
#'   candidate swap that conflicts with an already-committed swap belonging
#'   to a structurally collinear partner is queued and retried once that
#'   partner's own swap has been applied (rather than simply dropped).
#' @param n.af.passes Integer, default `1L`. Number of AF selection passes.
#'   When `> 1`, cells with an initial AF abundance at or above the
#'   `refine.af.quantile` quantile are re-scored against their residual in
#'   subsequent passes, accumulating additional AF abundance where the
#'   refinement score is < 1.0. The AF *index* used in the output is always
#'   the one selected on the first pass.
#' @param refine.af.quantile Numeric in `[0, 1]`, default `0.5`. Quantile
#'   (type 7) of the first-pass AF abundance used to select which cells are
#'   eligible for AF refinement passes.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Numeric matrix (cells x (n_fluorophores + 2)) with column names
#'   matching the row names of `spectra` plus `"AF"` (abundance of the selected
#'   AF component) and `"AF Index"` (1-based index into `af.spectra` of the
#'   selected AF spectrum).
#'
#' @seealso [unmix.autospectral()], [unmix.autospectral.rcpp()],
#'   [get.spectral.variants()], [create.parallel.lapply()]
#'
#' @export

unmix.autospectral.joint <- function(
    raw.data,
    spectra,
    af.spectra,
    asp,
    spectra.variants     = NULL,
    n.passes              = 1L,
    parallel              = TRUE,
    threads               = NULL,
    cell.weight           = FALSE,
    noise.floor           = NULL,
    alpha                 = 0.5,
    collinear.thresh      = 0.5,
    joint.pair.resolution = TRUE,
    n.af.passes           = 1L,
    refine.af.quantile    = 0.5,
    verbose               = TRUE
) {

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  if ( is.null( af.spectra ) || !is.matrix( af.spectra ) || nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided as a matrix.", call. = FALSE )

  raw.data.cols <- colnames( raw.data )
  spectra.cols  <- colnames( spectra )

  if ( !identical( raw.data.cols, spectra.cols ) ) {
    if ( all( spectra.cols %in% raw.data.cols ) &&
         length( spectra.cols ) == length( raw.data.cols ) ) {
      raw.data <- raw.data[ , spectra.cols ]
      message( "Columns reordered to match spectra." )
    } else {
      stop(
        "Column names in spectra and raw.data do not match perfectly;
        cannot reorder by name alone.",
        call. = FALSE
      )
    }
  }

  if ( n.af.passes < 1 )
    stop( "n.af.passes must be >= 1.", call. = FALSE )
  if ( refine.af.quantile < 0 || refine.af.quantile > 1 )
    stop( "refine.af.quantile must be between 0 and 1.", call. = FALSE )

  fluorophores <- rownames( spectra )
  F            <- nrow( spectra )
  D            <- ncol( spectra )
  N            <- nrow( raw.data )
  nAF          <- nrow( af.spectra )

  # noise floor: scalar fallback / expansion, matches C++ default of 125.0
  if ( is.null( noise.floor ) ) {
    noise.floor <- rep( 125.0, D )
  } else if ( length( noise.floor ) == 1L ) {
    noise.floor <- rep( noise.floor, D )
  }

  # -------------------------------------------------------------------------
  # Section 1: Global pre-computations  (mirrors C++ Section 1)
  # -------------------------------------------------------------------------

  # Global detector weight vector: 1 / max(mean detector signal, noise.floor)
  # when cell.weight is requested, all-ones otherwise.
  if ( cell.weight ) {
    col.means   <- colMeans( raw.data )
    w.global    <- 1 / pmax( col.means, noise.floor )
  } else {
    w.global    <- rep( 1, D )
  }
  sqrt.w.global <- sqrt( w.global )

  # Global weighted pseudoinverse P = (S_w S_w^T)^{-1} S_w, re-scaled back
  # into unweighted detector space: P[f,d] = P_w[f,d] * sqrt.w.global[d].
  spectra.w <- sweep( spectra, 2, sqrt.w.global, `*` )
  P.w       <- solve( spectra.w %*% t( spectra.w ), spectra.w )
  P         <- sweep( P.w, 2, sqrt.w.global, `*` )

  # AF helpers, mirroring the C++ weighted/unweighted split exactly.
  v_lib_af      <- P %*% t( af.spectra )                             # F x nAF
  r_lib_af      <- t( af.spectra ) - t( spectra ) %*% v_lib_af       # D x nAF
  r_dots_af     <- pmax( colSums( ( r_lib_af^2 ) * w.global ), 1e-10 )       # nAF, weighted
  r_lib_af_w2   <- sweep( r_lib_af, 1, w.global^2, `*` )             # D x nAF
  r_dots_af_raw <- colSums( r_lib_af^2 )                              # nAF, unweighted

  # Covariance-propagated AF endmember weights.
  af.cov.mat <- P %*% stats::cov( af.spectra ) %*% t( P )
  w_af       <- sqrt( abs( diag( af.cov.mat ) ) ) + 1e-8              # length F

  af.only <- is.null( spectra.variants ) || length( spectra.variants ) == 0

  # Per-cell (or per-residual) AF scoring closure — used both for the initial
  # AF pass and any refinement passes. Mirrors the `score_af` lambda in the
  # C++ pipeline exactly, including the rank-1 residual-norm update trick.
  score.af <- function( active.raw ) {
    init.f          <- as.numeric( P %*% active.raw )
    base.resid.af   <- active.raw - as.numeric( t( spectra ) %*% init.f )
    base.resid.sq   <- max( sum( base.resid.af^2 ), 1e-16 )
    base.resid.norm <- sqrt( base.resid.sq )
    base.fluor.l1   <- max( sum( w_af * abs( init.f ) ), 1e-8 )

    k.af.vec <- pmax( as.numeric( t( r_lib_af_w2 ) %*% active.raw ), 0 ) / r_dots_af

    cross.af    <- as.numeric( t( r_lib_af ) %*% base.resid.af )
    resid.sq.af <- base.resid.sq - 2 * ( k.af.vec * cross.af ) +
                   ( k.af.vec^2 * r_dots_af_raw )
    presid.af   <- sqrt( pmax( resid.sq.af, 0 ) ) / base.resid.norm

    diffs.af  <- sweep( v_lib_af, 2, k.af.vec, `*` )
    diffs.af  <- sweep( diffs.af, 1, init.f, `-` )
    pfluor.af <- as.numeric( t( abs( diffs.af ) ) %*% w_af ) / base.fluor.l1

    score.af.vec <- presid.af * pfluor.af
    best.j       <- which.min( score.af.vec )

    list( j = best.j, k = k.af.vec[ best.j ], score = score.af.vec[ best.j ] )
  }

  # -------------------------------------------------------------------------
  # Section 2: Per-fluorophore pre-computations  (mirrors C++ Section 2)
  # -------------------------------------------------------------------------
  precomp      <- list()
  active.names <- character( 0 )

  if ( !af.only ) {

    variants       <- spectra.variants$variants
    delta.list.in  <- spectra.variants$delta.list
    pos.thresholds <- rep( Inf, F )
    names( pos.thresholds ) <- fluorophores
    if ( !is.null( spectra.variants$thresholds ) )
      pos.thresholds[ names( spectra.variants$thresholds ) ] <- spectra.variants$thresholds

    opt.rec   <- spectra.variants$optimize.recommended
    opt.names <- if ( !is.null( opt.rec ) ) names( opt.rec )[ opt.rec ] else names( variants )
    opt.names <- intersect( opt.names, fluorophores )

    for ( fl in opt.names ) {
      if ( !fl %in% names( variants ) ) next

      v.mats <- variants[[ fl ]]     # n_variants x D
      n.var  <- nrow( v.mats )
      if ( n.var == 0 ) next

      master.row <- spectra[ fl, , drop = TRUE ]
      delta      <- v.mats - matrix( master.row, nrow = n.var, ncol = D, byrow = TRUE )

      other.fl <- fluorophores[ fluorophores != fl ]
      S_nof    <- spectra[ other.fl, , drop = FALSE ]
      U_nof    <- solve( S_nof %*% t( S_nof ), S_nof )

      v_lib <- U_nof %*% t( delta )   # (F-1) x n_variants

      # covariance-propagated leakage weights: ridge 1e-4 matches C++ Section 2
      delta.obs         <- if ( !is.null( delta.list.in[[ fl ]] ) )
                              delta.list.in[[ fl ]] else delta
      delta.cov         <- stats::cov( delta.obs )
      diag( delta.cov ) <- diag( delta.cov ) + 1e-4
      leakage.cov       <- U_nof %*% delta.cov %*% t( U_nof )
      w.leakage         <- sqrt( abs( diag( leakage.cov ) ) ) + 1e-8   # length F-1

      # rank-1 residual helpers: r_lib(:,v) = delta(v) - S^T P delta(v)
      r_lib    <- t( delta ) - t( spectra ) %*% ( P %*% t( delta ) )   # D x n_var
      r_lib_sq <- r_lib^2                                              # D x n_var
      r_dots   <- colSums( r_lib^2 )                                   # n_var, unweighted

      precomp[[ fl ]] <- list(
        master.idx = which( fluorophores == fl ),
        other.fl   = other.fl,
        v.mats     = v.mats,
        n.var      = n.var,
        v_lib      = v_lib,
        w.leakage  = w.leakage,
        r_lib      = r_lib,
        r_lib_sq   = r_lib_sq,
        r_dots     = r_dots,
        pos.thresh = pos.thresholds[ fl ]
      )
      active.names <- c( active.names, fl )
    }
  }

  # -------------------------------------------------------------------------
  # Section 2B: Structural collinearity precompute  (mirrors C++ Section 2B)
  # -------------------------------------------------------------------------
  is.collinear <- matrix( FALSE, F, F )
  if ( !af.only && length( active.names ) > 1 ) {
    idxs <- vapply( active.names, function( fl ) precomp[[ fl ]]$master.idx, integer( 1 ) )
    n.act <- length( idxs )
    for ( a in seq_len( n.act - 1 ) ) {
      for ( b in ( a + 1 ):n.act ) {
        fa <- idxs[ a ]; fb <- idxs[ b ]
        c.val <- abs( sum( P[ fa, ] * P[ fb, ] ) ) /
          ( sqrt( sum( P[ fa, ]^2 ) ) * sqrt( sum( P[ fb, ]^2 ) ) + 1e-12 )
        if ( c.val > collinear.thresh ) {
          is.collinear[ fa, fb ] <- TRUE
          is.collinear[ fb, fa ] <- TRUE
        }
      }
    }
  }

  if ( verbose ) {
    if ( af.only )
      message( "Running joint AutoSpectral pipeline (pure R): AF extraction only." )
    else
      message( paste0(
        "Running joint AutoSpectral pipeline (pure R): AF extraction + ",
        "variant optimisation for: ",
        paste( active.names, collapse = ", " ), "."
      ) )
  }

  # -------------------------------------------------------------------------
  # Section 3: parallel execution, staged to match the C++ pipeline
  # -------------------------------------------------------------------------

  # resolve threads
  if ( parallel ) {
    if ( is.null( threads ) ) threads <- asp$worker.process.n
    if ( !is.null( threads ) && threads == 0 )
      threads <- parallelly::availableCores()
  } else {
    threads <- 1L
  }

  # Only *static* (unchanging-throughout-the-pipeline) objects are exported to
  # the parallel backend. Anything that changes between stages (per-cell
  # residuals, AF abundance/index) is passed explicitly through the `lapply`
  # argument itself, so a PSOCK backend (which snapshots exports once, at
  # cluster creation) produces identical results to a forked backend.
  worker.exports <- c(
    "spectra", "af.spectra", "F", "D", "nAF", "fluorophores",
    "P", "v_lib_af", "r_lib_af", "r_lib_af_w2", "r_dots_af", "r_dots_af_raw", "w_af",
    "score.af",
    "af.only", "precomp", "active.names", "is.collinear",
    "cell.weight", "sqrt.w.global", "noise.floor",
    "n.passes", "alpha", "collinear.thresh", "joint.pair.resolution"
  )

  result.setup <- create.parallel.lapply(
    asp,
    exports    = worker.exports,
    parallel   = parallel,
    threads    = threads,
    export.env = environment()
  )
  lapply.function <- result.setup$lapply

  # -----------------------------------------------------------------------
  # Per-cell worker: fluorophore solve + joint variant selection, run once
  # AF selection/refinement (Stages A/B below) has already produced a final
  # residual, AF abundance, and AF index for every cell. Mirrors C++
  # Section 3's main `#pragma omp for` loop.
  # -----------------------------------------------------------------------
  process.cell <- function( input ) {

    cell.raw   <- input$raw
    cell.resid <- input$resid
    k.af       <- input$k.af
    best.j.af  <- input$j.af

    cell.S <- spectra   # F x D — AF already subtracted into cell.resid

    if ( cell.weight ) {
      S.F.w      <- sweep( cell.S, 2, sqrt.w.global, `*` )
      coeff.init <- as.numeric( solve(
        S.F.w %*% t( S.F.w ),
        S.F.w %*% ( cell.resid * sqrt.w.global )
      ) )
      y.hat <- as.numeric( t( cell.S ) %*% coeff.init ) + ( cell.raw - cell.resid )
      sw    <- 1 / sqrt( pmax( abs( y.hat ), noise.floor ) )
    } else {
      sw <- rep( 1, D )
    }

    cell.S.w      <- sweep( cell.S, 2, sw, `*` )
    fl.unm        <- as.numeric( solve(
      cell.S.w %*% t( cell.S.w ),
      cell.S.w %*% ( cell.resid * sw )
    ) )
    resid.raw     <- cell.resid - as.numeric( t( cell.S ) %*% pmax( fl.unm, 0 ) )
    resid         <- resid.raw * sw

    if ( af.only || length( active.names ) == 0 ) {
      return( c( fl.unm, k.af, best.j.af ) )
    }

    cell.resid.ss <- sum( cell.resid^2 )
    best.v        <- stats::setNames( rep( -1L, length( active.names ) ), active.names )
    y.vec         <- cell.resid * sw
    rss.accepted  <- NULL   # (re)initialised at the top of every pass below

    # ---------------------------------------------------------------------
    # try.commit: verifies a single candidate swap by re-solving the F x F
    # normal equations with that row of cell.S replaced, and only keeps it
    # if it reduces the currently accepted (weighted) RSS. Mathematically
    # equivalent to the C++ incremental Gram-matrix update — both simply
    # recompute A = S_w S_w^T for the current cell.S — just without the
    # rank-1 shortcut, which R doesn't need for correctness.
    # ---------------------------------------------------------------------
    try.commit <- function( fl, v ) {
      pc  <- precomp[[ fl ]]
      idx <- pc$master.idx

      prev.row      <- cell.S[ idx, ]
      cell.S[ idx, ] <<- pc$v.mats[ v, ]

      cell.S.w.trial <- sweep( cell.S, 2, sw, `*` )
      A.trial <- cell.S.w.trial %*% t( cell.S.w.trial )
      b.trial <- as.numeric( cell.S.w.trial %*% y.vec )

      trial.unmixed <- tryCatch( solve( A.trial, b.trial ), error = function( e ) NULL )
      if ( is.null( trial.unmixed ) ) {
        cell.S[ idx, ] <<- prev.row
        return( invisible( FALSE ) )
      }

      trial.resid.raw <- cell.resid - as.numeric( t( cell.S ) %*% pmax( trial.unmixed, 0 ) )
      trial.resid     <- trial.resid.raw * sw
      trial.rss       <- sum( trial.resid^2 )

      if ( trial.rss < rss.accepted ) {
        best.v[ fl ]  <<- v
        fl.unm        <<- trial.unmixed
        resid.raw     <<- trial.resid.raw
        resid         <<- trial.resid
        rss.accepted  <<- trial.rss
        invisible( TRUE )
      } else {
        cell.S[ idx, ] <<- prev.row
        invisible( FALSE )
      }
    }

    # ---------------------------------------------------------------------
    # C. JOINT VARIANT SELECTION  (n.passes passes)
    # ---------------------------------------------------------------------
    for ( pass in seq_len( n.passes ) ) {

      rss.curr        <- max( sum( resid^2 ), 1e-12 )
      rss.curr.sqrt   <- sqrt( rss.curr )
      ratio.thresh.sq <- 1.1025 * rss.curr   # (1.05^2) * rss.curr
      rss.accepted    <- sum( resid^2 )      # reset every pass, matches C++

      candidates    <- list()
      rsw           <- resid * sw
      w.eff         <- sw * sw
      q.ready       <- stats::setNames( rep( FALSE, length( active.names ) ), active.names )
      q.by.active   <- list()

      for ( fl in active.names ) {
        pc    <- precomp[[ fl ]]
        abund <- fl.unm[ pc$master.idx ]
        if ( abund < pc$pos.thresh ) next

        other.unm <- fl.unm[ -pc$master.idx ]
        base.leak <- max( sum( pc$w.leakage * abs( other.unm ) ), 1e-8 )
        cur.v     <- best.v[ fl ]

        if ( cell.weight ) {
          if ( !q.ready[ fl ] ) {
            q.by.active[[ fl ]] <- as.numeric( t( pc$r_lib_sq ) %*% w.eff )
            q.ready[ fl ]       <- TRUE
          }
          q.ref <- q.by.active[[ fl ]]
        } else {
          q.ref <- pc$r_dots
        }

        cross.v <- as.numeric( t( pc$r_lib ) %*% rsw )
        if ( cur.v < 0L ) {
          drsq.v <- q.ref
        } else {
          g.cur   <- if ( cell.weight )
            as.numeric( t( pc$r_lib ) %*% ( pc$r_lib[ , cur.v ] * w.eff ) )
          else
            as.numeric( t( pc$r_lib ) %*% pc$r_lib[ , cur.v ] )
          drsq.v  <- q.ref + q.ref[ cur.v ] - 2 * g.cur
          cross.v <- cross.v - cross.v[ cur.v ]
        }

        abund2 <- abund * abund

        for ( v in seq_len( pc$n.var ) ) {
          new.rss <- rss.curr - 2 * abund * cross.v[ v ] + abund2 * drsq.v[ v ]
          if ( new.rss > ratio.thresh.sq ) next   # fast reject: resid_ratio > 1.05

          vl <- pc$v_lib[ , v ]
          if ( cur.v < 0L ) {
            leak.num <- sum( pc$w.leakage * abs( other.unm - abund * vl ) )
          } else {
            vlc      <- pc$v_lib[ , cur.v ]
            leak.num <- sum( pc$w.leakage * abs( other.unm - abund * ( vl - vlc ) ) )
          }
          leakage.ratio <- leak.num / base.leak
          resid.ratio   <- sqrt( max( new.rss, 0 ) ) / rss.curr.sqrt

          joint.score <- max( resid.ratio,   1e-8 )^alpha *
                         max( leakage.ratio, 1e-8 )^( 1 - alpha )

          if ( joint.score < 1.0 )
            candidates[[ length( candidates ) + 1L ]] <- list(
              score = joint.score, fl = fl, v = v
            )
        }
      }

      if ( length( candidates ) == 0L ) break

      # sort candidates best-first
      scores.ord <- order( vapply( candidates, `[[`, numeric( 1 ), "score" ) )
      candidates <- candidates[ scores.ord ]

      # conflict resolution: commit non-overlapping swaps, queueing retries
      # for candidates that only conflict with a structurally collinear
      # committed partner.
      committed        <- stats::setNames( rep( FALSE, length( active.names ) ), active.names )
      committed.deltas <- list()
      commits          <- list()
      queued.retries   <- list()

      for ( cand in candidates ) {
        fl <- cand$fl
        if ( committed[ fl ] ) next

        pc    <- precomp[[ fl ]]
        abund <- fl.unm[ pc$master.idx ]
        cur.v <- best.v[ fl ]
        v     <- cand$v

        dr <- if ( cur.v < 0L ) pc$r_lib[ , v ] * abund
              else ( pc$r_lib[ , v ] - pc$r_lib[ , cur.v ] ) * abund
        dr.norm <- max( sqrt( sum( dr^2 ) ), 1e-12 )

        conflict <- FALSE
        for ( cd in committed.deltas ) {
          cosine <- abs( sum( dr * cd$dr ) ) / ( dr.norm * cd$norm )
          if ( cosine > 0.5 ) {
            conflict <- TRUE
            winner.master  <- precomp[[ cd$fl ]]$master.idx
            collinear.pair <- is.collinear[ pc$master.idx, winner.master ]
            if ( joint.pair.resolution && collinear.pair )
              queued.retries[[ length( queued.retries ) + 1L ]] <- list( fl = fl, v = v )
            break
          }
        }
        if ( conflict ) next

        committed[ fl ] <- TRUE
        committed.deltas[[ length( committed.deltas ) + 1L ]] <- list( dr = dr, norm = dr.norm, fl = fl )
        commits[[ length( commits ) + 1L ]] <- list( fl = fl, v = v )
      }

      if ( length( commits ) == 0L ) break

      for ( cm in commits ) try.commit( cm$fl, cm$v )

      # Joint-pair retry: replay candidates discarded for conflicting with a
      # structurally-collinear committed candidate, now that partner's own
      # swap (if any) has already been applied above.
      if ( joint.pair.resolution ) {
        for ( qr in queued.retries ) try.commit( qr$fl, qr$v )
      }

      if ( sum( resid.raw^2 ) < 1e-16 * cell.resid.ss ) break
    }

    return( c( fl.unm, k.af, best.j.af ) )
  }

  # -----------------------------------------------------------------------
  # Stage A: initial per-cell AF selection (parallel). raw.data is static
  # for the whole pipeline, so it's safe to reference it directly here even
  # under a PSOCK backend.
  # -----------------------------------------------------------------------
  stage.a <- tryCatch(
    expr    = lapply.function( seq_len( N ), function( i ) score.af( raw.data[ i, ] ) ),
    finally = NULL
  )
  af.index.vec     <- vapply( stage.a, `[[`, integer( 1 ), "j" )
  af.abundance.vec <- vapply( stage.a, `[[`, numeric( 1 ), "k" )
  resid.mat        <- raw.data - af.abundance.vec * af.spectra[ af.index.vec, , drop = FALSE ]

  # -----------------------------------------------------------------------
  # Stage B: optional AF refinement passes. Only cells at/above the
  # refine.af.quantile of first-pass abundance are refined. The AF *index*
  # recorded in the output is always the one chosen on the first pass.
  # -----------------------------------------------------------------------
  still.active <- rep( FALSE, N )
  if ( n.af.passes > 1 ) {
    af.refine.cutoff <- stats::quantile(
      af.abundance.vec, probs = refine.af.quantile, type = 7, names = FALSE
    )
    still.active <- af.abundance.vec >= af.refine.cutoff
  }

  for ( af.pass in seq_len( n.af.passes - 1L ) ) {
    active.idx <- which( still.active )
    if ( length( active.idx ) == 0L ) break

    resid.list  <- lapply( active.idx, function( i ) resid.mat[ i, ] )
    refine.res  <- lapply.function( resid.list, score.af )

    for ( k in seq_along( active.idx ) ) {
      i <- active.idx[ k ]
      r <- refine.res[[ k ]]
      if ( r$score < 1.0 ) {
        af.abundance.vec[ i ] <- af.abundance.vec[ i ] + r$k
        resid.mat[ i, ]       <- resid.mat[ i, ] - r$k * af.spectra[ r$j, ]
      } else {
        still.active[ i ] <- FALSE
      }
    }
  }

  # -----------------------------------------------------------------------
  # Stage C: main per-cell loop — fluorophore solve + joint variant
  # selection. Per-cell dynamic state (raw row, post-AF residual, AF
  # abundance/index) is bundled explicitly per call, as noted above.
  # -----------------------------------------------------------------------
  cell.results <- tryCatch(
    expr = {
      cell.inputs <- lapply( seq_len( N ), function( i ) list(
        raw   = raw.data[ i, ],
        resid = resid.mat[ i, ],
        k.af  = af.abundance.vec[ i ],
        j.af  = af.index.vec[ i ]
      ) )
      lapply.function( cell.inputs, process.cell )
    },
    finally = {
      if ( !is.null( result.setup$cleanup ) ) result.setup$cleanup()
    }
  )

  result <- do.call( rbind, cell.results )
  colnames( result ) <- c( fluorophores, "AF", "AF Index" )

  if ( verbose )
    message( paste0( "\033[32m", "Joint pipeline complete.", "\033[0m" ) )

  return( result )
}
