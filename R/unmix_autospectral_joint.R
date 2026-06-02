# unmix_autospectral_joint.R

#' @title Joint AutoSpectral Unmixing (Pure R)
#'
#' @description
#' Pure-R implementation of the joint covariance-weighted AutoSpectral pipeline.
#' Mirrors the algorithm in `unmix_autospectral_joint_cpp` without requiring
#' `AutoSpectralRcpp`. For each cell, AF variant selection and fluorophore
#' spectral variant optimisation are performed in a single jointly-scored pass
#' rather than sequentially.
#'
#' The AF selection score is a multiplicative product of the residual ratio and
#' a covariance-propagated fluorophore leakage ratio, matching Section A of the
#' C++ pipeline. Fluorophore variant selection uses the same composite score
#' across all optimisable fluorophores before any swaps are committed, with
#' conflict resolution by residual-delta cosine similarity (Section C).
#'
#' Parallelisation follows the same pattern as `unmix.autospectral()`, using
#' `create.parallel.lapply()` to handle platform differences between Windows
#' (PSOCK cluster) and Linux (forked `mclapply`).
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
#' @param n.passes Integer, default `2L`. Number of joint optimisation passes
#'   per cell. Matches `n_passes` in `unmix_autospectral_joint_cpp`. Higher
#'   values allow more variant swaps to be committed per cell at the cost of
#'   additional computation.
#' @param n.variants Integer, default `1L`. Maximum number of top-scoring
#'   variants to evaluate per fluorophore per pass. Corresponds to `k_opt`
#'   in the legacy pipeline. Set via `speed` in the calling function, or
#'   supply directly.
#' @param parallel Logical, default `TRUE`. Whether to use parallel processing
#'   across cells. Uses `create.parallel.lapply()` to handle platform
#'   differences.
#' @param threads Numeric, default `NULL`. Number of worker threads. When
#'   `NULL`, `asp$worker.process.n` is used. When `threads = 0`,
#'   `parallelly::availableCores()` is used. Ignored when `parallel = FALSE`.
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
    spectra.variants = NULL,
    n.passes         = 2L,
    n.variants       = 1L,
    parallel         = TRUE,
    threads          = NULL,
    verbose          = TRUE
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

  fluorophores <- rownames( spectra )
  F            <- nrow( spectra )
  D            <- ncol( spectra )
  N            <- nrow( raw.data )
  nAF          <- nrow( af.spectra )

  # -------------------------------------------------------------------------
  # Section 1: Global pre-computations  (mirrors C++ Section 1)
  # -------------------------------------------------------------------------

  # Global pseudoinverse P = (S S^T)^{-1} S,  shape F x D
  P <- solve( spectra %*% t( spectra ), spectra )

  # AF helpers: v_lib_af and r_lib_af for rank-1 update scoring
  v_lib_af <- P %*% t( af.spectra )                           # F   x nAF
  r_lib_af <- t( af.spectra ) - t( spectra ) %*% v_lib_af    # D   x nAF
  r_dots_af <- pmax( colSums( r_lib_af^2 ), 1e-10 )           # nAF

  # Covariance-propagated AF weights
  af.cov.mat <- P %*% stats::cov( af.spectra ) %*% t( P )
  w_af       <- sqrt( abs( diag( af.cov.mat ) ) ) + 1e-8     # length F

  # -------------------------------------------------------------------------
  # Section 2: Per-fluorophore pre-computations  (mirrors C++ Section 2)
  # -------------------------------------------------------------------------
  af.only <- is.null( spectra.variants ) || length( spectra.variants ) == 0

  precomp      <- list()
  active.names <- character( 0 )

  if ( !af.only ) {

    variants      <- spectra.variants$variants
    delta.list.in <- spectra.variants$delta.list
    pos.thresholds <- rep( Inf, F )
    names( pos.thresholds ) <- fluorophores
    if ( !is.null( spectra.variants$thresholds ) )
      pos.thresholds[ names( spectra.variants$thresholds ) ] <- spectra.variants$thresholds

    # respect optimize.recommended filter if present
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
      U_nof    <- tryCatch(
        solve( S_nof %*% t( S_nof ), S_nof ),
        error = function( e ) {
          # Moore-Penrose fallback if S_nof S_nof^T is singular
          sv    <- svd( S_nof )
          tol   <- max( dim( S_nof ) ) * .Machine$double.eps * sv$d[ 1 ]
          d.inv <- ifelse( sv$d > tol, 1 / sv$d, 0 )
          t( sv$v %*% diag( d.inv, length( d.inv ) ) %*% t( sv$u ) )
        }
      )

      v_lib <- U_nof %*% t( delta )   # (F-1) x n_variants

      # covariance-propagated leakage weights: ridge 1e-4 matches C++ Section 2
      delta.obs <- if ( !is.null( delta.list.in[[ fl ]] ) )
        delta.list.in[[ fl ]] else delta
      delta.cov        <- if ( nrow( delta.obs ) > 1 ) stats::cov( delta.obs )
                          else t( delta.obs ) %*% delta.obs / max( nrow( delta.obs ) - 1, 1 )
      diag( delta.cov ) <- diag( delta.cov ) + 1e-4
      leakage.cov      <- U_nof %*% delta.cov %*% t( U_nof )
      w.leakage        <- sqrt( abs( diag( leakage.cov ) ) ) + 1e-8  # length F-1

      # rank-1 residual helpers: r_lib(:,v) = delta(v) - S^T P delta(v)
      r_lib  <- t( delta ) - t( spectra ) %*% ( P %*% t( delta ) )   # D x n_var
      r_dots <- pmax( colSums( r_lib^2 ), 1e-12 )                     # n_var

      precomp[[ fl ]] <- list(
        master.idx = which( fluorophores == fl ),
        other.fl   = other.fl,
        v.mats     = v.mats,
        n.var      = n.var,
        v_lib      = v_lib,
        w.leakage  = w.leakage,
        r_lib      = r_lib,
        r_dots     = r_dots,
        pos.thresh = pos.thresholds[ fl ]
      )
      active.names <- c( active.names, fl )
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
  # Section 3: Per-cell parallel loop  (mirrors C++ Section 3)
  # -------------------------------------------------------------------------

  # resolve threads
  if ( parallel ) {
    if ( is.null( threads ) ) threads <- asp$worker.process.n
    if ( !is.null( threads ) && threads == 0 )
      threads <- parallelly::availableCores()
  } else {
    threads <- 1L
  }

  # objects the worker function closes over: listed explicitly for export
  worker.exports <- c(
    "raw.data", "spectra", "af.spectra",
    "F", "D", "nAF", "fluorophores",
    "P", "v_lib_af", "r_lib_af", "r_dots_af", "w_af",
    "af.only", "precomp", "active.names",
    "n.passes", "n.variants"
  )

  result.setup <- create.parallel.lapply(
    asp,
    exports    = worker.exports,
    parallel   = parallel,
    threads    = threads,
    export.env = environment()
  )
  lapply.function <- result.setup$lapply

  # per-cell worker: closed over all pre-computed objects
  process.cell <- function( i ) {

    cell.raw <- raw.data[ i, , drop = TRUE ]

    # -----------------------------------------------------------------
    # A. AF SELECTION
    # -----------------------------------------------------------------
    init.f        <- as.numeric( P %*% cell.raw )
    base.resid.af <- cell.raw - as.numeric( t( spectra ) %*% init.f )
    base.resid.n  <- max( sqrt( sum( base.resid.af^2 ) ), 1e-8 )
    base.fluor.l1 <- max( sum( w_af * abs( init.f ) ), 1e-8 )

    best.score.af <- Inf
    best.j.af     <- 1L

    for ( j in seq_len( nAF ) ) {
      k.j     <- max( 0, sum( cell.raw * r_lib_af[ , j ] ) / r_dots_af[ j ] )
      r.j     <- base.resid.af - k.j * r_lib_af[ , j ]
      p.resid <- sqrt( sum( r.j^2 ) ) / base.resid.n
      p.fluor <- sum( w_af * abs( init.f - k.j * v_lib_af[ , j ] ) ) / base.fluor.l1
      score   <- p.resid * p.fluor
      if ( score < best.score.af ) {
        best.score.af <- score
        best.j.af     <- j
      }
    }

    # -----------------------------------------------------------------
    # B. INITIAL JOINT SOLVE
    # -----------------------------------------------------------------
    cell.S     <- rbind( spectra, af.spectra[ best.j.af, , drop = FALSE ] )
    StS        <- t( cell.S ) %*% cell.S
    uf         <- as.numeric( solve( StS, t( cell.S ) %*% cell.raw ) )
    fl.unm     <- uf[ seq_len( F ) ]
    resid      <- cell.raw - as.numeric( t( cell.S ) %*% uf )

    if ( af.only || length( active.names ) == 0 ) {
      return( c( fl.unm, uf[ F + 1L ], best.j.af ) )
    }

    best.v <- stats::setNames( rep( -1L, length( active.names ) ), active.names )

    # -----------------------------------------------------------------
    # C. JOINT VARIANT SELECTION  (n.passes passes)
    # -----------------------------------------------------------------
    for ( pass in seq_len( n.passes ) ) {

      rss.curr      <- max( sum( resid^2 ), 1e-12 )
      rss.curr.sqrt <- sqrt( rss.curr )

      # score all candidate swaps across all active fluorophores
      candidates <- list()

      for ( fl in active.names ) {
        pc    <- precomp[[ fl ]]
        abund <- fl.unm[ pc$master.idx ]
        if ( abund < pc$pos.thresh ) next

        other.unm <- fl.unm[ -pc$master.idx ]
        base.leak <- max( sum( pc$w.leakage * abs( other.unm ) ), 1e-8 )
        cur.v     <- best.v[ fl ]

        # limit candidates to top n.variants by residual-update score
        resid.scores <- numeric( pc$n.var )
        for ( v in seq_len( pc$n.var ) ) {
          dr       <- if ( cur.v < 0L ) pc$r_lib[ , v ]
                      else pc$r_lib[ , v ] - pc$r_lib[ , cur.v ]
          cross    <- sum( resid * dr )
          dr.sq    <- sum( dr^2 )
          new.rss  <- max( rss.curr - 2 * abund * cross + abund^2 * dr.sq, 0 )
          resid.scores[ v ] <- sqrt( new.rss )
        }
        top.v <- order( resid.scores )[ seq_len( min( n.variants, pc$n.var ) ) ]

        for ( v in top.v ) {
          if ( cur.v < 0L ) {
            dr         <- pc$r_lib[ , v ]
            dv.leakage <- pc$v_lib[ , v ]
          } else {
            dr         <- pc$r_lib[ , v ] - pc$r_lib[ , cur.v ]
            dv.leakage <- pc$v_lib[ , v ] - pc$v_lib[ , cur.v ]
          }

          cross     <- sum( resid * dr )
          dr.sq     <- sum( dr^2 )
          new.rss   <- max( rss.curr - 2 * abund * cross + abund^2 * dr.sq, 0 )

          resid.ratio   <- sqrt( new.rss ) / rss.curr.sqrt
          new.other     <- other.unm - abund * dv.leakage
          leakage.ratio <- sum( pc$w.leakage * abs( new.other ) ) / base.leak
          composite     <- resid.ratio * leakage.ratio

          if ( composite < 1.0 )
            candidates[[ length( candidates ) + 1L ]] <- list(
              score = composite, fl = fl, v = v
            )
        }
      }

      if ( length( candidates ) == 0L ) break

      # sort candidates best-first
      scores.ord <- order( sapply( candidates, `[[`, "score" ) )
      candidates <- candidates[ scores.ord ]

      # conflict resolution: commit non-overlapping swaps
      # two swaps conflict when their residual deltas have cosine similarity > 0.5
      committed     <- stats::setNames( rep( FALSE, length( active.names ) ), active.names )
      committed.drs <- list()
      commits       <- list()

      for ( cand in candidates ) {
        fl <- cand$fl
        if ( committed[ fl ] ) next

        pc    <- precomp[[ fl ]]
        abund <- fl.unm[ pc$master.idx ]
        cur.v <- best.v[ fl ]
        v     <- cand$v

        dr      <- if ( cur.v < 0L ) pc$r_lib[ , v ] * abund
                   else ( pc$r_lib[ , v ] - pc$r_lib[ , cur.v ] ) * abund
        dr.norm <- max( sqrt( sum( dr^2 ) ), 1e-12 )

        conflict <- FALSE
        for ( cd in committed.drs ) {
          cosine <- abs( sum( dr * cd ) ) /
            ( dr.norm * max( sqrt( sum( cd^2 ) ), 1e-12 ) )
          if ( cosine > 0.5 ) { conflict <- TRUE; break }
        }
        if ( conflict ) next

        committed[ fl ]    <- TRUE
        committed.drs[[ length( committed.drs ) + 1L ]] <- dr
        commits[[ length( commits ) + 1L ]] <- list( fl = fl, v = v )
      }

      if ( length( commits ) == 0L ) break

      # apply committed swaps and re-solve once per pass
      for ( cm in commits ) {
        best.v[ cm$fl ]                       <- cm$v
        cell.S[ precomp[[ cm$fl ]]$master.idx, ] <- precomp[[ cm$fl ]]$v.mats[ cm$v, ]
      }

      StS    <- t( cell.S ) %*% cell.S
      uf     <- as.numeric( solve( StS, t( cell.S ) %*% cell.raw ) )
      fl.unm <- uf[ seq_len( F ) ]
      resid  <- cell.raw - as.numeric( t( cell.S ) %*% uf )

      # early exit if residual is negligible
      if ( sqrt( sum( resid^2 ) ) < 1e-8 * sqrt( sum( cell.raw^2 ) ) ) break
    }

    return( c( fl.unm, uf[ F + 1L ], best.j.af ) )
  }

  # run across all cells
  cell.results <- tryCatch(
    expr = {
      lapply.function( seq_len( N ), process.cell )
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
