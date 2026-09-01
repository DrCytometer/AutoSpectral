# correct_spectra_glasso.R

#' @title Correct Spectra by Graphical Lasso
#'
#' @description
#' Estimates the residual spillover matrix and re-measures fluorophore
#' signatures from a fully stained sample, the same task `fix.my.unmix()`
#' performs, but replaces its pair-at-a-time negative-envelope estimator with
#' a joint, L1-penalised regression: for every target fluorophore, its
#' apparent (compensated) abundance among events negative for it is regressed
#' on every other fluorophore's abundance at once, with an L1 penalty that
#' drives the coefficient of an uninvolved fluorophore to exactly zero. This
#' is neighbourhood selection (Meinshausen & Buhlmann 2006) applied one row at
#' a time to the residual spillover matrix rather than to a covariance matrix,
#' so the returned coefficients are directly the row of `spillover` that
#' `fix.my.unmix()` estimates pair by pair.
#'
#' The two functions share an identification strategy and diverge in how they
#' handle co-expression. Both restrict the fit to events negative for the
#' target, since a marker-negative population must read zero in every channel
#' it is negative for - whatever a negative population's apparent abundance
#' correlates with there is spillover, not real expression, because the cells
#' contributing it do not express the target at all. `fix.my.unmix()` then
#' estimates one source at a time and relies on `source.dominant` masking to
#' stop a second bright fluorophore's own spillover from being booked against
#' the source under test. This function instead regresses the target
#' jointly on every other fluorophore's abundance in the same fit, so a
#' source that only looks correlated with the target because it is itself
#' correlated with the true culprit has that correlation absorbed by the true
#' culprit's own coefficient, and the L1 penalty turns the "is this source
#' really contributing" question into automatic variable selection instead of
#' a hand-tuned dominance heuristic. The trade-off is that pairs which are
#' both spectrally similar and genuinely co-expressed remain as hard for this
#' estimator to separate as for the pairwise one; nothing about fitting them
#' jointly changes what two collinear predictors can and cannot identify.
#'
#' As in `fix.my.unmix()`, nothing is back-solved out of the spillover matrix.
#' The matrix is used only to refine which events count as negative for each
#' target; the spectra themselves are re-measured directly from background
#' subtracted raw data by `extract.raw.signature()`, exactly as
#' `fix.my.unmix()`'s second phase does, with one difference: the
#' co-fluorophores carried into that fit (`active`) are the ones this
#' function's lasso selected for that target, not the whole panel run through
#' a ridge penalty. A pair the lasso found no evidence of coupling for is left
#' out of the signature fit entirely, rather than being included and shrunk.
#'
#' A coefficient's sign carries the identification, and it is not symmetric
#' in target and source. A marker-negative population cannot read below zero
#' from real emission - nothing biological subtracts signal - so where a
#' positive coefficient is confounded with genuine co-expression (see
#' `CONTEXT_information_sources.md`'s hypernegativity argument), a negative
#' one is unambiguous, and it identifies a different row than a positive one
#' would. Corrupting fluorophore j's reference row with eps of fluorophore
#' k's shape (`S.wrong[j,] <- normalize(S.true[j,] + eps * S.true[k,])`)
#' leaves the pair fit the way one might expect - k's population regressed on
#' j - reading essentially zero, because unmixing is a change of basis and
#' the two rows actually coupled by that change are the other way round: j's
#' own population, regressed on k, reads close to `-eps`. The error is in row
#' j, but its fingerprint is a negative coefficient found while fitting row
#' k, not a positive one found while fitting row j's own. Every target's
#' per-row lasso fit already estimates every cell of this matrix correctly,
#' once every fluorophore has been run as a target in turn - the fitted
#' matrix itself is not wrong - but a target's active set, if it were read
#' only from its own row, would never see a negative coefficient discovered
#' while fitting someone else's row. `active.set` folds this in explicitly:
#' for every negative coefficient found anywhere in the fitted matrix, it is
#' the fluorophore whose row produced the fit (the source, not the target of
#' that particular fit) whose reference shape needs the fitted target added
#' to its own signature re-measurement.
#'
#' @section How the lasso row is estimated:
#' For target \eqn{j} and its target-negative event set \eqn{N_j}, let
#' \eqn{y} be \eqn{j}'s compensated abundance over \eqn{N_j} and \eqn{X} be
#' every other fluorophore's abundance over the same events. The fit solves
#' \deqn{ \hat\beta_j = \arg\min_\beta \tfrac{1}{2n}\lVert y - X\beta -
#'   \alpha \rVert_2^2 + \lambda \sum_{k \ne j} |\beta_k| }
#' by cyclic coordinate descent along a path of decreasing \eqn{\lambda}, with
#' the intercept \eqn{\alpha} unpenalised so it can absorb a constant
#' background offset the way the pairwise estimator's anchor bin does.
#' \eqn{\lambda} is chosen per target by a two-way split-half holdout: fit on
#' one half, score squared error on the other, and back; the two error curves
#' are summed over a shared \eqn{\lambda} grid, and the largest \eqn{\lambda}
#' (sparsest fit) within `margin.frac` of the minimum is kept, then refit on
#' the full target-negative set at that \eqn{\lambda}. Preferring the sparser
#' of two disagreeing halves is the same safe-direction bias
#' `fix.my.unmix()`'s target-negative truncation uses when co-expression
#' would otherwise shrink a slope toward zero.
#'
#' In lay terms: instead of asking "does fluorophore A leak into fluorophore
#' B's channel, checked one pair at a time," this asks fluorophore B's whole
#' negative population "which of everyone else in the panel actually explains
#' the small amount of signal you have here," with a rule that only keeps an
#' answer if it clearly earns its place - a fluorophore whose apparent
#' contribution could just as well be noise, or could just as well be
#' somebody else's spillover showing up correlated with it, gets a flat zero
#' rather than a small uncertain number. And once an answer does earn its
#' place, which row it corrects depends on its sign: a fluorophore that
#' looks brighter than it should, in step with something else in the panel,
#' has real spillover to remove from its own row; a fluorophore that looks
#' dimmer than it should, in step with something else, cannot really be
#' losing signal to biology, so what is wrong is not its own row but the
#' other fluorophore's reference shape, which has borrowed a slice of this
#' one's.
#'
#' @param spectra The spectral matrix, fluorophores x detectors, L-infinity
#'   normalised. Ideally the output of `correct.unmixing.signatures()`.
#' @param unstained.sample File path and name for a raw unstained sample,
#'   acquired the same day and matching the autofluorescence of the fully
#'   stained sample.
#' @param fully.stained.sample File path and name for a raw fully stained
#'   sample.
#' @param flow.control The flow.control list.
#' @param asp The AutoSpectral parameter list.
#' @param variants The variant list returned by `get.spectral.variants()`.
#'   Only `variants$spillover.spread` is used, to set abundance-dependent
#'   positivity boundaries. If `NULL`, boundaries are flat. Default `NULL`.
#' @param af.name Character or `NULL`, the name of an autofluorescence row in
#'   `spectra`. Never treated as a panel fluorophore and never corrected.
#'   Default `"AF"`.
#' @param af.basis Optional matrix (components x detectors) from
#'   `get.af.basis()`. When `NULL` and `bg.mode = "af.deconv"`, it is built
#'   from the unstained sample. Default `NULL`.
#' @param af.n.pc Integer or `"auto"`, passed to `get.af.basis()`. Default
#'   `"auto"`.
#' @param bg.mode Character. `"af.deconv"` (default) fits a multi-component
#'   autofluorescence basis jointly with the panel; `"af.row"` uses the
#'   single `af.name` row already in `spectra`; `"global.mean"` uses the mean
#'   unstained spectrum as a single background row; `"none"` fits no
#'   background.
#' @param large.gate Logical, whether to use a large scatter gate. Default
#'   `TRUE`.
#' @param downsample Logical or numeric. `FALSE` disables downsampling; a
#'   numeric gives the number of events to use, stratified by dominant
#'   fluorophore under the starting spectra. Default `20000`.
#' @param downsample.background.frac Numeric in (0, 1), the share of
#'   `downsample` reserved for events dominant for nothing. Default `0.3`.
#' @param downsample.min.stratum Integer, the floor below which a
#'   fluorophore's own positive population is kept whole. Default `2000`.
#' @param unstained.threshold Numeric in (0, 1), the percentile of the
#'   unstained control defining positivity. Default `0.99`.
#' @param unstained.margin Numeric, multiplier applied to that threshold.
#'   Default `1.3`.
#' @param spread.kappa Numeric, how many spillover-spread standard
#'   deviations above the flat threshold still count as negative. Default
#'   `2`.
#' @param min.negative.events Integer, the fewest target-negative events a
#'   fluorophore must have before its row is fit at all; below this the
#'   fluorophore is left uncorrected for this iteration. Default `200`.
#' @param max.truncated.events Integer, cap on the target-negative events
#'   used per lasso fit; above the cap the population is subsampled, since
#'   the fit's cost scales with it directly. Default `20000`.
#' @param max.mask.passes Integer, how many times the target-negative
#'   selection is recomputed with the current fitted row's joint spillover
#'   contribution removed. Default `3`.
#' @param mask.tolerance Numeric, the fraction of the target-negative set
#'   that may change between mask passes before the mask is treated as
#'   settled. Default `0.05`.
#' @param min.span Numeric, minimum abundance span a selected source must
#'   have across the target-negative population, in units of the source's
#'   own flat threshold, or its coefficient is discarded. Default `5`.
#' @param min.rise Numeric, the fitted rise a selected source's coefficient
#'   must produce across its own span, in standard deviations of the
#'   target's negative population, or the coefficient is discarded. Default
#'   `1`.
#' @param max.coefficient Numeric, the largest residual spillover
#'   coefficient accepted; a larger fitted value is discarded rather than
#'   kept, on the same small-error premise `fix.my.unmix()` uses. Default
#'   `0.5`.
#' @param coefficient.tolerance Numeric, fitted coefficients smaller than
#'   this in absolute value are treated as exactly zero. Default `1e-4`.
#' @param n.lambda Integer, points on the lasso regularisation path. Default
#'   `40`.
#' @param lambda.min.ratio Numeric in (0, 1), the smallest path value as a
#'   fraction of `lambda.max`, the value at which every coefficient is
#'   already zero. Default `1e-3`.
#' @param margin.frac Numeric, how far above the minimum two-way holdout
#'   error a larger, sparser `lambda` may sit and still be preferred.
#'   Default `0.05`.
#' @param max.iter Integer, maximum outer spillover-matrix iterations.
#'   Default `5`.
#' @param convergence.threshold Numeric, residual spillover coefficient at
#'   which iteration stops. Default `0.01`.
#' @param convergence.quantile Numeric, the quantile of the off-diagonal
#'   coefficients the convergence test uses. Default `0.95`.
#' @param step.spillover Numeric in (0, 1], the fraction of each outer
#'   iteration's fitted spillover update applied. The lasso already shrinks
#'   and selects each row, so unlike `fix.my.unmix()`'s per-coefficient
#'   `trust` weighting this is a single damping factor for the whole update,
#'   only useful if the outer iteration oscillates. Default `1`.
#' @param update.spectra Logical, whether to run the raw-space signature
#'   phase. When `FALSE`, only the spillover and compensation matrices are
#'   returned. Default `TRUE`.
#' @param step Numeric, the fraction of each accepted signature change
#'   applied. Default `1`.
#' @param null.fit Optional, the result of running this function on the
#'   control the reference spectra were extracted from, where it should
#'   return exactly what it was given. Its bias is subtracted from both the
#'   spillover matrix and the accepted signatures. Default `NULL`.
#' @param max.resid.ratio Numeric, multiple of the null run's median
#'   `resid.rel` above which a fit is refused. Ignored without `null.fit`.
#'   Default `3`.
#' @param max.intercept.ratio Numeric, the same for `intercept.rel`. Default
#'   `3`.
#' @param intercept Logical, whether the signature fit carries an intercept.
#'   Default `TRUE`.
#' @param min.explained Numeric, minimum `explained.total`. Default `0.8`.
#' @param max.explained Numeric, maximum `explained.total`. Default `1.2`.
#' @param max.resid Numeric, maximum relative fit residual. Overridden by
#'   `null.fit`. Default `0.03`.
#' @param max.intercept Numeric, maximum relative intercept. Overridden by
#'   `null.fit`. Default `0.03`.
#' @param min.bg.align Numeric, minimum cosine between the fitted intercept
#'   and the candidate signature. Default `-0.9`.
#' @param max.clamp.frac Numeric, maximum fraction of a candidate row's
#'   absolute mass non-negativity clamping may remove. Default `0.15`.
#' @param max.anchor Numeric, maximum unremoved background relative to the
#'   brightest fitted signal. Default `0.10`.
#' @param max.vif Numeric, maximum variance inflation factor for the
#'   fluorophore's own abundance within its population. Default `500`.
#' @param max.condition.increase Numeric, the factor by which a single
#'   accepted row may increase the condition number of the unmixing design.
#'   Default `1.05`.
#' @param peak.shift.min.rel Numeric, see `extract.raw.signature()`. Default
#'   `0.7`.
#' @param max.angle Numeric, maximum angular change of a row in degrees.
#'   Default `15`.
#' @param max.hotspot Numeric, hotspot scale above which a fluorophore is
#'   frozen as inseparable from the autofluorescence basis. Default `5`.
#' @param leakage.margin Numeric, the fraction by which held-out leakage
#'   must increase before a candidate is refused for it. Default `0.05`.
#' @param n.levels Integer, maximum abundance bins for the signature fit,
#'   passed to `extract.raw.signature()`. Default `60`.
#' @param min.bin.events Integer, fewest events per abundance bin, passed to
#'   `extract.raw.signature()`. Default `50`.
#' @param multivariate Logical, passed to `extract.raw.signature()`. Default
#'   `TRUE`.
#' @param ridge Numeric, ridge penalty for that joint fit. Default `1e-6`.
#' @param output.suffix Character, appended to the csv and figure filenames
#'   so a run of this function does not overwrite `fix.my.unmix()`'s output
#'   in the same directory. Default `"_glasso"`.
#' @param figures Logical, whether to write the spillover heatmap. Default
#'   `TRUE`.
#' @param save Logical, whether to write the csv outputs. Default `TRUE`.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`spectra`}{The corrected spectra, fluorophores x detectors. Rows
#'     that failed a gate are unchanged.}
#'   \item{`spectra.backsolved`}{`spillover %*% spectra`, the algebraic
#'     back-solve. Diagnostic only.}
#'   \item{`spillover`}{The estimated residual spillover matrix, fluorophores
#'     x fluorophores. Read `spillover[a, b]` as the estimated contribution
#'     of fluorophore a's true abundance to fluorophore b's apparent
#'     abundance under the starting spectra; it is filled cell by cell while
#'     each fluorophore in turn is the fit's target, and it is not
#'     symmetric.}
#'   \item{`compensation`}{Its inverse.}
#'   \item{`coefficient.log`}{Per-(source, target) pair, one row per pair,
#'     from the final iteration: the fitted coefficient (`0` if the lasso
#'     dropped it or a gate discarded it) and whether it was selected. Raw
#'     fit output; a negative `beta` here means the *source* row is the one
#'     implicated, not the target - see `active.set`.}
#'   \item{`lambda.log`}{Per-target diagnostics from the final iteration:
#'     events used, candidate sources, sources selected, and the chosen
#'     `lambda`.}
#'   \item{`signature.log`}{Per-fluorophore signature statistics and gate
#'     outcomes, from `extract.raw.signature()`.}
#'   \item{`convergence.log`}{Per-iteration delta history.}
#'   \item{`active.set`}{Named list, the fluorophores carried into each
#'     target's signature fit: itself, whichever sources its own row
#'     selected with a positive coefficient (genuine spillover into it,
#'     which must be modelled to isolate its true signal), and whichever
#'     other fluorophores' own fits named it as a negative-coefficient
#'     source (evidence its own reference shape needs that fluorophore
#'     modelled out too). See `reversed.donors` for the second route on its
#'     own.}
#'   \item{`reversed.donors`}{Named list, per fluorophore, the donors folded
#'     into its `active.set` by the negative-coefficient route: every other
#'     fluorophore whose own row-fit found this one as a negative-
#'     coefficient candidate source.}
#'   \item{`af.basis`, `af.hotspot`, `af.frozen`}{The autofluorescence basis,
#'     its coupling to the panel, and the fluorophores frozen because of it.}
#' }
#'
#' @importFrom sp point.in.polygon
#' @importFrom stats approx mad median quantile sd setNames
#'
#' @export

correct.spectra.glasso <- function(
    spectra,
    unstained.sample,
    fully.stained.sample,
    flow.control,
    asp,
    variants                    = NULL,
    af.name                     = "AF",
    af.basis                    = NULL,
    af.n.pc                     = "auto",
    bg.mode                     = c( "af.deconv", "af.row", "global.mean", "none" ),
    large.gate                  = TRUE,
    downsample                  = 20000,
    downsample.background.frac  = 0.3,
    downsample.min.stratum      = 2000L,
    unstained.threshold         = 0.99,
    unstained.margin            = 1.3,
    spread.kappa                = 2,
    min.negative.events         = 200L,
    max.truncated.events        = 20000L,
    max.mask.passes             = 3L,
    mask.tolerance               = 0.05,
    min.span                     = 5,
    min.rise                     = 1,
    max.coefficient               = 0.5,
    coefficient.tolerance         = 1e-4,
    n.lambda                      = 40L,
    lambda.min.ratio               = 1e-3,
    margin.frac                    = 0.05,
    max.iter                       = 5L,
    convergence.threshold          = 0.01,
    convergence.quantile           = 0.95,
    step.spillover                 = 1,
    update.spectra                 = TRUE,
    step                           = 1,
    null.fit                       = NULL,
    max.resid.ratio                = 3,
    max.intercept.ratio            = 3,
    intercept                      = TRUE,
    min.explained                  = 0.8,
    max.explained                  = 1.2,
    max.resid                      = 0.03,
    max.intercept                  = 0.03,
    min.bg.align                   = -0.9,
    max.clamp.frac                 = 0.15,
    max.anchor                     = 0.10,
    max.vif                        = 500,
    max.condition.increase         = 1.05,
    peak.shift.min.rel             = 0.7,
    max.angle                      = 15,
    max.hotspot                    = 5,
    leakage.margin                 = 0.05,
    n.levels                       = 60L,
    min.bin.events                 = 50L,
    multivariate                   = TRUE,
    ridge                          = 1e-6,
    output.suffix                  = "_glasso",
    figures                        = TRUE,
    save                           = TRUE,
    verbose                        = TRUE
) {

  bg.mode <- match.arg( bg.mode )

  # Set once for the whole run, matching fix.my.unmix(): the lasso path's
  # own subsampling and the downsample below both draw on it.
  if ( !is.null( asp$bird.seed ) ) set.seed( asp$bird.seed )

  null.spillover <- null.fit$spillover
  null.spectra   <- null.fit$spectra

  if ( !is.null( null.fit$signature.log ) ) {

    max.resid <- max.resid.ratio * stats::median(
      null.fit$signature.log$resid.rel, na.rm = TRUE )

    max.intercept <- max.intercept.ratio * stats::median(
      null.fit$signature.log$intercept.rel, na.rm = TRUE )
  }

  spectra <- as.matrix( spectra )

  if ( is.null( rownames( spectra ) ) )
    stop( "`spectra` must have fluorophore row names.", call. = FALSE )

  fluorophores  <- setdiff( rownames( spectra ), af.name )
  fluorophore.n <- length( fluorophores )

  if ( fluorophore.n < 2 )
    stop( "At least two panel fluorophores are required.", call. = FALSE )

  if ( save && ! dir.exists( asp$fix.unmixing.dir ) )
    dir.create( asp$fix.unmixing.dir, recursive = TRUE )

  spillover.spread <- variants$spillover.spread

  if ( is.null( spillover.spread ) )
    warning( paste0( "No `spillover.spread` in `variants`; positivity ",
                     "boundaries will be flat." ), call. = FALSE )

  # ---------------------------------------------------------------------------
  # Read and gate - identical to fix.my.unmix()
  # ---------------------------------------------------------------------------

  read.gated <- function( file.name, gate.polygon = NULL, label ) {

    if ( verbose )
      message( sprintf( "\033[34mReading %s.\033[0m", label ) )

    expr.data <- readFCS( file.name, columns = flow.control$scatter.and.channel.spectral )
    gate.data <- expr.data[ , flow.control$scatter.parameter ]

    if ( is.null( gate.polygon ) )
      gate.polygon <- do.gate(
        gate.data, viability.gate = FALSE, large.gate = large.gate,
        samp = label,
        scatter.and.channel.label = flow.control$scatter.and.channel.label,
        control.type = "cells", asp )

    keep <- which( sp::point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ],
                                         gate.polygon$x, gate.polygon$y ) != 0 )

    list( data = expr.data[ keep, flow.control$spectral.channel, drop = FALSE ],
          gate = gate.polygon )
  }

  unstained.in <- read.gated( unstained.sample, NULL, "unstained raw" )
  stained.in   <- read.gated( fully.stained.sample, unstained.in$gate,
                              "fully stained raw" )

  unstained.raw <- unstained.in$data
  stained.raw   <- stained.in$data

  if ( ncol( stained.raw ) != ncol( spectra ) )
    stop( "`spectra` and the spectral channels of the data do not match.",
          call. = FALSE )

  # ---------------------------------------------------------------------------
  # Background basis and the projection it defines - identical to fix.my.unmix()
  # ---------------------------------------------------------------------------

  spectra.fluor <- spectra[ fluorophores, , drop = FALSE ]

  background.basis <- switch(
    bg.mode,
    af.deconv = {
      if ( is.null( af.basis ) )
        af.basis <- get.af.basis( unstained.raw, n.pc = af.n.pc,
                                  verbose = verbose )
      as.matrix( af.basis )
    },
    af.row = {
      if ( is.null( af.name ) || ! af.name %in% rownames( spectra ) )
        stop( "`bg.mode = \"af.row\"` requires an `af.name` row in `spectra`.",
              call. = FALSE )
      spectra[ af.name, , drop = FALSE ]
    },
    global.mean = {
      bg <- colMeans( unstained.raw )
      matrix( bg / max( sqrt( sum( bg^2 ) ), .Machine$double.eps ), nrow = 1,
              dimnames = list( "BG", colnames( unstained.raw ) ) )
    },
    none = NULL
  )

  design <- if ( is.null( background.basis ) ) spectra.fluor else
    rbind( background.basis, spectra.fluor )

  if ( nrow( design ) > ncol( design ) )
    stop( "The background basis plus the panel exceeds the detector count.",
          call. = FALSE )

  af.hotspot <- calculate.hotspot.matrix( design )

  af.frozen <- character( 0 )
  if ( !is.null( background.basis ) ) {
    coupling  <- apply(
      af.hotspot[ fluorophores, rownames( background.basis ), drop = FALSE ],
      1, max )
    af.frozen <- names( coupling )[ coupling > max.hotspot ]
  }

  if ( verbose && length( af.frozen ) > 0 )
    message( sprintf(
      "\033[31mFrozen, inseparable from the background basis: %s.\033[0m",
      paste( af.frozen, collapse = ", " ) ) )

  project <- function( y ) {

    coefs <- unmix.ols.fast( y, design )
    colnames( coefs ) <- rownames( design )

    background <- 0
    if ( !is.null( background.basis ) ) {
      b <- coefs[ , rownames( background.basis ), drop = FALSE ]
      b[ b < 0 ] <- 0
      background <- b %*% background.basis
    }

    list( abundance = coefs[ , fluorophores, drop = FALSE ],
          residual  = y - background )
  }

  stained.fit   <- project( stained.raw )
  unstained.fit <- project( unstained.raw )

  if ( !isFALSE( downsample ) && !is.null( downsample ) ) {

    n.keep <- min( as.integer( downsample ), nrow( stained.raw ) )

    if ( n.keep < nrow( stained.raw ) ) {

      thr.rough <- unstained.margin * apply(
        unstained.fit$abundance, 2, stats::quantile,
        probs = unstained.threshold, names = FALSE )
      names( thr.rough ) <- fluorophores

      thr.mat.rough <- get.spread.thresholds(
        unmixed          = stained.fit$abundance,
        thresholds       = thr.rough,
        spillover.spread = spillover.spread,
        spread.kappa     = spread.kappa,
        verbose          = FALSE )

      keep.idx <- .fix.stratified.sample(
        abundance        = stained.fit$abundance,
        thresholds       = thr.rough,
        threshold.matrix = thr.mat.rough,
        n.total          = n.keep,
        background.frac  = downsample.background.frac,
        min.stratum      = downsample.min.stratum )

      stained.raw <- stained.raw[ keep.idx, , drop = FALSE ]
      stained.fit <- list(
        abundance = stained.fit$abundance[ keep.idx, , drop = FALSE ],
        residual  = stained.fit$residual[  keep.idx, , drop = FALSE ] )
    }
  }

  # ---------------------------------------------------------------------------
  # Phase one: residual spillover matrix by per-target neighbourhood lasso
  # ---------------------------------------------------------------------------

  spillover.curr <- diag( fluorophore.n )
  dimnames( spillover.curr ) <- list( fluorophores, fluorophores )

  beta.by.target <- vector( "list", fluorophore.n )
  names( beta.by.target ) <- fluorophores

  convergence.log <- data.frame(
    iter = integer(), delta = numeric(), delta.quantile = numeric(),
    delta.max = numeric(), n.fitted = integer(), stringsAsFactors = FALSE )

  delta.history <- rep( NA_real_, 3L )
  pair.log      <- NULL
  lambda.log    <- NULL

  for ( iter in seq_len( as.integer( max.iter ) ) ) {

    compensation.curr <- tryCatch( solve( spillover.curr ),
                                   error = function( e ) NULL )

    if ( is.null( compensation.curr ) ) {
      warning( "Spillover matrix became singular; stopping iteration.",
               call. = FALSE )
      break
    }

    unmixed.comp   <- stained.fit$abundance   %*% compensation.curr
    unstained.comp <- unstained.fit$abundance %*% compensation.curr

    thresholds <- unstained.margin * apply(
      unstained.comp, 2, stats::quantile, probs = unstained.threshold,
      names = FALSE )
    names( thresholds ) <- fluorophores

    neg.var <- apply( unstained.comp, 2, stats::mad )^2
    names( neg.var ) <- fluorophores

    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed.comp,
      thresholds       = thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      verbose          = FALSE )

    marker.spillover <- diag( fluorophore.n )
    dimnames( marker.spillover ) <- list( fluorophores, fluorophores )

    pair.log   <- list()
    lambda.log <- list()

    for ( target in fluorophores ) {

      sources <- setdiff( fluorophores, target )

      neg.idx <- which( unmixed.comp[ , target ] <= threshold.matrix[ , target ] )

      row.fit  <- NULL
      beta.row <- NULL

      if ( length( neg.idx ) >= min.negative.events ) {

        for ( pass in seq_len( as.integer( max.mask.passes ) ) ) {

          row.fit <- .glasso.neighborhood.row(
            abundance         = unmixed.comp,
            target            = target,
            sources           = sources,
            neg.idx           = neg.idx,
            n.lambda          = n.lambda,
            lambda.min.ratio  = lambda.min.ratio,
            margin.frac       = margin.frac,
            max.events        = max.truncated.events )

          beta.pass <- row.fit$beta
          beta.pass[ abs( beta.pass ) < coefficient.tolerance ] <- 0
          beta.pass[ abs( beta.pass ) > max.coefficient ]       <- 0

          # Kept as soon as it is fitted, not only once the mask below is
          # judged settled: a pass whose *next* mask would fall below
          # `min.negative.events` still produced a legitimate fit on the
          # population it was given, and breaking out below must not throw
          # that fit away.
          beta.row <- beta.pass

          fitted.spillover <- as.numeric(
            unmixed.comp[ , sources, drop = FALSE ] %*% beta.pass )

          neg.idx.next <- which(
            ( unmixed.comp[ , target ] - fitted.spillover ) <=
              threshold.matrix[ , target ] )

          if ( length( neg.idx.next ) < min.negative.events ) break

          moved <- length( union( neg.idx, neg.idx.next ) ) -
            length( intersect( neg.idx, neg.idx.next ) )

          settled <- ( moved / max( length( neg.idx ), 1L ) ) <= mask.tolerance

          neg.idx <- neg.idx.next

          if ( settled ) break
        }
      }

      # Coefficients that survive the lasso and the coefficient cap still
      # have to earn their place the way a pairwise-estimated coefficient
      # does: enough abundance span in the source to be identifiable, and
      # enough fitted rise across that span to be distinguishable from the
      # target's own negative-population noise. This gate applies the same
      # bar to both signs - identifiability does not depend on which row a
      # coefficient will eventually be attributed to.
      if ( !is.null( beta.row ) ) {

        selected <- names( beta.row )[ beta.row != 0 ]

        for ( src in selected ) {

          span.src  <- diff( range( unmixed.comp[ neg.idx, src ] ) )
          noise.tgt <- sqrt( max( neg.var[ target ], 0 ) )

          identifiable <- span.src > min.span * abs( thresholds[ src ] ) &&
            abs( beta.row[ src ] ) * span.src >= min.rise * noise.tgt

          if ( !identifiable ) beta.row[ src ] <- 0
        }
      }

      beta.by.target[[ target ]] <- beta.row

      if ( !is.null( beta.row ) )
        marker.spillover[ sources, target ] <- beta.row[ sources ]

      n.selected <- if ( is.null( beta.row ) ) 0L else sum( beta.row != 0 )

      lambda.log[[ length( lambda.log ) + 1L ]] <- data.frame(
        target       = target,
        n.events     = length( neg.idx ),
        n.candidates = length( sources ),
        n.selected   = n.selected,
        lambda       = if ( is.null( row.fit ) ) NA_real_ else row.fit$lambda,
        row.names    = NULL, stringsAsFactors = FALSE )

      for ( src in sources )
        pair.log[[ length( pair.log ) + 1L ]] <- data.frame(
          source    = src,
          target    = target,
          beta      = if ( is.null( beta.row ) ) 0 else unname( beta.row[ src ] ),
          selected  = if ( is.null( beta.row ) ) FALSE else
            unname( beta.row[ src ] != 0 ),
          row.names = NULL, stringsAsFactors = FALSE )
    }

    slope.error <- marker.spillover - diag( fluorophore.n )

    off.diagonal <- abs( slope.error[ row( slope.error ) != col( slope.error ) ] )

    delta          <- stats::sd( slope.error )
    delta.quantile <- unname( stats::quantile(
      off.diagonal, probs = convergence.quantile, names = FALSE ) )
    delta.max      <- max( off.diagonal )

    n.fitted <- sum( vapply( beta.by.target, function( b )
      if ( is.null( b ) ) 0L else sum( b != 0 ), integer( 1 ) ) )

    convergence.log[ nrow( convergence.log ) + 1L, ] <-
      list( iter, delta, delta.quantile, delta.max, n.fitted )

    if ( verbose )
      message( sprintf(
        "\033[34miter %d: delta %.5f, q%.0f %.5f, max %.5f, %d coefficient(s) selected\033[0m",
        iter, delta, 100 * convergence.quantile, delta.quantile, delta.max,
        n.fitted ) )

    spillover.next <- spillover.curr + step.spillover * slope.error %*% spillover.curr

    diag.next <- diag( spillover.next )

    if ( any( !is.finite( diag.next ) ) || any( diag.next <= 0 ) ) {
      warning( "Spillover diagonal left the positive range; stopping iteration.",
               call. = FALSE )
      break
    }

    spillover.curr <- sweep( spillover.next, 1, diag.next, "/" )

    if ( delta.quantile < convergence.threshold ) {
      if ( verbose ) message( "\033[32mConverged.\033[0m" )
      break
    }

    delta.history <- c( delta.history[ -1L ], delta )

    if ( all( is.finite( delta.history ) ) &&
         mean( diff( delta.history ) ) >= 0 ) {
      if ( verbose ) message( "\033[33mRefinement stalled, stopping.\033[0m" )
      break
    }

    if ( iter == max.iter && verbose )
      message( "\033[33mReached iteration limit.\033[0m" )
  }

  coefficient.log <- do.call( rbind, pair.log )
  lambda.log      <- do.call( rbind, lambda.log )

  if ( !is.null( null.spillover ) ) {

    shared <- intersect( fluorophores, rownames( null.spillover ) )

    if ( length( shared ) > 1 ) {

      bias <- null.spillover[ shared, shared, drop = FALSE ] -
        diag( length( shared ) )

      spillover.curr[ shared, shared ] <-
        spillover.curr[ shared, shared ] - bias

      diag( spillover.curr ) <- 1
    }
  }

  compensation.curr <- solve( spillover.curr )

  # Each fluorophore's own row-fit only reports the coefficients that
  # explain ITS negative population; whether that population's own row is
  # itself the wrong one is invisible to its own fit. That evidence turns
  # up as a negative coefficient in whichever OTHER fluorophore's fit sees
  # this one as a candidate source: a marker-negative population cannot
  # read below zero from real emission, so a negative coefficient is never
  # evidence that the candidate spills into the target being fit - it is
  # evidence that the candidate's own reference row carries some of the
  # target's shape. That correction belongs on the candidate's row, not the
  # target's, so it is routed there before active sets are built.
  #
  # own.set[[j]]: j itself, plus every source j's own fit found a genuine
  #   (positive) spillover contribution from. A target the lasso could not
  #   fit at all (too few negative events throughout) falls back to the
  #   whole panel, extract.raw.signature()'s own default, rather than being
  #   left with no nuisance removal at all.
  # reversed.donors[[j]]: every fluorophore whose own fit found j as a
  #   negative-coefficient source - i.e. every fluorophore whose reference
  #   row is implicated as the true donor of whatever is wrong with j's.

  own.set <- lapply( fluorophores, function( j ) {
    b <- beta.by.target[[ j ]]
    if ( is.null( b ) ) fluorophores else c( j, names( b )[ b > 0 ] )
  } )
  names( own.set ) <- fluorophores

  reversed.donors <- stats::setNames(
    vector( "list", fluorophore.n ), fluorophores )

  for ( q in fluorophores ) {

    b <- beta.by.target[[ q ]]
    if ( is.null( b ) ) next

    corrupted.by.q <- names( b )[ b < 0 ]

    for ( j in corrupted.by.q )
      reversed.donors[[ j ]] <- c( reversed.donors[[ j ]], q )
  }

  # The active set each target's signature fit will carry is the union of
  # what its own row of the spillover matrix says it needs (genuine
  # spillover into it, from own.set) and what every other row's negative
  # coefficients say about it (evidence its own reference shape is
  # contaminated, from reversed.donors).
  active.set <- lapply( fluorophores, function( j )
    unique( c( own.set[[ j ]], reversed.donors[[ j ]] ) ) )
  names( active.set ) <- fluorophores

  # ---------------------------------------------------------------------------
  # Phase two: signatures re-measured in raw space
  # ---------------------------------------------------------------------------

  spectra.new   <- spectra
  signature.log <- NULL

  if ( update.spectra ) {

    unmixed.comp <- stained.fit$abundance %*% compensation.curr
    residual     <- stained.fit$residual

    thresholds <- unstained.margin * apply(
      unstained.fit$abundance %*% compensation.curr, 2, stats::quantile,
      probs = unstained.threshold, names = FALSE )
    names( thresholds ) <- fluorophores

    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed.comp,
      thresholds       = thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      verbose          = FALSE )

    dyn.range <- apply( unmixed.comp, 2, stats::quantile, probs = 0.999,
                        names = FALSE ) - thresholds
    dyn.range <- pmax( dyn.range, .Machine$double.eps )

    score    <- sweep( pmax( unmixed.comp - threshold.matrix, 0 ), 2,
                       dyn.range, "/" )
    dominant <- max.col( score, ties.method = "first" )
    top      <- score[ cbind( seq_len( nrow( score ) ), dominant ) ]
    dominant[ top <= 0 ] <- 0L

    background.idx <- which( dominant == 0L )

    condition.curr <- calculate.condition.number( design )

    log.rows <- list()

    for ( f in seq_along( fluorophores ) ) {

      j   <- fluorophores[ f ]
      idx <- which( dominant == f )

      reject <- NA_character_

      if ( j %in% af.frozen ) reject <- "af.coupled"

      candidate <- NULL

      if ( is.na( reject ) ) {

        candidate <- extract.raw.signature(
          raw.data       = residual[ idx, , drop = FALSE ],
          spectra        = spectra.new[ fluorophores, , drop = FALSE ],
          abundance      = unmixed.comp[ idx, , drop = FALSE ],
          target         = j,
          active         = active.set[[ j ]],
          intercept      = intercept,
          multivariate   = multivariate,
          ridge          = ridge,
          n.levels       = n.levels,
          min.bin.events = min.bin.events,
          min.events     = min.negative.events,
          background.raw = if ( length( background.idx ) > 0 )
            residual[ background.idx, , drop = FALSE ] else NULL )

        if ( is.null( candidate ) ) reject <- "no.fit"
      }

      if ( is.na( reject ) ) {

        st <- candidate$stats

        reject <- NA_character_

        if ( st$x.span <= min.span * abs( thresholds[ j ] ) )
          reject <- "span" else
            if ( !is.finite( st$explained.total ) ||
                 st$explained.total < min.explained ||
                 st$explained.total > max.explained )
              reject <- "explained" else
                if ( st$resid.rel > max.resid )
                  reject <- "fit" else
                    if ( st$intercept.rel > max.intercept )
                      reject <- "offset" else
                        if ( is.finite( st$bg.align ) && st$bg.align < min.bg.align )
                          reject <- "background.confound" else
                            if ( st$clamp.frac > max.clamp.frac )
                              reject <- "negative.mass" else
                                if ( st$deg.change > max.angle )
                                  reject <- "angle" else
                                    if ( is.finite( st$anchor.rel ) && st$anchor.rel > max.anchor )
                                      reject <- "background" else
                                        if ( st$vif.target > max.vif )
                                          reject <- "collinear" else
                                            if ( st$peak.new != st$peak.curr &&
                                                 st$peak.new.rel < peak.shift.min.rel )
                                              reject <- "peak.shift"
      }

      leak.before     <- NA_real_
      leak.after      <- NA_real_
      condition.after <- NA_real_

      if ( is.na( reject ) ) {

        proposed <- candidate$signature

        if ( !is.null( null.spectra ) && j %in% rownames( null.spectra ) )
          proposed <- proposed -
            ( null.spectra[ j, colnames( spectra ) ] - spectra[ j, ] )

        stepped <- pmax( spectra.new[ j, ] +
                           step * ( proposed - spectra.new[ j, ] ), 0 )

        if ( max( stepped ) <= 0 ) {
          reject <- "no.fit"
        } else {

          trial <- spectra.new
          trial[ j, ] <- stepped / max( stepped )
        }
      }

      if ( is.na( reject ) ) {

        design.trial <- if ( is.null( background.basis ) )
          trial[ fluorophores, , drop = FALSE ] else
            rbind( background.basis, trial[ fluorophores, , drop = FALSE ] )

        condition.after <- calculate.condition.number( design.trial )

        leak.seed <- ( if ( is.null( asp$bird.seed ) ) 0L else asp$bird.seed ) +
          f * 1009L

        set.seed( leak.seed )
        leak.before <- .fix.leakage( residual[ idx, , drop = FALSE ],
                                     design, j, fluorophores )
        set.seed( leak.seed )
        leak.after  <- .fix.leakage( residual[ idx, , drop = FALSE ],
                                     design.trial, j, fluorophores )

        reject <- if ( condition.after >
                       max.condition.increase * condition.curr )
          "conditioning" else
            if ( is.finite( leak.before ) && is.finite( leak.after ) &&
                 leak.after >= ( 1 + leakage.margin ) * leak.before )
              "leakage" else NA_character_

        if ( is.na( reject ) ) {
          spectra.new    <- trial
          design         <- design.trial
          condition.curr <- condition.after
        }
      }

      stats.block <- if ( is.null( candidate ) )
        .fix.empty.signature.stats() else
          candidate$stats[ , setdiff( names( candidate$stats ),
                                      c( "fluorophore", "n.events" ) ),
                           drop = FALSE ]

      # Diagnostic only, never used for acceptance: the same fit with every
      # fluorophore active, fix.my.unmix()'s convention, run purely to show
      # whether a rejected row's residual comes down when nothing is left
      # out of the nuisance removal. A gap here (full clears max.resid, the
      # active-set fit does not) points at real coupling outside the lasso's
      # selection; no gap points at bin noise from a small population.
      resid.rel.full <- NA_real_

      if ( is.na( reject ) || reject == "fit" ) {

        candidate.full <- extract.raw.signature(
          raw.data       = residual[ idx, , drop = FALSE ],
          spectra        = spectra.new[ fluorophores, , drop = FALSE ],
          abundance      = unmixed.comp[ idx, , drop = FALSE ],
          target         = j,
          active         = fluorophores,
          intercept      = intercept,
          multivariate   = multivariate,
          ridge          = ridge,
          n.levels       = n.levels,
          min.bin.events = min.bin.events,
          min.events     = min.negative.events,
          background.raw = if ( length( background.idx ) > 0 )
            residual[ background.idx, , drop = FALSE ] else NULL )

        if ( !is.null( candidate.full ) )
          resid.rel.full <- candidate.full$stats$resid.rel
      }

      log.rows[[ length( log.rows ) + 1L ]] <- data.frame(
        fluorophore = j,
        n.events    = length( idx ),
        n.active    = length( active.set[[ j ]] ),
        stats.block,
        resid.rel.full  = resid.rel.full,
        leak.before     = leak.before,
        leak.after      = leak.after,
        condition.after = condition.after,
        accepted        = is.na( reject ),
        reason          = if ( is.na( reject ) ) "accepted" else reject,
        row.names       = NULL,
        stringsAsFactors = FALSE
      )
    }

    signature.log <- do.call( rbind, log.rows )

    if ( verbose )
      message( sprintf(
        "\033[32mAccepted signature updates for %d of %d fluorophore(s)%s.\033[0m",
        sum( signature.log$accepted ), fluorophore.n,
        if ( any( signature.log$accepted ) ) paste0(
          ": ", paste( signature.log$fluorophore[ signature.log$accepted ],
                       collapse = ", " ) ) else "" ) )
  }

  # ---------------------------------------------------------------------------
  # Outputs
  # ---------------------------------------------------------------------------

  spectra.backsolved <- spillover.curr %*% spectra[ fluorophores, , drop = FALSE ]
  row.max <- apply( spectra.backsolved, 1, max )
  row.max[ !is.finite( row.max ) | row.max <= 0 ] <- 1
  spectra.backsolved <- spectra.backsolved / row.max

  if ( figures )
    spectral.heatmap(
      spectra = spillover.curr,
      title = paste0( "CorrectSpectraGlasso_spillover_heatmap", output.suffix ),
      plot.dir = asp$fix.unmixing.dir,
      legend.label = "Spillover",
      save = TRUE
    )

  if ( save ) {
    utils::write.csv( spillover.curr,
                      file = file.path( asp$fix.unmixing.dir,
                                        paste0( "spillover", output.suffix, ".csv" ) ) )
    utils::write.csv( compensation.curr,
                      file = file.path( asp$fix.unmixing.dir,
                                        paste0( "compensation", output.suffix, ".csv" ) ) )
    utils::write.csv( spectra.new,
                      file = file.path( asp$fix.unmixing.dir,
                                        paste0( "spectra", output.suffix, ".csv" ) ) )
  }

  list(
    spectra            = spectra.new,
    spectra.backsolved = spectra.backsolved,
    spillover          = spillover.curr,
    compensation       = compensation.curr,
    coefficient.log    = coefficient.log,
    lambda.log         = lambda.log,
    signature.log      = signature.log,
    convergence.log    = convergence.log,
    active.set         = active.set,
    reversed.donors    = reversed.donors,
    af.basis           = background.basis,
    af.hotspot         = af.hotspot,
    af.frozen          = af.frozen
  )
}


#' One row of the residual spillover matrix by neighbourhood-selection lasso:
#' target fluorophore's compensated abundance, regressed on every other
#' fluorophore's abundance at once, over the target's own negative-population
#' events.
#'
#' @param abundance Numeric matrix (events x fluorophores), currently
#'   compensated abundances.
#' @param target Character, the fluorophore whose row is being estimated.
#' @param sources Character vector, every other fluorophore.
#' @param neg.idx Integer vector, row indices of `abundance` currently
#'   treated as negative for `target`.
#' @param n.lambda,lambda.min.ratio,margin.frac See
#'   `correct.spectra.glasso()`.
#' @param max.events Integer, cap on events used; above it `neg.idx` is
#'   subsampled.
#'
#' @return A named list: `beta` (named numeric, one coefficient per source),
#'   `intercept`, `lambda` (the chosen penalty), `n.events`, `path` (the
#'   two-way holdout error curve, for diagnostics).
#' @noRd
.glasso.neighborhood.row <- function( abundance, target, sources, neg.idx,
                                      n.lambda         = 40L,
                                      lambda.min.ratio = 1e-3,
                                      margin.frac       = 0.05,
                                      max.events        = 20000L ) {

  if ( length( neg.idx ) > max.events )
    neg.idx <- sample( neg.idx, max.events )

  X <- abundance[ neg.idx, sources, drop = FALSE ]
  y <- abundance[ neg.idx, target ]

  sel <- .glasso.select.lambda(
    X, y, n.lambda = n.lambda, lambda.min.ratio = lambda.min.ratio,
    margin.frac = margin.frac )

  list(
    beta      = sel$beta,
    intercept = sel$intercept,
    lambda    = sel$lambda,
    n.events  = length( neg.idx ),
    path      = sel$path
  )
}


#' Chooses the lasso penalty for one neighbourhood row by a two-way
#' split-half holdout, then refits at the chosen penalty on the full
#' (target-negative) data supplied.
#'
#' Events are split into two halves; a path is fit on each half and scored on
#' the other, over a single lambda grid shared by both fits so the two error
#' curves are comparable point for point. Their sum is the two-way holdout
#' error curve. Rather than take its minimiser outright, the largest lambda
#' (sparsest fit) whose error is within `margin.frac` of the minimum is
#' preferred - the same bias toward the more conservative of two disagreeing
#' estimates `correct.unmixing.signatures()`'s held-out step search already
#' uses, taking the smaller of two step sizes when the two orderings
#' disagree.
#'
#' @param X Numeric matrix (events x sources).
#' @param y Numeric vector, length `nrow(X)`.
#' @param weights Optional numeric vector of per-event weights. Default
#'   `NULL` (uniform).
#' @param n.lambda,lambda.min.ratio,margin.frac See
#'   `correct.spectra.glasso()`.
#'
#' @return A named list: `lambda`, `beta` (named numeric over
#'   `colnames(X)`), `intercept`, `path` (data frame of the grid and both
#'   halves' errors).
#' @noRd
.glasso.select.lambda <- function( X, y, weights = NULL, n.lambda = 40L,
                                   lambda.min.ratio = 1e-3,
                                   margin.frac       = 0.05 ) {

  n <- nrow( X )
  if ( is.null( weights ) ) weights <- rep( 1, n )

  fit.full <- .glasso.lasso.path( X, y, weights = weights,
                                  n.lambda = n.lambda,
                                  lambda.min.ratio = lambda.min.ratio )
  grid <- fit.full$lambda

  # A deterministic split rather than an extra call to sample(): stable
  # under repeated calls at the same seed for the same event count, and
  # needs no separate `set.seed()` bookkeeping the way `sample()` would to
  # stay reproducible across mask passes.
  half <- rep_len( c( 1L, 2L ), n )

  eval.holdout <- function( train, test ) {

    fit <- .glasso.lasso.path( X[ train, , drop = FALSE ], y[ train ],
                               weights = weights[ train ],
                               n.lambda = n.lambda,
                               lambda.min.ratio = lambda.min.ratio )

    beta.grid <- apply( fit$beta, 1, function( b )
      stats::approx( fit$lambda, b, xout = grid, rule = 2 )$y )
    if ( is.null( dim( beta.grid ) ) )
      beta.grid <- matrix( beta.grid, ncol = 1,
                           dimnames = list( NULL, colnames( X ) ) )

    int.grid <- stats::approx( fit$lambda, fit$intercept, xout = grid,
                               rule = 2 )$y

    pred <- sweep( X[ test, , drop = FALSE ] %*% t( beta.grid ), 2,
                   int.grid, "+" )

    colSums( weights[ test ] * ( y[ test ] - pred )^2 ) /
      sum( weights[ test ] )
  }

  err.2 <- eval.holdout( which( half == 1L ), which( half == 2L ) )
  err.1 <- eval.holdout( which( half == 2L ), which( half == 1L ) )

  total <- err.1 + err.2

  best      <- which.min( total )
  tolerated <- total <= ( 1 + margin.frac ) * total[ best ]

  chosen <- max( grid[ tolerated ] )
  li     <- which( grid == chosen )[ 1 ]

  list(
    lambda    = chosen,
    beta      = fit.full$beta[ , li ],
    intercept = fit.full$intercept[ li ],
    path      = data.frame( lambda = grid, error.1 = err.1, error.2 = err.2,
                            total = total )
  )
}


#' Pathwise coordinate-descent lasso, base R only. Predictors and response
#' are weighted-centred and predictors weighted-standardised so a single
#' `lambda` path is comparable across sources in different abundance units;
#' the intercept is recovered afterward in closed form rather than treated
#' as a penalised coordinate. Coordinates are updated by the standard
#' running-residual trick (add a coordinate's own contribution back into the
#' residual before its update, remove the new contribution after), which
#' keeps a full sweep at O(events x sources) rather than O(events x
#' sources^2); this is the part a future Rcpp port, following the pattern
#' already set by `fix_huber_slope_rcpp()`, would target first.
#'
#' @param X Numeric matrix (events x predictors).
#' @param y Numeric vector, length `nrow(X)`.
#' @param weights Optional numeric vector of per-event weights. Default
#'   `NULL` (uniform).
#' @param n.lambda Integer, points on the path. Default `40`.
#' @param lambda.min.ratio Numeric in (0, 1), smallest path value as a
#'   fraction of `lambda.max`. Default `1e-3`.
#' @param tol Numeric, coordinate-descent convergence tolerance on the
#'   largest coefficient change within a sweep. Default `1e-7`.
#' @param max.iter Integer, maximum sweeps per `lambda`. Default `1000`.
#'
#' @return A named list: `lambda` (length `n.lambda`), `beta` (predictors x
#'   `n.lambda`), `intercept` (length `n.lambda`), and the centring/scaling
#'   constants used.
#' @noRd
.glasso.lasso.path <- function( X, y, weights = NULL, n.lambda = 40L,
                                lambda.min.ratio = 1e-3, tol = 1e-7,
                                max.iter = 1000L ) {

  X <- as.matrix( X )
  n <- nrow( X )
  p <- ncol( X )
  fluor.names <- colnames( X )

  if ( is.null( weights ) ) weights <- rep( 1, n )
  w.sum <- sum( weights )

  x.center <- colSums( weights * X ) / w.sum
  y.mean   <- sum( weights * y ) / w.sum

  Xc <- sweep( X, 2, x.center, "-" )
  yc <- y - y.mean

  x.scale <- sqrt( colSums( weights * Xc^2 ) / w.sum )
  x.scale[ x.scale <= 0 ] <- 1

  Z   <- sweep( Xc, 2, x.scale, "/" )
  wz  <- weights * Z
  wz2 <- colSums( wz * Z ) / w.sum
  wz2[ wz2 <= 0 ] <- 1

  grad0      <- as.numeric( crossprod( Z, weights * yc ) ) / w.sum
  lambda.max <- max( abs( grad0 ) )

  if ( !is.finite( lambda.max ) || lambda.max <= 0 ) {

    beta.zero <- matrix( 0, p, n.lambda, dimnames = list( fluor.names, NULL ) )

    return( list(
      lambda    = rep( 0, n.lambda ),
      beta      = beta.zero,
      intercept = rep( y.mean, n.lambda ),
      x.center  = x.center, x.scale = x.scale, y.mean = y.mean ) )
  }

  lambda.seq <- lambda.max *
    lambda.min.ratio ^ ( seq( 0, 1, length.out = n.lambda ) )

  beta.path <- matrix( 0, p, n.lambda, dimnames = list( fluor.names, NULL ) )
  beta.curr <- rep( 0, p )
  resid     <- yc

  for ( li in seq_len( n.lambda ) ) {

    lam <- lambda.seq[ li ]

    for ( sweep.i in seq_len( max.iter ) ) {

      max.change <- 0

      for ( k in seq_len( p ) ) {

        beta.old.k <- beta.curr[ k ]

        r.k <- resid + Z[ , k ] * beta.old.k
        z.k <- sum( wz[ , k ] * r.k ) / w.sum

        beta.new.k <- sign( z.k ) * max( abs( z.k ) - lam, 0 ) / wz2[ k ]

        if ( beta.new.k != beta.old.k ) {
          resid          <- resid - Z[ , k ] * ( beta.new.k - beta.old.k )
          beta.curr[ k ] <- beta.new.k
        }

        max.change <- max( max.change, abs( beta.new.k - beta.old.k ) )
      }

      if ( max.change < tol ) break
    }

    beta.path[ , li ] <- beta.curr
  }

  beta.orig      <- beta.path / x.scale
  intercept.path <- y.mean - as.numeric( x.center %*% beta.orig )

  list(
    lambda    = lambda.seq,
    beta      = beta.orig,
    intercept = intercept.path,
    x.center  = x.center,
    x.scale   = x.scale,
    y.mean    = y.mean
  )
}
