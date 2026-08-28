# fix_my_unmix.R

#' @title Fix My Unmix
#'
#' @description
#' Corrects the component of a reference spectra error that lives inside the row
#' space of the spectra themselves, using a fully stained sample and the fact
#' that a marker-negative population must read zero in every channel it is
#' negative for.
#'
#' This is the complement of `correct.unmixing.signatures()`, which corrects the
#' shape of individual rows from detector-space residuals and is structurally
#' blind to in-span error. Run this second, on the spectra that function
#' returns, because in unmixed space a row-shape error and a genuine spillover
#' error are indistinguishable.
#'
#' The correction runs in two phases. The first estimates the residual spillover
#' matrix pair by pair, from a robust line through the events that are negative
#' for the target channel: their apparent abundance there is spillover from the
#' source, so its slope against source abundance is the residual coefficient.
#' Truncating the target from above shrinks that slope toward zero when the two
#' markers are co-expressed, which is the safe direction. The second phase does
#' not back-solve a spectrum out of the matrix. It uses the matrix only to
#' define clean populations, then re-measures each fluorophore's signature
#' directly from background subtracted raw data with `extract.raw.signature()`,
#' so the result is a physical spectrum rather than whichever combination of the
#' existing rows happens to reproduce the required compensation.
#'
#' Both phases have a bias that can be measured without any ground truth. Run on
#' the control the reference spectra came from, this function should return
#' exactly what it was given; whatever it returns instead is its own artefact.
#' Supplying that run as `null.fit` subtracts the artefact from the coefficients
#' and from the signatures, which recovers substantially more than refusing the
#' same quantities on the same information.
#'
#' Autofluorescence is fitted jointly with the panel from a multi-component
#' basis (`bg.mode = "af.deconv"`), not represented by a single row. In a real
#' sample the background a single row leaves behind tracks cell type, cell type
#' determines marker expression, and the leftover is then indistinguishable from
#' spillover by any channel-against-channel fit.
#'
#' Every row update must pass an acceptance stack; a fluorophore that fails any
#' gate keeps its starting spectrum.
#'
#' @importFrom sp point.in.polygon
#' @importFrom stats mad median qnorm quantile sd setNames lm.wfit
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
#' @param variants The variant list returned by `get.spectral.variants()`. Only
#'   `variants$spillover.spread` is used, to set abundance-dependent positivity
#'   boundaries and to correct the envelope for the widening of the negative
#'   population. If `NULL`, boundaries are flat and the correction is biased at
#'   the bright end.
#' @param af.name Character or `NULL`, the name of an autofluorescence row in
#'   `spectra`. It is never treated as a panel fluorophore and never corrected.
#'   Default `"AF"`.
#' @param af.basis Optional matrix (components x detectors) from
#'   `get.af.basis()`. When `NULL` and `bg.mode = "af.deconv"`, it is built from
#'   the unstained sample. Default `NULL`.
#' @param af.n.pc Integer or `"auto"`, passed to `get.af.basis()`. Default
#'   `"auto"`.
#' @param bg.mode Character. `"af.deconv"` (default) fits a multi-component
#'   autofluorescence basis jointly with the panel; `"af.row"` uses the single
#'   `af.name` row already in `spectra`; `"global.mean"` uses the mean unstained
#'   spectrum as a single background row; `"none"` fits no background.
#' @param large.gate Logical, whether to use a large scatter gate. Default
#'   `TRUE`.
#' @param max.iter Integer, maximum spillover-matrix iterations. Default `20`.
#' @param downsample Logical or numeric. `FALSE` disables downsampling; a
#'   numeric gives the number of events to use. Values above the event count
#'   are reduced to it by a stratified sample over each event's dominant
#'   fluorophore under the starting spectra, so a dim or rare dye's own
#'   positive population is not thinned at the same rate as the background
#'   bulk. Default `20000`.
#' @param downsample.background.frac Numeric in (0, 1), the share of
#'   `downsample` reserved for events dominant for nothing. Default `0.3`.
#' @param downsample.min.stratum Integer, the floor below which a
#'   fluorophore's own positive population is kept whole by the stratified
#'   downsample rather than thinned further. Default `2000`.
#' @param unstained.threshold Numeric in (0, 1), the percentile of the unstained
#'   control defining positivity. Default `0.99`.
#' @param unstained.margin Numeric, multiplier applied to that threshold.
#'   Default `1.3`.
#' @param spread.kappa Numeric, how many spillover-spread standard deviations
#'   above the flat threshold still count as negative. Default `2`.
#' @param envelope.quantiles Numeric pair. The first is the envelope the slope
#'   is measured on; the second is a comparison quantile whose disagreement with
#'   the first reports co-expression. Default `c( 0.05, 0.5 )`.
#' @param min.negative.frac Numeric, the smallest fraction of target-negative
#'   events the brightest source bins may contain before the pair is declared
#'   unidentifiable and left alone. Default `0.10`.
#' @param max.disagreement Numeric, envelope-versus-median slope disagreement
#'   above which the coefficient's trust weight is reduced. Default `0.5`.
#' @param leakage.prior Logical, whether to weight coefficients by the spillover
#'   a plausible spectral error could produce, from `get.variant.leakage.prior()`.
#'   When `FALSE`, or when no variants are supplied, the estimate supplies its
#'   own prior variance, which lets a large spurious coefficient justify itself.
#'   Default `TRUE`.
#' @param span.fraction Numeric in (0, 1], passed to
#'   `get.variant.leakage.prior()`. The fraction of real spectral error the
#'   variant family is expected to span; lower values give corrections more
#'   freedom outside the observed variant directions. Default `0.6`.
#' @param min.negative.events Integer, minimum events for a pair fit. Default
#'   `200`.
#' @param min.bin.negative Integer, minimum target-negative events an abundance
#'   bin must contain before its envelope quantile is used. Default `25`.
#' @param min.span Numeric, minimum source abundance span in units of that
#'   fluorophore's own threshold. Default `5`.
#' @param min.rise Numeric, the fitted rise across the source abundance span,
#'   in standard deviations of the target's negative population, below which a
#'   coefficient is treated as unresolvable and set to zero. Default `1`.
#' @param n.levels Integer, the maximum number of abundance bins for the
#'   signature fit. Must exceed the panel size by at least three when
#'   `multivariate` is `TRUE`. Default `60`.
#' @param min.bin.events Integer, the fewest events an abundance bin may contain,
#'   passed to `extract.raw.signature()`. Default `50`.
#' @param n.levels.pair Integer, abundance bins for the pair estimator. Default
#'   `10`.
#' @param multivariate Logical, passed to `extract.raw.signature()`, whether
#'   every active abundance is fitted jointly rather than the nuisance rows
#'   being subtracted beforehand. Default `TRUE`.
#' @param ridge Numeric, ridge penalty for that joint fit. Default `1e-6`.
#' @param max.truncated.events Integer, cap on the events used for the robust
#'   pair fit. Everything above the source threshold is kept and the negative
#'   bulk below it is subsampled, since the bulk sits at the origin and carries
#'   no leverage. Default `20000`.
#' @param max.mask.passes Integer, how many times the target-negative selection
#'   is recomputed with the current spillover estimate removed. Deciding
#'   negativity on the raw target caps the fit's source range at the target
#'   threshold divided by the coefficient, so the estimator loses its leverage
#'   precisely as the coefficient it is measuring grows. Default `3`.
#' @param estimator Character, which pair estimator supplies the coefficients.
#'   `"truncated"` (default) is a robust line through the target-negative
#'   events; `"envelope"` is the binned lower envelope.
#' @param source.dominant Logical, whether source-positive events dominated by
#'   another fluorophore are dropped from the pair fit. Their apparent source
#'   abundance is that other fluorophore's brightness read through a spectral
#'   error, so binning on it returns that error as this pair's spillover.
#'   Default `TRUE`.
#' @param spread.addback Logical, whether the measured spillover spread is added
#'   back to the envelope before its slope is taken. Default `FALSE`.
#' @param anchor.weight Numeric, the weight ceiling for the zero-abundance
#'   anchor bin, in multiples of the best positive bin. Default `1`.
#' @param max.coefficient Numeric, the largest residual spillover coefficient
#'   accepted. Rows are L-infinity normalised, so a coefficient this large means
#'   one row is wrong by that fraction of another row's whole spectrum, which
#'   violates the small-error premise the linearisation rests on. Default `0.2`.
#' @param convergence.threshold Numeric, residual spillover coefficient at which
#'   iteration stops. Default `0.01`.
#' @param convergence.quantile Numeric, the quantile of the off-diagonal
#'   coefficients the convergence test uses, so that one pathological pair
#'   cannot define convergence. Default `0.95`.
#' @param update.spectra Logical, whether to run the raw-space signature phase.
#'   When `FALSE`, only the spillover and compensation matrices are returned.
#'   Default `TRUE`.
#' @param step Numeric, the fraction of each accepted signature change applied.
#'   Default `1`.
#' @param null.fit Optional, the result of running this function on the control
#'   the reference spectra were extracted from, where it should return exactly
#'   what it was given. What it returns instead is its own artefact, measured
#'   with no ground truth, and it is subtracted from both the spillover matrix
#'   and the signatures. Its `signature.log` also sets `max.resid` and
#'   `max.intercept` from the panel's own distribution. Default `NULL`.
#' @param max.resid.ratio Numeric, multiple of the null run's median
#'   `resid.rel` above which a fit is refused. Ignored without `null.fit`.
#'   Default `3`.
#' @param max.intercept.ratio Numeric, the same for `intercept.rel`. Default
#'   `3`.
#' @param intercept Logical, whether the signature fit carries an intercept,
#'   which absorbs any constant offset the removal of the other panel rows left
#'   behind. Default `TRUE`.
#' @param min.explained Numeric, minimum `explained.total`, the fraction of the
#'   brightest bin's signal the whole signature fit reproduces. Default `0.8`.
#' @param max.explained Numeric, the maximum of the same quantity. A value above
#'   one means the nuisance removal over-subtracted. Default `1.2`.
#' @param max.resid Numeric, maximum relative fit residual. Overridden by
#'   `null.fit` when supplied. Default `0.03`.
#' @param max.intercept Numeric, maximum relative intercept. Overridden by
#'   `null.fit` when supplied. Default `0.03`.
#' @param min.bg.align Numeric, minimum cosine between the fitted intercept and
#'   the candidate signature. A value near minus one means the fit is trading a
#'   background floor against the slope, the event-level signature of an
#'   unmodelled background confound. Default `-0.9`.
#' @param gate.on.bias Logical, whether `min.impact.ratio` actually blocks a
#'   candidate, versus being computed and logged only. Default `FALSE`.
#'   Evidence from a counterfactual audit (score whether accepting each
#'   rejected candidate would have moved the row toward or away from ground
#'   truth): of every candidate this gate rejected, 79.9% would have helped
#'   if accepted (n=144, two substrates) - 90.7% on one substrate, a
#'   coin-flip 47.2% on the other. The gate is net-harmful on the substrate
#'   where it discriminates worst and only weakly useful on the other, so it
#'   defaults off rather than removed - `bias.impact` and the ratio are still
#'   computed and logged either way, so the evidence for re-enabling it on a
#'   new substrate remains visible without code changes.
#' @param min.impact.ratio Numeric, how many times a proposed step's abundance
#'   effect must exceed the effect of the same phase's bias on its own control,
#'   measured through the current unmixing operator. Degrees weight every
#'   detector alike and every dye alike; this weights each row by the abundances
#'   it actually produces, so a large rotation of a dim row is cheap and a small
#'   rotation of a bright collinear row is not. Requires `null.spectra` and
#'   `gate.on.bias = TRUE`. Default `2`.
#' @param max.angle Numeric, maximum angular change of a row in degrees.
#' @param max.clamp.frac Numeric, maximum fraction of a candidate row's
#'   absolute mass that non-negativity clamping may remove. Default `0.15`.
#' @param max.anchor Numeric, maximum unremoved background relative to the
#'   brightest fitted signal. Default `0.10`.
#' @param max.vif Numeric, maximum variance inflation factor for the
#'   fluorophore's own abundance within its population. Default `500`.
#' @param max.condition.increase Numeric, the factor by which a single accepted
#'   row may increase the condition number of the unmixing design. Default
#'   `1.05`.
#' @param leakage.margin Numeric, the fraction by which held-out leakage must
#'   increase, not merely be numerically larger, before a candidate is
#'   refused for it. `.fix.leakage()` is a held-out statistic on a
#'   dominance population that can be small, and "any increase blocks" was
#'   itself producing false positives - candidates blocked on multiple
#'   consecutive fits that would demonstrably have helped if accepted.
#'   Default `0.05`.
#' @param peak.shift.min.rel Numeric, a candidate row whose peak detector
#'   differs from the current row is only accepted if the current signature's
#'   own value at that new peak channel, relative to the current signature's
#'   own peak, is at least this large - a close secondary peak taking over is
#'   plausible, a shift to a channel the current signature barely uses is not.
#'   Use `Inf` to never allow a peak shift and `0` to always allow one.
#'   Default `0.7`.
#' @param max.hotspot Numeric, hotspot scale above which a fluorophore is
#'   considered inseparable from the autofluorescence basis and is frozen.
#'   Default `5`.
#' @param n.threads Integer, OpenMP threads for the batched pair estimator
#'   (`fix_envelope_truncated_batch_rcpp()`). Keep at the default unless this
#'   call is not itself running inside another parallel context (e.g. one
#'   sample of several under `mclapply()`), since the two layers of
#'   parallelism would otherwise compete for the same cores. Default `1L`.
#' @param figures Logical, whether to write the spillover heatmap. Default
#'   `TRUE`.
#' @param save Logical, whether to write the csv outputs. Default `TRUE`.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`spectra`}{The corrected spectra, fluorophores x detectors. Rows that
#'     failed a gate are unchanged.}
#'   \item{`spectra.backsolved`}{`spillover %*% spectra`, the algebraic
#'     back-solve. Diagnostic only: it reproduces the compensation exactly but is
#'     not constrained to be a physical spectrum.}
#'   \item{`spillover`}{The estimated residual spillover matrix, fluorophores x
#'     fluorophores.}
#'   \item{`compensation`}{Its inverse.}
#'   \item{`trust`}{Per-coefficient trust weights used to damp the update.}
#'   \item{`coefficient.log`}{Per-pair fit diagnostics from the final
#'     iteration.}
#'   \item{`signature.log`}{Per-fluorophore signature statistics and gate
#'     outcomes.}
#'   \item{`convergence.log`}{Per-iteration delta history.}
#'   \item{`af.basis`, `af.hotspot`, `af.frozen`}{The autofluorescence basis, its
#'     coupling to the panel, and the fluorophores frozen because of it.}
#' }
#'
#' @export

fix.my.unmix <- function(
    spectra,
    unstained.sample,
    fully.stained.sample,
    flow.control,
    asp,
    variants               = NULL,
    af.name                = "AF",
    af.basis               = NULL,
    af.n.pc                = "auto",
    bg.mode                = c( "af.deconv", "af.row", "global.mean", "none" ),
    large.gate             = TRUE,
    max.iter               = 20L,
    downsample             = 20000,
    downsample.background.frac = 0.3,
    downsample.min.stratum     = 2000L,
    unstained.threshold    = 0.99,
    unstained.margin       = 1.3,
    spread.kappa           = 2,
    envelope.quantiles     = c( 0.05, 0.5 ),
    min.negative.frac      = 0.10,
    max.disagreement       = 0.5,
    leakage.prior          = TRUE,
    span.fraction          = 0.6,
    min.negative.events    = 200L,
    min.bin.negative       = 25L,
    min.span               = 5,
    min.rise               = 1,
    n.levels               = 60L,
    min.bin.events         = 50L,
    n.levels.pair          = 10L,
    multivariate           = TRUE,
    ridge                  = 1e-6,
    max.truncated.events   = 20000L,
    max.mask.passes        = 3L,
    estimator              = c( "truncated", "envelope" ),
    source.dominant        = TRUE,
    spread.addback         = FALSE,
    anchor.weight          = 1,
    max.coefficient        = 0.5,
    convergence.threshold  = 0.01,
    convergence.quantile   = 0.95,
    update.spectra         = TRUE,
    step                   = 1,
    null.fit               = NULL,
    max.resid.ratio        = 3,
    max.intercept.ratio    = 3,
    intercept              = TRUE,
    min.explained          = 0.8,
    max.explained          = 1.2,
    max.resid              = 0.03,
    max.intercept          = 0.03,
    min.bg.align           = -0.9,
    gate.on.bias           = FALSE,
    min.impact.ratio       = 2,
    max.angle              = 15,
    max.clamp.frac         = 0.15,
    max.anchor             = 0.10,
    max.vif                = 500,
    max.condition.increase = 1.05,
    peak.shift.min.rel     = 0.7,
    max.hotspot            = 5,
    leakage.margin         = 0.05,
    n.threads              = 1L,
    figures                = TRUE,
    save                   = TRUE,
    verbose                = TRUE
) {

  bg.mode   <- match.arg( bg.mode )
  estimator <- match.arg( estimator )

  # Set once, not only around the downsample, because the pair estimator also
  # subsamples the negative bulk and the whole run should be reproducible.
  if ( !is.null( asp$bird.seed ) ) set.seed( asp$bird.seed )

  # The same estimator, run on the control the reference spectra were extracted
  # from, should return exactly what it was given. What it returns instead is
  # its own artefact, measured with no ground truth, and it is removed from both
  # the spillover matrix and the signatures rather than used to refuse them.
  null.spillover <- null.fit$spillover
  null.spectra   <- null.fit$spectra

  # `resid.rel` and `intercept.rel` describe the row's population and its
  # background, not whether the spectra are wrong, so they reproduce almost
  # exactly on the null run and their absolute scale is a property of the
  # particle type. The thresholds are therefore set against the panel's own
  # distribution.
  if ( !is.null( null.fit$signature.log ) ) {

    max.resid <- max.resid.ratio * stats::median(
      null.fit$signature.log$resid.rel, na.rm = TRUE )

    max.intercept <- max.intercept.ratio * stats::median(
      null.fit$signature.log$intercept.rel, na.rm = TRUE )
  }

  spectra <- as.matrix( spectra )

  if ( is.null( rownames( spectra ) ) )
    stop( "`spectra` must have fluorophore row names.", call. = FALSE )

  if ( length( envelope.quantiles ) != 2 )
    stop( "`envelope.quantiles` must be a pair.", call. = FALSE )

  fluorophores  <- setdiff( rownames( spectra ), af.name )
  fluorophore.n <- length( fluorophores )

  if ( fluorophore.n < 2 )
    stop( "At least two panel fluorophores are required.", call. = FALSE )

  if ( save && ! dir.exists( asp$fix.unmixing.dir ) )
    dir.create( asp$fix.unmixing.dir, recursive = TRUE )

  spillover.spread <- variants$spillover.spread

  if ( is.null( spillover.spread ) )
    warning( paste0( "No `spillover.spread` in `variants`; positivity ",
                     "boundaries will be flat and the envelope will not be ",
                     "corrected for the widening of the negative population, ",
                     "biasing coefficients at the bright end." ),
             call. = FALSE )

  # ---------------------------------------------------------------------------
  # Read and gate
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
  # Background basis and the projection it defines
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

  prior <- NULL

  if ( leakage.prior ) {

    if ( is.null( variants ) ) {
      warning( paste0( "`leakage.prior` requires `variants`; falling back to ",
                       "the estimate's own variance." ), call. = FALSE )
    } else {
      prior <- get.variant.leakage.prior(
        spectra       = spectra,
        variants      = variants,
        extra.rows    = background.basis,
        af.name       = af.name,
        span.fraction = span.fraction,
        verbose       = verbose )
    }
  }

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

  # A uniform downsample keeps every population in the same proportion the
  # raw file already had it in, so it thins a dim or rare fluorophore's own
  # positive population at the same rate as a background bulk that was
  # already more than large enough to set its own threshold. Stratifying by
  # dominant fluorophore under the starting spectra fixes that: every
  # fluorophore's stratum gets its own quota, with a floor that protects a
  # dim dye's whole positive population and a separate cap on the events
  # dominant for nothing.
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
  # Phase one: residual spillover matrix from the negative envelope
  # ---------------------------------------------------------------------------

  spillover.curr <- diag( fluorophore.n )
  dimnames( spillover.curr ) <- list( fluorophores, fluorophores )

  trust <- spillover.curr
  pair.log <- NULL

  # Warm-starts the truncated estimator's IRLS from the previous outer
  # iteration's fitted slope for the same pair, since spillover.curr moves by
  # less each iteration as it converges and each pair's true slope moves with
  # it. Zero for every pair on the first iteration, which is the same OLS
  # start the estimator has always used.
  slope.warm <- matrix( 0, fluorophore.n, fluorophore.n,
                        dimnames = list( fluorophores, fluorophores ) )

  convergence.log <- data.frame(
    iter = integer(), delta = numeric(), delta.quantile = numeric(),
    delta.max = numeric(), n.fitted = integer(), stringsAsFactors = FALSE )

  delta.history <- rep( NA_real_, 3L )

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

    # Thresholds and negative spread are recomputed in the current compensated
    # space. Holding them at their starting values lets the boundary drift away
    # from the data as the matrix is refined.
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

    dominant <- NULL

    if ( source.dominant ) {

      dyn.range <- pmax(
        apply( unmixed.comp, 2, stats::quantile, probs = 0.999,
               names = FALSE ) - thresholds, .Machine$double.eps )

      # Each abundance is scored against its own event's positivity boundary,
      # not the flat cut. A channel a bright dye spills into has a wider
      # negative population there, and scoring it flat while dividing by that
      # channel's own dynamic range hands the bright dye's mid-range events to
      # whichever narrower neighbour it spreads into.
      dominance.score <- sweep(
        pmax( unmixed.comp - threshold.matrix, 0 ), 2, dyn.range, "/" )

      dominant <- max.col( dominance.score, ties.method = "first" )
      dominant[ dominance.score[ cbind( seq_along( dominant ),
                                        dominant ) ] <= 0 ] <- 0L
    }

    marker.spillover <- diag( fluorophore.n )
    dimnames( marker.spillover ) <- list( fluorophores, fluorophores )

    trust <- diag( fluorophore.n )
    dimnames( trust ) <- list( fluorophores, fluorophores )

    pair.log <- list()

    # One batched call per source instead of one scalar call per pair -
    # everything that depends only on the source (the dominance mask, the
    # source-positivity split, the abundance bins) is computed once and
    # shared across every channel it might spill into. Per-channel
    # acceptance logic below is unchanged from the scalar loop; only how
    # `est` is obtained differs.
    for ( source in fluorophores ) {

      channel.set <- setdiff( fluorophores, source )

      spread.var.vec <- rep( 0, length( channel.set ) )
      if ( !is.null( spillover.spread ) &&
           source %in% rownames( spillover.spread ) ) {
        avail <- intersect( channel.set, colnames( spillover.spread ) )
        spread.var.vec[ match( avail, channel.set ) ] <-
          spillover.spread[ source, avail ]
        spread.var.vec[ !is.finite( spread.var.vec ) ] <- 0
        spread.var.vec <- pmax( spread.var.vec, 0 )
      }

      est.batch <- .fix.envelope.slope.batch(
        x.source          = unmixed.comp[ , source ],
        X.target          = unmixed.comp[ , channel.set, drop = FALSE ],
        threshold.source  = threshold.matrix[ , source ],
        Threshold.target  = threshold.matrix[ , channel.set, drop = FALSE ],
        spread.var        = spread.var.vec,
        neg.var           = neg.var[ channel.set ],
        source.mask       = if ( is.null( dominant ) ) NULL else
          dominant == match( source, fluorophores ),
        quantiles             = envelope.quantiles,
        n.levels              = n.levels.pair,
        min.events            = min.negative.events,
        min.bin.negative      = min.bin.negative,
        spread.addback        = spread.addback,
        anchor.weight         = anchor.weight,
        max.truncated.events  = max.truncated.events,
        max.coefficient       = max.coefficient,
        max.mask.passes       = max.mask.passes,
        start.slope           = slope.warm[ source, channel.set ],
        n.threads             = n.threads )

      for ( ci in seq_along( channel.set ) ) {

        channel <- channel.set[ ci ]
        est <- if ( is.null( est.batch ) ) NULL else
          as.list( est.batch[ ci, , drop = FALSE ] )

        if ( !is.null( est ) && is.finite( est$slope.truncated ) )
          slope.warm[ source, channel ] <- est$slope.truncated

        w <- 0

        slope.use <- if ( is.null( est ) ) NA_real_ else
          if ( estimator == "truncated" ) est$slope.truncated else est$slope

        if ( !is.null( est ) && is.finite( slope.use ) &&
             est$coverage >= min.negative.frac &&
             abs( slope.use ) <= max.coefficient &&
             abs( slope.use ) * est$span >= min.rise * est$noise &&
             est$span > min.span * abs( thresholds[ source ] ) ) {

          prior.var <- if ( is.null( prior ) ) slope.use^2 else
            prior$variance[ source, channel ]

          if ( !is.finite( prior.var ) || prior.var <= 0 )
            prior.var <- slope.use^2

          w <- if ( is.finite( est$se ) && est$se > 0 )
            prior.var / ( prior.var + est$se^2 ) else 1

          if ( is.finite( est$disagreement ) &&
               est$disagreement > max.disagreement )
            w <- w * max.disagreement / est$disagreement

          marker.spillover[ source, channel ] <- slope.use
        }

        trust[ source, channel ] <- w

        pair.log[[ length( pair.log ) + 1L ]] <- data.frame(
          source       = source,
          channel      = channel,
          slope           = if ( is.null( est ) ) NA_real_ else est$slope,
          slope.truncated = if ( is.null( est ) ) NA_real_ else
            est$slope.truncated,
          slope.median    = if ( is.null( est ) ) NA_real_ else est$slope.alt,
          se           = if ( is.null( est ) ) NA_real_ else est$se,
          disagreement = if ( is.null( est ) ) NA_real_ else est$disagreement,
          coverage     = if ( is.null( est ) ) NA_real_ else est$coverage,
          span         = if ( is.null( est ) ) NA_real_ else est$span,
          trust        = w,
          row.names    = NULL,
          stringsAsFactors = FALSE
        )
      }
    }

    slope.error <- marker.spillover - diag( fluorophore.n )

    off.diagonal <- abs( slope.error[ row( slope.error ) != col( slope.error ) ] )

    delta          <- stats::sd( slope.error )
    delta.quantile <- unname( stats::quantile(
      off.diagonal, probs = convergence.quantile, names = FALSE ) )
    delta.max      <- max( off.diagonal )

    convergence.log[ nrow( convergence.log ) + 1L, ] <-
      list( iter, delta, delta.quantile, delta.max,
            sum( trust > 0 ) - fluorophore.n )

    if ( verbose )
      message( sprintf(
        "\033[34miter %d: delta %.5f, q%.0f %.5f, max %.5f, %d coefficient(s) fitted\033[0m",
        iter, delta, 100 * convergence.quantile, delta.quantile, delta.max,
        sum( trust > 0 ) - fluorophore.n ) )

    # The trust weight damps the step; it must not touch the estimate or the
    # convergence test, or the fixed point becomes trust-weighted identity
    # rather than identity.
    spillover.next <- spillover.curr + ( trust * slope.error ) %*% spillover.curr

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

  # Run on the control the reference spectra were extracted from, this estimator
  # should return the identity. Whatever it returns instead is its own bias,
  # measured without any ground truth, and it is subtracted here. The correction
  # is worth more than refusing the same coefficients on the same information.
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

    # The abundance each row actually produces, used to weight a proposed step
    # by the error it would move rather than by its angle.
    abundance.high <- apply( unmixed.comp, 2, stats::quantile, probs = 0.999,
                             names = FALSE )
    names( abundance.high ) <- fluorophores

    score    <- sweep( pmax( unmixed.comp - threshold.matrix, 0 ), 2,
                       dyn.range, "/" )
    dominant <- max.col( score, ties.method = "first" )
    top      <- score[ cbind( seq_len( nrow( score ) ), dominant ) ]
    dominant[ top <= 0 ] <- 0L

    background.idx <- which( dominant == 0L )

    condition.curr <- calculate.condition.number( design )

    # How far the signature phase moves each row on a control where there is
    # nothing to correct. A proposed change that is not several times this is
    # the estimator moving the row rather than the data.
    deg.bias <- stats::setNames( rep( NA_real_, fluorophore.n ), fluorophores )

    if ( !is.null( null.spectra ) ) {

      shared <- intersect( fluorophores, rownames( null.spectra ) )

      if ( length( shared ) > 0 ) {

        a <- null.spectra[ shared, colnames( spectra ), drop = FALSE ]
        b <- spectra[ shared, , drop = FALSE ]

        cos.bias <- rowSums( a * b ) /
          pmax( sqrt( rowSums( a^2 ) ) * sqrt( rowSums( b^2 ) ),
                .Machine$double.eps )

        deg.bias[ shared ] <-
          180 / pi * acos( pmin( 1, pmax( -1, cos.bias ) ) )
      }
    }

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
          active         = fluorophores,
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

      step.impact <- NA_real_
      bias.impact <- NA_real_

      if ( is.na( reject ) ) {

        unmixing.curr <- solve( tcrossprod( design ), design )[
          fluorophores, , drop = FALSE ]

        step.impact <- abundance.high[ j ] * sqrt( sum(
          ( ( candidate$signature - spectra.new[ j, ] ) %*%
              t( unmixing.curr ) )^2 ) )

        if ( !is.null( null.spectra ) && j %in% rownames( null.spectra ) )
          bias.impact <- abundance.high[ j ] * sqrt( sum(
            ( ( null.spectra[ j, colnames( spectra ) ] - spectra[ j, ] ) %*%
                t( unmixing.curr ) )^2 ) )

        st <- candidate$stats

        # Gates are ordered cheapest and most decisive first, and the first one
        # to fire names the reason. Under the joint fit the fluorophore's own
        # term is a partial slope, so a collinear neighbour can take a share of
        # `explained` without the decomposition being wrong; `explained.total`
        # is the one that must stay near one.
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
                        if ( gate.on.bias &&
                             is.finite( step.impact ) && is.finite( bias.impact ) &&
                             step.impact < min.impact.ratio * bias.impact )
                          reject <- "bias" else
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

      leak.before <- NA_real_
      leak.after  <- NA_real_
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

        # .fix.leakage() subsamples internally above its own max.events;
        # matching the seed for both calls means "before" and "after" score
        # the same held-out events, so a difference reflects the design
        # change and not which random subset happened to be drawn each call.
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

      # The statistics block is taken from `extract.raw.signature()` when there
      # is one and built from the same column set when there is not, so a change
      # to what that function reports cannot silently break the row bind.
      stats.block <- if ( is.null( candidate ) )
        .fix.empty.signature.stats() else
          candidate$stats[ , setdiff( names( candidate$stats ),
                                      c( "fluorophore", "n.events" ) ),
                           drop = FALSE ]

      log.rows[[ length( log.rows ) + 1L ]] <- data.frame(
        fluorophore = j,
        n.events    = length( idx ),
        stats.block,
        deg.bias        = unname( deg.bias[ j ] ),
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
      title = "FixMyUnmix_spillover_heatmap",
      plot.dir = asp$fix.unmixing.dir,
      legend.label = "Spillover",
      save = TRUE
    )

  if ( save ) {
    utils::write.csv( spillover.curr,
                      file = file.path( asp$fix.unmixing.dir,
                                        asp$fix.spillover.filename ) )
    utils::write.csv( compensation.curr,
                      file = file.path( asp$fix.unmixing.dir,
                                        asp$fix.compensation.filename ) )
    utils::write.csv( spectra.new,
                      file = file.path( asp$fix.unmixing.dir,
                                        asp$fix.spectra.filename ) )
  }

  list(
    spectra            = spectra.new,
    spectra.backsolved = spectra.backsolved,
    spillover          = spillover.curr,
    compensation       = compensation.curr,
    trust              = trust,
    coefficient.log    = coefficient.log,
    signature.log      = signature.log,
    convergence.log    = convergence.log,
    af.basis           = background.basis,
    af.hotspot         = af.hotspot,
    af.frozen          = af.frozen,
    leakage.prior      = prior
  )
}


#' Residual spillover slope for one source-target pair, by two estimators.
#'
#' `slope.truncated` is a robust line through the target-negative events, with
#' no binning and no zero-abundance anchor. Negativity is decided on the target
#' with this pair's current spillover estimate already removed, and the two are
#' iterated together. Deciding it on the raw target caps the fit's source range
#' at threshold over coefficient: the larger the residual spillover, the more of
#' the bright source population it pushes across the target's boundary, until
#' the events carrying the answer are the ones excluded and the pair reports
#' nothing at all. The estimator is then least able to measure exactly the
#' coefficients that matter most.
#'
#' `slope` is the binned envelope. Source-positive events are binned by source
#' abundance, everything below the source threshold is collapsed into a single
#' zero-abundance anchor, and in each bin the envelope is a fixed quantile of
#' the target-negative events only. It is retained as an alternative and for the
#' statistics the binning produces. When the population cannot support the
#' binning, the truncated estimate is still returned; the binned quantities come
#' back as `NA`.
#'
#' `coverage` is the fraction of target-negative events in the brightest third
#' of the source-positive population, measured on the corrected target, and is
#' the identifiability test for co-expression. A target positive only because
#' the source spills into it returns to the negative set once the estimate is
#' removed; a target the cell genuinely expresses does not, so when `coverage`
#' collapses the pair carries no information about spillover whichever slope is
#' taken. `slope.alt` is fitted to the plain median of the whole bin, so
#' disagreement between the two reports co-expression rather than instability.
#' @noRd
.fix.envelope.slope <- function( x.source, x.target,
                                 threshold.source, threshold.target,
                                 spread.var, neg.var,
                                 source.mask          = NULL,
                                 quantiles            = c( 0.05, 0.5 ),
                                 n.levels             = 10L,
                                 min.events           = 200L,
                                 min.bin.negative     = 25L,
                                 spread.addback       = TRUE,
                                 anchor.weight        = 1,
                                 max.truncated.events = 20000L,
                                 max.coefficient      = 0.2,
                                 max.mask.passes      = 3L,
                                 mask.tolerance       = 0.05,
                                 start.slope          = 0 ) {

  n <- length( x.source )
  if ( n < min.events ) return( NULL )

  if ( length( threshold.source ) == 1L )
    threshold.source <- rep( threshold.source, n )
  if ( length( threshold.target ) == 1L )
    threshold.target <- rep( threshold.target, n )

  # An event that clears the source threshold while another fluorophore is
  # brighter carries that other fluorophore's abundance read through a spectral
  # error, not this source's. It corrupts both estimators the same way, so the
  # filter is applied before either: for the robust line it removes exactly the
  # high-leverage points a wrong dye supplies, which is how a dim source
  # acquires a coefficient of order one.
  if ( !is.null( source.mask ) ) {

    keep.event <- source.mask | x.source <= threshold.source

    x.source         <- x.source[ keep.event ]
    x.target         <- x.target[ keep.event ]
    threshold.source <- threshold.source[ keep.event ]
    threshold.target <- threshold.target[ keep.event ]

    n <- length( x.source )
    if ( n < min.events ) return( NULL )
  }

  # Everything above the source threshold is kept and the negative bulk below it
  # is subsampled. The bulk sits at the origin and carries no leverage, so the
  # slope is unchanged while the robust fit's cost stops scaling with file size.
  select.negative <- function( slope ) {

    index <- which( x.target - slope * x.source < threshold.target )

    if ( length( index ) <= max.truncated.events ) return( index )

    bright <- index[ x.source[ index ] > threshold.source[ index ] ]
    bulk   <- setdiff( index, bright )
    n.bulk <- max( max.truncated.events - length( bright ), min.events )

    if ( length( bulk ) > n.bulk ) bulk <- sample( bulk, n.bulk )

    c( bright, bulk )
  }

  fit.truncated <- function( index, start = NULL ) {

    if ( length( index ) < min.events ) return( NA_real_ )

    .fix.huber.slope( x.source[ index ], x.target[ index ], start = start )[ 2 ]
  }

  start.slope <- if ( is.finite( start.slope ) &&
                      abs( start.slope ) <= max.coefficient ) start.slope else 0

  truncated.index <- select.negative( start.slope )
  slope.truncated <- fit.truncated( truncated.index, start = c( 0, start.slope ) )

  for ( pass in seq_len( as.integer( max.mask.passes ) ) ) {

    if ( !is.finite( slope.truncated ) ||
         abs( slope.truncated ) > max.coefficient ) break

    index.next <- select.negative( slope.truncated )
    slope.next <- fit.truncated( index.next, start = c( 0, slope.truncated ) )

    if ( !is.finite( slope.next ) ) break

    settled <- abs( slope.next - slope.truncated ) <=
      mask.tolerance * abs( slope.next )

    slope.truncated <- slope.next
    truncated.index <- index.next

    if ( settled ) break
  }

  slope.mask <- if ( is.finite( slope.truncated ) &&
                     abs( slope.truncated ) <= max.coefficient )
    slope.truncated else 0

  target.negative <- x.target - slope.mask * x.source < threshold.target

  source.positive <- x.source > threshold.source

  coverage <- 1

  if ( sum( source.positive ) >= min.bin.negative ) {

    bright.index <- which( source.positive )
    bright.cut   <- stats::quantile( x.source[ bright.index ], probs = 2 / 3,
                                     names = FALSE )
    top.index    <- bright.index[ x.source[ bright.index ] >= bright.cut ]

    if ( length( top.index ) >= min.bin.negative )
      coverage <- mean( target.negative[ top.index ] )
  }

  span.truncated <- if ( length( truncated.index ) > 1L )
    diff( range( x.source[ truncated.index ] ) ) else 0

  # The truncated estimator needs none of the binning, so a population that
  # cannot support the bins must not take the whole pair down with it.
  truncated.only <- function() {

    if ( !is.finite( slope.truncated ) ) return( NULL )

    list( slope           = NA_real_,
          slope.truncated = slope.truncated,
          se              = NA_real_,
          slope.alt       = NA_real_,
          disagreement    = NA_real_,
          coverage        = coverage,
          n               = n,
          span            = span.truncated,
          noise           = sqrt( max( neg.var, 0 ) ) )
  }

  if ( sum( source.positive ) < min.events ) return( truncated.only() )
  if ( sum( !source.positive ) < min.bin.negative ) return( truncated.only() )

  brk <- unique( stats::quantile(
    x.source[ source.positive ],
    probs = seq( 0, 1, length.out = n.levels + 1 ), names = FALSE ) )
  if ( length( brk ) < 4 ) return( truncated.only() )

  bin <- integer( n )
  bin[ source.positive ] <- as.integer( cut(
    x.source[ source.positive ], breaks = brk, include.lowest = TRUE ) )

  bins <- sort( unique( bin ) )
  if ( length( bins ) < 5 ) return( truncated.only() )

  # Each bin's row indices are split out once and reused for every per-bin
  # statistic below, rather than re-scanning `bin == b` over all n events
  # once per bin for each of centre / neg.n / bin.n / envelope.value /
  # median.value in turn.
  idx.by.bin <- split( seq_len( n ), bin )
  idx.by.bin <- idx.by.bin[ as.character( bins ) ]

  centre <- unname( vapply( idx.by.bin, function( idx )
    stats::median( x.source[ idx ] ), numeric( 1 ) ) )
  neg.n <- unname( vapply( idx.by.bin, function( idx )
    sum( target.negative[ idx ] ), integer( 1 ) ) )
  bin.n <- unname( lengths( idx.by.bin ) )

  usable <- neg.n >= min.bin.negative
  if ( sum( usable ) < 4L ) return( truncated.only() )

  sd.bin <- sqrt( pmax( neg.var + spread.var * pmax( centre, 0 ), 0 ) )
  z      <- if ( spread.addback ) -stats::qnorm( quantiles ) else c( 0, 0 )

  # Bins are weighted by the precision of their own quantile, events over
  # variance, and the zero-abundance anchor is capped at the best positive bin.
  # Uncapped, the anchor holds most of the file and sets the slope on its own,
  # which turns the fit into a comparison between two different populations
  # rather than a trend within one.
  anchor.cap <- function( w ) {
    if ( !any( bins > 0L ) || !any( bins == 0L ) ) return( w )
    w[ bins == 0L ] <- pmin( w[ bins == 0L ],
                             anchor.weight * max( w[ bins > 0L ] ) )
    w
  }

  weight.envelope <- anchor.cap( neg.n / pmax( sd.bin^2, .Machine$double.eps ) )
  weight.compare  <- anchor.cap( bin.n / pmax( sd.bin^2, .Machine$double.eps ) )

  fit.trace <- function( value, weight ) {

    keep <- usable & is.finite( value )
    if ( sum( keep ) < 4L ) return( list( slope = NA_real_, se = NA_real_ ) )

    x.bin <- centre[ keep ]
    y.bin <- value[ keep ]
    w.bin <- weight[ keep ]

    # Closed-form weighted least squares for the two-parameter fit, in place
    # of stats::lm.wfit(). At n.levels.pair bins a QR decomposition buys
    # nothing here, and dropping it removes this function's only call into
    # compiled linear algebra, which matters both for its own overhead, run
    # tens of thousands of times per correction, and for safety under
    # mclapply on macOS, where forking a process that has touched
    # Accelerate's threaded BLAS is unreliable.
    x.mean <- sum( w.bin * x.bin ) / sum( w.bin )
    y.mean <- sum( w.bin * y.bin ) / sum( w.bin )
    sxx    <- sum( w.bin * ( x.bin - x.mean )^2 )
    sxy    <- sum( w.bin * ( x.bin - x.mean ) * ( y.bin - y.mean ) )

    slope     <- if ( sxx > 0 ) sxy / sxx else NA_real_
    intercept <- y.mean - slope * x.mean
    residuals <- y.bin - ( intercept + slope * x.bin )

    dof <- sum( keep ) - 2L

    se <- if ( is.finite( slope ) && dof > 0 && sxx > 0 )
      sqrt( sum( w.bin * residuals^2 ) / dof / sxx ) else NA_real_

    list( slope = slope, se = se )
  }

  envelope.value <- unname( vapply( idx.by.bin, function( idx ) {
    v <- x.target[ idx ][ target.negative[ idx ] ]
    if ( length( v ) < min.bin.negative ) return( NA_real_ )
    stats::quantile( v, probs = quantiles[ 1 ], names = FALSE )
  }, numeric( 1 ) ) )

  median.value <- unname( vapply( idx.by.bin, function( idx )
    stats::quantile( x.target[ idx ], probs = quantiles[ 2 ], names = FALSE ),
    numeric( 1 ) ) )

  envelope <- fit.trace( envelope.value + z[ 1 ] * sd.bin, weight.envelope )
  compare  <- fit.trace( median.value   + z[ 2 ] * sd.bin, weight.compare )

  if ( !is.finite( envelope$slope ) && !is.finite( slope.truncated ) )
    return( NULL )

  disagreement <- abs( envelope$slope - compare$slope ) /
    ( abs( envelope$slope ) + abs( compare$slope ) + .Machine$double.eps )

  list( slope           = envelope$slope,
        slope.truncated = slope.truncated,
        se              = envelope$se,
        slope.alt       = compare$slope,
        disagreement    = disagreement,
        coverage        = coverage,
        n               = n,
        span            = max( centre ) - min( centre ),
        noise           = sqrt( max( neg.var, 0 ) ) )
}


#' An empty statistics block matching what `extract.raw.signature()` reports,
#' used for fluorophores whose population could not support a fit so that every
#' row of the signature log carries the same columns.
#' @noRd
.fix.empty.signature.stats <- function() {
  data.frame(
    n.bins           = NA_integer_,
    n.active         = NA_integer_,
    joint            = NA,
    x.span           = NA_real_,
    explained        = NA_real_,
    explained.total  = NA_real_,
    resid.rel        = NA_real_,
    intercept.rel    = NA_real_,
    bg.align         = NA_real_,
    clamp.frac       = NA_real_,
    anchor.rel       = NA_real_,
    vif.target       = NA_real_,
    vif.max          = NA_real_,
    vif.partner      = NA_character_,
    deg.change       = NA_real_,
    peak.curr        = NA_character_,
    peak.new         = NA_character_,
    row.names        = NULL,
    stringsAsFactors = FALSE )
}


#' Median off-target abundance within a fluorophore's own population, as a
#' fraction of its own signal. Needs no ground truth, so it is the acceptance
#' test that survives into production.
#' @noRd
.fix.leakage <- function( raw.data, design, target, panel, max.events = 20000L ) {

  if ( nrow( raw.data ) > max.events )
    raw.data <- raw.data[ sample.int( nrow( raw.data ), max.events ), ,
                          drop = FALSE ]

  unmixed <- unmix.ols.fast( raw.data, design )
  colnames( unmixed ) <- rownames( design )

  on <- stats::median( unmixed[ , target ] )
  if ( !is.finite( on ) || on <= 0 ) return( NA_real_ )

  off <- setdiff( panel, target )
  if ( length( off ) == 0 ) return( NA_real_ )

  sum( abs( apply( unmixed[ , off, drop = FALSE ], 2, stats::median ) ) ) / on
}


#' Stratified subsample of raw events by dominant fluorophore, for
#' downsampling before the correction loops.
#'
#' A uniform subsample keeps every population in the same proportion the raw
#' file already had it in. In a fully stained sample that means the biggest
#' population, usually events with most of the panel silent, sets how hard
#' everything else gets thinned; a dim or rare fluorophore's own positive
#' population is cut at the same rate as a background bulk that was already
#' more than large enough to define its threshold, and it is exactly that
#' population the pair estimator and the signature fit need the most events
#' from.
#'
#' Membership is by positivity, not by a single winner-take-all dominance
#' label: an event above threshold for two fluorophores at once belongs to
#' both strata, not whichever scores higher. That matters specifically for a
#' fully stained sample, where co-expression is real and is exactly the
#' population downstream identifiability checks rely on; a mutually-exclusive
#' assignment would credit a co-positive event to only one of the two markers
#' and let a large stratum's subsampling discard it from the other's
#' population by chance. On a concatenated single-stain pool the distinction
#' is close to moot, since the spread-scaled boundary already keeps one
#' tube's events out of another's positive population except right at the
#' edge, so this reduces to the same partition either way there.
#'
#' Every fluorophore's stratum is sampled against its own quota: a floor that
#' keeps a dim dye's whole positive population intact whenever it is smaller
#' than the floor, and a proportional share of what is left of the budget for
#' strata that can use more. The chosen events are the union across
#' fluorophores, so a co-positive event drawn by more than one stratum is
#' kept once and effectively has a higher retention chance, which is the
#' right direction. Events positive for nothing are capped separately as a
#' flat fraction of the total, since additional events there mostly refine a
#' threshold that is already well determined.
#'
#' @param abundance Numeric matrix (events x fluorophores), a first-pass
#'   unmixed abundance under the spectra downsampling will feed into.
#' @param thresholds Named numeric vector, flat per-fluorophore positivity
#'   thresholds in the same units as `abundance`.
#' @param threshold.matrix Numeric matrix, same shape as `abundance`, the
#'   spread-scaled positivity boundary from `get.spread.thresholds()`.
#' @param n.total Integer, the total event budget.
#' @param background.frac Numeric in (0, 1), the share of `n.total` reserved
#'   for events positive for nothing. Default `0.3`.
#' @param min.stratum Integer, the floor below which a fluorophore's positive
#'   population is kept whole rather than thinned. Default `2000`.
#'
#' @return Integer vector of row indices to keep.
#'
#' @noRd
.fix.stratified.sample <- function( abundance, thresholds, threshold.matrix,
                                    n.total, background.frac = 0.3,
                                    min.stratum = 2000L ) {

  n <- nrow( abundance )
  if ( n <= n.total ) return( seq_len( n ) )

  positive <- ( abundance - threshold.matrix ) > 0

  background.idx <- which( rowSums( positive ) == 0L )
  positive.idx    <- lapply( seq_len( ncol( abundance ) ),
                             function( j ) which( positive[ , j ] ) )

  n.background <- min( length( background.idx ),
                       round( background.frac * n.total ) )
  budget.positive <- n.total - n.background

  sizes   <- vapply( positive.idx, length, integer( 1 ) )
  floor.n <- pmin( sizes, as.integer( min.stratum ) )

  # A panel with many fluorophores can demand more floor than the budget
  # holds; when it does, every stratum's floor is scaled down by the same
  # factor rather than starving whichever strata come later in column order.
  quota <- if ( sum( floor.n ) > budget.positive ) {
    pmin( sizes, floor( floor.n * budget.positive / max( sum( floor.n ), 1 ) ) )
  } else {
    spare           <- pmax( sizes - floor.n, 0 )
    leftover.budget <- budget.positive - sum( floor.n )
    extra <- if ( sum( spare ) > 0 )
      leftover.budget * spare / sum( spare ) else rep( 0, length( sizes ) )
    pmin( sizes, floor.n + round( extra ) )
  }

  chosen.positive <- unique( unlist( Map( function( idx, k )
    if ( length( idx ) <= k ) idx else sample( idx, k ),
    positive.idx, quota ), use.names = FALSE ) )

  chosen.background <- if ( n.background >= length( background.idx ) )
    background.idx else sample( background.idx, n.background )

  c( chosen.positive, chosen.background )
}

#' Closed-form Huber M-estimator for a two-parameter (intercept + slope) fit,
#' in place of fit.robust.linear.model()'s call to MASS::rlm(). Same default
#' behaviour: weights from psi.huber with tuning constant k = 1.345, scale
#' re-estimated each iteration from the MAD of the current residuals,
#' iterating to a coefficient-change tolerance, initial weights of 1 so the
#' first fit is OLS. On non-convergence it returns a zero slope rather than
#' falling back to OLS, matching fit.robust.linear.model()'s fix.unmix = TRUE
#' behaviour, since a coefficient this function cannot pin down is exactly
#' the case the caller needs to see as untrustworthy rather than as some
#' other estimator's best guess.
#'
#' Every step is the same closed-form weighted least squares fit.trace() in
#' `.fix.envelope.slope()` already uses, so this makes no call to
#' stats::lm.wfit(), solve(), or any other compiled linear-algebra routine,
#' which is the second and larger BLAS dependency this pair loop had.
#'
#' @param x,y Numeric vectors, the predictor and response.
#' @param k Numeric, the Huber tuning constant. Default `1.345`.
#' @param max.iter Integer, maximum IRLS iterations. Default `100`.
#' @param tol Numeric, relative coefficient-change convergence tolerance.
#'   Default `1e-4`.
#' @param start Optional numeric `c(intercept, slope)`, a warm start. When
#'   supplied, the first iteration's weights come from its residuals instead
#'   of OLS weights, so a caller iterating the same pair across outer passes
#'   can resume near its last estimate. Default `NULL`.
#'
#' @return Numeric vector `c(intercept, slope)`.
#'
#' @noRd
.fix.huber.slope <- function( x, y, k = 1.345, max.iter = 100L, tol = 1e-4,
                              start = NULL ) {

  if (
    requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    "fix_huber_slope_rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
  ) {
    AutoSpectralRcpp::fix_huber_slope_rcpp(
      x, y, k = k, max_iter = as.integer( max.iter ),
      tol = tol, start = start
    )
  } else {
    stop( "Install AutoSpectralRcpp" )
  }

}
