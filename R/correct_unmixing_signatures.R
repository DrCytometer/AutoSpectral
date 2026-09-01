# correct_unmixing_signatures.R

#' @title Correct Unmixing Signatures
#'
#' @description
#' Detects and corrects shape errors in a reference `spectra` matrix directly
#' from stained sample data, for the situation where the available reference
#' spectra are subtly wrong for the sample being analysed - typically because
#' bead-derived spectra are applied to cells, or spectra from another day,
#' lot, or instrument state are in use.
#'
#' Each event is assigned to the fluorophore it is most strongly positive
#' for, as a fraction of that fluorophore's own dynamic range above its
#' positivity threshold. Within each such dominance population, the data are
#' re-unmixed against a restricted design (the dominant fluorophore plus any
#' genuinely co-active ones), binned by abundance, anchored with the
#' background population, and the detector-space residual is regressed on the
#' restricted abundances. The dominant fluorophore's regression slope
#' estimates the error in its spectrum, and the row is updated by a step
#' chosen through a held-out search.
#'
#' A correction is only applied when it passes every acceptance gate:
#' \itemize{
#'   \item the population spans enough abundance to identify a slope;
#'   \item the fluorophore's own term explains most of the
#'     background-subtracted signal at its brightest (`min.explained`), so a
#'     control too dim to characterise is left untouched;
#'   \item the held-out step search finds a step that reduces the restricted
#'     residual (`min.gain`);
#'   \item the correction is small relative to the row it corrects
#'     (`max.step`) and does not inflate the fluorophore's own apparent
#'     abundance (`max.span.drift`);
#'   \item the correction is not a background confound: when the apparent
#'     abundance trend is a common-mode background residual rather than a
#'     spectral error, the event-level regression splits that one physical
#'     direction into an intercept and a slope that are anti-collinear, and
#'     the correction is rejected (`max.bg.alignment`).
#' }
#' Rejected fluorophores keep their starting spectra; the method is designed
#' to leave a row untouched rather than risk making it worse.
#'
#' Autofluorescence is handled by background removal before fitting rather
#' than as a design row, so a fluorophore's abundance and the background
#' amount are not estimated from the same spectral vector. For cells,
#' scatter-matched per-event subtraction (`bg.mode = "scatter.knn"`) uses an
#' unstained control matched on scatter, an independent measurement channel.
#' For beads and other uniform particles, whose scatter is uninformative,
#' the mean of the sample's own background events is subtracted
#' (`bg.mode = "global.mean"`).
#'
#' When `scatter` is supplied and `gate.main` is `TRUE`, events are first
#' gated to the main scatter population by density, removing noise, debris
#' and most aggregates, whose distinct autofluorescence otherwise
#' contaminates the dominance populations and the background estimate.
#'
#' @param raw.data Numeric matrix (events x detectors), raw detector-space
#'   data. Pooled or concatenated single-stained controls, or a fully
#'   stained sample with well-separated populations. Columns must match the
#'   columns of `spectra`.
#' @param spectra Numeric matrix (fluorophores x detectors), the starting
#'   reference spectra to be corrected, L-infinity normalised.
#' @param unmixed.thresholds Named numeric vector covering every fluorophore
#'   in `spectra` (the autofluorescence row may be omitted), giving the
#'   positivity threshold in unmixed space, typically the 99.5th percentile
#'   of an unstained control unmixed against `spectra`.
#' @param asp Optional AutoSpectral parameter list from
#'   `get.autospectral.param()`. Used only to seed the random number
#'   generator (`asp$bird.seed`) for reproducible subsampling. Default
#'   `NULL`.
#' @param af.name Character, the name of the autofluorescence row in
#'   `spectra`, or `NULL` if there is none. The AF row is never treated as a
#'   panel fluorophore and is never corrected. Default `"AF"`.
#' @param scatter Optional numeric matrix (events x scatter parameters),
#'   row-matched to `raw.data`. Required for `gate.main` and for
#'   `bg.mode = "scatter.knn"`. Default `NULL`.
#' @param gate.main Logical, whether to gate events to the main scatter
#'   population before fitting. Requires `scatter`. Default `TRUE`.
#' @param gate.level Numeric in (0, 1). Events are kept when their 2D
#'   scatter density exceeds this fraction of the modal density. Default
#'   `0.1`.
#' @param spillover.spread Optional matrix from
#'   `get.spectral.variants()$spillover.spread`. When supplied, co-activity
#'   is judged against per-event thresholds widened by the spillover spread
#'   each bright fluorophore contributes, so spillover from the dominant dye
#'   is not mistaken for a co-active fluorophore. Default `NULL` (flat
#'   thresholds).
#' @param spread.kappa Numeric, how many spillover-spread standard
#'   deviations above the flat threshold still count as negative. Default
#'   `2`.
#' @param bg.mode Character, background removal mode: `"global.mean"`
#'   (default), `"scatter.knn"`, or `"none"`. See Description.
#' @param unstained Optional numeric matrix (events x detectors), raw
#'   unstained control, required for `bg.mode = "scatter.knn"`.
#' @param unstained.scatter Optional numeric matrix (events x scatter
#'   parameters), row-matched to `unstained`, required for
#'   `bg.mode = "scatter.knn"`. When `gate.main` is `TRUE` the unstained
#'   control is gated with the same density rule.
#' @param k.neighbors Integer, neighbours for scatter-matched background
#'   subtraction. Default `20`.
#' @param n.levels Integer, abundance bins per dominance population.
#'   Default `10`.
#' @param n.iter Integer, maximum correction iterations. Iteration stops
#'   early when no fluorophore accepts a step. Default `6`.
#' @param min.events Integer, minimum events in a dominance population.
#'   Default `200`.
#' @param min.span Numeric, minimum abundance span in units of the
#'   fluorophore's own threshold magnitude. Default `5`.
#' @param min.explained Numeric in (0, 1), minimum fraction of the
#'   background-subtracted signal the fluorophore's own term must account
#'   for at its brightest bin. Default `0.5`.
#' @param min.gain Numeric, minimum held-out relative residual reduction.
#'   Default `0.002`.
#' @param step.grid Numeric vector of candidate step sizes for the held-out
#'   search. Default `c( 0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1 )`; the
#'   small steps matter for tandem dyes, whose abundance-dependent variant
#'   mixture makes the residual objective a narrow valley.
#' @param n.split.trials Integer, number of independent random 50/50 splits
#'   used by the held-out step search. Default `1`, which reproduces the
#'   original fixed, un-reseeded split exactly. Raising this trades a
#'   single split's noise for a vote across `n.split.trials` splits (see
#'   `min.split.frac`); useful for fluorophores whose true correction is
#'   real but small relative to per-event noise, where a lone 50/50 split
#'   can land on the unlucky side.
#' @param min.split.frac Numeric in (0, 1], the minimum fraction of
#'   `n.split.trials` splits that must independently find a beneficial step
#'   before one is accepted. Ignored when `n.split.trials = 1`. Default
#'   `0.6`.
#' @param max.step Numeric, maximum norm of the correction relative to the
#'   norm of the row it corrects. Larger proposed corrections are rejected
#'   outright rather than scaled down. Default `0.15`.
#' @param max.span.drift Numeric, maximum allowed growth of the
#'   fluorophore's apparent abundance span across iterations. Default
#'   `1.10`.
#' @param max.bg.alignment Numeric in `[-1, 0]`. The correction is rejected
#'   when the cosine between the event-level regression intercept and slope
#'   falls at or below this value. Default `-0.9`.
#' @param nuisance.frac Numeric in (0, 1), fraction of a dominance
#'   population that must be co-active for another fluorophore before it is
#'   carried as a nuisance term in the restricted design. Default `0.5`.
#' @param footprint.frac Numeric in `[0, 1)`. Restricts the slope fit and
#'   the held-out step search to detectors where the dominant dye's own
#'   current spectrum exceeds `footprint.frac` of its own peak. A dye
#'   cannot carry real shape-error signal in a detector it does not
#'   meaningfully emit into; for a narrow-emission dye read out on a wide
#'   detector array, those channels only dilute the held-out residual
#'   objective with noise from channels that are pure background for that
#'   dye. Abundance estimation and the `explained`/`bg.align` gates are
#'   unaffected; only the slope fit's response and the held-out
#'   objective's norm are restricted. `0` reproduces the previous,
#'   unrestricted behaviour exactly. Default `0.02`.
#' @param footprint.min.channels Integer, minimum detectors the
#'   `footprint.frac` mask must keep; below this the mask is dropped and
#'   every detector is used, so a pathologically narrow spectrum cannot
#'   leave too few channels to fit. Default `3L`.
#' @param background.n Integer, maximum background events used for the
#'   zero-abundance anchor. Default `5000`.
#' @param true.spectra Optional numeric matrix (fluorophores x detectors),
#'   independently-known ground truth with row names matching `spectra`.
#'   Purely diagnostic: when supplied, the returned `recovery` table reports
#'   the angular error before and after correction per fluorophore.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`spectra`}{The corrected reference spectra matrix. Rows that
#'     failed any acceptance gate are unchanged.}
#'   \item{`fit.log`}{Data frame, one row per fluorophore per iteration,
#'     with the fit statistics and every gate quantity, including
#'     `bg.align`.}
#'   \item{`accepted`}{Named logical vector, whether each panel fluorophore
#'     accepted at least one correction step.}
#'   \item{`dominant`}{Integer vector over the (gated) events: the index
#'     into `panel` of each event's dominant fluorophore, `0` for
#'     background.}
#'   \item{`panel`}{Character vector, the fluorophores eligible for
#'     correction (`spectra` rows minus `af.name`).}
#'   \item{`gate.keep`}{Logical vector over the input rows of `raw.data`,
#'     `TRUE` for events inside the main-population gate, or `NULL` if no
#'     gating was applied.}
#'   \item{`recovery`}{Data frame of angular errors against `true.spectra`,
#'     or `NULL` if `true.spectra` was not supplied.}
#' }
#'
#' @importFrom MASS bandwidth.nrd ginv kde2d
#' @importFrom stats coef lm.fit mad quantile sd
#'
#' @export

correct.unmixing.signatures <- function(
    raw.data,
    spectra,
    unmixed.thresholds,
    asp                = NULL,
    af.name            = "AF",
    scatter            = NULL,
    gate.main          = TRUE,
    gate.level         = 0.1,
    spillover.spread   = NULL,
    spread.kappa       = 2,
    bg.mode            = c( "global.mean", "scatter.knn", "none" ),
    unstained          = NULL,
    unstained.scatter  = NULL,
    k.neighbors        = 20L,
    n.levels           = 10L,
    n.iter             = 6L,
    min.events         = 200L,
    min.span           = 5,
    min.explained      = 0.5,
    min.gain           = 0.002,
    step.grid          = c( 0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1 ),
    n.split.trials     = 1L,
    min.split.frac     = 0.6,
    max.step           = 0.15,
    max.span.drift     = 1.10,
    max.bg.alignment   = -0.9,
    nuisance.frac      = 0.5,
    footprint.frac        = 0.02,
    footprint.min.channels = 3L,
    background.n       = 5000L,
    true.spectra       = NULL,
    verbose            = TRUE
) {

  bg.mode <- match.arg( bg.mode )

  spectra  <- as.matrix( spectra )
  raw.data <- as.matrix( raw.data )

  if ( is.null( rownames( spectra ) ) )
    stop( "`spectra` must have fluorophore row names.", call. = FALSE )

  if ( ncol( raw.data ) != ncol( spectra ) )
    stop( "`raw.data` and `spectra` must have the same detector columns.",
          call. = FALSE )

  panel <- setdiff( rownames( spectra ), af.name )

  if ( !all( panel %in% names( unmixed.thresholds ) ) )
    stop( "`unmixed.thresholds` must be named and cover every panel fluorophore.",
          call. = FALSE )

  if ( bg.mode == "scatter.knn" &&
       ( is.null( scatter ) || is.null( unstained ) ||
         is.null( unstained.scatter ) ) )
    stop( paste0( "`bg.mode = \"scatter.knn\"` requires `scatter`, ",
                  "`unstained` and `unstained.scatter`." ), call. = FALSE )

  if ( !is.null( asp ) && !is.null( asp$bird.seed ) )
    set.seed( asp$bird.seed )

  # ---------------------------------------------------------------------------
  # Main-population gate
  # ---------------------------------------------------------------------------

  gate.keep <- NULL

  if ( gate.main && !is.null( scatter ) ) {

    gate.keep <- .signature.main.gate( scatter, gate.level = gate.level )

    if ( verbose )
      message( sprintf(
        "Main-population gate: kept %d of %d events (%.1f%%).",
        sum( gate.keep ), length( gate.keep ), 100 * mean( gate.keep ) ) )

    raw.data <- raw.data[ gate.keep, , drop = FALSE ]
    scatter  <- scatter[ gate.keep, , drop = FALSE ]

    if ( !is.null( unstained ) && !is.null( unstained.scatter ) ) {

      unst.keep <- .signature.main.gate( unstained.scatter,
                                         gate.level = gate.level )

      if ( verbose )
        message( sprintf(
          "Main-population gate, unstained: kept %d of %d events (%.1f%%).",
          sum( unst.keep ), length( unst.keep ), 100 * mean( unst.keep ) ) )

      unstained         <- unstained[ unst.keep, , drop = FALSE ]
      unstained.scatter <- unstained.scatter[ unst.keep, , drop = FALSE ]
    }
  }

  # ---------------------------------------------------------------------------
  # Dominance assignment and co-activity
  # ---------------------------------------------------------------------------

  unmixed <- unmix.ols.fast( raw.data, spectra )

  if ( is.null( spillover.spread ) ) {

    if ( verbose )
      message( paste0(
        "`spillover.spread` not supplied: co-activity thresholds are flat. ",
        "This measurably weakens nuisance-set assignment and gain detection ",
        "for some fluorophores; supplying get.spectral.variants()$spillover.spread ",
        "is recommended whenever available for this particle type." ) )

    above <- sweep( unmixed[ , panel, drop = FALSE ],
                    2, unmixed.thresholds[ panel ], ">" )
  } else {
    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed,
      thresholds       = unmixed.thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      verbose          = FALSE
    )
    above <- unmixed[ , panel, drop = FALSE ] >
      threshold.matrix[ , panel, drop = FALSE ]
  }

  # Dominance is scored as a fraction of each dye's own dynamic range above
  # its threshold, not as a multiple of the threshold: thresholds vary by
  # hundreds of fold across a panel and can be negative when taken against
  # the wrong spectra, so dividing by them hands every dim or ambiguous
  # event to whichever dye happens to have the smallest threshold.
  excess <- pmax( sweep( unmixed[ , panel, drop = FALSE ],
                         2, unmixed.thresholds[ panel ], "-" ), 0 )

  dyn.range <- apply( unmixed[ , panel, drop = FALSE ], 2, stats::quantile,
                      probs = 0.999 ) - unmixed.thresholds[ panel ]
  dyn.range <- pmax( dyn.range, .Machine$double.eps )

  score    <- sweep( excess, 2, dyn.range, "/" )
  dominant <- max.col( score, ties.method = "first" )
  top      <- score[ cbind( seq_len( nrow( score ) ), dominant ) ]
  dominant[ top <= 0 ] <- 0L

  background.idx <- which( dominant == 0L )

  # ---------------------------------------------------------------------------
  # Background removal
  # ---------------------------------------------------------------------------

  if ( bg.mode == "scatter.knn" ) {

    raw.data <- .signature.knn.subtract(
      raw.data, scatter, unstained, unstained.scatter,
      k.neighbors = k.neighbors )

  } else if ( bg.mode == "global.mean" &&
              length( background.idx ) >= min.events ) {

    raw.data <- sweep( raw.data, 2,
                       colMeans( raw.data[ background.idx, , drop = FALSE ] ),
                       "-" )
  }

  if ( length( background.idx ) > background.n )
    background.idx <- sample( background.idx, background.n )

  # ---------------------------------------------------------------------------
  # Correction loop
  # ---------------------------------------------------------------------------

  spectra.new <- spectra
  fit.log     <- list()
  span.first  <- stats::setNames( rep( NA_real_, length( panel ) ), panel )

  for ( iter in seq_len( n.iter ) ) {

    iter.accepted <- FALSE

    for ( f in seq_along( panel ) ) {

      j   <- panel[ f ]
      idx <- which( dominant == f )
      if ( length( idx ) < min.events ) next

      # Co-active dyes are carried as nuisance predictors, not discarded, so
      # a dye whose spectrally collinear partner is always co-positive still
      # gets a usable population.
      co.frac  <- colMeans( above[ idx, , drop = FALSE ] )
      nuisance <- setdiff( panel[ co.frac > nuisance.frac ], j )
      active   <- c( j, nuisance )

      # A dye cannot carry real shape-error signal in a detector it does not
      # meaningfully emit into; restricting the fit and the held-out
      # objective to its own footprint keeps those channels' noise from
      # diluting the signal for narrow-emission dyes. Falls back to every
      # detector if the mask would leave too few to fit.
      footprint.mask <- spectra.new[ j, ] > footprint.frac * max( spectra.new[ j, ] )
      if ( sum( footprint.mask ) < footprint.min.channels )
        footprint.mask <- rep( TRUE, ncol( spectra.new ) )

      # Slope from a chosen subset of the group's events: bin by abundance,
      # anchor with the background population, regress the restricted
      # residual on all restricted abundances, and take the dominant dye's
      # coefficients. The multivariate fit prevents a co-varying nuisance
      # dye's own spectral error from being booked against the dominant dye.
      fit.slope <- function( use.idx ) {

        if ( length( use.idx ) < min.events %/% 2 ) return( NULL )

        y.use <- raw.data[ use.idx, , drop = FALSE ]
        x.use <- .signature.restricted.unmix( y.use, spectra.new, active )[ , 1 ]

        brk <- unique( stats::quantile(
          x.use, probs = seq( 0, 1, length.out = n.levels + 1 ) ) )
        if ( length( brk ) < 3 ) return( NULL )

        bin        <- as.integer( cut( x.use, breaks = brk, include.lowest = TRUE ) )
        bins.use   <- sort( unique( bin ) )
        idx.by.bin <- split( seq_along( bin ), bin )
        idx.by.bin <- idx.by.bin[ as.character( bins.use ) ]
        y.bin <- t( vapply( idx.by.bin, function( idx )
          colMeans( y.use[ idx, , drop = FALSE ] ),
          numeric( ncol( y.use ) ) ) )

        if ( length( background.idx ) >= min.events )
          y.bin <- rbind(
            colMeans( raw.data[ background.idx, , drop = FALSE ] ), y.bin )

        if ( nrow( y.bin ) < length( active ) + 3 ) return( NULL )

        x.bin <- .signature.restricted.unmix( y.bin, spectra.new, active )
        r.bin <- y.bin - x.bin %*% spectra.new[ active, , drop = FALSE ]

        # The slope is always fit on every detector: footprint.frac sharpens
        # which channels decide whether and how far to step (the held-out
        # objective in residual.gain, below), not what the eventual
        # correction touches. Restricting the fit's response to the footprint
        # would permanently cap that row's achievable correction at whatever
        # fraction of its true error happens to fall inside the mask - a real
        # cost for a dye whose error is broad but small, and one that bites
        # hardest on already-well-corrected rows, where even a few percent of
        # uncorrected error left outside the mask is a large relative hit to
        # a small residual angle.
        fit.m <- stats::lm.fit( x = cbind( 1, x.bin ), y = r.bin )
        slope <- stats::coef( fit.m )[ 2, ]
        slope[ !is.finite( slope ) ] <- 0

        resid.m <- fit.m$residuals[ , footprint.mask, drop = FALSE ]
        ss.res  <- sum( resid.m^2 )
        ss.tot  <- sum( sweep( r.bin[ , footprint.mask, drop = FALSE ], 2,
                               colMeans( r.bin[ , footprint.mask, drop = FALSE ] ) )^2 )

        # At its brightest, how much of the background-subtracted signal
        # does this dye's own term account for? A control where the answer
        # is small is not really a control for this dye.
        top.bin   <- which.max( x.bin[ , 1 ] )
        top.norm  <- sqrt( sum( y.bin[ top.bin, ]^2 ) )
        explained <- if ( top.norm > 0 )
          sqrt( sum( ( x.bin[ top.bin, 1 ] * spectra.new[ j, ] )^2 ) ) /
          top.norm else 0

        list(
          slope     = slope,
          r.sq      = if ( ss.tot > 0 ) 1 - ss.res / ss.tot else 0,
          x.span    = max( x.bin[ , 1 ] ) - min( x.bin[ , 1 ] ),
          explained = explained,
          n.bins    = nrow( r.bin )
        )
      }

      full <- fit.slope( idx )
      if ( is.null( full ) ) next

      # Background-confound statistic. A genuine spectral error gives the
      # event-level regression an intercept and a slope with independent
      # physical origins. A common-mode background residual whose magnitude
      # tracks brightness is one physical direction split in two by the
      # regression, so the intercept and slope come out anti-collinear.
      # Strong anti-alignment marks a correction that would walk the row
      # into the background direction rather than fix its shape.
      y.evt <- raw.data[ idx, , drop = FALSE ]
      x.evt <- .signature.restricted.unmix( y.evt, spectra.new, active )
      r.evt <- y.evt - x.evt %*% spectra.new[ active, , drop = FALSE ]

      fit.evt   <- stats::lm.fit( x = cbind( 1, x.evt[ , 1 ] ), y = r.evt )
      alpha.evt <- stats::coef( fit.evt )[ 1, ]
      beta.evt  <- stats::coef( fit.evt )[ 2, ]
      alpha.evt[ !is.finite( alpha.evt ) ] <- 0
      beta.evt[ !is.finite( beta.evt ) ]   <- 0

      alpha.norm <- sqrt( sum( alpha.evt^2 ) )
      beta.norm  <- sqrt( sum( beta.evt^2 ) )
      bg.align   <- if ( alpha.norm > 0 && beta.norm > 0 )
        sum( alpha.evt * beta.evt ) / ( alpha.norm * beta.norm ) else NA_real_

      # Held-out step search: fit on one half, pick the step size that most
      # reduces the restricted residual on the other. This bounds the step
      # but cannot on its own tell a corrected row from a repurposed one -
      # hence the magnitude, drift and confound gates below.
      #
      # A single 50/50 split is a noisy instrument for a fluorophore whose
      # true correction is real but small relative to per-event noise: the
      # default n.split.trials = 1 reproduces the original fixed,
      # un-reseeded rep_len(c(1L, 2L), length(idx)) split exactly, for
      # backward compatibility. Raising n.split.trials instead draws that
      # many independent random 50/50 splits and requires only
      # min.split.frac of them to agree a step helps, taking the median of
      # the agreeing splits' step and gain - trading one noisy verdict for
      # a vote.
      residual.gain <- function( fit.idx, test.idx ) {

        fs <- fit.slope( fit.idx )
        if ( is.null( fs ) ) return( NULL )

        y.test <- raw.data[ test.idx, , drop = FALSE ]
        base   <- sqrt( sum( y.test[ , footprint.mask, drop = FALSE ]^2 ) )
        if ( base <= 0 ) return( NULL )

        obj <- vapply( step.grid, function( t ) {
          s.try <- spectra.new
          s.try[ j, ] <- pmax( s.try[ j, ] + t * fs$slope, 0 )
          if ( max( s.try[ j, ] ) <= 0 ) return( Inf )
          s.try <- .signature.renorm( s.try )
          x.try <- .signature.restricted.unmix( y.test, s.try, active )
          resid <- y.test - x.try %*% s.try[ active, , drop = FALSE ]
          sqrt( sum( resid[ , footprint.mask, drop = FALSE ]^2 ) ) / base
        }, numeric( 1 ) )

        best <- which.min( obj )
        list( t = step.grid[ best ], gain = obj[ 1 ] - obj[ best ] )
      }

      split.trial <- function( half ) {

        idx.a <- idx[ half == 1L ]
        idx.b <- idx[ half == 2L ]

        gain.ab <- residual.gain( idx.a, idx.b )
        gain.ba <- residual.gain( idx.b, idx.a )

        if ( is.null( gain.ab ) || is.null( gain.ba ) )
          return( c( t = 0, gain = NA_real_ ) )

        c( t = min( gain.ab$t, gain.ba$t ),
           gain = min( gain.ab$gain, gain.ba$gain ) )
      }

      if ( n.split.trials <= 1L ) {

        one   <- split.trial( rep_len( c( 1L, 2L ), length( idx ) ) )
        t.hat <- one[ "t" ]
        gain  <- one[ "gain" ]

      } else {

        trials <- t( vapply( seq_len( n.split.trials ), function( trial )
          split.trial( sample( rep_len( c( 1L, 2L ), length( idx ) ) ) ),
          numeric( 2 ) ) )

        passed <- trials[ , "t" ] > 0

        if ( any( passed ) && mean( passed ) >= min.split.frac ) {
          t.hat <- stats::median( trials[ passed, "t" ] )
          gain  <- stats::median( trials[ passed, "gain" ] )
        } else {
          t.hat <- 0
          gain  <- NA_real_
        }
      }

      slope     <- full$slope
      step.norm <- sqrt( sum( slope^2 ) )
      row.norm  <- sqrt( sum( spectra.new[ j, ]^2 ) )
      rel.step  <- step.norm / row.norm

      if ( is.na( span.first[ j ] ) ) span.first[ j ] <- full$x.span
      span.drift <- full$x.span / span.first[ j ]

      # A correction larger than the row it corrects is wrong on its face,
      # whatever it does to the residual, so a step that would need scaling
      # down to fit is rejected outright rather than applied in miniature.
      # Span drift catches a row rotating toward whatever else the sample
      # contains, which makes its own apparent abundance grow.
      accepted <- full$x.span > min.span * abs( unmixed.thresholds[ j ] ) &&
        is.finite( full$explained ) && full$explained > min.explained &&
        t.hat > 0 && is.finite( gain ) && gain > min.gain &&
        is.finite( rel.step ) && rel.step > 0 && rel.step <= max.step &&
        is.finite( span.drift ) && span.drift <= max.span.drift &&
        ( !is.finite( bg.align ) || bg.align > max.bg.alignment )

      if ( accepted ) {
        spectra.new[ j, ] <- pmax( spectra.new[ j, ] + t.hat * slope, 0 )
        iter.accepted     <- TRUE
      }

      fit.log[[ length( fit.log ) + 1L ]] <- data.frame(
        iter        = iter,
        fluorophore = j,
        n.events    = length( idx ),
        n.nuisance  = length( nuisance ),
        n.bins      = full$n.bins,
        x.span      = full$x.span,
        explained   = full$explained,
        r.squared   = full$r.sq,
        bg.align    = bg.align,
        t.hat       = t.hat,
        gain        = gain,
        rel.step    = rel.step,
        span.drift  = span.drift,
        accepted    = accepted,
        row.names   = NULL
      )
    }

    spectra.new <- .signature.renorm( spectra.new )

    if ( !iter.accepted ) break
  }

  fit.log <- if ( length( fit.log ) > 0 ) do.call( rbind, fit.log ) else NULL

  accepted <- vapply( panel, function( j ) {
    if ( is.null( fit.log ) ) return( FALSE )
    any( fit.log$accepted[ fit.log$fluorophore == j ] )
  }, logical( 1 ) )

  if ( verbose )
    message( sprintf(
      "Accepted corrections for %d of %d fluorophore(s)%s.",
      sum( accepted ), length( panel ),
      if ( sum( accepted ) > 0 )
        paste0( ": ", paste( panel[ accepted ], collapse = ", " ) ) else "" ) )

  # ---------------------------------------------------------------------------
  # Optional recovery diagnostics against ground truth
  # ---------------------------------------------------------------------------

  recovery <- NULL

  if ( !is.null( true.spectra ) ) {

    common <- intersect( panel, rownames( true.spectra ) )

    if ( length( common ) > 0 ) {

      deg.start <- .signature.row.angle(
        spectra[ common, , drop = FALSE ],
        true.spectra[ common, colnames( spectra ), drop = FALSE ] )
      deg.after <- .signature.row.angle(
        spectra.new[ common, , drop = FALSE ],
        true.spectra[ common, colnames( spectra ), drop = FALSE ] )

      recovery <- data.frame(
        fluorophore = common,
        deg.start   = deg.start,
        deg.after   = deg.after,
        recovered   = ( deg.start - deg.after ) / deg.start,
        accepted    = accepted[ common ],
        row.names   = NULL
      )
    }
  }

  list(
    spectra   = spectra.new,
    fit.log   = fit.log,
    accepted  = accepted,
    dominant  = dominant,
    panel     = panel,
    gate.keep = gate.keep,
    recovery  = recovery
  )
}


#' L-infinity renormalisation of spectra rows.
#' @noRd
.signature.renorm <- function( s ) {
  row.max <- apply( s, 1, max )
  row.max[ row.max <= 0 ] <- 1
  s / row.max
}


#' Row angles in degrees between two matched spectra matrices.
#' @noRd
.signature.row.angle <- function( a, b ) {
  a  <- as.matrix( a )
  b  <- as.matrix( b )
  cs <- rowSums( a * b ) / ( sqrt( rowSums( a^2 ) ) * sqrt( rowSums( b^2 ) ) )
  180 / pi * acos( pmin( 1, pmax( -1, cs ) ) )
}


#' Restricted unmix of a data matrix against a chosen fluorophore set,
#' falling back to a pseudoinverse when the restricted design is singular.
#' @noRd
.signature.restricted.unmix <- function( y, spectra, active ) {

  sub <- spectra[ active, , drop = FALSE ]

  sub.plus <- tryCatch(
    solve.default( tcrossprod( sub ), sub ),
    error = function( e ) MASS::ginv( tcrossprod( sub ) ) %*% sub
  )

  y %*% t( sub.plus )
}


#' Main-population scatter gate: keep events whose 2D scatter density
#' exceeds a fraction of the modal density. Density-relative rather than
#' polygonal, so it needs no per-run landmarks.
#' @noRd
.signature.main.gate <- function( scatter, gate.level = 0.1, grid.n = 128L,
                                  max.events = 100000L ) {

  sc <- as.matrix( scatter[ , 1:2, drop = FALSE ] )

  fit.idx <- seq_len( nrow( sc ) )
  if ( length( fit.idx ) > max.events )
    fit.idx <- sample( fit.idx, max.events )

  bw <- .safe.bandwidth( sc[ fit.idx, , drop = FALSE ] )
  bw <- pmax( bw, .Machine$double.eps )

  kde <- MASS::kde2d( sc[ fit.idx, 1 ], sc[ fit.idx, 2 ],
                      h = bw, n = grid.n )

  ix   <- findInterval( sc[ , 1 ], kde$x, all.inside = TRUE )
  iy   <- findInterval( sc[ , 2 ], kde$y, all.inside = TRUE )
  dens <- kde$z[ cbind( ix, iy ) ]

  dens >= gate.level * max( kde$z )
}


#' Scatter-matched per-event background subtraction. Each event has the mean
#' spectrum of its k nearest unstained neighbours in scatter space
#' subtracted. Matching on scatter rather than on the spectral profile is
#' the point: the amount subtracted is set by an independent measurement
#' channel, so it cannot be absorbed back into the fluorophore abundance the
#' way a spectrally-derived background estimate can.
#' @noRd
.signature.knn.subtract <- function( raw.data, scatter, unstained,
                                     unstained.scatter,
                                     k.neighbors = 20L,
                                     max.reference = 50000L ) {

  ref.idx <- seq_len( nrow( unstained ) )
  if ( length( ref.idx ) > max.reference )
    ref.idx <- sample( ref.idx, max.reference )

  knn.idx <- FNN::knnx.index(
    data  = as.matrix( unstained.scatter[ ref.idx, , drop = FALSE ] ),
    query = as.matrix( scatter ),
    k     = k.neighbors
  )

  bg <- matrix( 0, nrow( raw.data ), ncol( raw.data ),
                dimnames = dimnames( raw.data ) )

  for ( ki in seq_len( k.neighbors ) )
    bg <- bg +
      as.matrix( unstained[ ref.idx, , drop = FALSE ] )[
        knn.idx[ , ki ], , drop = FALSE ]

  raw.data - bg / k.neighbors
}
