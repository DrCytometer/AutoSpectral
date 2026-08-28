# estimate_noise_model.R

#' @title Estimate Detector Noise Model
#'
#' @description
#' Estimates the two parameters of the detector noise model used by
#' `unmix.gls()`: a per-detector additive variance floor (read noise plus
#' any residual dark-offset contribution) and a per-detector photon
#' conversion factor `counts.per.unit`, such that
#' \deqn{Var(y_d) = read.var_d + \mu_d / \kappa_d}
#'
#' Estimation uses the mean-variance relationship of the *residual* after
#' least-squares projection onto `spectra` (optionally augmented with
#' `af.spectra`). Residuals are orthogonal to the fitted values by
#' construction, so binning on the fitted value does not bias the variance
#' estimate. A degrees-of-freedom correction accounts for the variance
#' removed by the projection.
#'
#' Best estimated on bead controls (no autofluorescence, no biological
#' spread) or on an unstained cell control together with its AF dictionary.
#' On fully stained samples, unmodelled spectral variation inflates both
#' parameters; that is not necessarily wrong (it reflects true predictive
#' uncertainty) but it is no longer a pure instrument measurement.
#'
#' @param raw.data Numeric matrix (events x detectors), or a path to an FCS
#'   file. Detector columns must include all columns of `spectra`.
#' @param spectra Numeric matrix (fluorophores x detectors).
#' @param af.spectra Optional AF dictionary (AF components x detectors), as
#'   returned by `get.af.spectra()`. Strongly recommended for cell controls.
#' @param n.bins Integer, number of quantile bins per detector. Default `40`.
#' @param min.bin.n Integer, minimum events per bin. Default `50`.
#' @param trim.quantile Numeric, upper quantile of each detector discarded
#'   before binning, to exclude saturation. Default `0.999`.
#' @param verbose Logical. Default `TRUE`.
#' @param af.pc.n Integer, number of principal components to use to describe
#' `af.spectra`. Default is `5`.
#' @param af.raw.data Optional numeric matrix, or a path to an FCS file, of
#'   individual autofluorescence-only events -- typically the same unstained
#'   sample `af.spectra` was built from. When supplied, the AF shape basis is
#'   fit by PCA on these per-cell events (L-infinity normalised per event,
#'   matching `get.af.spectra()`'s own convention) instead of on the rows of
#'   `af.spectra`. `af.spectra` rows are SOM node centroids, already
#'   compressed from the underlying cells; if within-node AF spread is
#'   substantially real rather than noise (check `shape.frac` in
#'   `attr(af.spectra, "af.model")`), the centroids alone understate true
#'   per-cell AF shape diversity, and the basis fit here from raw events will
#'   need more components to reach the same variance-explained target.
#'   Ignored when `af.spectra` is `NULL`. Default `NULL`, which fits the
#'   basis on `af.spectra` as before.
#' @param af.basis.n.cells Integer, default `20000`. Maximum events from
#'   `af.raw.data` used for the PCA; downsampled if more are supplied. Only
#'   used when `af.raw.data` is supplied.
#' @param read.var.floor Optional named numeric vector, per-detector additive
#' @param file.id Optional vector (character or factor), length
#'   `nrow(raw.data)`, identifying which source file/acquisition each row
#'   came from. When supplied, per-file residual means are subtracted before
#'   the mean-variance regression, so between-file offsets (gain drift,
#'   voltage differences across acquisition days) are not counted as
#'   photon-driven variance. Recommended whenever `raw.data` was pooled
#'   across multiple control files, e.g. from `pool.scc.for.noise.model()`.
#'   Default `NULL`.
#' @param unstained.data Optional numeric matrix (events x detectors) of an
#'   unstained/AF-only control, used only to anchor the lower tranche of the
#'   linearity check at `dark.quantile` of the true-dark distribution. When
#'   `NULL`, falls back to `dark.quantile` of `raw.data` itself, which is a
#'   weaker proxy since it is not a genuinely dark reference. Default `NULL`.
#' @param dark.quantile Numeric in (0, 1), the quantile of `unstained.data`
#'   (or of `raw.data` when `unstained.data` is `NULL`) used as the upper
#'   bound of the lowest linearity tranche. Default `0.95`.
#' @param n.tranche Integer, number of intensity tranches for the linearity
#'   check. Default `10L`.
#'
#' @return A named list:
#' \describe{
#'   \item{`read.var`}{Named numeric vector, additive variance floor per
#'     detector, in squared instrument units.}
#'   \item{`counts.per.unit`}{Named numeric vector, \eqn{\kappa_d}.}
#'   \item{`kappa.pooled`}{Numeric scalar, median of `counts.per.unit`.}
#'   \item{`fit.table`}{Data frame of per-bin means and variances, for QC
#'     plotting.}
#'   \item{`r.squared`}{Named numeric vector, per-detector regression fit.}
#'   \item{`curvature.coef`}{Named numeric vector, the quadratic coefficient
#'     of `variance ~ mean + mean^2` fitted across `n.tranche` intensity
#'     tranches per detector. Negative indicates variance growing more
#'     slowly than linear at high signal (compression, e.g. approaching
#'     saturation); positive indicates faster-than-linear growth.}
#'   \item{`curvature.p`}{Named numeric vector, the p-value of
#'     `curvature.coef`. Values below roughly `0.01` indicate the linear
#'     mean-variance model is a poor fit for that detector's full range.}
#'   \item{`curvature.table`}{Data frame of per-tranche means, variances,
#'     and counts, one block per detector, for QC plotting -- the same
#'     shape as `fit.table` but at tranche rather than bin resolution.}
#'   \item{`noise.floor.signal`}{`sqrt(read.var)`, in *signal* units, for
#'     passing to the existing `noise.floor` argument of the C++ pipeline.
#'     Note the unit difference: `read.var` is a variance, the C++
#'     `noise.floor` is a signal level.}
#'   \item{`read.var.source`}{Named character vector, `"external"` or
#'     `"regression"` per detector, indicating whether `read.var` came from
#'     `read.var.floor` or the fitted intercept.}
#'   \item{`kappa.source`}{Named character vector, `"regression"` or
#'     `"pooled.median"` per detector, flagging detectors where the fitted
#'     slope was non-positive and `counts.per.unit` was filled from the
#'     pooled median instead of measured directly. Treat these detectors'
#'     `counts.per.unit`, and any ratio computed from it, with caution.}
#' }
#'
#' @export

estimate.noise.model <- function(
    raw.data,
    spectra,
    af.spectra    = NULL,
    n.bins        = 40L,
    min.bin.n     = 50L,
    trim.quantile = 0.999,
    verbose       = TRUE,
    af.pc.n       = 5L,
    af.raw.data      = NULL,
    af.basis.n.cells = 20000L,
    read.var.floor = NULL,
    file.id       = NULL,
    unstained.data = NULL,
    dark.quantile  = 0.95,
    n.tranche      = 10L

) {

  if ( is.character( raw.data ) ) {
    if ( verbose ) message( "Reading FCS file: ", raw.data )
    probe.cols <- colnames( readFCS( raw.data, start.row = 1, end.row = 1 ) )
    raw.data   <- readFCS( raw.data, columns = intersect( colnames( spectra ), probe.cols ) )
  }
  raw.data <- as.matrix( raw.data )

  det.names <- colnames( spectra )
  if ( !all( det.names %in% colnames( raw.data ) ) )
    stop( "`raw.data` does not contain all detectors present in `spectra`.",
          call. = FALSE )
  raw.data <- raw.data[ , det.names, drop = FALSE ]

  if ( !is.null( unstained.data ) ) {
    unstained.data <- as.matrix( unstained.data )
    if ( !all( det.names %in% colnames( unstained.data ) ) )
      stop( "`unstained.data` does not contain all detectors present in `spectra`.",
            call. = FALSE )
    unstained.data <- unstained.data[ , det.names, drop = FALSE ]
  }

  # ---------------------------------------------------------------------------
  # Projection basis
  # ---------------------------------------------------------------------------

  # The AF dictionary can hold 50-100+ node spectra which, with the
  # fluorophores, meets or exceeds the detector count and leaves no residual
  # to measure. AF lives in a low-dimensional subspace, so a handful of
  # principal components spans nearly all of it.

  basis <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  if ( !is.null( af.spectra ) ) {

    af.mat <- as.matrix( af.spectra )[ , det.names, drop = FALSE ]

    if ( !is.null( af.raw.data ) ) {

      # Fit the shape basis from raw per-cell AF events rather than the
      # dictionary centroids. af.spectra rows are already L-infinity
      # normalised by get.af.spectra(); raw events are not, so
      # normalisation happens here, per event.

      if ( is.character( af.raw.data ) ) {
        if ( verbose ) message( "Reading FCS file for AF basis: ", af.raw.data )
        af.raw.data <- readFCS( af.raw.data, columns = det.names )
      }
      af.events <- as.matrix( af.raw.data )[ , det.names, drop = FALSE ]

      if ( nrow( af.events ) > af.basis.n.cells )
        af.events <- af.events[ sample( nrow( af.events ), af.basis.n.cells ), , drop = FALSE ]

      af.peak <- apply( abs( af.events ), 1, max )
      af.norm <- af.events[ af.peak > 0, , drop = FALSE ] / af.peak[ af.peak > 0 ]

      af.svd   <- svd( af.norm, nu = 0L,
                       nv = min( af.pc.n, ncol( af.norm ) ) )
      af.basis <- t( af.svd$v )
      rownames( af.basis ) <- paste0( "AFPC", seq_len( nrow( af.basis ) ) )
      colnames( af.basis ) <- det.names

      if ( verbose )
        message( sprintf(
          "AF shape basis from %d raw events: %d PCs, %.1f%% of variance",
          nrow( af.norm ), nrow( af.basis ),
          100 * sum( af.svd$d[ seq_len( nrow( af.basis ) ) ]^2 ) /
            sum( af.svd$d^2 ) ) )

    } else if ( nrow( af.mat ) <= af.pc.n ) {
      af.basis <- af.mat
    } else {
      af.svd   <- svd( af.mat, nu = 0L,
                       nv = min( af.pc.n, ncol( af.mat ) ) )
      af.basis <- t( af.svd$v )
      rownames( af.basis ) <- paste0( "AFPC", seq_len( nrow( af.basis ) ) )
      colnames( af.basis ) <- det.names

      if ( verbose )
        message( sprintf(
          "AF dictionary (%d spectra) reduced to %d PCs, %.1f%% of variance",
          nrow( af.mat ), nrow( af.basis ),
          100 * sum( af.svd$d[ seq_len( nrow( af.basis ) ) ]^2 ) /
            sum( af.svd$d^2 ) ) )
    }

    basis <- rbind( basis, af.basis )
  }

  det.n   <- ncol( basis )
  basis.n <- nrow( basis )

  if ( basis.n >= det.n - 4L )
    stop( "Projection basis (", basis.n, " rows) leaves too few residual ",
          "degrees of freedom for ", det.n, " detectors. Reduce `af.pc.n`.",
          call. = FALSE )

  x.hat <- unmix.ols.fast( raw.data, basis )
  y.hat <- x.hat %*% basis
  resid <- raw.data - y.hat

  # Between-file offsets (gain drift, voltage differences across acquisition
  # days) show up as a mean shift in the residual that has nothing to do
  # with photon statistics. Left in, it inflates whichever of read.var/kappa
  # absorbs a shift correlated with brightness. Centering residuals within
  # each file removes it while leaving the within-file photon-driven spread
  # -- the thing kappa actually measures -- untouched.
  if ( !is.null( file.id ) ) {
    if ( length( file.id ) != nrow( raw.data ) )
      stop( "`file.id` must have one entry per row of `raw.data`.", call. = FALSE )
    file.id <- as.factor( file.id )
    resid   <- apply( resid, 2, function( col ) col - stats::ave( col, file.id ) )
    colnames( resid ) <- det.names
  }

  # residual variance is deflated by the projection: E[r_d^2] = sigma_d^2 (1 - h_d),
  # where h_d is the d-th diagonal of the detector-space hat matrix
  # H = t(basis) (basis t(basis))^{-1} basis. Leverage is far from uniform:
  # detectors under emission peaks lose much more residual variance than empty
  # detectors, so a single average correction biases per-detector kappa low at
  # exactly the peak detectors. Correct per detector instead.
  hat.diag <- colSums( ( MASS::ginv( basis %*% t( basis ) ) %*% basis ) * basis )
  hat.diag <- pmin( pmax( hat.diag, 0 ), 1 - 1e-3 )
  dof.correction <- 1 / ( 1 - hat.diag )

  # ---------------------------------------------------------------------------
  # Per-detector mean-variance regression
  # ---------------------------------------------------------------------------

  read.var       <- stats::setNames( rep( NA_real_, det.n ), det.names )
  kappa          <- stats::setNames( rep( NA_real_, det.n ), det.names )
  r.squared      <- stats::setNames( rep( NA_real_, det.n ), det.names )
  extrap         <- stats::setNames( rep( NA_real_, det.n ), det.names )
  curvature.coef <- stats::setNames( rep( NA_real_, det.n ), det.names )
  curvature.p    <- stats::setNames( rep( NA_real_, det.n ), det.names )
  fit.rows       <- vector( "list", det.n )
  tranche.rows   <- vector( "list", det.n )

  for ( d in seq_len( det.n ) ) {

    d.name <- det.names[ d ]
    fit.d <- y.hat[ , d ]
    res.d <- resid[ , d ]

    keep <- fit.d <= stats::quantile( fit.d, trim.quantile, na.rm = TRUE ) &
      is.finite( fit.d ) & is.finite( res.d )
    fit.d <- fit.d[ keep ]
    res.d <- res.d[ keep ]

    if ( length( fit.d ) < min.bin.n * 4L ) next

    breaks <- unique( stats::quantile(
      fit.d, probs = seq( 0, 1, length.out = n.bins + 1L ), na.rm = TRUE ) )
    if ( length( breaks ) < 4L ) next

    bin <- cut( fit.d, breaks = breaks, include.lowest = TRUE, labels = FALSE )

    bin.n    <- tabulate( bin, nbins = length( breaks ) - 1L )
    good.bin <- which( bin.n >= min.bin.n )
    if ( length( good.bin ) < 4L ) next

    bin.mean <- vapply( good.bin, function( b )
      mean( fit.d[ bin == b ] ), numeric( 1 ) )
    bin.var  <- vapply( good.bin, function( b )
      stats::mad( res.d[ bin == b ] )^2 * dof.correction[ d ], numeric( 1 ) )

    ok <- is.finite( bin.mean ) & is.finite( bin.var ) & bin.var > 0
    if ( sum( ok ) < 4L ) next

    fit.rows[[ d ]] <- data.frame(
      detector = det.names[ d ],
      mean     = bin.mean[ ok ],
      variance = bin.var[ ok ],
      n        = bin.n[ good.bin ][ ok ],
      stringsAsFactors = FALSE
    )

    bin.df <- fit.rows[[ d ]]

    has.floor <- !is.null( read.var.floor ) &&
      d.name %in% names( read.var.floor ) &&
      is.finite( read.var.floor[ d.name ] ) &&
      read.var.floor[ d.name ] > 0

    if ( has.floor ) {

      # Floor supplied externally (e.g. spectral.variants$noise.floor, from
      # SCC negative-event quantiles). The intercept is no longer free: only
      # the slope (1 / kappa) is fit, through the known floor. This avoids
      # the non-identifiability of a freely-fit intercept on cell data, where
      # unmodelled AF shape variance inflates the low-signal bins.

      floor.d <- read.var.floor[ d.name ]
      bin.df$resid.var <- pmax( bin.df$variance - floor.d, 0 )

      w.bin  <- bin.df$n / pmax( bin.df$variance, 1e-8 )^2
      lm.fit <- stats::lm( resid.var ~ mean - 1, data = bin.df, weights = w.bin )

      for ( pass in seq_len( 2L ) ) {
        v.fit  <- pmax( stats::fitted( lm.fit ) + floor.d, 1e-8 )
        w.bin  <- bin.df$n / v.fit^2
        lm.fit <- stats::lm( resid.var ~ mean - 1, data = bin.df, weights = w.bin )
      }

      slope <- stats::coef( lm.fit )[ 1L ]

      tss <- sum( ( bin.df$resid.var -
                      stats::weighted.mean( bin.df$resid.var, w.bin ) )^2 * w.bin )
      rss <- sum( stats::residuals( lm.fit )^2 * w.bin )
      r.squared[ d ] <- if ( tss > 0 ) 1 - rss / tss else NA_real_
      extrap[ d ]    <- NA_real_

      if ( is.finite( slope ) && slope > 0 ) kappa[ d ] <- 1 / slope
      read.var[ d ] <- floor.d

    } else {

      # Var( s^2 ) ~ 2 sigma^4 / n, so weight by n / sigma^4. Weighting by n
      # alone lets the brightest bins dictate the intercept, which is a long
      # extrapolation from them. Two IRLS passes using the fitted variance.
      w.bin  <- bin.df$n / pmax( bin.df$variance, 1e-8 )^2
      lm.fit <- stats::lm( variance ~ mean, data = bin.df, weights = w.bin )

      for ( pass in seq_len( 2L ) ) {
        v.fit  <- pmax( stats::fitted( lm.fit ), 1e-8 )
        w.bin  <- bin.df$n / v.fit^2
        lm.fit <- stats::lm( variance ~ mean, data = bin.df, weights = w.bin )
      }

      cf <- stats::coef( lm.fit )
      r.squared[ d ] <- summary( lm.fit )$r.squared

      # how far the intercept is extrapolated: variance at the lowest bin
      # divided by the fitted intercept. Values >> 1 mean the floor is not
      # identifiable from this file.
      extrap[ d ] <- min( bin.df$variance ) / max( cf[ 1L ], 1e-8 )

      if ( is.finite( cf[ 2L ] ) && cf[ 2L ] > 0 ) kappa[ d ] <- 1 / cf[ 2L ]
      if ( is.finite( cf[ 1L ] ) ) read.var[ d ] <- max( cf[ 1L ], 1e-8 )

      # ---------------------------------------------------------------------
      # Linearity check: ten tranches, not two
      # ---------------------------------------------------------------------
      # The single-slope fit above assumes one kappa across the whole
      # detector range. Rather than compare two arbitrary halves, split the
      # range into n.tranche levels and fit a quadratic term across their
      # (mean, variance) points -- the t-test on that term is a direct lack-
      # of-fit test for the linear model above. Tranche 1 is anchored to
      # dark.quantile of the unstained population (when supplied) rather than
      # an internal quantile of whatever got pooled into fit.d, so it isolates
      # the regime a genuinely dark event sits in -- exactly where baseline
      # restoration would show up.

      dark.cut <- if ( !is.null( unstained.data ) )
        stats::quantile( unstained.data[ , det.names[ d ] ], dark.quantile, na.rm = TRUE ) else
          stats::quantile( fit.d, dark.quantile, na.rm = TRUE )

      above <- fit.d > dark.cut

      if ( sum( !above ) >= min.bin.n && sum( above ) >= min.bin.n * ( n.tranche - 1L ) ) {

        tr.breaks <- unique( c(
          -Inf, dark.cut,
          stats::quantile( fit.d[ above ],
                           probs = seq( 0, 1, length.out = n.tranche )[ -1L ],
                           na.rm = TRUE )
        ) )

        if ( length( tr.breaks ) >= 4L ) {

          tranche <- cut( fit.d, breaks = tr.breaks, include.lowest = TRUE, labels = FALSE )
          tr.n    <- tabulate( tranche, nbins = length( tr.breaks ) - 1L )
          tr.good <- which( tr.n >= min.bin.n )

          if ( length( tr.good ) >= 4L ) {

            tr.mean <- vapply( tr.good, function( b )
              mean( fit.d[ tranche == b ] ), numeric( 1 ) )
            tr.var  <- vapply( tr.good, function( b )
              stats::mad( res.d[ tranche == b ] )^2 * dof.correction[ d ], numeric( 1 ) )

            tr.ok <- is.finite( tr.mean ) & is.finite( tr.var ) & tr.var > 0

            if ( sum( tr.ok ) >= 4L ) {

              tr.df <- data.frame(
                detector = det.names[ d ],
                tranche  = tr.good[ tr.ok ],
                mean     = tr.mean[ tr.ok ],
                variance = tr.var[ tr.ok ],
                n        = tr.n[ tr.good ][ tr.ok ],
                stringsAsFactors = FALSE
              )

              tranche.rows[[ d ]] <- tr.df

              quad.fit <- tryCatch(
                stats::lm( variance ~ mean + I( mean^2 ), data = tr.df, weights = n ),
                error = function( e ) NULL )

              if ( !is.null( quad.fit ) ) {
                quad.cf <- summary( quad.fit )$coefficients
                if ( "I(mean^2)" %in% rownames( quad.cf ) ) {
                  curvature.coef[ d ] <- quad.cf[ "I(mean^2)", "Estimate" ]
                  curvature.p[ d ]    <- quad.cf[ "I(mean^2)", "Pr(>|t|)" ]
                }
              }
            }
          }
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Fill failures from the pooled estimate
  # ---------------------------------------------------------------------------

  read.var.source <- stats::setNames( rep( "regression", det.n ), det.names )
  if ( !is.null( read.var.floor ) ) {
    ext <- det.names %in% names( read.var.floor ) &
      is.finite( read.var.floor[ det.names ] ) &
      read.var.floor[ det.names ] > 0
    ext[ is.na( ext ) ] <- FALSE
    read.var.source[ ext ] <- "external"
  }

  kappa.source <- stats::setNames( rep( "regression", det.n ), det.names )

  bad.kappa <- !is.finite( kappa ) | kappa <= 0
  if ( any( bad.kappa ) ) {
    if ( all( bad.kappa ) )
      stop( "Mean-variance regression failed at every detector. Check that ",
            "`raw.data` spans a range of intensities and that `spectra` ",
            "matches the file.", call. = FALSE )
    warning( "Non-positive variance slope at ", sum( bad.kappa ), " detector(s): ",
             paste( det.names[ bad.kappa ], collapse = ", " ),
             "; filled from the pooled median.", call. = FALSE )
    kappa.source[ bad.kappa ] <- "pooled.median"
    kappa[ bad.kappa ] <- stats::median( kappa[ !bad.kappa ] )
  }

  bad.floor <- !is.finite( read.var ) | read.var <= 0
  if ( all( bad.floor ) )
    stop( "The variance floor could not be estimated at any detector. ",
          "Supply `read.var.floor` or a control spanning low intensities.",
          call. = FALSE )
  if ( any( bad.floor ) )
    read.var[ bad.floor ] <- stats::median( read.var[ !bad.floor ] )

  n.tested  <- sum( !is.na( curvature.p ) )
  n.curved  <- sum( curvature.p < 0.01, na.rm = TRUE )

  if ( verbose ) {
    message( sprintf(
      "Noise model: median read SD %.1f, median counts.per.unit %.3g, median R2 %.3f",
      stats::median( sqrt( read.var ) ),
      stats::median( kappa ),
      stats::median( r.squared, na.rm = TRUE )
    ) )
    message( sprintf(
      "  %d / %d detectors show significant curvature (quadratic term p < 0.01)",
      n.curved, n.tested ) )
  }

  if ( stats::median( extrap, na.rm = TRUE ) > 20 )
    warning( "The variance floor is being extrapolated from data far above ",
             "zero signal (median ratio ", round( stats::median( extrap, na.rm = TRUE ) ),
             "). `read.var` is poorly identified. Estimate on beads or on a ",
             "blank acquisition if the floor matters for your application.",
             call. = FALSE )

  if ( n.tested > 0 && n.curved / n.tested > 0.25 )
    warning( "Significant curvature (quadratic term p < 0.01) detected at ",
             n.curved, " of ", n.tested, " detectors (",
             round( 100 * n.curved / n.tested ), "%). A single kappa per ",
             "detector is an approximation there; treat unmix.gls() ",
             "weighting as approximate at the extremes of the intensity ",
             "range. See `curvature.table` for the tranche-level fits.",
             call. = FALSE )

  if ( any( read.var.source == "regression" ) &&
       stats::median( extrap[ read.var.source == "regression" ], na.rm = TRUE ) > 20 )
    warning( "The variance floor is being extrapolated from data far above ",
             "zero signal at detector(s) without an external floor (median ratio ",
             round( stats::median( extrap[ read.var.source == "regression" ], na.rm = TRUE ) ),
             "). `read.var` is poorly identified from this file. Pass ",
             "`read.var.floor` (e.g. `spectral.variants$noise.floor` from ",
             "single-stained control negatives) or an unstained bead ",
             "population instead; a blank/buffer acquisition is not possible ",
             "on this class of instrument.",
             call. = FALSE )

  list(
    read.var           = read.var,
    counts.per.unit    = kappa,
    kappa.pooled       = stats::median( kappa ),
    fit.table          = do.call( rbind, fit.rows ),
    r.squared          = r.squared,
    extrap             = extrap,
    curvature.coef     = curvature.coef,
    curvature.p        = curvature.p,
    curvature.table    = do.call( rbind, tranche.rows ),
    noise.floor.signal = sqrt( read.var )
  )
}
