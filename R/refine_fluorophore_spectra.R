# refine_fluorophore_spectra.R

#' @title Refine Fluorophore Spectra
#'
#' @description
#' Optional second pass over `get.fluorophore.spectra()`'s or
#' `get.spectra.automated()`'s output, re-measuring each row directly from its
#' own single-color control with `extract.raw.signature()` -- the same
#' raw-space engine `fix.my.unmix()` uses for its own phase two.
#'
#' Checking a spectrum against the same gated, downsampled population it was
#' fit from cannot surface bias that population selection itself introduced
#' -- exactly the effect documented for `fix.my.unmix()`. A fresh read is the
#' only way this check is independent of whatever selection produced the
#' first-pass spectra, whether that was `define.flow.control()`'s scatter
#' gate, `clean.controls()`'s downsampling, or `get.spectra.automated()`'s
#' cosine-similarity filter.
#'
#' That fresh read is done one control at a time, not by building a second,
#' whole-panel `flow.control` in memory: `.read.fcs.clean()` (the private FCS
#' reader `get.spectra.automated()` already uses) reads one file, drops
#' saturating and doublet events, and returns a plain matrix; nothing else is
#' retained. Where a control needs autofluorescence removed, its matching
#' negative is read the same way and the two are handed to `remove.af()` in a
#' throwaway two-entry list -- `remove.af()` only ever touches
#' `clean.expr[[samp]]` and `clean.expr[[matching.negative]]`, so it does not
#' need `clean.controls()`'s whole-panel structure to run correctly on a
#' single pair. Every control's data is discarded once that control's fit is
#' done. The trade is disk I/O: every one of `n.iter` passes re-reads every
#' control file, and its matching negative, from scratch. Nothing is cached
#' across controls or across iterations, so peak memory is bounded by roughly
#' one control file and one negative file, regardless of panel size.
#'
#' Every fluorophore's control is unmixed against its own row alone (`active
#' = target` in `extract.raw.signature()`), not the whole panel: a
#' single-color control has only one fluorophore truly present, so the
#' population every other row should read zero in is known, not inferred, and
#' no nuisance subtraction is needed. A companion diagnostic,
#' `estimate.residual.spillover()`, checks the same read's abundance in every
#' *other* fluorophore's column, with its negative mask set to "every event",
#' unconditionally, for the same reason, reusing the same read rather than a
#' second pass over the file. That check can only ever see the component of a
#' row's error that lies in the span of the current library -- unmixing
#' removes a spectral error by projecting it through the compensation
#' operator, so only the projection is visible to any abundance-space
#' residual -- so it is logged and optionally warned on, not used to update a
#' row; `extract.raw.signature()`, which regresses the raw detector trace
#' directly rather than inverting anything, is what recovers the full error
#' and is what actually updates `spectra`. Because it reuses this iteration's
#' read rather than the just-updated matrix, it reports cross-talk under the
#' spectra this iteration started from, not the spectra it ends with.
#'
#' Every row update passes the same acceptance stack `fix.my.unmix()` phase
#' two uses (`max.angle`, `min.explained`/`max.explained`, `max.resid`,
#' `max.intercept`, `min.bg.align`); a row that fails any gate keeps its
#' starting spectrum for that iteration. Every fluorophore is refit against
#' the same starting `spectra` within one iteration and all accepted updates
#' are applied together at the end of it, so results do not depend on the
#' order controls happen to be listed in.
#'
#' @param marker.spectra Numeric matrix (samples x detectors), the first-pass
#'   spectra. Rownames must be the control `sample` identifiers (not
#'   necessarily the dye identity -- a dye can have more than one control);
#'   `attr(marker.spectra, "fluorophore")`, if present, is used to keep
#'   replicate controls of the same dye out of each other's cross-talk check.
#' @param control.dir,control.def.file As in `define.flow.control()`. Should
#'   normally be the same files used to build the spectra `marker.spectra`
#'   was fit from.
#' @param asp The AutoSpectral parameter list.
#' @param af.remove Logical, default `TRUE`. Whether to read each cell-type
#'   control's matching negative and run `remove.af()`'s intrusive-AF gate on
#'   it. Bead controls are never AF-removed, matching `clean.controls()`.
#'   Set `FALSE` for a bead-only panel.
#' @param af.figures Logical, default `FALSE`. Whether `remove.af()` writes
#'   its own AF-removal diagnostic figures for each control.
#' @param singlet.quantiles,remove.doublets As in `get.spectra.automated()`,
#'   passed to `.read.fcs.clean()` for every file read here.
#' @param allow.duplicate.controls Logical, default `FALSE`. As in
#'   `define.flow.control()`/`get.spectra.automated()`; set `TRUE` if the
#'   control file the first pass used permits multiple controls per dye.
#' @param n.iter Integer, maximum refine iterations. Default `3L`.
#' @param intercept,multivariate,ridge,n.levels,min.bin.events As in
#'   `extract.raw.signature()`.
#' @param min.events Integer, minimum events for both the signature fit and
#'   the cross-talk pair estimator. Default `200L`.
#' @param max.angle Numeric, degrees. A candidate row is rejected if its
#'   cosine distance from the current row exceeds this. Default `5`.
#' @param min.explained,max.explained,max.resid,max.intercept,min.bg.align As
#'   in `fix.my.unmix()`'s phase-two acceptance gates.
#' @param unstained.threshold,unstained.margin Numeric, used to set the
#'   cross-talk check's source threshold from the control's matching
#'   negative, when one is defined. Same convention as `fix.my.unmix()`.
#' @param crosstalk.check Logical, whether to run the diagnostic
#'   `estimate.residual.spillover()` pass each iteration. Default `TRUE`.
#' @param max.crosstalk Numeric, the abundance-space coefficient above which
#'   a console warning is printed. Diagnostic only; does not gate a row.
#'   Default `0.1`.
#' @param n.levels.pair Integer, abundance bins for the cross-talk pair
#'   estimator. Default `10L`.
#' @param convergence.threshold Numeric, degrees. Iteration stops early once
#'   the largest accepted `deg.change` in a pass falls below this. Default
#'   `0.5`.
#' @param step Numeric in (0, 1], the fraction of each accepted change
#'   applied. Default `1`.
#' @param n.threads Integer, threads for the cross-talk pair estimator.
#'   Default `1L`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`spectra`}{The refined matrix, same shape and rownames as
#'     `marker.spectra`, with `attr(spectra, "fluorophore")` preserved.}
#'   \item{`log`}{Data frame, one row per sample per iteration, with the
#'     acceptance decision and fit diagnostics.}
#'   \item{`crosstalk`}{Data frame of `estimate.residual.spillover()` output
#'     for every sample and iteration, or `NULL` if `crosstalk.check =
#'     FALSE`.}
#' }
#'
#' @importFrom stats quantile setNames
#'
#' @export

refine.fluorophore.spectra <- function(
    marker.spectra,
    control.dir,
    control.def.file,
    asp,
    af.remove = TRUE,
    af.figures = FALSE,
    singlet.quantiles = c( 0.85, 0.975 ),
    remove.doublets = TRUE,
    allow.duplicate.controls = FALSE,
    n.iter = 3L,
    intercept = TRUE,
    multivariate = TRUE,
    ridge = 1e-6,
    n.levels = 60L,
    min.bin.events = 50L,
    min.events = 200L,
    max.angle = 5,
    min.explained = 0.8,
    max.explained = 1.2,
    max.resid = 0.05,
    max.intercept = 0.05,
    min.bg.align = -0.9,
    unstained.threshold = 0.99,
    unstained.margin = 1.3,
    crosstalk.check = TRUE,
    max.crosstalk = 0.1,
    n.levels.pair = 10L,
    convergence.threshold = 0.5,
    step = 1,
    n.threads = 1L,
    verbose = TRUE
) {

  marker.spectra <- as.matrix( marker.spectra )

  if ( is.null( rownames( marker.spectra ) ) )
    stop( "`marker.spectra` must have sample row names.", call. = FALSE )

  fluorophores <- rownames( marker.spectra )
  detectors    <- colnames( marker.spectra )

  if ( length( fluorophores ) < 2 )
    stop( "At least two rows are required.", call. = FALSE )

  fluor.attr <- attr( marker.spectra, "fluorophore" )
  if ( is.null( fluor.attr ) ) fluor.attr <- fluorophores
  names( fluor.attr ) <- fluorophores

  if ( verbose )
    message( paste0( "\033[34m", "Checking control file for refine step",
                     "\033[0m" ) )

  check.control.file(
    control.dir, control.def.file, asp, legacy = FALSE,
    allow.duplicate.controls = allow.duplicate.controls
  )

  ctrl.path <- if ( file.exists( control.def.file ) )
    control.def.file else file.path( control.dir, control.def.file )

  ctrl.tbl <- utils::read.csv(
    ctrl.path, stringsAsFactors = FALSE, strip.white = TRUE )

  ctrl.tbl$sample <- .build.control.sample.names(
    ctrl.tbl$fluorophore, ctrl.tbl$control.type, ctrl.tbl$marker )

  refine.rows <- which( ctrl.tbl$sample %in% fluorophores )

  if ( length( refine.rows ) < 2 )
    stop( paste0( "Fewer than two `marker.spectra` rows match samples in ",
                  "`control.def.file`." ), call. = FALSE )

  # resolve each refine row's negative source exactly as
  # get.spectra.automated() does, so a `universal.negative` column that is
  # blank (internal-negative), a filename, or `TRUE` (use the row(s) marked
  # TRUE) all resolve the same way in both pipelines
  uneg.bool.all <- suppressWarnings( as.logical( ctrl.tbl$universal.negative ) )

  unstained.sources <- lapply( refine.rows, function( i ) {
    src <- .parse.unstained.source( ctrl.tbl$universal.negative[ i ] )
    if ( src$type == "global.true" ) {
      true.rows <- which( !is.na( uneg.bool.all ) & uneg.bool.all )
      if ( length( true.rows ) > 0 )
        return( list( type = "file", file = ctrl.tbl$filename[ true.rows[ 1L ] ] ) )
      return( list( type = "internal", file = NULL ) )
    }
    src
  } )
  names( unstained.sources ) <- ctrl.tbl$sample[ refine.rows ]

  sat.value        <- if ( !is.null( asp$expr.data.max ) ) asp$expr.data.max else Inf
  scatter.channels <- read.scatter.parameter( asp )

  first.fcs         <- file.path( control.dir, ctrl.tbl$filename[ refine.rows[ 1L ] ] )
  spectral.channels <- .derive.spectral.channels( first.fcs, asp )

  if ( !identical( spectral.channels, detectors ) )
    stop( paste0( "Detectors read from `control.def.file` do not match ",
                  "`marker.spectra`; is this the same control file the ",
                  "first pass was built from?" ), call. = FALSE )

  if ( af.remove && af.figures ) {
    if ( !dir.exists( asp$figure.clean.control.dir ) )
      dir.create( asp$figure.clean.control.dir, recursive = TRUE )
    if ( !dir.exists( asp$figure.spectral.ribbon.dir ) )
      dir.create( asp$figure.spectral.ribbon.dir, recursive = TRUE )
    if ( !dir.exists( asp$figure.scatter.dir.base ) )
      dir.create( asp$figure.scatter.dir.base, recursive = TRUE )
  }

  read.one <- function( filename, label )
    .read.fcs.clean(
      file.path( control.dir, filename ), label,
      spectral.channels, scatter.channels, sat.value,
      singlet.quantiles, remove.doublets, asp, verbose )

  # unmix against the target's own row alone -- no nuisance subtraction
  project.one <- function( y, target.row )
    as.numeric( unmix.ols.fast( y[ , detectors, drop = FALSE ], target.row ) )

  # unmix against the whole current panel, for the cross-talk diagnostic only
  project.panel <- function( y, spectra ) {
    coefs <- unmix.ols.fast( y[ , detectors, drop = FALSE ], spectra )
    colnames( coefs ) <- rownames( spectra )
    coefs[ , fluorophores, drop = FALSE ]
  }

  log.rows       <- list()
  crosstalk.rows <- list()

  for ( iter in seq_len( n.iter ) ) {

    updated      <- marker.spectra
    accepted.deg <- numeric( 0 )

    for ( i in refine.rows ) {

      samp    <- ctrl.tbl$sample[ i ]
      is.cell <- identical( ctrl.tbl$control.type[ i ], "cells" )
      src     <- unstained.sources[[ samp ]]

      if ( verbose )
        message( sprintf( "\033[32mRefining %s (iteration %d)\033[0m", samp, iter ) )

      pos.mat <- tryCatch(
        read.one( ctrl.tbl$filename[ i ], samp ),
        error = function( e ) {
          if ( verbose )
            message( sprintf( "\033[31m%s: failed to read (%s)\033[0m",
                              samp, e$message ) )
          NULL
        } )

      if ( is.null( pos.mat ) || nrow( pos.mat ) < min.events ) {
        log.rows[[ length( log.rows ) + 1 ]] <- data.frame(
          iter = iter, sample = samp, accepted = FALSE,
          reject = "read.fail.or.min.events", deg.change = NA_real_,
          explained.total = NA_real_, resid.rel = NA_real_,
          stringsAsFactors = FALSE )
        next
      }

      neg.mat <- NULL
      if ( af.remove && is.cell && src$type == "file" ) {
        neg.mat <- tryCatch(
          read.one( src$file, paste0( samp, " (negative)" ) ),
          error = function( e ) {
            if ( verbose )
              message( sprintf( "\033[33m%s: negative failed to read (%s)\033[0m",
                                samp, e$message ) )
            NULL
          } )
      }

      y <- pos.mat

      if ( !is.null( neg.mat ) ) {

        clean.pair <- list()
        clean.pair[[ samp ]]  <- pos.mat
        clean.pair[[ "NEG" ]] <- neg.mat

        y <- tryCatch(
          remove.af(
            samp                 = samp,
            clean.expr           = clean.pair,
            spectral.channel     = spectral.channels,
            peak.channel         = stats::setNames( ctrl.tbl$channel[ i ], samp ),
            universal.negative   = stats::setNames( "NEG", samp ),
            asp                  = asp,
            scatter.param        = scatter.channels,
            scatter.match        = FALSE,
            main.figures         = af.figures,
            intermediate.figures = FALSE,
            verbose              = verbose ),
          error = function( e ) {
            if ( verbose )
              message( sprintf(
                "\033[33m%s: AF removal failed (%s), using ungated data\033[0m",
                samp, e$message ) )
            pos.mat
          } )
      }

      abundance <- project.one( y, marker.spectra[ samp, , drop = FALSE ] )

      candidate <- extract.raw.signature(
        raw.data       = y[ , detectors, drop = FALSE ],
        spectra        = marker.spectra,
        abundance      = matrix( abundance, ncol = 1,
                                 dimnames = list( NULL, samp ) ),
        target         = samp,
        active         = samp,
        intercept      = intercept,
        multivariate   = multivariate,
        ridge          = ridge,
        n.levels       = n.levels,
        min.bin.events = min.bin.events,
        min.events     = min.events )

      reject <- NA_character_
      if ( is.null( candidate ) ) reject <- "no.fit"

      if ( is.na( reject ) ) {
        st <- candidate$stats
        if ( st$deg.change > max.angle )
          reject <- "angle"
        else if ( !is.finite( st$explained.total ) ||
                  st$explained.total < min.explained ||
                  st$explained.total > max.explained )
          reject <- "explained"
        else if ( st$resid.rel > max.resid )
          reject <- "fit"
        else if ( st$intercept.rel > max.intercept )
          reject <- "offset"
        else if ( is.finite( st$bg.align ) && st$bg.align < min.bg.align )
          reject <- "bg.align"
      }

      accepted <- is.na( reject )

      if ( accepted ) {
        new.row <- step * candidate$signature +
          ( 1 - step ) * marker.spectra[ samp, ]
        new.row <- new.row / max( new.row )
        updated[ samp, ] <- new.row
        accepted.deg <- c( accepted.deg, candidate$stats$deg.change )
      }

      log.rows[[ length( log.rows ) + 1 ]] <- data.frame(
        iter = iter, sample = samp, accepted = accepted, reject = reject,
        deg.change = if ( is.null( candidate ) ) NA_real_ else candidate$stats$deg.change,
        explained.total = if ( is.null( candidate ) ) NA_real_ else candidate$stats$explained.total,
        resid.rel = if ( is.null( candidate ) ) NA_real_ else candidate$stats$resid.rel,
        stringsAsFactors = FALSE )

      if ( crosstalk.check ) {

        panel.abundance <- project.panel( y, marker.spectra )

        threshold.source <- if ( !is.null( neg.mat ) && nrow( neg.mat ) >= min.events ) {
          neg.abundance <- project.panel( neg.mat, marker.spectra )
          unstained.margin * stats::quantile(
            neg.abundance[ , samp ], unstained.threshold, names = FALSE )
        } else {
          stats::quantile( panel.abundance[ , samp ], 0.10, names = FALSE )
        }

        targets <- fluorophores[
          fluor.attr != fluor.attr[ samp ] & fluorophores != samp ]

        if ( length( targets ) > 0 ) {

          neg.mask <- matrix( TRUE, nrow = nrow( panel.abundance ), ncol = length( targets ) )

          est <- estimate.residual.spillover(
            unmixed          = panel.abundance,
            source            = samp,
            targets           = targets,
            negative.mask      = neg.mask,
            threshold.source   = threshold.source,
            n.levels          = n.levels.pair,
            min.events        = min.events,
            n.threads         = n.threads )

          if ( !is.null( est ) ) {

            est$iter   <- iter
            est$sample <- samp
            crosstalk.rows[[ length( crosstalk.rows ) + 1 ]] <- est

            bad <- which( is.finite( est$slope.truncated ) &
                           abs( est$slope.truncated ) > max.crosstalk )

            if ( length( bad ) > 0 && verbose )
              message( sprintf(
                "\033[33m%s: residual cross-talk into %s (max %.3f)\033[0m",
                samp, paste( est$target[ bad ], collapse = ", " ),
                max( abs( est$slope.truncated[ bad ] ) ) ) )
          }
        }
      }

      rm( pos.mat, neg.mat, y )
    }

    marker.spectra <- updated

    max.deg <- if ( length( accepted.deg ) > 0 ) max( accepted.deg ) else NA_real_

    if ( verbose )
      message( sprintf(
        "Refine iteration %d: %d/%d rows accepted, max angle change %s",
        iter, length( accepted.deg ), length( refine.rows ),
        if ( is.na( max.deg ) ) "NA" else sprintf( "%.2f deg", max.deg ) ) )

    if ( !is.na( max.deg ) && max.deg < convergence.threshold ) break
  }

  attr( marker.spectra, "fluorophore" ) <- fluor.attr

  list(
    spectra   = marker.spectra,
    log       = do.call( rbind, log.rows ),
    crosstalk = if ( length( crosstalk.rows ) > 0 )
      do.call( rbind, crosstalk.rows ) else NULL
  )
}
