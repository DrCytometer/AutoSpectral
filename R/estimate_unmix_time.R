# estimate_unmix_time.r

#' @title Estimate Unmixing Time
#'
#' @description
#' Benchmarks a small, representative sample of events from an FCS file
#' through the exact unmixing dispatch `unmix.fcs()` would use, then scales
#' the measured per-event rate up to the full file to give a fast,
#' approximate estimate of total unmixing time. Most useful for the
#' `AutoSpectral` method with `spectra.variants` supplied, where per-cell
#' variant selection makes unmixing the slowest step and its cost is hard
#' to predict analytically: it depends on the number of events, detectors,
#' fluorophores, the size of `spectra.variants`, and how many cells in the
#' file are above positivity thresholds and therefore trigger variant
#' search. Because the probe runs the real dispatch code on real data, all
#' of these factors are captured automatically without being modelled
#' explicitly.
#'
#' @details
#' For the range of event counts typical of flow cytometry files (tens of
#' thousands to tens of millions), per-event unmixing cost dominates any
#' fixed, event-count-independent overhead (thread pool setup, small object
#' allocation), so total dispatch time is well approximated as directly
#' proportional to the number of events. A per-event rate is measured once
#' from a small probe of `sample.events` drawn from the middle of the file
#' (not the first events, to avoid acquisition start-up artefacts) and then
#' used to extrapolate to the full event count. This is expected to be
#' accurate to within roughly 4-5x under normal conditions. It will be less
#' accurate if the fraction of multi-positive events, which trigger more
#' per-cell variant search, varies substantially over the course of
#' acquisition, since the probe's local complexity may then not represent
#' the file as a whole.
#'
#' File reading is timed and extrapolated the same way; file writing
#' (`writeFCS()`) is not estimated and is excluded from the total, since it
#' is typically small compared to `AutoSpectral` per-cell unmixing when
#' `spectra.variants` is supplied.
#'
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectra A matrix containing the spectral data, as passed to
#'   `unmix.fcs()`.
#' @param asp The AutoSpectral parameter list, as passed to `unmix.fcs()`.
#' @param method A character string specifying the unmixing method, as in
#'   `unmix.fcs()`. Default `"AutoSpectral"`.
#' @param af.spectra,spectra.variants,use.dist0,speed AutoSpectral-specific
#'   arguments, passed through unchanged to the same unmixing dispatch
#'   `unmix.fcs()` uses; see `?unmix.fcs` for details on each.
#' @param parallel,threads,n.variants,pipeline Execution and pipeline
#'   selection arguments, passed through unchanged; see `?unmix.fcs`.
#' @param n.passes,n.af.passes,cell.weight,noise.floor Joint-pipeline
#'   tuning arguments, passed through unchanged; see `?unmix.fcs`.
#' @param alpha,collinear.threshold,joint.pair.resolution,refine.af.quantile
#'   Further joint-pipeline tuning arguments, passed through unchanged;
#'   see `?unmix.fcs`.
#' @param weights,divergence.threshold,divergence.handling,balance.weight
#'   Weighting and IRLS-divergence arguments (`WLS`/`Poisson`/
#'   `FastPoisson`), passed through unchanged; see `?unmix.fcs`.
#' @param sample.events Numeric, number of events to draw for the timing
#'   probe. Default `5000`. Larger samples give a more stable estimate at
#'   the cost of a longer probe run.
#' @param chunk.size Numeric, as passed to `unmix.fcs()`; used only to
#'   report how many chunks the full run will use. Default `2e6`.
#' @param verbose Logical, whether to print the estimate. Default `TRUE`.
#'
#' @return Invisibly, a list with `total.events`, `sample.events`,
#'   `events.per.second`, `estimated.unmix.seconds`,
#'   `estimated.read.seconds`, `estimated.total.seconds`, and `chunk.n`;
#'   or `NULL` (invisibly) if `method` is not `"AutoSpectral"` or
#'   `spectra.variants` is `NULL`, in which case the timing probe is
#'   skipped as unnecessary and, if `verbose`, a message explains why.
#'
#' @export

estimate.unmix.time <- function(
    fcs.file,
    spectra,
    asp,
    method                 = c( "AutoSpectral", "OLS", "WLS", "Poisson", "FastPoisson" ),
    af.spectra             = NULL,
    spectra.variants       = NULL,
    use.dist0              = TRUE,
    speed                  = c( "fast", "medium", "slow" ),
    parallel               = TRUE,
    threads                = if ( parallel ) 0 else 1,
    n.variants             = NULL,
    pipeline               = c( "joint", "legacy" ),
    n.passes               = 1L,
    n.af.passes            = 1L,
    cell.weight            = if ( asp$cytometer == "ID7000" ) TRUE else FALSE,
    noise.floor            = 125,
    alpha                  = 0.5,
    collinear.threshold    = 0.5,
    joint.pair.resolution  = TRUE,
    refine.af.quantile     = 0.5,
    weights                = NULL,
    divergence.threshold   = 1e4,
    divergence.handling    = "Balance",
    balance.weight         = 0.5,
    sample.events          = 5000,
    chunk.size             = 2e6,
    verbose                = TRUE
) {

  method       <- match.arg( method )
  pipeline.arg <- match.arg( pipeline )

  if ( !.unmix.method.is.slow( method, spectra.variants ) )
    return( invisible( NULL ) )

  if ( method == "AutoSpectral" && is.null( n.variants ) )
    n.variants <- switch( match.arg( speed ), "fast" = 1L, "medium" = 3L, "slow" = 10L )

  if ( is.null( threads ) ) threads <- asp$worker.process.n
  if ( parallel && threads == 0 ) threads <- parallelly::availableCores()

  # header only, to get total event count cheaply
  import.meta  <- readFCS( fcs.file, return.keywords = TRUE, start.row = 1, end.row = 1 )
  total.events <- as.numeric( import.meta$keywords[[ "$TOT" ]] )

  spectral.channel <- colnames( spectra )

  sample.n <- min( sample.events, total.events )

  # draw the probe from the middle of the file, not the start, to reduce
  # the chance of catching acquisition start-up debris unrepresentative of
  # the file as a whole
  mid.start <- max( 1, floor( ( total.events - sample.n ) / 2 ) + 1 )
  mid.end   <- mid.start + sample.n - 1

  read.start <- Sys.time()
  probe.data <- readFCS(
    fcs.file,
    return.keywords = FALSE,
    start.row = mid.start,
    end.row = mid.end,
    columns = spectral.channel
  )
  read.elapsed <- as.numeric( Sys.time() - read.start, units = "secs" )

  probe.spectral <- probe.data[ , spectral.channel, drop = FALSE ]

  probe <- .probe.unmix.rate(
    chunk.spectral         = probe.spectral,
    spectra                = spectra,
    method                 = method,
    sample.n               = sample.n,
    af.spectra             = af.spectra,
    spectra.variants       = spectra.variants,
    use.dist0              = use.dist0,
    speed                  = speed,
    parallel               = parallel,
    threads                = threads,
    n.variants             = n.variants,
    pipeline.arg           = pipeline.arg,
    n.af.passes            = n.af.passes,
    n.passes               = n.passes,
    cell.weight            = cell.weight,
    noise.floor            = noise.floor,
    alpha                  = alpha,
    collinear.threshold    = collinear.threshold,
    joint.pair.resolution  = joint.pair.resolution,
    refine.af.quantile     = refine.af.quantile,
    asp                    = asp,
    weights                = weights,
    divergence.threshold   = divergence.threshold,
    divergence.handling    = divergence.handling,
    balance.weight         = balance.weight
  )

  events.per.second <- probe$events.per.second
  read.per.second   <- sample.n / read.elapsed

  estimated.unmix.seconds <- total.events / events.per.second
  estimated.read.seconds  <- total.events / read.per.second
  estimated.total.seconds <- estimated.unmix.seconds + estimated.read.seconds

  chunk.n <- ceiling( total.events / chunk.size )

  result <- list(
    total.events            = total.events,
    sample.events           = sample.n,
    events.per.second       = events.per.second,
    estimated.unmix.seconds = estimated.unmix.seconds,
    estimated.read.seconds  = estimated.read.seconds,
    estimated.total.seconds = estimated.total.seconds,
    chunk.n                 = chunk.n
  )

  if ( verbose ) {
    message( sprintf(
      "Probe: %d events in %.2fs unmix (%.0f events/sec), %.2fs read.",
      sample.n, probe$elapsed.seconds, events.per.second, read.elapsed
    ) )
    message( sprintf(
      "Estimated total: %s for %d events across %d chunk(s) (unmix ~%s, read ~%s). Treat as accurate to within roughly 4-5x.",
      .format.duration( estimated.total.seconds ),
      total.events,
      chunk.n,
      .format.duration( estimated.unmix.seconds ),
      .format.duration( estimated.read.seconds )
    ) )
  }

  invisible( result )
}


# ---------------------------------------------------------------------------
# Internal: format a duration in seconds as a human-readable string
# ---------------------------------------------------------------------------

.format.duration <- function( seconds ) {
  if ( !is.finite( seconds ) ) return( "unknown" )
  if ( seconds < 60 ) return( sprintf( "%.0fs", seconds ) )
  if ( seconds < 3600 ) return( sprintf( "%.1fmin", seconds / 60 ) )
  sprintf( "%.1fhr", seconds / 3600 )
}


# ---------------------------------------------------------------------------
# Internal: does this method/argument combination unmix slowly enough that
# a timing probe is worth its own cost? Shared by estimate.unmix.time() and
# the automatic check inside unmix.fcs().
# ---------------------------------------------------------------------------

.unmix.method.is.slow <- function( method, spectra.variants ) {
  method %in% c( "Poisson", "FastPoisson" ) ||
    ( method == "AutoSpectral" && !is.null( spectra.variants ) )
}


# ---------------------------------------------------------------------------
# Internal: time a dispatch call on an already-loaded block of spectral
# data and return the measured per-event rate. Shared by
# estimate.unmix.time() (which reads its own probe block from disk) and
# unmix.fcs() (which reuses its first chunk's already-loaded data).
# ---------------------------------------------------------------------------

.probe.unmix.rate <- function(
    chunk.spectral,
    spectra,
    method,
    sample.n                = min( 5000, nrow( chunk.spectral ) ),
    af.spectra               = NULL,
    spectra.variants         = NULL,
    use.dist0                = TRUE,
    speed                    = c( "fast", "medium", "slow" ),
    parallel                 = TRUE,
    threads                  = 1,
    n.variants                = NULL,
    pipeline.arg              = c( "joint", "legacy" ),
    n.af.passes                = 1L,
    n.passes                   = 1L,
    cell.weight                = FALSE,
    noise.floor                = 125,
    alpha                       = 0.5,
    collinear.threshold          = 0.5,
    joint.pair.resolution         = TRUE,
    refine.af.quantile             = 0.5,
    asp,
    weights                         = NULL,
    divergence.threshold             = 1e4,
    divergence.handling               = "Balance",
    balance.weight                     = 0.5
) {
  n.avail  <- nrow( chunk.spectral )
  sample.n <- min( sample.n, n.avail )

  # spread the probe indices across the available block rather than just
  # taking the first rows, to reduce sensitivity to any local run of
  # atypical events (e.g. debris at the very start of acquisition)
  probe.idx <- if ( sample.n >= n.avail ) seq_len( n.avail ) else
    round( seq( 1, n.avail, length.out = sample.n ) )

  probe.spectral <- chunk.spectral[ probe.idx, , drop = FALSE ]

  if ( ( method %in% c( "WLS", "Poisson", "FastPoisson" ) ) && is.null( weights ) )
    weights <- 1 / pmax( abs( colMeans( probe.spectral ) ), noise.floor )

  probe.start <- Sys.time()
  invisible( .dispatch.unmix.chunk(
    chunk.spectral         = probe.spectral,
    spectra                = spectra,
    method                 = method,
    af.spectra             = af.spectra,
    spectra.variants       = spectra.variants,
    use.dist0              = use.dist0,
    speed                  = speed,
    parallel               = parallel,
    threads                = threads,
    n.variants             = n.variants,
    pipeline.arg           = pipeline.arg,
    n.af.passes            = n.af.passes,
    n.passes               = n.passes,
    cell.weight            = cell.weight,
    noise.floor            = noise.floor,
    alpha                  = alpha,
    collinear.threshold    = collinear.threshold,
    joint.pair.resolution  = joint.pair.resolution,
    refine.af.quantile     = refine.af.quantile,
    asp                    = asp,
    weights                = weights,
    divergence.threshold   = divergence.threshold,
    divergence.handling    = divergence.handling,
    balance.weight         = balance.weight,
    verbose                = FALSE
  ) )
  probe.elapsed <- as.numeric( Sys.time() - probe.start, units = "secs" )

  list(
    sample.n          = sample.n,
    elapsed.seconds    = probe.elapsed,
    events.per.second = sample.n / probe.elapsed
  )
}


