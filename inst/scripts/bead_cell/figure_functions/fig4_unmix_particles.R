# fig4_unmix_particles.R
#
# Panel e of Figure 4: one single-stained cell control unmixed several times
# in the same full panel, each time with the stained fluorophore represented
# by the spectrum from a different particle type (cells, or one of the bead
# types). The same positive events are followed through every unmixing.
#
# Reuses AutoSpectral readFCS(), apply.gate(), unmix.ols(), calculate.ssi(),
# calculate.condition.number() and cosine.similarity().

#' Locate a saved landmark gate (newest file) or define it again
#'
#' define.gate.landmarks() saves each gate as
#' <dir>/landmark_gate_definition_<gate.name>_<timestamp>.rds. The newest file
#' is used; if none exists the gate is defined again from the control file.
#' @return list(gate = polygon list(x, y, ...), source = file path or "defined")
get.landmark.gate <- function( gate.dir, gate.name, control.file, control.dir, asp,
                               output.dir = tempdir() ) {
  pat <- paste0( "^landmark_gate_definition_", gate.name, "_[0-9_]+\\.rds$" )
  files <- if ( dir.exists( gate.dir ) ) list.files( gate.dir, pattern = pat, full.names = TRUE ) else character( 0 )
  if ( length( files ) > 0 ) {
    f <- sort( files, decreasing = TRUE )[ 1 ]
    message( "Using saved gate: ", f )
    return( list( gate = readRDS( f ), source = f ) )
  }
  message( "No saved '", gate.name, "' gate in ", gate.dir, "; defining it with define.gate.landmarks()" )
  g <- define.gate.landmarks( control.file = control.file, control.dir = control.dir, asp = asp,
                              gate.name = gate.name, output.dir = output.dir )
  list( gate = g, source = "defined" )
}

#' Unmix one stained control with spectra from several particle types
#'
#' @param stained.file,unstained.file FCS paths (stained cells, unstained cells).
#' @param panel Full-panel spectral matrix (fluorophores x detectors) built
#'   from the reference particle's spectra.
#' @param particle.spectra Named list (particle type -> named numeric vector or
#'   1-row matrix over the detectors of `panel`): the stained fluorophore's
#'   spectrum from each particle type. The first entry should be the reference.
#' @param fluorophore Stained fluorophore (row of `panel`).
#' @param asp AutoSpectral parameter list.
#' @param gate Scatter gate polygon (list with x, y).
#' @param reference Particle type that defines positivity (same events for
#'   every unmixing) and the displayed unstained population.
#' @param positive.quantile,positive.multiplier Stained events are kept if
#'   their value in `fluorophore`, unmixed with the reference matrix, exceeds
#'   positive.multiplier x the positive.quantile of the unmixed unstained cells.
#' @param add.af Add an autofluorescence row "AF" to every matrix: the
#'   per-detector median of the gated unstained cells, divided by its maximum
#'   (replaces an "AF" row already in `panel`).
#' @param y.exclude Unmixed channels never chosen automatically as the y axis
#'   (default "AF").
#' @param y.channel Unmixed channel for the y axis. NULL or NA (default) =
#'   automatic: the channel with the largest |SSI| over the non-reference
#'   (bead) unmixings. A channel given here is used for the SSI table, and the
#'   caller uses the returned y.channel for the spread boundaries and axes.
#' @param unmix.fun Unmixing function (default unmix.ols).
#' @return list(raw, unmixed (per particle: stained, unstained), keep,
#'   threshold, y.channel, y.channel.auto, ssi (per particle at y.channel), ssi.all (channel x
#'   bead particle), spectra.list, af.spectrum, condition.number,
#'   spectrum.comparison)
unmix.particle.spectra <- function(
    stained.file, unstained.file, panel, particle.spectra, fluorophore, asp, gate,
    reference = names( particle.spectra )[ 1 ],
    positive.quantile = 0.995, positive.multiplier = 2,
    add.af = TRUE, y.exclude = "AF", y.channel = NULL,
    unmix.fun = unmix.ols
) {
  spec.ch <- colnames( panel )
  scatter <- asp$default.scatter.parameter
  if ( !fluorophore %in% rownames( panel ) ) stop( fluorophore, " is not in the panel", call. = FALSE )
  if ( !reference %in% names( particle.spectra ) ) stop( "reference not in particle.spectra", call. = FALSE )

  # read, drop saturated events, gate
  read.gated <- function( f ) {
    d <- readFCS( f, columns = union( scatter, spec.ch ) )
    d <- d[ rowSums( d[ , spec.ch, drop = FALSE ] >= asp$expr.data.max ) == 0, , drop = FALSE ]
    apply.gate( d, gate, scatter.param = scatter, asp = asp )
  }
  raw <- list( stained = read.gated( stained.file ), unstained = read.gated( unstained.file ) )
  message( sprintf( "Gated events: stained %d, unstained %d", nrow( raw$stained ), nrow( raw$unstained ) ) )

  # autofluorescence: median of the gated unstained cells, peak = 1
  af.spectrum <- NULL
  if ( add.af ) {
    af.spectrum <- apply( raw$unstained[ , spec.ch, drop = FALSE ], 2, stats::median )
    af.spectrum <- af.spectrum / max( af.spectrum )
    panel <- rbind( panel[ setdiff( rownames( panel ), "AF" ), , drop = FALSE ],
                    AF = af.spectrum[ spec.ch ] )
  }

  # one full-panel matrix per particle type
  spectra.list <- lapply( particle.spectra, function( s ) {
    s <- if ( is.matrix( s ) ) s[ 1, spec.ch ] else s[ spec.ch ]
    if ( any( is.na( s ) ) ) stop( "Particle spectrum lacks some panel detectors", call. = FALSE )
    m <- panel
    m[ fluorophore, ] <- as.numeric( s )
    m
  } )

  unmixed <- lapply( spectra.list, function( m ) list(
    stained   = unmix.fun( raw$stained[ , spec.ch, drop = FALSE ], m ),
    unstained = unmix.fun( raw$unstained[ , spec.ch, drop = FALSE ], m )
  ) )

  # positive events, from the reference unmixing (same events everywhere)
  thr <- unname( stats::quantile( unmixed[[ reference ]]$unstained[ , fluorophore ], positive.quantile ) ) *
    positive.multiplier
  keep <- which( unmixed[[ reference ]]$stained[ , fluorophore ] > thr )
  message( sprintf( "Positivity threshold (%s unmixing): %.1f; %d positive events", reference, thr, length( keep ) ) )
  if ( length( keep ) < 10 ) warning( "Fewer than 10 positive events", call. = FALSE )

  ssi.of <- function( p ) suppressWarnings( calculate.ssi(
    unmixed[[ p ]]$unstained, unmixed[[ p ]]$stained[ keep, , drop = FALSE ], fluorophore ) )

  # y axis: channel with the largest |SSI| over all non-reference (bead) unmixings
  beads <- setdiff( names( spectra.list ), reference )
  ssi.all <- sapply( beads, function( p ) ssi.of( p )$SSI )
  ssi.all <- matrix( ssi.all, ncol = length( beads ), dimnames = list( rownames( panel ), beads ) )
  ssi.all[ c( fluorophore, intersect( y.exclude, rownames( ssi.all ) ) ), ] <- NA
  best <- which( abs( ssi.all ) == max( abs( ssi.all ), na.rm = TRUE ), arr.ind = TRUE )[ 1, ]
  y.channel.auto <- rownames( ssi.all )[ best[ 1 ] ]
  message( "Largest |SSI| over the bead unmixings: ", y.channel.auto, " (", colnames( ssi.all )[ best[ 2 ] ],
           ", SSI ", round( ssi.all[ best[ 1 ], best[ 2 ] ], 2 ), ")" )
  if ( is.null( y.channel ) || is.na( y.channel ) || y.channel == "" ) {
    y.channel <- y.channel.auto
  } else {
    if ( !y.channel %in% rownames( panel ) )
      stop( "y.channel '", y.channel, "' is not an unmixed channel. Available: ",
            paste( rownames( panel ), collapse = ", " ), call. = FALSE )
    if ( y.channel == fluorophore ) stop( "y.channel must differ from the stained fluorophore", call. = FALSE )
    message( "y axis set by the user: ", y.channel, " (automatic choice: ", y.channel.auto, ")" )
  }

  ssi <- do.call( rbind, lapply( names( spectra.list ), function( p ) {
    r <- ssi.of( p )
    data.frame( Particle = p, Fluorophore = fluorophore, Y.channel = y.channel,
                N.events = length( keep ), SSI = r[ y.channel, "SSI" ],
                Spillover = r[ y.channel, "Spillover" ], stringsAsFactors = FALSE )
  } ) )

  ref.s <- spectra.list[[ reference ]][ fluorophore, ]
  spectrum.comparison <- do.call( rbind, lapply( names( spectra.list ), function( p ) {
    s <- spectra.list[[ p ]][ fluorophore, ]
    cs <- sum( ref.s * s ) / sqrt( sum( ref.s^2 ) * sum( s^2 ) )
    data.frame( Particle = p, Cosine.vs.reference = cs,
                Angle.vs.reference.deg = acos( min( 1, max( -1, cs ) ) ) * 180 / pi,
                Condition.number = calculate.condition.number( spectra.list[[ p ]] ),
                stringsAsFactors = FALSE )
  } ) )

  list(
    raw = raw, unmixed = unmixed, keep = keep, threshold = thr,
    y.channel = y.channel, y.channel.auto = y.channel.auto, ssi = ssi, ssi.all = ssi.all,
    spectra.list = spectra.list, af.spectrum = af.spectrum,
    condition.number = stats::setNames( spectrum.comparison$Condition.number, spectrum.comparison$Particle ),
    spectrum.comparison = spectrum.comparison
  )
}

# suggest.panel.swaps() is an exported AutoSpectral function
# (R/suggest_panel_swaps.R).
