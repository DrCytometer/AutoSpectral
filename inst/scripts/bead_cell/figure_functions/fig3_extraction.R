# fig3_extraction.R
#
# Event-level helpers for Figure 3 (streptavidin-biotin experiment,
# 2025-01-10, Aurora 5L). They reproduce, on explicitly chosen sets of
# events, the steps the AutoSpectral legacy pipeline uses to turn a gated
# single-stained control plus a negative into a reference spectrum:
#
#   * reading                                read.fcs.clean() (saturated events
#                                            and doublets removed)
#   * landmark scatter gates per cell type   define.gate.landmarks() procedure
#                                            via gate.scatter.match()
#   * density scatter gate for beads         define.gate.density()
#   * gating                                 apply.gate()
#   * positive events: above the asp$positivity.threshold quantile of the
#     negative in the peak channel, brightest asp$positive.n
#                                            (get.universal.negative())
#   * negative events, one of
#       "internal"  lowest asp$negative.n events of the stained file itself
#                   (downsample.control(), i.e. clean.controls(
#                   universal.negative = FALSE))
#       "random"    asp$negative.n random events of a separate unstained
#                   population (clean.controls(scatter.match = FALSE))
#       "knn"       the k nearest unstained neighbours in FSC/SSC of every
#                   selected positive event, unique union
#                   (clean.controls(scatter.match = TRUE) /
#                   get.universal.negative())
#   * spectrum: robust (rlm) slope of every detector on the peak detector
#     over the pooled positive + negative events, then divided by its
#     maximum (get.fluorophore.spectra(use.clean.expr = TRUE))
#
# Doing this outside define.flow.control() is what allows the positive
# events of one cell type to be paired with the background of another
# (panel c), which the control-file machinery does not allow.

# ---- files -------------------------------------------------------------------

#' Find the single-stained SA files of one particle type
#'
#' File names look like
#'   "A5 Splenocytes_Ly-6G-biotin_SA-FITC_Splenocytes SA.fcs"
#'   "C1 BD_CD4-biotin_SA-PCPCy55_BD Comp SA.fcs"
#'
#' @param dir Folder with the files.
#' @param markers Marker names as written in the files ("CD4", "Ly-6G", ...).
#' @param fluor.tags Named character vector: names = fluorophore names used in
#'   the spectral library, values = tag used in the file names.
#' @return data.frame(marker, fluor, file), one row per marker x fluorophore.
find.sa.files <- function( dir, markers, fluor.tags ) {

  if ( !dir.exists( dir ) ) stop( "Folder not found: ", dir, call. = FALSE )
  all.files <- list.files( dir, pattern = "\\.fcs$" )

  out <- do.call( rbind, lapply( markers, function( mk ) {
    do.call( rbind, lapply( names( fluor.tags ), function( fl ) {
      pattern <- paste0( "_", mk, "-biotin_", fluor.tags[[ fl ]], "_" )
      hit <- all.files[ grepl( pattern, all.files, fixed = TRUE ) ]
      if ( length( hit ) != 1 )
        stop( "Expected one file matching '", pattern, "' in ", dir, ", found ",
              length( hit ), call. = FALSE )
      data.frame( marker = mk, fluor = fl, file = file.path( dir, hit ),
                  stringsAsFactors = FALSE )
    } ) )
  } ) )
  rownames( out ) <- paste( out$fluor, out$marker, sep = "|" )
  out
}

#' Read one or more FCS files (pooled) with AutoSpectral's read.fcs.clean()
#' (as in get.spectra.automated()): events saturated in any spectral detector
#' or in FSC-A/SSC-A are removed, then doublets by the two-pass scatter-ratio
#' rule (FSC-A/FSC-H below its `singlet.quantiles[1]` quantile, then
#' SSC-A/SSC-H below its `singlet.quantiles[2]` quantile), per file.
#'
#' Removing scatter-saturated events matters here: in these files the
#' brightest events in the peak detector are mostly aggregates at the FSC/SSC
#' ceiling, which would otherwise pull the landmark gates into that corner.
#' @return matrix with the scatter and spectral columns.
read.events <- function( files, spectral.channels, asp,
                         singlet.quantiles = c( 0.85, 0.975 ), remove.doublets = TRUE,
                         verbose = FALSE ) {
  sc <- asp$default.scatter.parameter
  m <- do.call( rbind, lapply( files, function( f ) read.fcs.clean(
    path = f, label = basename( f ), spectral.channels = spectral.channels,
    scatter.channels = sc, sat.value = asp$expr.data.max,
    singlet.quantiles = singlet.quantiles, remove.doublets = remove.doublets,
    asp = asp, verbose = verbose ) ) )
  m[ , union( sc, spectral.channels ), drop = FALSE ]
}

#' Logical index of events inside a scatter gate (same test as apply.gate())
in.gate <- function( m, gate, asp ) {
  pip <- utils::getFromNamespace( ".point.in.polygon", "AutoSpectral" )
  sc <- asp$default.scatter.parameter
  pip( m[ , sc[ 1 ] ], m[ , sc[ 2 ] ], gate$x, gate$y ) != 0
}

# a gate name that passes AutoSpectral's control-file checks
.gate.label <- function( x ) gsub( "[^A-Za-z0-9]", ".", x )

# ---- gates -------------------------------------------------------------------

#' One landmark FSC/SSC gate per marker
#'
#' Same procedure as define.gate.landmarks(): for each marker, the brightest
#' `n.cells` events in the peak detector of each of its stained files (one per
#' fluorophore) are pooled and the gate is the convex hull of the densest
#' `percentile` % of them (2-D kernel density, `grid.n` grid; the package's
#' gate.scatter.match() does exactly this). The difference is that the events
#' come from read.fcs.clean() (see read.events()), whereas
#' define.gate.landmarks() ranks the raw file, where the brightest events are
#' scatter-saturated aggregates in these data.
#'
#' @param events Named list of cleaned event matrices (names = rownames(files)).
#' @param files data.frame from find.sa.files().
#' @param peak.channels Named character vector, fluorophore -> peak detector.
#' @param background Optional matrix (e.g. unstained cells) drawn in grey on
#'   the gate plots.
#' @return Named list of gates (names = markers); a PNG of each gate is
#'   written to `output.dir`.
define.marker.landmark.gates <- function(
    events, files, peak.channels, asp, output.dir,
    n.cells = 2000, percentile = 70, grid.n = 100, bandwidth.factor = 1,
    background = NULL
) {
  if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
  sc <- asp$default.scatter.parameter
  gate.fun <- utils::getFromNamespace( "gate.scatter.match", "AutoSpectral" )

  gates <- lapply( unique( files$marker ), function( mk ) {
    ids <- rownames( files )[ files$marker == mk ]
    top <- do.call( rbind, lapply( ids, function( id ) {
      m <- events[[ id ]]
      pk <- peak.channels[[ files[ id, "fluor" ] ]]
      m[ utils::head( order( m[ , pk ], decreasing = TRUE ), n.cells ), sc, drop = FALSE ]
    } ) )
    gate <- gate.fun( top, percentile = percentile / 100, bw.factor = bandwidth.factor, grid.n = grid.n )

    # record of the gate
    p <- ggplot2::ggplot() +
      { if ( !is.null( background ) )
          scattermore::geom_scattermore(
            data = data.frame( x = background[ , sc[ 1 ] ], y = background[ , sc[ 2 ] ] ),
            ggplot2::aes( x, y ), colour = "grey75", pointsize = 1 ) } +
      scattermore::geom_scattermore( data = data.frame( x = top[ , 1 ], y = top[ , 2 ] ),
                                     ggplot2::aes( x, y ), colour = "#21908C", pointsize = 1.5 ) +
      ggplot2::geom_polygon( data = data.frame( x = gate$x, y = gate$y ), ggplot2::aes( x, y ),
                             fill = NA, colour = "black" ) +
      ggplot2::coord_cartesian( xlim = c( 0, asp$scatter.data.max.x ),
                                ylim = c( 0, asp$scatter.data.max.y ) ) +
      ggplot2::labs( x = sc[ 1 ], y = sc[ 2 ],
                     title = sprintf( "%s landmark gate (%d brightest events x %d files)",
                                      mk, n.cells, length( ids ) ),
                     subtitle = "teal: landmark events; grey: unstained cells" ) +
      ggplot2::theme_bw( base_size = 9 )
    ggplot2::ggsave( file.path( output.dir, sprintf( "Fig3_landmark_gate_%s.png", .gate.label( mk ) ) ),
                     p, width = 4.5, height = 4.5, dpi = 150, bg = "white" )
    gate
  } )
  names( gates ) <- unique( files$marker )
  gates
}

#' Density FSC/SSC gate for the beads (pooled files, staged in one folder)
define.bead.density.gate <- function( bead.files, stage.dir, asp, output.dir,
                                      gate.name = "BD.beads", verbose = TRUE ) {
  if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
  staged.names <- sprintf( "bead_%02d.fcs", seq_along( bead.files ) )
  .stage.fcs( bead.files, staged.names, stage.dir )   # from extract_reference_spectrum.R
  gate.table <- data.frame(
    filename           = staged.names,
    fluorophore        = NA_character_,
    marker             = NA_character_,
    channel            = NA_character_,
    control.type       = "beads",
    universal.negative = NA_character_,
    large.gate         = NA,
    gate.name          = gate.name,
    gate.define        = TRUE,
    is.viability       = NA,
    stringsAsFactors   = FALSE
  )
  define.gate.density(
    control.file  = NULL,
    control.dir   = stage.dir,
    asp           = asp,
    gate.name     = gate.name,
    output.dir    = output.dir,
    verbose       = verbose,
    control.table = gate.table
  )
}

# ---- event selection ---------------------------------------------------------

#' Positive and negative events for one spectrum
#'
#' @param stained Gated stained events (matrix: scatter + spectral).
#' @param peak.channel Peak detector of the fluorophore.
#' @param mode "internal", "random" or "knn" (see top of file).
#' @param background Matrix of unstained events (ignored for "internal"). For
#'   "random" this is already restricted to the chosen gate; for "knn" it is
#'   the pool searched for scatter neighbours.
#' @param positive.n,negative.n,threshold As asp$positive.n, asp$negative.n,
#'   asp$positivity.threshold.
#' @param k.neighbors Neighbours per positive event ("knn").
#' @param seed Random seed (asp$bird.seed).
#' @param pos Optional matrix of positive events to use as they are (the
#'   selection step is skipped and only the negative is chosen). Used so that
#'   two spectra compared in panel c share exactly the same positive events.
#' @return list(pos, neg) (matrices) and the threshold used.
select.events <- function(
    stained, peak.channel, mode = c( "internal", "random", "knn" ),
    background = NULL, asp,
    positive.n = asp$positive.n, negative.n = asp$negative.n,
    threshold = asp$positivity.threshold, k.neighbors = 3L,
    seed = asp$bird.seed, pos = NULL
) {
  mode <- match.arg( mode )
  peak <- stained[ , peak.channel ]

  if ( mode == "internal" ) {
    # downsample.control(): brightest positive.n and dimmest negative.n
    o <- order( peak, decreasing = TRUE )
    if ( is.null( pos ) ) pos <- stained[ utils::head( o, positive.n ), , drop = FALSE ]
    neg <- stained[ utils::head( rev( o ), negative.n ), , drop = FALSE ]
    return( list( pos = pos, neg = neg, threshold = NA_real_ ) )
  }

  if ( is.null( background ) || nrow( background ) < 20 )
    stop( "Too few background events (", NROW( background ), ").", call. = FALSE )

  # get.universal.negative(): positives above the threshold quantile of the negative
  thr <- stats::quantile( background[ , peak.channel ], threshold )
  if ( is.null( pos ) ) {
    above <- which( peak > thr )
    if ( length( above ) < asp$min.cell.stop.n )
      stop( "Only ", length( above ), " positive events above the threshold.", call. = FALSE )
    if ( length( above ) < asp$min.cell.warning.n )
      warning( "Only ", length( above ), " positive events above the threshold.", call. = FALSE )
    above <- above[ order( peak[ above ], decreasing = TRUE ) ]
    pos <- stained[ utils::head( above, positive.n ), , drop = FALSE ]
  }

  if ( mode == "random" ) {
    set.seed( seed )
    idx <- if ( nrow( background ) > negative.n ) sample( nrow( background ), negative.n ) else
      seq_len( nrow( background ) )
    neg <- background[ idx, , drop = FALSE ]
  } else {
    sc <- asp$default.scatter.parameter
    knn <- FNN::knnx.index( data = background[ , sc, drop = FALSE ],
                            query = pos[ , sc, drop = FALSE ],
                            k = min( k.neighbors, nrow( background ) ) )
    neg <- background[ unique( as.vector( knn ) ), , drop = FALSE ]
  }
  list( pos = pos, neg = neg, threshold = unname( thr ) )
}

# ---- spectrum ----------------------------------------------------------------

#' Reference spectrum from positive + negative events (get.fluorophore.spectra)
#'
#' @param method "rlm" (package default: robust slope per detector, as in
#'   get.fluorophore.spectra()) or "ols" (ordinary least squares slope, used for
#'   speed in the bootstrap only if requested).
#' @return Named numeric vector, maximum = 1.
spectrum.from.events <- function( pos, neg, peak.channel, spectral.channels, asp,
                                  method = c( "rlm", "ols" ) ) {
  method <- match.arg( method )
  d <- rbind( pos[ , spectral.channels, drop = FALSE ], neg[ , spectral.channels, drop = FALSE ] )
  x <- d[ , peak.channel ]
  coef <- vapply( spectral.channels, function( ch ) {
    if ( ch == peak.channel ) return( 1 )
    if ( method == "rlm" )
      fit.robust.linear.model( x, d[ , ch ], peak.channel, ch, asp$rlm.iter.max )[ 2 ]
    else
      stats::cov( x, d[ , ch ] ) / stats::var( x )
  }, numeric( 1 ) )
  coef / max( coef )
}

#' One spectrum, end to end
extract.arm <- function( stained, peak.channel, spectral.channels, asp, mode,
                         background = NULL, k.neighbors = 3L, method = "rlm", pos = NULL ) {
  ev <- select.events( stained, peak.channel, mode, background, asp,
                       k.neighbors = k.neighbors, pos = pos )
  list(
    spectrum = spectrum.from.events( ev$pos, ev$neg, peak.channel, spectral.channels, asp, method ),
    pos = ev$pos, neg = ev$neg, threshold = ev$threshold,
    n.pos = nrow( ev$pos ), n.neg = nrow( ev$neg )
  )
}

#' Bootstrap of the difference between two spectra of the SAME positive events
#' with two different negatives (delta = spectrum(pos, neg.b) - spectrum(pos, neg.a))
#'
#' Positive events and both negative sets are resampled with replacement.
#' @return matrix n.boot x detectors.
bootstrap.delta <- function( pos, neg.a, neg.b, peak.channel, spectral.channels, asp,
                             n.boot = 200, method = "rlm", cores = 1L, seed = 42 ) {
  one <- function( b ) {
    set.seed( seed + b )
    p  <- pos[ sample( nrow( pos ), replace = TRUE ), , drop = FALSE ]
    na <- neg.a[ sample( nrow( neg.a ), replace = TRUE ), , drop = FALSE ]
    nb <- neg.b[ sample( nrow( neg.b ), replace = TRUE ), , drop = FALSE ]
    spectrum.from.events( p, nb, peak.channel, spectral.channels, asp, method ) -
      spectrum.from.events( p, na, peak.channel, spectral.channels, asp, method )
  }
  res <- if ( cores > 1 && .Platform$OS.type == "unix" )
    parallel::mclapply( seq_len( n.boot ), one, mc.cores = cores ) else
    lapply( seq_len( n.boot ), one )
  do.call( rbind, res )
}
