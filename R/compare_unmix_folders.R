# compare_unmix_folders.R

## Internal helper. Selects the positive and negative event groups shared by
## every metric for one (folder, fluorophore) combination: positive events are
## single-stained cells whose on-target channel exceeds the
## `positivity.quantile` quantile of the unstained control in that same
## channel. `pos.all` keeps every such event (needed for the spillover ratio
## and false-positive rate, both of which need the true positive population,
## not a brightest-only subset); `pos.capped` keeps only the brightest
## `n.positive` of them (matching the convention used by `calculate.ssi()`,
## for robustness and speed on the SSI/delta-MFI/Mahalanobis metrics) — if
## `n.positive` exceeds the number of events actually above threshold, every
## positive event is used. `neg.capped` is a random subsample of the
## unstained control of at most `max.negative` events. Returns `NULL` (with a
## warning) if fewer than 10 positive events are found.
##
## @keywords internal
.get.positive.negative.groups <- function(
    unstained.mat,
    single.stained.mat,
    on.channel,
    positivity.quantile,
    n.positive,
    max.negative,
    seed,
    label,
    fluorophore
) {

  threshold <- stats::quantile( unstained.mat[ , on.channel ], positivity.quantile )
  pos.mask  <- single.stained.mat[ , on.channel ] > threshold
  pos.all   <- single.stained.mat[ pos.mask, , drop = FALSE ]

  if ( nrow( pos.all ) < 10 ) {
    warning(
      paste0(
        "Folder '", label, "', fluorophore '", fluorophore,
        "': fewer than 10 positive events found; skipping."
      ),
      call. = FALSE
    )
    return( NULL )
  }

  if ( nrow( pos.all ) > n.positive ) {
    top.idx <- order(
      pos.all[ , on.channel ], decreasing = TRUE
    )[ seq_len( n.positive ) ]
    pos.capped <- pos.all[ top.idx, , drop = FALSE ]
  } else {
    pos.capped <- pos.all
  }

  if ( nrow( unstained.mat ) > max.negative ) {
    set.seed( seed )
    neg.idx    <- sample( nrow( unstained.mat ), max.negative )
    neg.capped <- unstained.mat[ neg.idx, , drop = FALSE ]
  } else {
    neg.capped <- unstained.mat
  }

  list( pos.all = pos.all, pos.capped = pos.capped, neg.capped = neg.capped )
}


## Internal helper. Computes SSI, |delta MFI|, spillover ratio, and
## false-positive rate for every off-target channel of one (folder,
## fluorophore) combination. `off.channels` is a named vector (names =
## off-target fluorophore identity, values = matched detector channel);
## `off.thresholds` is a named vector (names = fluorophore identity) of
## precomputed unstained-derived positivity thresholds for every channel in
## the panel, so the same quantile is not recomputed for every fluorophore
## pairing.
##
## The spillover ratio is deliberately a plain median-based calculation
## (`(median(pos) - median(neg)) / median(neg)`, on every positive event, not
## a robust/outlier-downweighted regression): a robust fit would clean up
## exactly the AF- and outlier-driven spillover this metric is meant to
## surface, so plain medians are used instead.
##
## @keywords internal
.calculate.channel.metrics <- function(
    pos.capped,
    pos.all,
    neg.capped,
    off.channels,
    off.thresholds
) {

  off.fluors <- names( off.channels )
  n.ch       <- length( off.channels )

  ssi             <- numeric( n.ch )
  delta.mfi       <- numeric( n.ch )
  spillover.ratio <- numeric( n.ch )
  fpr             <- numeric( n.ch )

  for ( j in seq_len( n.ch ) ) {

    ch <- off.channels[[ j ]]
    fl <- off.fluors[ j ]

    median.neg.capped <- stats::median( neg.capped[ , ch ] )
    median.pos.capped <- stats::median( pos.capped[ , ch ] )
    mad.neg           <- stats::mad(    neg.capped[ , ch ] )

    delta.mfi[ j ] <- abs( median.pos.capped - median.neg.capped )
    ssi[ j ]       <- abs( median.pos.capped - median.neg.capped ) / ( 2 * mad.neg )

    # spillover ratio uses every positive event (not just the brightest
    # `n.positive`), against the same unstained reference used throughout
    median.pos.all <- stats::median( pos.all[ , ch ] )

    spillover.ratio[ j ] <- abs( ( median.pos.all - median.neg.capped ) / median.neg.capped )

    fpr[ j ] <- mean( pos.all[ , ch ] > off.thresholds[[ fl ]] )
  }

  data.frame(
    off.target.fluorophore = off.fluors,
    ssi             = ssi,
    delta.mfi       = delta.mfi,
    spillover.ratio = spillover.ratio,
    fpr             = fpr,
    stringsAsFactors = FALSE
  )
}


## Internal helper. Median Mahalanobis distance of the positive events from
## the unstained distribution, computed jointly over every off-target channel
## (the on-target channel is excluded beforehand by the caller). The
## covariance matrix is estimated from `neg.capped` and lightly ridge-
## regularized before pseudo-inversion (`MASS::ginv()`) to guard against
## near-singular covariance in large, collinear panels; `dimnames` are
## reattached after `ginv()`, which otherwise drops them.
##
## @keywords internal
.calculate.mahalanobis.summary <- function(
    pos.capped,
    neg.capped,
    off.channels,
    ridge.factor
) {

  off.ch <- unname( off.channels )

  neg.sub <- neg.capped[ , off.ch, drop = FALSE ]
  pos.sub <- pos.capped[ , off.ch, drop = FALSE ]

  mu    <- colMeans( neg.sub )
  sigma <- stats::cov( neg.sub )

  if ( ridge.factor > 0 ) {
    sigma <- sigma + diag( ridge.factor * mean( diag( sigma ) ), nrow( sigma ) )
  }

  sigma.inv           <- MASS::ginv( sigma )
  dimnames( sigma.inv ) <- dimnames( sigma )

  centered <- sweep( pos.sub, 2, mu, "-" )
  d2       <- rowSums( ( centered %*% sigma.inv ) * centered )
  dist     <- sqrt( pmax( d2, 0 ) )

  stats::median( dist )
}


## Internal helper. Converts an HSL color (hue in degrees [0, 360); saturation
## and lightness in [0, 1]) to a hex color string. Implemented locally
## (standard HSL -> RGB algorithm) rather than pulling in a colour-space
## dependency for something this small.
##
## @keywords internal
.hsl.to.hex <- function( h, s, l ) {

  h <- ( h %% 360 ) / 360

  if ( s == 0 ) {
    r <- g <- b <- l
  } else {

    q <- if ( l < 0.5 ) l * ( 1 + s ) else l + s - l * s
    p <- 2 * l - q

    hue.to.rgb <- function( p, q, t ) {
      if ( t < 0 ) t <- t + 1
      if ( t > 1 ) t <- t - 1
      if ( t < 1 / 6 ) return( p + ( q - p ) * 6 * t )
      if ( t < 1 / 2 ) return( q )
      if ( t < 2 / 3 ) return( p + ( q - p ) * ( 2 / 3 - t ) * 6 )
      p
    }

    r <- hue.to.rgb( p, q, h + 1 / 3 )
    g <- hue.to.rgb( p, q, h )
    b <- hue.to.rgb( p, q, h - 1 / 3 )
  }

  grDevices::rgb( r, g, b )
}


## Internal helper. Approximate perceptual luminance of a hex color (simple
## gamma-uncorrected weighted sum, not full WCAG relative luminance). Used
## only as a heuristic to keep shaded palette colors visible against a white
## plot background and distinguishable from the pure-black unstained points,
## not for any accessibility guarantee.
##
## @keywords internal
.approx.luminance <- function( hex ) {
  rgb.vals <- grDevices::col2rgb( hex ) / 255
  0.2126 * rgb.vals[ 1, ] + 0.7152 * rgb.vals[ 2, ] + 0.0722 * rgb.vals[ 3, ]
}


## Internal helper. For a fixed hue/saturation, finds the lightest and
## darkest HSL lightness values whose approximate luminance stays within
## [luminance.floor, luminance.cap]. Falls back to a single mid-range
## lightness if no value in the search grid satisfies both bounds (only
## possible for extreme saturation values).
##
## @keywords internal
.find.lightness.range <- function( h, s, luminance.floor, luminance.cap ) {

  l.grid   <- seq( 0.05, 0.95, by = 0.01 )
  lum.grid <- vapply(
    l.grid, function( l ) .approx.luminance( .hsl.to.hex( h, s, l ) ), numeric( 1 )
  )

  ok <- lum.grid >= luminance.floor & lum.grid <= luminance.cap

  if ( !any( ok ) ) {
    target <- mean( c( luminance.floor, luminance.cap ) )
    best.l <- l.grid[ which.min( abs( lum.grid - target ) ) ]
    return( c( l.min = best.l, l.max = best.l ) )
  }

  c( l.min = min( l.grid[ ok ] ), l.max = max( l.grid[ ok ] ) )
}


## Internal helper. Maps every fluorophore in `fluorophore.names` to a hex
## color for consistent use across every plot. Hue is set by excitation
## laser (Red -> red, YellowGreen -> yellow-green, Blue -> blue, Violet ->
## indigo, UV -> violet), matching how those laser lines are conventionally
## drawn on a spectral cytometer. Within a laser group, color is then shaded
## from pale to dark by rank of emission (`nominal.wavelength`): the
## shortest-emission fluorophore on that laser gets the palest shade, the
## longest gets the darkest, evenly spread across the group's available
## lightness range to maximise contrast for however many fluorophores share
## that laser in this panel. The available lightness range is itself bounded
## by an approximate-luminance heuristic so the palest shade stays visible
## against a white plot background and the darkest stays visually distinct
## from the pure-black unstained points. `excitation.laser` values of "YG"
## are treated as a synonym of "YellowGreen" (both appear in
## `fluorophore_database.csv`). Fluorophores that cannot be matched to
## `fluorophore.database`, or that lack a recorded laser or emission
## wavelength, fall back to mid-grey with a warning.
##
## @keywords internal
.build.fluorophore.palette <- function(
    fluorophore.names,
    fluorophore.database,
    luminance.floor = 0.18,
    luminance.cap   = 0.78
) {

  fluorophore.names <- unique( fluorophore.names[ !is.na( fluorophore.names ) ] )

  base.hue <- c(
    Red         = 355,
    YellowGreen = 75,
    Blue        = 215,
    Violet      = 260,
    UV          = 290
  )
  base.sat <- c(
    Red         = 0.75,
    YellowGreen = 0.75,
    Blue        = 0.70,
    Violet      = 0.65,
    UV          = 0.60
  )

  db <- fluorophore.database[
    , c( "fluorophore", "excitation.laser", "nominal.wavelength" ), drop = FALSE
  ]
  db$excitation.laser[ db$excitation.laser %in% "YG" ] <- "YellowGreen"

  lookup             <- db[ match( fluorophore.names, db$fluorophore ), ]
  rownames( lookup ) <- fluorophore.names

  unmatched <- fluorophore.names[
    is.na( lookup$excitation.laser ) | is.na( lookup$nominal.wavelength ) |
      !lookup$excitation.laser %in% names( base.hue )
  ]
  if ( length( unmatched ) > 0 ) {
    warning(
      paste0(
        "Could not resolve excitation laser / emission wavelength for: ",
        paste( unmatched, collapse = ", " ), "; plotted in grey."
      ),
      call. = FALSE
    )
  }

  color.map <- stats::setNames(
    rep( "grey50", length( fluorophore.names ) ), fluorophore.names
  )

  for ( laser in names( base.hue ) ) {

    in.group <- fluorophore.names[
      !fluorophore.names %in% unmatched & lookup$excitation.laser == laser
    ]
    if ( length( in.group ) == 0 ) next

    l.range <- .find.lightness.range(
      h = base.hue[[ laser ]], s = base.sat[[ laser ]],
      luminance.floor = luminance.floor, luminance.cap = luminance.cap
    )

    wl  <- lookup[ in.group, "nominal.wavelength" ]
    ord <- order( wl )
    n   <- length( in.group )

    l.vals <- if ( n == 1 ) {
      mean( c( l.range[[ "l.min" ]], l.range[[ "l.max" ]] ) )
    } else {
      seq( l.range[[ "l.max" ]], l.range[[ "l.min" ]], length.out = n )
    }

    hex <- vapply(
      l.vals,
      function( l ) .hsl.to.hex( base.hue[[ laser ]], base.sat[[ laser ]], l ),
      character( 1 )
    )

    color.map[ in.group[ ord ] ] <- hex
  }

  color.map
}


## are coloured by the `color.col` column of `df`, mapped through
## `color.map` (named vector, value -> hex color), when both are supplied;
## otherwise every point is plotted black. Set `log.scale = TRUE` for a
## log10 y-axis with intermediate (2-9x) tick marks on the left edge, via
## `annotation_logticks()`, so the axis reads as clearly logarithmic rather
## than looking like an ordinary axis with oddly-spaced labels; every
## metric plotted here is non-negative, but exact-zero values have no log
## and are dropped from the plot, with a warning.
##
## @keywords internal
.plot.metric.boxplot <- function(
    df,
    title,
    y.label,
    file.path.out,
    color.col           = NULL,
    color.map           = NULL,
    log.scale           = FALSE,
    plot.width          = 7,
    plot.height         = 5,
    base.font.size      = 11,
    title.size          = NULL,
    point.size          = 1.5,
    point.alpha         = 0.6,
    text.angle          = 45,
    legend.max.rows     = 25,
    legend.ncol         = NULL,
    legend.font.size    = NULL,
    legend.key.size     = 0.8,
    legend.width.per.col = 1.1
) {

  use.color <- !is.null( color.col ) && !is.null( color.map ) &&
    color.col %in% colnames( df )

  if ( use.color ) {
    df$.color.group <- df[[ color.col ]]

    n.legend.items <- length( unique( stats::na.omit( df$.color.group ) ) )

    legend.ncol.final <- if ( !is.null( legend.ncol ) ) {
      legend.ncol
    } else {
      max( 1L, ceiling( n.legend.items / legend.max.rows ) )
    }

    legend.font.size.final <- if ( !is.null( legend.font.size ) ) {
      legend.font.size
    } else {
      base.font.size * 0.7
    }
  }

  if ( log.scale ) {
    n.non.positive <- sum( df$value <= 0, na.rm = TRUE )
    if ( n.non.positive > 0 ) {
      warning(
        paste0(
          "'", title, "': ", n.non.positive,
          " point(s) with value <= 0 omitted from the log10 y-axis."
        ),
        call. = FALSE
      )
    }
  }

  p <- ggplot2::ggplot( df, ggplot2::aes( x = folder, y = value ) ) +
    ggplot2::geom_boxplot( outlier.shape = NA, fill = NA )

  p <- if ( use.color ) {
    p +
      ggplot2::geom_jitter(
        ggplot2::aes( color = .color.group ),
        width = 0.15, alpha = point.alpha, size = point.size
      ) +
      ggplot2::scale_color_manual(
        values = color.map,
        name   = "Fluorophore",
        guide  = ggplot2::guide_legend(
          ncol         = legend.ncol.final,
          override.aes = list( size = 3, alpha = 1 )
        )
      )
  } else {
    p +
      ggplot2::geom_jitter(
        width = 0.15, alpha = point.alpha, size = point.size, color = "black"
      )
  }

  if ( log.scale ) {
    p <- p +
      ggplot2::scale_y_log10() +
      ggplot2::annotation_logticks( sides = "l" )
  }

  p <- p +
    ggplot2::labs( title = title, x = "Unmixing", y = y.label ) +
    ggplot2::theme_classic( base_size = base.font.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = text.angle, hjust = 1 )
    )

  if ( use.color ) {
    p <- p +
      ggplot2::theme(
        legend.text     = ggplot2::element_text( size = legend.font.size.final ),
        legend.title    = ggplot2::element_text( size = legend.font.size.final + 1 ),
        legend.key.size = grid::unit( legend.key.size, "lines" ),
        legend.position = "right"
      )
  }

  if ( !is.null( title.size ) ) {
    p <- p + ggplot2::theme( plot.title = ggplot2::element_text( size = title.size ) )
  }

  effective.width <- if ( use.color ) {
    plot.width + max( 0L, legend.ncol.final - 1L ) * legend.width.per.col
  } else {
    plot.width
  }

  ggplot2::ggsave(
    file.path.out, p, device = ragg::agg_jpeg, width = effective.width, height = plot.height
  )
  print( p )

  invisible( p )
}


#' @title Compare Unmixing Quality Across Operators or Methods
#'
#' @description
#' Computes six unmixing-quality metrics from folders of already-unmixed FCS
#' files (for example, one folder per operator plus one for AutoSpectral),
#' using the setup tables written by `setup.unmix.comparison()`. By default,
#' produces one overview figure per metric — one point per fluorophore, one
#' column per folder, box-and-whisker overlay — plus the full per-off-target-
#' channel results and a per-fluorophore summary, both as CSVs. Set
#' `plot.per.fluorophore = TRUE` to additionally get one figure per
#' fluorophore per metric, for closer inspection.
#'
#' Every value is keyed by fluorophore identity rather than raw detector name,
#' so results are directly comparable across folders even when different
#' unmixing software uses different channel-naming conventions.
#'
#' Metrics:
#' \enumerate{
#'   \item \strong{Unstained rSD (MAD)}: robust standard deviation of the
#'     unstained control in every mapped fluorescence channel. One figure
#'     total (not per fluorophore), since it does not depend on any
#'     single-stained control.
#'   \item \strong{Secondary Stain Index (SSI)}: `abs(median(pos) - median(neg))
#'     / (2 * rSD(neg))`, in every off-target channel. Always non-negative,
#'     since only the magnitude of separation from the unstained control is
#'     of interest, not its sign.
#'   \item \strong{Delta MFI}: `abs(median(pos) - median(neg))` in every
#'     off-target channel.
#'   \item \strong{Spillover ratio}: `(median(pos) - median(neg)) /
#'     median(neg)`, fitted to every positive event. Deliberately a plain
#'     median-based ratio rather than a robust regression: a robust fit would
#'     downweight exactly the AF- and outlier-driven spillover this metric is
#'     meant to be sensitive to.
#'   \item \strong{False-positive rate}: fraction of positive events whose
#'     off-target channel value exceeds the unstained-derived positivity
#'     threshold for that channel.
#'   \item \strong{Mahalanobis distance}: median multivariate distance of
#'     positive events from the unstained distribution, computed over all
#'     off-target channels jointly. This is a single summary value per
#'     folder/fluorophore rather than a per-channel metric, so it is plotted
#'     as one point per folder rather than a box-and-whisker.
#' }
#'
#' Positive events for a fluorophore are defined, in every folder, as events
#' in that fluorophore's single-stained control whose on-target channel value
#' exceeds the `positivity.quantile` quantile of the unstained control in that
#' same channel. Metrics 2, 3, and 6 use up to `n.positive` of the brightest
#' such events (matching the convention used elsewhere in AutoSpectral by
#' `calculate.ssi()`) — if fewer than `n.positive` events clear the threshold,
#' every event that does is used. Metrics 4 and 5 always use every positive
#' event, since the spillover ratio and false-positive rate both need the
#' true positive population rather than a brightest-only subset.
#'
#' @param folders Named character vector of directory paths, one per unmixed
#' result set (must use the same names as when calling
#' `setup.unmix.comparison()`).
#' @param setup.files Named character vector of paths to the (optionally
#' hand-edited) setup CSVs produced by `setup.unmix.comparison()`. Names must
#' match `folders`.
#' @param positivity.quantile Numeric, default `0.99`. Quantile of the
#' unstained control, in the relevant channel, used to threshold positive
#' events and (per off-target channel) false-positive events.
#' @param n.positive Numeric, default `2000`. Maximum number of brightest
#' positive events used for metrics 2, 3, and 6. If more than `n.positive`
#' events clear the positivity threshold, only the brightest `n.positive` are
#' used; if fewer clear it, every positive event is used.
#' @param max.negative Numeric, default `2000`. Maximum number of unstained
#' events randomly sampled as the negative reference.
#' @param ridge.factor Numeric, default `1e-4`. Relative diagonal shrinkage
#' applied to the off-target covariance matrix before inversion for metric 6,
#' to guard against near-singular covariance in large panels. Set to `0` to
#' disable.
#' @param seed Numeric, default `42`. Random seed used for all downsampling.
#' @param fluorophore.database Data frame of fluorophore names and metadata
#' (must include `fluorophore`, `excitation.laser`, and `nominal.wavelength`
#' columns), used only to build the fluorophore color palette for the
#' figures. Default `NULL` loads the bundled `fluorophore_database.csv`.
#' @param plot.dir Character. Directory for output figures. Created if
#' absent. Default `"./figure_unmix_comparison"`.
#' @param output.csv Character. Path for the long-format, per-off-target-
#' channel results CSV. Default `"unmix_comparison_results.csv"`.
#' @param summary.csv Character. Path for the per-fluorophore summary CSV that
#' the overview figures are built from: for `SSI`, `Delta.MFI`,
#' `Spillover.ratio`, and `FPR` this is the sum of that
#' fluorophore's off-target-channel values; for `Mahalanobis` it is the same
#' per-fluorophore median distance already computed (no further summarizing
#' possible, since it has no per-channel spread). Default
#' `"unmix_comparison_summary.csv"`.
#' @param plot.per.fluorophore Logical, default `FALSE`. The overview figures
#' (one point per fluorophore, one column per folder) are always produced.
#' Set `TRUE` to additionally produce the more granular one-figure-per-
#' fluorophore plots (one point per off-target channel) for closer
#' inspection of a specific fluorophore.
#' @param log.scale Logical, default `FALSE`. Plot the y-axis on a log10
#' scale. Passed straight through to `unmix.comparison.plot()`; see there
#' for the caveat on exact-zero values.
#' @param plot.width,plot.height Numeric, defaults `7` and `5` (inches).
#' Dimensions passed to `ggplot2::ggsave()` for every figure.
#' @param base.font.size Numeric, default `11`. Base font size (points)
#' passed to `ggplot2::theme_classic()`; every other text element (axis
#' text, legend) scales from this.
#' @param title.size Numeric, default `NULL`. Explicit plot-title font size
#' (points). `NULL` leaves the `theme_classic()` default (scaled from
#' `base.font.size`) in place.
#' @param point.size Numeric, default `1.5`. Size of the jittered points.
#' @param point.alpha Numeric, default `0.6`. Opacity of the jittered
#' points, in `[0, 1]`.
#' @param text.angle Numeric, default `45`. Rotation angle (degrees) of the
#' x-axis (folder) labels.
#' @param legend.max.rows Numeric, default `25`. When a figure's legend
#' (fluorophore color key) would need more than this many entries in a
#' single column, it wraps into additional columns rather than
#' overflowing the top of the plot. Ignored when `legend.ncol` is set.
#' @param legend.ncol Numeric, default `NULL`. Fixes the number of legend
#' columns explicitly, overriding `legend.max.rows`. Useful for keeping a
#' consistent legend layout across a batch of figures with varying
#' fluorophore counts.
#' @param legend.font.size Numeric, default `NULL`. Font size (points) for
#' legend text. `NULL` scales it from `base.font.size` (70%).
#' @param legend.key.size Numeric, default `0.8`. Size (lines) of the
#' legend color swatches; smaller values fit more entries per column.
#' @param legend.width.per.col Numeric, default `1.1`. Extra figure width
#' (inches) added for each legend column beyond the first, so a wrapped
#' legend never crowds the plot panel.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list with two elements: `results`, the full
#' long-format data frame (columns `folder`, `on.target.fluorophore`,
#' `off.target.fluorophore`, `metric`, `value`) also written to
#' `output.csv`; and `summary`, the per-fluorophore roll-up (columns
#' `folder`, `fluorophore`, `metric`, `value`) also written to
#' `summary.csv`. Figures are written to `plot.dir` and printed to the active
#' graphics device.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_point labs
#' @importFrom ggplot2 theme_classic theme element_text ggsave scale_color_manual
#' @importFrom ggplot2 scale_y_log10 annotation_logticks
#' @importFrom ragg agg_jpeg
#' @importFrom MASS ginv
#' @importFrom grDevices rgb col2rgb
#'
#' @export

compare.unmix.folders <- function(
    folders,
    setup.files,
    positivity.quantile   = 0.99,
    n.positive            = 2000,
    max.negative          = 2000,
    ridge.factor          = 1e-4,
    seed                  = 42,
    fluorophore.database  = NULL,
    plot.dir              = "./figure_unmix_comparison",
    output.csv            = "unmix_comparison_results.csv",
    summary.csv           = "unmix_comparison_summary.csv",
    plot.per.fluorophore  = FALSE,
    log.scale             = FALSE,
    plot.width            = 7,
    plot.height           = 5,
    base.font.size        = 11,
    title.size            = NULL,
    point.size            = 1.5,
    point.alpha           = 0.6,
    text.angle            = 45,
    legend.max.rows       = 25,
    legend.ncol           = NULL,
    legend.font.size      = NULL,
    legend.key.size       = 0.8,
    legend.width.per.col  = 1.1,
    verbose               = TRUE
) {

  # --- input validation -------------------------------------------------

  if ( !setequal( names( folders ), names( setup.files ) ) ) {
    stop( "`folders` and `setup.files` must have identical names.", call. = FALSE )
  }

  missing.setup <- setup.files[ !file.exists( setup.files ) ]
  if ( length( missing.setup ) > 0 ) {
    stop(
      paste0( "Setup CSV(s) not found: ", paste( missing.setup, collapse = ", " ) ),
      call. = FALSE
    )
  }

  labels <- names( folders )

  # --- load setup tables ---------------------------------------------------

  setup.tables <- lapply( setup.files, function( f ) {
    tb <- utils::read.csv( f, stringsAsFactors = FALSE, strip.white = TRUE )
    tb[ tb == "" ] <- NA
    tb
  } )
  names( setup.tables ) <- labels

  bad.rows <- vapply(
    setup.tables, function( tb ) any( tb$flag != "OK", na.rm = TRUE ), logical( 1 )
  )
  if ( any( bad.rows ) ) {
    warning(
      paste0(
        "\033[31mUnresolved 'flag' entries remain in the setup table(s) for: ",
        paste( labels[ bad.rows ], collapse = ", " ),
        ". Affected fluorophores will be skipped.\033[0m"
      ),
      call. = FALSE
    )
  }

  all.results <- list()

  # --- process each folder --------------------------------------------------

  for ( label in labels ) {

    if ( verbose ) message( sprintf( "\033[34mProcessing folder: %s\033[0m", label ) )

    folder.path <- folders[[ label ]]
    tb          <- setup.tables[[ label ]]

    unstained.row <- tb[ tb$fluorophore == "Unstained", ]
    if ( nrow( unstained.row ) != 1 || is.na( unstained.row$filename[ 1 ] ) ) {
      warning(
        paste0( "Folder '", label, "': no usable unstained row; skipping folder." ),
        call. = FALSE
      )
      next
    }

    unstained.mat <- readFCS( file.path( folder.path, unstained.row$filename[ 1 ] ) )

    ok.rows <- tb[
      tb$fluorophore != "Unstained" &
        tb$flag == "OK" &
        !is.na( tb$channel )  & tb$channel  != "No match" &
        !is.na( tb$filename ) & tb$filename != "No match",
      , drop = FALSE
    ]

    if ( nrow( ok.rows ) == 0 ) {
      warning(
        paste0( "Folder '", label, "': no usable single-stained rows; skipping folder." ),
        call. = FALSE
      )
      next
    }

    fluor.channel.map  <- stats::setNames( ok.rows$channel, ok.rows$fluorophore )
    all.fluor.channels <- unname( fluor.channel.map )

    missing.channels <- setdiff( all.fluor.channels, colnames( unstained.mat ) )
    if ( length( missing.channels ) > 0 ) {
      stop(
        paste0(
          "Folder '", label, "': mapped channel(s) not found in the unstained ",
          "FCS file: ", paste( missing.channels, collapse = ", " )
        ),
        call. = FALSE
      )
    }

    # --- metric 1: unstained rSD (MAD), one row per fluorophore/channel ---

    mad.vals <- vapply(
      fluor.channel.map, function( ch ) stats::mad( unstained.mat[ , ch ] ), numeric( 1 )
    )

    all.results[[ length( all.results ) + 1 ]] <- data.frame(
      folder                  = label,
      on.target.fluorophore   = NA_character_,
      off.target.fluorophore  = names( mad.vals ),
      metric                  = "Unstained.rSD",
      value                   = unname( mad.vals ),
      stringsAsFactors        = FALSE
    )

    # --- precompute per-channel positivity thresholds once per folder -----

    off.thresholds <- vapply(
      fluor.channel.map,
      function( ch ) stats::quantile( unstained.mat[ , ch ], positivity.quantile ),
      numeric( 1 )
    )

    # --- per-fluorophore metrics -------------------------------------------

    for ( i in seq_len( nrow( ok.rows ) ) ) {

      fl         <- ok.rows$fluorophore[ i ]
      on.channel <- ok.rows$channel[ i ]
      ss.file    <- ok.rows$filename[ i ]

      off.fluors   <- setdiff( names( fluor.channel.map ), fl )
      off.channels <- fluor.channel.map[ off.fluors ]

      if ( length( off.channels ) == 0 ) {
        warning(
          paste0(
            "Folder '", label, "', fluorophore '", fl,
            "': no off-target channels available; skipping."
          ),
          call. = FALSE
        )
        next
      }

      single.mat <- readFCS( file.path( folder.path, ss.file ) )

      needed.channels <- unique( c( on.channel, off.channels ) )
      missing.ss       <- setdiff( needed.channels, colnames( single.mat ) )
      if ( length( missing.ss ) > 0 ) {
        warning(
          paste0(
            "Folder '", label, "', fluorophore '", fl,
            "': missing channel(s) in single-stained file (",
            paste( missing.ss, collapse = ", " ), "); skipping."
          ),
          call. = FALSE
        )
        next
      }

      groups <- .get.positive.negative.groups(
        unstained.mat       = unstained.mat,
        single.stained.mat  = single.mat,
        on.channel          = on.channel,
        positivity.quantile = positivity.quantile,
        n.positive          = n.positive,
        max.negative        = max.negative,
        seed                = seed,
        label               = label,
        fluorophore         = fl
      )

      if ( is.null( groups ) ) next

      # --- metrics 2, 3, 4, 5: per off-target channel ---

      channel.metrics <- .calculate.channel.metrics(
        pos.capped     = groups$pos.capped,
        pos.all        = groups$pos.all,
        neg.capped     = groups$neg.capped,
        off.channels   = off.channels,
        off.thresholds = off.thresholds
      )

      n.ch <- nrow( channel.metrics )

      metric.names <- c( "SSI", "Delta.MFI", "Spillover.ratio", "FPR" )

      value.vec <- c(
        channel.metrics$ssi, channel.metrics$delta.mfi,
        channel.metrics$spillover.ratio, channel.metrics$fpr
      )

      all.results[[ length( all.results ) + 1 ]] <- data.frame(
        folder                  = label,
        on.target.fluorophore   = fl,
        off.target.fluorophore  = rep( channel.metrics$off.target.fluorophore, length( metric.names ) ),
        metric                  = rep( metric.names, each = n.ch ),
        value                   = value.vec,
        stringsAsFactors        = FALSE
      )

      # --- metric 6: single Mahalanobis summary ---

      mahal.value <- .calculate.mahalanobis.summary(
        pos.capped   = groups$pos.capped,
        neg.capped   = groups$neg.capped,
        off.channels = off.channels,
        ridge.factor = ridge.factor
      )

      all.results[[ length( all.results ) + 1 ]] <- data.frame(
        folder                  = label,
        on.target.fluorophore   = fl,
        off.target.fluorophore  = NA_character_,
        metric                  = "Mahalanobis",
        value                   = mahal.value,
        stringsAsFactors        = FALSE
      )
    }
  }

  if ( length( all.results ) == 0 ) {
    stop( "No usable results were produced from any folder.", call. = FALSE )
  }

  results.df <- do.call( rbind, all.results )
  rownames( results.df ) <- NULL

  utils::write.csv( results.df, output.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote results: %s\033[0m", output.csv ) )

  # --- per-fluorophore summary: sum across off-target channels for the
  # four per-channel metrics, plus the Mahalanobis median as-is (it is already
  # a single value per folder/fluorophore, with no per-channel spread to
  # summarize). This is the table the overview figures are built from.

  summary.rows <- list()

  sum.metrics <- c( "SSI", "Delta.MFI", "Spillover.ratio", "FPR" )

  for ( m in sum.metrics ) {

    metric.rows <- results.df[ results.df$metric == m, ]
    if ( nrow( metric.rows ) == 0 ) next

    agg <- stats::aggregate(
      value ~ folder + on.target.fluorophore,
      data = metric.rows,
      FUN  = sum
    )

    summary.rows[[ length( summary.rows ) + 1 ]] <- data.frame(
      folder      = agg$folder,
      fluorophore = agg$on.target.fluorophore,
      metric      = m,
      value       = agg$value,
      stringsAsFactors = FALSE
    )
  }

  mahal.rows <- results.df[ results.df$metric == "Mahalanobis", ]
  if ( nrow( mahal.rows ) > 0 ) {
    summary.rows[[ length( summary.rows ) + 1 ]] <- data.frame(
      folder      = mahal.rows$folder,
      fluorophore = mahal.rows$on.target.fluorophore,
      metric      = "Mahalanobis",
      value       = mahal.rows$value,
      stringsAsFactors = FALSE
    )
  }

  summary.df           <- do.call( rbind, summary.rows )
  rownames( summary.df ) <- NULL

  utils::write.csv( summary.df, summary.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote summary: %s\033[0m", summary.csv ) )

  # --- plots ---------------------------------------------------------------
  # delegated to unmix.comparison.plot() so the same plotting code can be
  # re-run standalone later, without repeating the unmixing above

  unmix.comparison.plot(
    results              = results.df,
    summary              = summary.df,
    setup.files          = setup.files,
    fluorophore.database = fluorophore.database,
    plot.dir             = plot.dir,
    plot.per.fluorophore = plot.per.fluorophore,
    log.scale            = log.scale,
    plot.width           = plot.width,
    plot.height          = plot.height,
    base.font.size       = base.font.size,
    title.size           = title.size,
    point.size           = point.size,
    point.alpha          = point.alpha,
    text.angle           = text.angle,
    legend.max.rows      = legend.max.rows,
    legend.ncol          = legend.ncol,
    legend.font.size     = legend.font.size,
    legend.key.size      = legend.key.size,
    legend.width.per.col = legend.width.per.col,
    verbose              = verbose
  )

  return( invisible( list( results = results.df, summary = summary.df ) ) )
}
