# plot_spectra_automated_steps.R
#
# Manuscript figure helper: illustrates the internal steps of
# get.spectra.automated() for one or more single-stained controls, using the
# same biplot / trace / scatter-match building blocks as the rest of
# AutoSpectral (create.biplot(), spectral.trace(), scatter.match.plot()) so
# the panels are stylistically consistent with the package's other output.
#
# Panels produced per fluorophore:
#   A. -A vs -H singlet gating (FSC and SSC, linear scale), black = retained,
#      grey = excluded, red = the final "clean positive" events selected by
#      the cosine-similarity filter (panel D) and carried into panel F
#   B. Peak-finding traces: unstained AF, pre-orthogonalization and
#      post-orthogonalization single-stained control (spectral.trace()),
#      with an arrow marking the identified empirical peak detector
#   C. Brightest-event candidate selection: KDE-smoothed 1D histogram of all
#      singlet-gated single-stained control events on the empirical peak
#      channel, with a bracket marking the range of the top n.candidates
#      events selected before cosine filtering
#   D. Cosine-similarity filter biplot: empirical peak channel (x) vs. the
#      non-colliding peak channel of the unstained AF vector (y), single-
#      stained control events coloured by a continuous cosine-similarity-
#      to-AF gradient; events selected for the final profile are plotted at
#      selected.size.mult times the size of the rest (legend split into its
#      own column to keep both biplots the same width)
#   E. kNN scatter-matched AF subtraction -- generated via the existing
#      scatter.match.plot() and embedded (or referenced) in the composite.
#   F. Final spectral profile comparison: the per-cell background-subtracted
#      signature ("<Fluorophore> (Cells)"), the reference profile for the
#      same fluorophore ("<Fluorophore> (Beads)" -- from a paired bead
#      control in control.def.file if provided, otherwise the cytometer's
#      static spectral reference library), and the matched-unstained AF
#      trace ("Unstained (Autofluorescence)") (spectral.trace()).

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

## Biexponential transform matching the package's "legacy" AF/gate plots
## (gate.af.sample.plot(), remove.af()), built from asp$default.transformation.param.
.biexp.transform.legacy <- function( asp ) {
  flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )
}

## Shared ggplot theme matching create.biplot() / gate.af.sample.plot().
.theme.biplot <- function( asp, legend.position = "none" ) {
  ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      legend.position   = legend.position,
      axis.ticks        = ggplot2::element_line( linewidth = asp$figure.panel.line.size ),
      axis.text         = ggplot2::element_text( size = asp$figure.axis.text.size ),
      axis.title        = ggplot2::element_text( size = asp$figure.axis.title.size ),
      panel.border      = ggplot2::element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank()
    )
}

## Shared x/y scale builder (biexp-transformed, asp$ribbon.breaks ticks),
## matching create.biplot()'s axis presentation.
.biplot.scales <- function( biexp, asp, x.lab, y.lab ) {
  breaks <- asp$ribbon.breaks
  limits <- c( asp$ribbon.plot.min, asp$expr.data.max )
  axis.labels <- sapply( breaks, function( v ) {
    if ( v == 0 ) "0" else parse( text = paste0( "10^", log10( abs( v ) ) ) )
  } )

  list(
    ggplot2::scale_x_continuous(
      name = x.lab, breaks = biexp( breaks ), limits = biexp( limits ), labels = axis.labels
    ),
    ggplot2::scale_y_continuous(
      name = y.lab, breaks = biexp( breaks ), limits = biexp( limits ), labels = axis.labels
    )
  )
}

## Read a single-stained control FCS file, apply saturation removal, and flag
## (rather than drop) events failing the two-stage singlet discrimination, so
## both retained and excluded events can be plotted together.
.read.fcs.raw.singlet.flags <- function(
    path, spectral.channels, scatter.channels, sat.value, singlet.quantiles
) {
  fsc.a <- scatter.channels[ 1L ]
  fsc.h <- sub( "-A$", "-H", fsc.a )
  ssc.a <- if ( length( scatter.channels ) >= 2L ) scatter.channels[ 2L ] else NA_character_
  ssc.h <- if ( !is.na( ssc.a ) ) sub( "-A$", "-H", ssc.a ) else NA_character_

  height.channels <- sub( "-A$", "-H", scatter.channels )
  cols.keep       <- c( scatter.channels, height.channels, spectral.channels )

  probe.cols <- colnames( readFCS( path, start.row = 1, end.row = 1 ) )
  present    <- intersect( cols.keep, probe.cols )
  mat        <- readFCS( path, columns = present )

  spec.present <- intersect( spectral.channels, colnames( mat ) )
  if ( length( spec.present ) > 0 ) {
    sat.keep <- rowSums( mat[ , spec.present, drop = FALSE ] >= sat.value ) == 0
    mat      <- mat[ sat.keep, , drop = FALSE ]
  }

  fsc.keep <- rep( TRUE, nrow( mat ) )
  if ( all( c( fsc.a, fsc.h ) %in% colnames( mat ) ) ) {
    fsc.ratio <- mat[ , fsc.a ] / ( mat[ , fsc.h ] + 1e-9 )
    fsc.cut   <- stats::quantile( fsc.ratio, probs = singlet.quantiles[ 1L ] )
    fsc.keep  <- fsc.ratio < fsc.cut
  }

  ssc.keep <- rep( TRUE, nrow( mat ) )
  if ( !is.na( ssc.a ) && all( c( ssc.a, ssc.h ) %in% colnames( mat ) ) &&
       length( singlet.quantiles ) >= 2L ) {
    ssc.ratio <- mat[ , ssc.a ] / ( mat[ , ssc.h ] + 1e-9 )
    # mirrors .read.fcs.clean: the second-stage quantile is computed on the
    # subset that already passed the first-stage (FSC) gate
    ssc.cut  <- stats::quantile( ssc.ratio[ fsc.keep ], probs = singlet.quantiles[ 2L ] )
    ssc.keep <- ssc.ratio < ssc.cut
  }

  list(
    mat = mat, retained = fsc.keep & ssc.keep,
    fsc.a = fsc.a, fsc.h = fsc.h, ssc.a = ssc.a, ssc.h = ssc.h
  )
}

## Single -A vs -H biplot panel coloured by singlet-gate retention. Scatter
## parameters (FSC/SSC) are conventionally shown on a linear scale, unlike the
## biexponential fluorescence channels used elsewhere -- axis breaks/labels
## mirror scatter.match.plot()'s "Ne6" convention.
.singlet.gate.panel <- function(
    mat, x.col, y.col, retained, asp,
    retained.color = "black", excluded.color = "grey75", max.points = 5e4
) {
  df <- data.frame(
    x = mat[ , x.col ], y = mat[ , y.col ],
    Gate = factor(
      ifelse( retained, "Retained", "Excluded" ), levels = c( "Excluded", "Retained" )
    )
  )

  if ( nrow( df ) > max.points ) {
    set.seed( asp$bird.seed )
    df <- df[ sample( seq_len( nrow( df ) ), max.points ), , drop = FALSE ]
  }
  # draw retained events on top of excluded events
  df <- df[ order( df$Gate == "Retained" ), ]

  x.limits <- c( asp$scatter.data.min.x, asp$scatter.data.max.x )
  y.limits <- c( asp$scatter.data.min.y, asp$scatter.data.max.y )
  x.breaks <- seq( asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step )
  y.breaks <- seq( asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step )
  x.labels <- paste0( round( x.breaks / 1e6, 1 ), "e6" )
  y.labels <- paste0( round( y.breaks / 1e6, 1 ), "e6" )

  ggplot2::ggplot( df, ggplot2::aes( x, y, color = Gate ) ) +
    scattermore::geom_scattermore(
      pointsize = asp$figure.gate.point.size, alpha = 1, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c( Retained = retained.color, Excluded = excluded.color ), drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      name = x.col, breaks = x.breaks, labels = x.labels, limits = x.limits,
      expand = ggplot2::expansion( asp$figure.gate.scale.expand )
    ) +
    ggplot2::scale_y_continuous(
      name = y.col, breaks = y.breaks, labels = y.labels, limits = y.limits,
      expand = ggplot2::expansion( asp$figure.gate.scale.expand )
    ) +
    .theme.biplot( asp )
}

## Continuous colour gradient for cosine similarity to AF, in the same blue
## (least AF-like) -> red (most AF-like) direction as the legacy RdYlBu ramp
## used in .cosine.filter.single.plot(), but as a single smooth gradient bar
## rather than discrete tranches.
.cosine.gradient.scale <- function( name = "Cosine similarity\nto AF" ) {
  rdylbu <- c( "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090",
               "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695" )
  ggplot2::scale_color_gradientn( colors = rev( rdylbu ), name = name )
}

## Cowplot triple: unstained (black, left), single-stained control coloured by
## cosine similarity to AF (middle), and the colour-bar legend split out into
## its own column (right) so both biplots keep equal width. Mirrors the
## left/right layout of gate.af.sample.plot()'s legacy AF-gating figures.
.make.cosine.filter.biplot <- function(
    unstained.mat, stained.mat, x.dim, y.dim,
    cs.vals, selected, asp,
    unstained.point.color = "black", max.points = 5e4, y.lab.suffix = "",
    cosine.point.size = NULL, selected.size.mult = 5
) {
  biexp <- .biexp.transform.legacy( asp )
  y.lab <- if ( nzchar( y.lab.suffix ) ) paste( y.dim, y.lab.suffix ) else y.dim
  stained.pointsize <- if ( !is.null( cosine.point.size ) )
    cosine.point.size else asp$figure.gate.point.size * 1.3

  # -- unstained panel (black points, no colour mapping)
  df.u <- data.frame( x = unstained.mat[ , x.dim ], y = unstained.mat[ , y.dim ] )
  if ( nrow( df.u ) > max.points ) {
    set.seed( asp$bird.seed )
    df.u <- df.u[ sample( seq_len( nrow( df.u ) ), max.points ), , drop = FALSE ]
  }
  df.u$x.trans <- biexp( df.u$x )
  df.u$y.trans <- biexp( df.u$y )

  p.unstained <- ggplot2::ggplot( df.u, ggplot2::aes( x.trans, y.trans ) ) +
    scattermore::geom_scattermore(
      pointsize = asp$figure.gate.point.size, color = unstained.point.color,
      alpha = 1, na.rm = TRUE
    ) +
    .biplot.scales( biexp, asp, x.dim, y.lab ) +
    ggplot2::labs( title = "Unstained" ) +
    .theme.biplot( asp ) +
    ggplot2::theme( aspect.ratio = 1 )

  # -- single-stained panel (coloured by cosine similarity, continuous);
  # events selected for the final profile are drawn selected.size.mult times
  # larger than the rest, in place of a gate (nothing is actually excluded
  # here -- the selection feeds forward as an average, not a subset drawn
  # from a polygon)
  df.s <- data.frame(
    x = stained.mat[ , x.dim ], y = stained.mat[ , y.dim ],
    CosineSim = cs.vals, Selected = selected
  )
  if ( nrow( df.s ) > max.points ) {
    set.seed( asp$bird.seed )
    idx  <- sample( seq_len( nrow( df.s ) ), max.points )
    df.s <- df.s[ idx, , drop = FALSE ]
  }
  df.s$x.trans <- biexp( df.s$x )
  df.s$y.trans <- biexp( df.s$y )

  p.stained <- ggplot2::ggplot(
    df.s, ggplot2::aes( x.trans, y.trans, color = CosineSim )
  ) +
    scattermore::geom_scattermore(
      data = df.s[ !df.s$Selected, , drop = FALSE ],
      pointsize = stained.pointsize, alpha = 0.3, na.rm = TRUE
    ) +
    scattermore::geom_scattermore(
      data = df.s[ df.s$Selected, , drop = FALSE ],
      pointsize = stained.pointsize * selected.size.mult, alpha = 1, na.rm = TRUE
    ) +
    .cosine.gradient.scale() +
    .biplot.scales( biexp, asp, x.dim, y.lab ) +
    ggplot2::labs( title = "Single-stained control" ) +
    .theme.biplot( asp, legend.position = "right" ) +
    ggplot2::theme( aspect.ratio = 1 )

  # extract the colour-bar legend into its own column so the two biplots
  # above stay the same width as each other
  legend.grob <- cowplot::get_legend( p.stained )
  p.stained   <- p.stained + ggplot2::theme( legend.position = "none" )

  cowplot::plot_grid(
    p.unstained, p.stained, legend.grob,
    ncol = 3, rel_widths = c( 1, 1, 0.35 )
  )
}

## KDE-smoothed 1D histogram of all singlet-gated single-stained control
## events on the empirical peak channel, with a bracket above the curve
## marking the range of the top-expressing "brightest events" candidate pool
## selected before cosine filtering.
.brightest.selection.panel <- function(
    values, selected.range, asp, x.lab,
    fill.color = "steelblue", line.color = "black"
) {
  biexp   <- .biexp.transform.legacy( asp )
  x.trans <- biexp( values )

  dens     <- stats::density( x.trans, na.rm = TRUE )
  dens.df  <- data.frame( x = dens$x, y = dens$y )
  y.max    <- max( dens.df$y, na.rm = TRUE )
  bar.y    <- y.max * 1.15
  bar.x    <- biexp( selected.range )

  bracket.df <- data.frame(
    x    = c( bar.x[ 1L ], bar.x[ 1L ], bar.x[ 2L ], bar.x[ 2L ] ),
    y    = c( bar.y * 0.95, bar.y, bar.y, bar.y * 0.95 ),
    group  = c( 1L, 1L, 1L, 1L )
  )

  breaks <- asp$ribbon.breaks
  limits <- c( asp$ribbon.plot.min, asp$expr.data.max )
  axis.labels <- sapply( breaks, function( v ) {
    if ( v == 0 ) "0" else parse( text = paste0( "10^", log10( abs( v ) ) ) )
  } )

  ggplot2::ggplot( dens.df, ggplot2::aes( x, y ) ) +
    ggplot2::geom_area( fill = fill.color, color = NA, alpha = 0.7, na.rm = TRUE ) +
    ggplot2::geom_line( color = line.color, linewidth = asp$figure.spectra.line.size, na.rm = TRUE ) +
    ggplot2::geom_path(
      data = bracket.df, ggplot2::aes( x, y, group = group ),
      color = line.color, linewidth = 0.8, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text", x = mean( bar.x ), y = bar.y * 1.12,
      label = "Brightest events selected", size = 3.2, color = line.color
    ) +
    ggplot2::scale_x_continuous(
      name = x.lab, breaks = biexp( breaks ), limits = biexp( limits ), labels = axis.labels
    ) +
    ggplot2::scale_y_continuous( name = "Density" ) +
    ggplot2::coord_cartesian( ylim = c( 0, bar.y * 1.3 ) ) +
    .theme.biplot( asp )
}


.embed.or.placeholder <- function(
    img.path, label = "kNN scatter-matched AF subtraction"
) {
  if ( !is.null( img.path ) && file.exists( img.path ) &&
       requireNamespace( "magick", quietly = TRUE ) ) {
    img <- magick::image_read( img.path )
    cowplot::ggdraw() + cowplot::draw_image( img )
  } else {
    msg <- if ( is.null( img.path ) ) {
      paste0( label, ": insufficient scatter channels available; step skipped." )
    } else {
      paste0(
        label, " was saved separately at:\n", img.path,
        "\n(install the 'magick' package to embed it inline)"
      )
    }
    cowplot::ggdraw() + cowplot::draw_label( msg, size = 10, fontface = "italic" )
  }
}

## Identify rows of a control table flagged as bead controls via an optional
## `control.type` column (values "beads"/"cells", case-insensitive), mirroring
## the column already used by define.flow.control()/check.control.file(). If
## `control.type` is absent (older-style control files), every row reports
## FALSE and behaviour is unchanged from before this column existed.
.is.beads.row <- function( ctrl.tbl ) {
  if ( !"control.type" %in% colnames( ctrl.tbl ) )
    return( rep( FALSE, nrow( ctrl.tbl ) ) )
  tolower( trimws( ctrl.tbl$control.type ) ) %in% "beads"
}

## Resolve the comparator profile shown as "<fluor> (Beads)" in the final
## spectral-trace panel. If `ctrl.tbl` has a `control.type == "beads"` row for
## `fluor` with a file-type `universal.negative` (a paired unstained bead
## sample), the profile is calculated directly from that pair of FCS files:
## mean of the positive bead events minus mean of the matched negative bead
## events, L-infinity normalised. Otherwise falls back to the cytometer's
## static spectral reference library via .load.ref.library(). Returns a
## named numeric vector over `spectral.channels` (NA where neither source has
## data for a channel), or NULL if no source is available at all.
.get.reference.profile <- function(
    fluor, ctrl.tbl, control.dir, spectral.channels, scatter.channels,
    sat.value, singlet.quantiles, asp, db.col, verbose = TRUE
) {
  beads.rows <- which( .is.beads.row( ctrl.tbl ) & ctrl.tbl$fluorophore == fluor )

  if ( length( beads.rows ) > 0 ) {
    bead.row  <- beads.rows[ 1L ]
    bead.file <- ctrl.tbl$filename[ bead.row ]
    bead.path <- file.path( control.dir, bead.file )
    bead.neg  <- .parse.unstained.source( ctrl.tbl$universal.negative[ bead.row ] )

    if ( bead.neg$type != "file" ) {
      warning(
        "Bead control for '", fluor, "' has no paired negative in 'universal.negative'; ",
        "falling back to the spectral reference library for the Beads trace.",
        call. = FALSE
      )
    } else if ( !file.exists( bead.path ) ||
                !file.exists( file.path( control.dir, bead.neg$file ) ) ) {
      warning(
        "Bead control or negative file not found for '", fluor, "'; ",
        "falling back to the spectral reference library for the Beads trace.",
        call. = FALSE
      )
    } else {
      bead.pos <- .read.fcs.clean(
        bead.path, paste0( fluor, " (beads, positive)" ),
        spectral.channels, scatter.channels, sat.value, singlet.quantiles,
        asp = asp, verbose = verbose
      )
      bead.neg.mat <- .read.fcs.clean(
        file.path( control.dir, bead.neg$file ), paste0( fluor, " (beads, negative)" ),
        spectral.channels, scatter.channels, sat.value, singlet.quantiles,
        asp = asp, verbose = verbose
      )
      spec.cols <- intersect(
        spectral.channels, intersect( colnames( bead.pos ), colnames( bead.neg.mat ) )
      )

      if ( nrow( bead.pos ) >= 2L && nrow( bead.neg.mat ) >= 2L && length( spec.cols ) > 0L ) {
        bead.profile <- colMeans( bead.pos[ , spec.cols, drop = FALSE ] ) -
          colMeans( bead.neg.mat[ , spec.cols, drop = FALSE ] )
        bead.max <- max( abs( bead.profile ) )
        if ( is.finite( bead.max ) && bead.max > 0 ) {
          out <- stats::setNames( rep( NA_real_, length( spectral.channels ) ), spectral.channels )
          out[ spec.cols ] <- bead.profile / bead.max
          return( out )
        }
      }
      warning(
        "Too few clean events in the bead control for '", fluor, "'; ",
        "falling back to the spectral reference library for the Beads trace.",
        call. = FALSE
      )
    }
  }

  ref.lib <- .load.ref.library( db.col, spectral.channels )
  if ( !is.null( ref.lib ) && fluor %in% rownames( ref.lib ) )
    return( ref.lib[ fluor, ] )

  NULL
}

## Downward-pointing arrow and label marking a named Detector factor level
## (e.g. the empirical peak channel) on a spectral.trace() plot. `y.max`
## should be the largest value actually plotted -- spectral.trace() input
## here is normalised to a maximum magnitude of 1 per row, so the default
## traces peak around 1 -- so the annotation clears them regardless of scale.
.add.peak.arrow <- function(
    plot, peak.channel, y.max, label = "Peak detector", arrow.color = "black"
) {
  if ( !is.finite( y.max ) || y.max <= 0 ) y.max <- 1
  arrow.base <- y.max * 1.35
  arrow.tip  <- y.max * 1.08

  plot +
    ggplot2::annotate(
      "segment", x = peak.channel, xend = peak.channel,
      y = arrow.base, yend = arrow.tip,
      arrow = grid::arrow( length = grid::unit( 0.12, "inches" ), type = "closed" ),
      color = arrow.color, linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text", x = peak.channel, y = arrow.base * 1.08,
      label = label, size = 3.2, color = arrow.color
    ) +
    ggplot2::expand_limits( y = arrow.base * 1.15 )
}

## Adds a coloured highlight layer of "clean positive" events on top of an
## existing biplot panel (e.g. panel A's -A vs -H singlet gate), so it's
## visible where the events ultimately taken forward for profile calculation
## sit relative to the full singlet-gated population.
.add.highlight.layer <- function( plot, x.vals, y.vals, color = "red", pointsize = 3 ) {
  plot +
    scattermore::geom_scattermore(
      data = data.frame( x = x.vals, y = y.vals ),
      mapping = ggplot2::aes( x, y ), inherit.aes = FALSE,
      pointsize = pointsize, color = color, alpha = 1, na.rm = TRUE
    )
}


# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' @title Plot Automated Spectra Extraction Pipeline Steps
#'
#' @description
#' Builds a manuscript-ready, multi-panel figure illustrating the internal
#' steps of [get.spectra.automated()] for one or more single-stained control
#' samples: singlet gating (linear scale), AF-orthogonalisation peak-finding,
#' brightest-event candidate selection (KDE-smoothed histogram), the
#' cosine-similarity filter (biexponential biplot, continuous colour
#' gradient), and the kNN scatter-matched AF subtraction. Uses the same
#' building blocks as the rest of AutoSpectral ([create.biplot()]-style
#' biexponential biplots, [spectral.trace()], [scatter.match.plot()]) so
#' panel styling matches the package's other figures.
#'
#' @param control.dir Character. Path to the directory containing the
#'   single-stained control FCS files.
#' @param control.def.file Character. Path to (or filename of) the control
#'   definition CSV, as used by [get.spectra.automated()]. If the file has a
#'   `control.type` column (as used by [define.flow.control()] /
#'   [check.control.file()]), a row with `control.type == "beads"` sharing a
#'   fluorophore name with a "cells" row is treated as a paired bead control:
#'   its `filename` is the positive bead sample, and its `universal.negative`
#'   must point at the matching unstained bead file. When present, this pair
#'   supplies the "<Fluorophore> (Beads)" trace in panel F instead of the
#'   static spectral reference library.
#' @param asp The AutoSpectral parameter list from `get.autospectral.param()`.
#' @param fluorophores Character vector of fluorophore name(s) (as they
#'   appear in the `fluorophore` column of the control file) to illustrate.
#'   Default `NULL` illustrates the first external-negative fluorophore found.
#'   Fluorophores using internal-negative mode (no `universal.negative` entry)
#'   are skipped with a warning, since orthogonalisation / cosine-filter / kNN
#'   steps do not apply to them.
#' @param n.candidates Integer, default `1000`. Number of top-expressing
#'   candidate events shown as the "brightest events" selection (panel C) and
#'   carried into the cosine-filter biplot (panel D), matching the
#'   `n.candidates` argument of [get.spectra.automated()].
#' @param n.spectral Integer, default `200`. Number of the above candidates
#'   that would be retained after cosine filtering; these are drawn at full
#'   opacity, the remainder at reduced opacity, matching
#'   [get.spectra.automated()]'s `n.spectral` argument.
#' @param k.neighbors Integer, default `2`. Number of nearest neighbours used
#'   for the kNN scatter-match panel (panel E), matching
#'   [get.spectra.automated()]'s `k.neighbors` argument.
#' @param singlet.quantiles Numeric, default `c(0.85, 0.975)`. Quantile
#'   thresholds for the two-stage FSC/SSC singlet discrimination shown in
#'   panel A, matching [get.spectra.automated()].
#' @param trace.colors Named character vector of colours for the three traces
#'   in panel B. Names must be `"AF"`, `"Pre-orthogonalization"` and
#'   `"Post-orthogonalization"`.
#' @param cells.trace.color,beads.trace.color,af.trace.color Colours for the
#'   three traces in panel F: the calculated per-cell profile ("Cells"), the
#'   reference profile ("Beads"), and the matched-unstained AF trace.
#'   Defaults `"#D95F02"` / `"#377EB8"` / `"grey40"`.
#' @param singlet.retained.color,singlet.excluded.color Colours for retained
#'   and excluded events in panel A. Defaults `"black"` / `"grey75"`.
#' @param clean.positive.color Colour for the highlighted "clean positive"
#'   events overlaid on panel A. Default `"red"`.
#' @param clean.positive.point.size Numeric or `NULL` (default). Point size
#'   for the panel A highlight. If `NULL`, defaults to
#'   `asp$figure.gate.point.size * 1.5`.
#' @param brightest.fill.color,brightest.line.color Fill and line colours for
#'   the KDE-smoothed histogram in panel C. Defaults `"steelblue"` / `"black"`.
#' @param brightest.panel.width,brightest.panel.height Numeric or `NULL`
#'   (default). Width/height (inches) of panel C. If `NULL`, defaults to
#'   `panel.width * 2` / `panel.height`. When narrower than the overall
#'   composite width, the panel is centred with blank padding on either side.
#' @param unstained.point.color Colour for unstained events in panel D.
#'   Default `"black"`.
#' @param cosine.point.size Numeric or `NULL` (default). Point size for the
#'   non-selected single-stained control events in panel D (these are often a
#'   small, pre-selected "brightest events" pool, so the default point size
#'   can be too small to see clearly). If `NULL`, defaults to
#'   `asp$figure.gate.point.size * 1.3`.
#' @param selected.size.mult Numeric, default `5`. Events selected for the
#'   final profile are drawn in panel D at `cosine.point.size` times this
#'   multiple, in place of a gate polygon (nothing is actually excluded from
#'   the plot -- the selection feeds forward as an average).
#' @param max.points Integer. Maximum events plotted per panel (randomly
#'   downsampled beyond this for speed). Default `5e4`.
#' @param panel.width,panel.height Numeric. Width/height (inches) used per
#'   sub-panel when sizing the saved composite figure. Defaults `4` and `4`.
#' @param composite.width,composite.height Numeric or `NULL` (default).
#'   Override the overall saved figure dimensions (inches); if `NULL`, these
#'   are computed from `panel.width` / `panel.height` (and
#'   `brightest.panel.height`, if given).
#' @param output.dir Character or `NULL` (default). Directory to save the
#'   composite figure(s) and the kNN scatter-match panel. Defaults to the
#'   current working directory.
#' @param save Logical, default `TRUE`. Whether to save the composite figure
#'   for each fluorophore to `output.dir`.
#' @param file.type Character string, one of `"jpg"` (default), `"tiff"`,
#'   `"png"`, or `"pdf"`.
#' @param verbose Logical, default `TRUE`. Print progress messages.
#'
#' @return Invisibly, a named list (one entry per fluorophore) each
#'   containing the individual panel ggplot objects, the assembled
#'   `composite` cowplot object, and the empirical peak / AF-collision
#'   channels used.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_manual scale_alpha_manual theme_bw theme
#' @importFrom ggplot2 element_line element_text element_rect element_blank
#' @importFrom ggplot2 margin labs ggsave geom_area geom_line geom_path
#' @importFrom ggplot2 annotate coord_cartesian scale_color_gradientn
#' @importFrom scattermore geom_scattermore
#' @importFrom cowplot plot_grid ggdraw draw_image draw_label get_legend
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom ragg agg_jpeg agg_tiff agg_png
#' @importFrom stats density
#'
#' @seealso [get.spectra.automated()], [create.biplot()], [spectral.trace()],
#'   [scatter.match.plot()]
#'
#' @export

spectra.automated.steps.plot <- function(
    control.dir,
    control.def.file,
    asp,
    fluorophores           = NULL,
    n.candidates           = 1000L,
    n.spectral             = 200L,
    k.neighbors             = 2L,
    singlet.quantiles       = c( 0.85, 0.975 ),
    trace.colors            = c(
      "AF"                     = "grey40",
      "Pre-orthogonalization"  = "#1B9E77",
      "Post-orthogonalization" = "#D95F02"
    ),
    cells.trace.color       = "#D95F02",
    beads.trace.color       = "#377EB8",
    af.trace.color          = "grey40",
    singlet.retained.color  = "black",
    singlet.excluded.color  = "grey75",
    clean.positive.color      = "red",
    clean.positive.point.size = NULL,
    brightest.fill.color    = "steelblue",
    brightest.line.color    = "black",
    brightest.panel.width   = NULL,
    brightest.panel.height  = NULL,
    unstained.point.color   = "black",
    cosine.point.size       = NULL,
    selected.size.mult      = 5,
    max.points              = 5e4,
    panel.width             = 4,
    panel.height            = 4,
    composite.width         = NULL,
    composite.height        = NULL,
    output.dir              = NULL,
    save                    = TRUE,
    file.type               = "jpg",
    verbose                 = TRUE
) {

  # -- 0. Validate inputs
  if ( !dir.exists( control.dir ) )
    stop( "control.dir does not exist: ", control.dir, call. = FALSE )

  if ( is.null( output.dir ) ) output.dir <- getwd()
  if ( save && !dir.exists( output.dir ) )
    dir.create( output.dir, recursive = TRUE )

  file.type <- tolower( file.type )
  if ( file.type == "jpeg" ) file.type <- "jpg"
  file.type <- match.arg( file.type, c( "jpg", "tiff", "png", "pdf" ) )
  plot.device <- switch( file.type,
                         jpg  = ragg::agg_jpeg,
                         tiff = ragg::agg_tiff,
                         png  = ragg::agg_png,
                         pdf  = grDevices::pdf
  )

  # -- 1. Read control table (mirrors get.spectra.automated())
  ctrl.path <- if ( file.exists( control.def.file ) ) {
    control.def.file
  } else {
    file.path( control.dir, control.def.file )
  }
  ctrl.tbl <- utils::read.csv(
    ctrl.path, stringsAsFactors = FALSE, strip.white = TRUE
  )

  uneg.bool.all <- suppressWarnings( as.logical( ctrl.tbl$universal.negative ) )
  is.neg.row    <- !is.na( uneg.bool.all ) & uneg.bool.all
  is.beads.row  <- .is.beads.row( ctrl.tbl )

  fluor.rows.all <- which(
    !is.neg.row & !is.beads.row &
      !grepl( "negative", ctrl.tbl$fluorophore, ignore.case = TRUE ) &
      !grepl( "^AF$",     ctrl.tbl$fluorophore, ignore.case = TRUE )
  )

  if ( length( fluor.rows.all ) == 0 )
    stop( "No fluorophore rows found in ", control.def.file, call. = FALSE )

  all.fluor.names <- ctrl.tbl$fluorophore[ fluor.rows.all ]

  if ( is.null( fluorophores ) ) {
    fluorophores <- all.fluor.names[ 1L ]
    if ( verbose )
      message( "No `fluorophores` specified; illustrating: ", fluorophores )
  }

  missing.fluor <- setdiff( fluorophores, all.fluor.names )
  if ( length( missing.fluor ) > 0 )
    stop(
      "Fluorophore(s) not found in control file: ",
      paste( missing.fluor, collapse = ", " ), call. = FALSE
    )

  unstained.sources.all <- lapply( fluor.rows.all, function( i ) {
    src <- .parse.unstained.source( ctrl.tbl$universal.negative[ i ] )
    if ( src$type == "global.true" ) {
      true.rows <- which( !is.na( uneg.bool.all ) & uneg.bool.all )
      if ( length( true.rows ) > 0 )
        return( list( type = "file", file = ctrl.tbl$filename[ true.rows[ 1L ] ] ) )
      return( list( type = "internal", file = NULL ) )
    }
    src
  } )
  names( unstained.sources.all ) <- all.fluor.names

  # -- 2. Resolve channels
  first.fcs.path    <- file.path( control.dir, ctrl.tbl$filename[ fluor.rows.all[ 1L ] ] )
  spectral.channels <- .derive.spectral.channels( first.fcs.path, asp )
  scatter.channels  <- read.scatter.parameter( asp )
  sat.value         <- if ( !is.null( asp$expr.data.max ) ) asp$expr.data.max else Inf

  db.col <- .cytometer.to.db.col( asp$cytometer )
  clean.positive.point.size.use <- if ( !is.null( clean.positive.point.size ) )
    clean.positive.point.size else asp$figure.gate.point.size * 1.5

  # -- overall composite sizing (shared across fluorophores)
  composite.width.use  <- if ( !is.null( composite.width ) )  composite.width  else panel.width * 2
  brightest.width.use  <- if ( !is.null( brightest.panel.width ) )  brightest.panel.width  else composite.width.use
  brightest.height.use <- if ( !is.null( brightest.panel.height ) ) brightest.panel.height else panel.height

  row.heights <- c( panel.height, panel.height, brightest.height.use, panel.height, panel.height, panel.height )
  composite.height.use <- if ( !is.null( composite.height ) ) composite.height else sum( row.heights )

  results <- list()

  # -- 3. Loop over requested fluorophores
  for ( fluor in fluorophores ) {

    if ( verbose )
      message( "\033[34m-- Building pipeline-step figure for ", fluor, " --\033[0m" )

    row.i      <- fluor.rows.all[ match( fluor, all.fluor.names ) ]
    fluor.file <- ctrl.tbl$filename[ row.i ]
    src.i      <- unstained.sources.all[[ fluor ]]

    if ( src.i$type != "file" ) {
      warning(
        "'", fluor, "' uses internal-negative mode; orthogonalization / ",
        "cosine-filter / kNN steps do not apply. Skipping.", call. = FALSE
      )
      next
    }

    unstained.file <- src.i$file
    unstained.path <- file.path( control.dir, unstained.file )
    if ( !file.exists( unstained.path ) ) {
      warning( "Unstained file not found for '", fluor, "': ", unstained.file, call. = FALSE )
      next
    }

    ref.profile.i <- .get.reference.profile(
      fluor, ctrl.tbl, control.dir, spectral.channels, scatter.channels,
      sat.value, singlet.quantiles, asp, db.col, verbose = FALSE
    )

    # -- A. Read stained control raw (A + H columns) and flag singlets
    fcs.path.i <- file.path( control.dir, fluor.file )
    raw.i <- .read.fcs.raw.singlet.flags(
      fcs.path.i, spectral.channels, scatter.channels, sat.value, singlet.quantiles
    )

    singlet.panels <- list()
    if ( all( c( raw.i$fsc.a, raw.i$fsc.h ) %in% colnames( raw.i$mat ) ) )
      singlet.panels$fsc <- .singlet.gate.panel(
        raw.i$mat, raw.i$fsc.a, raw.i$fsc.h, raw.i$retained, asp,
        singlet.retained.color, singlet.excluded.color, max.points
      )
    if ( !is.na( raw.i$ssc.a ) &&
         all( c( raw.i$ssc.a, raw.i$ssc.h ) %in% colnames( raw.i$mat ) ) )
      singlet.panels$ssc <- .singlet.gate.panel(
        raw.i$mat, raw.i$ssc.a, raw.i$ssc.h, raw.i$retained, asp,
        singlet.retained.color, singlet.excluded.color, max.points
      )

    # singlet.plot is assembled later, once the panel D cosine filter has
    # identified the final "clean positive" events, so they can be
    # highlighted here in panel A

    spec.data.i <- raw.i$mat[ raw.i$retained, spectral.channels, drop = FALSE ]

    if ( nrow( spec.data.i ) < 2L ) {
      warning( "Too few retained singlet events for '", fluor, "'. Skipping.", call. = FALSE )
      next
    }

    # -- B. Unstained AF reference (fully cleaned, as in the main pipeline)
    unstained.mat <- .read.fcs.clean(
      unstained.path, paste0( "Unstained (", unstained.file, ")" ),
      spectral.channels, scatter.channels, sat.value, singlet.quantiles,
      asp = asp, verbose = FALSE
    )
    af.spectral.i <- unstained.mat[ , intersect( spectral.channels, colnames( unstained.mat ) ), drop = FALSE ]
    af.mean.i     <- colMeans( af.spectral.i )
    af.median.i   <- apply( af.spectral.i, 2, stats::median )

    # -- B continued: AF orthogonalisation -> empirical peak
    v.unit.i   <- af.mean.i / ( sqrt( sum( af.mean.i^2 ) ) + 1e-9 )
    proj.i     <- spec.data.i %*% v.unit.i
    mat.orth.i <- spec.data.i - proj.i %*% t( v.unit.i )
    empirical.peak <- names( which.max( colMeans( mat.orth.i ) ) )

    pre.mean  <- colMeans( spec.data.i )
    post.mean <- colMeans( mat.orth.i )

    trace.mat <- rbind(
      "AF"                     = af.mean.i / max( abs( af.mean.i ),  1e-9 ),
      "Pre-orthogonalization"  = pre.mean  / max( abs( pre.mean ),   1e-9 ),
      "Post-orthogonalization" = post.mean / max( abs( post.mean ),  1e-9 )
    )
    colnames( trace.mat ) <- spectral.channels

    trace.plot <- suppressMessages(
      spectral.trace(
        trace.mat, asp,
        title        = paste0( fluor, "_peak_finding_traces" ),
        split.lasers = FALSE, save = FALSE,
        figure.spectra.line.size  = asp$figure.spectra.line.size,
        figure.spectra.point.size = asp$figure.spectra.point.size
      ) +
        ggplot2::scale_color_manual( values = trace.colors, name = NULL ) +
        ggplot2::labs( title = paste( "Peak finding:", fluor ) )
    )
    trace.plot <- .add.peak.arrow( trace.plot, empirical.peak, max( trace.mat, na.rm = TRUE ) )

    # -- C. Brightest-event candidate selection
    n.cand.actual <- min( n.candidates, nrow( spec.data.i ) )
    peak.col      <- if ( empirical.peak %in% colnames( spec.data.i ) )
      empirical.peak else colnames( spec.data.i )[ 1L ]
    i.top   <- order( spec.data.i[ , peak.col ], decreasing = TRUE )[ seq_len( n.cand.actual ) ]
    top.mat <- spec.data.i[ i.top, , drop = FALSE ]

    brightest.panel <- .brightest.selection.panel(
      values         = spec.data.i[ , peak.col ],
      selected.range = range( top.mat[ , peak.col ] ),
      asp            = asp, x.lab = peak.col,
      fill.color     = brightest.fill.color,
      line.color     = brightest.line.color
    )

    if ( brightest.width.use < composite.width.use ) {
      pad.frac  <- ( composite.width.use - brightest.width.use ) / 2 / composite.width.use
      main.frac <- brightest.width.use / composite.width.use
      brightest.row <- cowplot::plot_grid(
        NULL, brightest.panel, NULL, ncol = 3,
        rel_widths = c( pad.frac, main.frac, pad.frac )
      )
    } else {
      brightest.row <- brightest.panel
    }

    # -- D. Cosine-similarity filter
    cs.vals <- .cosine.sim.rows( top.mat, af.median.i[ colnames( top.mat ) ] )

    n.spec.actual  <- min( n.spectral, length( cs.vals ) )
    order.cs       <- order( cs.vals )
    selected.local <- rep( FALSE, length( cs.vals ) )
    selected.local[ order.cs[ seq_len( n.spec.actual ) ] ] <- TRUE

    candidate.channels <- setdiff( colnames( af.spectral.i ), empirical.peak )
    y.channel.peak <- names( which.max( af.median.i[ candidate.channels ] ) )

    cosine.panel <- .make.cosine.filter.biplot(
      unstained.mat = af.spectral.i, stained.mat = top.mat,
      x.dim = empirical.peak, y.dim = y.channel.peak,
      cs.vals = cs.vals, selected = selected.local, asp = asp,
      unstained.point.color = unstained.point.color, max.points = max.points,
      cosine.point.size = cosine.point.size, selected.size.mult = selected.size.mult,
      y.lab.suffix = "(unstained peak channel)"
    )

    # -- E. kNN scatter-matching (reuse scatter.match.plot() as-is)
    scatter.match.file <- NULL
    subtraction.plot    <- NULL
    common.scat <- intersect( scatter.channels, colnames( unstained.mat ) )

    if ( length( common.scat ) >= 2L &&
         all( common.scat %in% colnames( raw.i$mat ) ) ) {

      i.spectral <- i.top[ order.cs[ seq_len( n.spec.actual ) ] ]
      clean.mat  <- raw.i$mat[ raw.i$retained, , drop = FALSE ][ i.spectral, , drop = FALSE ]

      if ( !is.null( singlet.panels$fsc ) )
        singlet.panels$fsc <- .add.highlight.layer(
          singlet.panels$fsc, clean.mat[ , raw.i$fsc.a ], clean.mat[ , raw.i$fsc.h ],
          color = clean.positive.color, pointsize = clean.positive.point.size.use
        )
      if ( !is.null( singlet.panels$ssc ) && !is.na( raw.i$ssc.a ) )
        singlet.panels$ssc <- .add.highlight.layer(
          singlet.panels$ssc, clean.mat[ , raw.i$ssc.a ], clean.mat[ , raw.i$ssc.h ],
          color = clean.positive.color, pointsize = clean.positive.point.size.use
        )

      fluor.scatter <- raw.i$mat[ raw.i$retained, , drop = FALSE ][ i.spectral, common.scat, drop = FALSE ]
      unstained.scat.mat <- unstained.mat[ , common.scat, drop = FALSE ]

      knn.idx <- tryCatch(
        FNN::knnx.index(
          data  = as.matrix( unstained.scat.mat ),
          query = as.matrix( fluor.scatter ),
          k     = k.neighbors
        ),
        error = function( e ) NULL
      )

      if ( !is.null( knn.idx ) ) {
        matched.scatter <- unstained.scat.mat[ unique( as.vector( knn.idx ) ), , drop = FALSE ]

        scatter.match.filename <- paste( fluor, asp$scatter.match.plot.filename, sep = "_" )
        scatter.match.dir <- if ( !is.null( asp$figure.scatter.dir.base ) )
          asp$figure.scatter.dir.base else output.dir
        scatter.match.file <- file.path( scatter.match.dir, scatter.match.filename )

        tryCatch(
          scatter.match.plot(
            pos.expr.data = fluor.scatter,
            neg.expr.data = matched.scatter,
            fluor.name    = fluor,
            scatter.param = common.scat,
            asp           = asp
          ),
          error = function( e ) {
            message( "\033[31m  scatter.match.plot error [", fluor, "]: ", e$message, "\033[0m" )
          }
        )

        # -- F. Per-event background subtraction using the kNN-matched
        # unstained events (unstained.scat.mat and af.spectral.i share row
        # indexing, both subset from the same cleaned unstained.mat)
        positive.spectral   <- spec.data.i[ i.spectral, , drop = FALSE ]
        background.spectral <- t( sapply( seq_len( nrow( knn.idx ) ), function( j ) {
          colMeans( af.spectral.i[ knn.idx[ j, ], , drop = FALSE ] )
        } ) )
        colnames( background.spectral ) <- colnames( af.spectral.i )
        subtracted.spectral <- positive.spectral - background.spectral

        bg.mean       <- colMeans( background.spectral )
        post.sub.mean <- colMeans( subtracted.spectral )

        # -- Final spectral profile comparison: the calculated per-cell
        # profile ("Cells"), the reference profile for the same fluorophore
        # ("Beads" -- a paired bead control if provided, otherwise the
        # static spectral reference library), and the matched-unstained AF
        cells.label <- paste0( fluor, " (Cells)" )
        beads.label <- paste0( fluor, " (Beads)" )
        af.label    <- "Unstained (Autofluorescence)"

        beads.trace <- if ( !is.null( ref.profile.i ) )
          ref.profile.i[ spectral.channels ] else rep( 0, length( spectral.channels ) )
        beads.trace[ !is.finite( beads.trace ) ] <- 0
        if ( is.null( ref.profile.i ) )
          warning(
            "No reference profile (bead control or spectral library) found for '",
            fluor, "'; plotting a flat zero line for '", beads.label, "'.", call. = FALSE
          )

        final.spectrum <- post.sub.mean / max( abs( post.sub.mean ), 1e-9 )
        final.spectrum[ final.spectrum < 0 ] <- 0

        subtraction.mat <- rbind(
          final.spectrum,
          beads.trace,
          bg.mean / max( abs( bg.mean ), 1e-9 )
        )
        rownames( subtraction.mat ) <- c( cells.label, beads.label, af.label )
        colnames( subtraction.mat ) <- spectral.channels

        final.trace.colors <- stats::setNames(
          c( cells.trace.color, beads.trace.color, af.trace.color ),
          c( cells.label, beads.label, af.label )
        )

        subtraction.plot <- suppressMessages(
          spectral.trace(
            subtraction.mat, asp,
            title        = paste0( fluor, "_final_spectral_profile" ),
            split.lasers = FALSE, save = FALSE,
            figure.spectra.line.size  = asp$figure.spectra.line.size,
            figure.spectra.point.size = asp$figure.spectra.point.size
          ) +
            ggplot2::scale_color_manual( values = final.trace.colors, name = NULL ) +
            ggplot2::labs( title = paste( "Final spectral profile:", fluor ) )
        )
      }
    }

    scatter.match.panel <- .embed.or.placeholder( scatter.match.file )
    if ( is.null( subtraction.plot ) )
      subtraction.plot <- .embed.or.placeholder( NULL, label = "Final spectral profile" )

    singlet.plot <- cowplot::plot_grid( plotlist = singlet.panels, nrow = 1 )

    # -- Assemble composite figure
    composite <- cowplot::plot_grid(
      cowplot::plot_grid( singlet.plot,         labels = "A" ),
      cowplot::plot_grid( trace.plot,           labels = "B" ),
      cowplot::plot_grid( brightest.row,        labels = "C" ),
      cowplot::plot_grid( cosine.panel,         labels = "D" ),
      cowplot::plot_grid( scatter.match.panel,  labels = "E" ),
      cowplot::plot_grid( subtraction.plot,     labels = "F" ),
      ncol = 1, rel_heights = row.heights
    )

    if ( save ) {
      out.file <- file.path(
        output.dir, sprintf( "%s_pipeline_steps.%s", fluor, file.type )
      )
      ggplot2::ggsave(
        out.file, plot = composite, device = plot.device,
        width = composite.width.use, height = composite.height.use, limitsize = FALSE
      )
      if ( verbose ) message( "\033[32m  Saved: ", out.file, "\033[0m" )
    }

    results[[ fluor ]] <- list(
      singlet.plot         = singlet.plot,
      trace.plot           = trace.plot,
      brightest.panel      = brightest.panel,
      cosine.panel         = cosine.panel,
      scatter.match.file   = scatter.match.file,
      subtraction.plot     = subtraction.plot,
      composite            = composite,
      empirical.peak       = empirical.peak,
      y.channel.peak       = y.channel.peak,
      reference.profile    = ref.profile.i
    )
  }

  invisible( results )
}
