# plot_spectra_standard_workflow.R
#
# Manuscript figure helper: illustrates a "standard" (e.g. SpectroFlo-style)
# manual single-stained-control workflow, for comparison against
# spectra.automated.steps.plot(). Depends on private helpers defined in
# plot_spectra_automated_steps.R (.biexp.transform.legacy(), .theme.biplot(),
# .biplot.scales(), .cosine.gradient.scale(), .cosine.sim.rows(),
# .derive.spectral.channels(), .embed.or.placeholder(), .add.highlight.layer(),
# .read.fcs.raw.singlet.flags()) -- both files must be loaded into the
# package together (e.g. both in R/, or both source()'d).
#
# Panels produced per fluorophore:
#   A. Octagon gate on FSC-A vs SSC-A (linear scale), pseudocolour density as
#      in create.biplot(), positioned on the highest-density region (typically
#      lymphocytes) via a 2D KDE peak search; red = the "true positive" events
#      for this control, identified independently of the standard workflow's
#      own gating/selection via `ground.truth.method` ("automated": replicates
#      get.spectra.automated()'s candidate + cosine-to-AF filter; "legacy":
#      the AF-removal gate.population.idx from clean.controls(); "none": the
#      standard workflow's own top n.highlight cosine-similarity selection,
#      kept only for reference since it is not an independent ground truth).
#      Highlighted events are plotted at their true FSC/SSC coordinates even
#      if they fall outside the octagon gate boundary -- that mismatch is
#      itself diagnostic of a bad gate.
#   B. Brightest/negative event selection: KDE-smoothed 1D histogram on the
#      fluorophore's peak channel (within the octagon gate), with brackets
#      marking the top n.bright.events "positive" events and the
#      negative.quantile.min-to-negative.quantile "negative" events; the
#      x-axis lower limit is set from the x.min.quantile quantile
#   C. Biplots for the negative fraction (left, plain black points) and
#      positive fraction (right, coloured by cosine similarity to the
#      negative fraction's median profile, boxed to indicate that all
#      positive-fraction events are carried forward): peak channel (x) vs.
#      peak AF channel (y)
#   D. Final spectral profile comparison: the population-level background-
#      subtracted signature ("<Fluorophore> (Cells)"), the reference profile
#      for the same fluorophore ("<Fluorophore> (Beads)" -- from a paired
#      bead control in control.def.file if provided, otherwise the
#      cytometer's static spectral reference library), and the negative-
#      fraction trace ("Unstained (Autofluorescence)") (spectral.trace())

# ---------------------------------------------------------------------------
# Private helpers specific to the standard workflow
# ---------------------------------------------------------------------------

## Look up a fluorophore's nominal peak channel from fluorophore_database.csv
## for the current cytometer, mirroring the column-selection logic in
## create.control.file(). CytoStellar channels are stored bare in the
## database (no -A/-H suffix); full suffix resolution lives in
## .resolve.cytostellar.suffix() (check_channels.R) and is not replicated
## here -- a plain "-A" is appended if the looked-up value lacks one.
.lookup.fluorophore.channel <- function( fluor.name, asp, fluorophore.database ) {
  col <- if ( asp$cytometer == "Aurora" ) {
    if ( !is.null( asp$cytometer.version ) && identical( asp$cytometer.version, "NL" ) )
      "channel.NL" else "channel.Aurora"
  } else if ( asp$cytometer == "ID7000" ) {
    "channel.ID7000"
  } else if ( asp$cytometer %in% c( "FACSDiscover A8", "FACSDiscover S8", "FACSDiscover" ) ) {
    "channel.s8"
  } else if ( asp$cytometer == "Opteon" ) {
    "channel.opteon"
  } else if ( asp$cytometer == "Mosaic" ) {
    "channel.mosaic"
  } else if ( asp$cytometer == "Xenith" ) {
    "channel.xenith"
  } else if ( asp$cytometer == "Symphony" ) {
    "channel.A5SE"
  } else if ( asp$cytometer == "CytoStellar" ) {
    "channel.cytostellar"
  } else {
    stop( "Unsupported cytometer: ", asp$cytometer, call. = FALSE )
  }

  if ( !col %in% colnames( fluorophore.database ) )
    stop( "Column '", col, "' not found in fluorophore.database.", call. = FALSE )

  detectors <- stats::setNames(
    fluorophore.database[[ col ]], fluorophore.database$fluorophore
  )
  channel <- detectors[[ fluor.name ]]

  if ( is.null( channel ) || is.na( channel ) || !nzchar( channel ) )
    stop(
      "No peak channel found for '", fluor.name, "' on ", asp$cytometer, ".",
      call. = FALSE
    )

  if ( asp$cytometer == "CytoStellar" && !grepl( "-[AH]$", channel ) )
    channel <- paste0( channel, "-A" )

  channel
}

## Simplified 2D KDE density-peak search, in the spirit of do.gate()'s
## master-density step but without the bounding/exclusion/Voronoi machinery --
## intended for finding a single dominant population (e.g. lymphocytes) to
## centre an octagon gate on.
.find.density.peak.2d <- function( x, y, n.grid = 100L ) {
  bw <- c( MASS::bandwidth.nrd( x ), MASS::bandwidth.nrd( y ) )
  bw[ bw <= 0 | !is.finite( bw ) ] <- 0.1
  dens <- MASS::kde2d( x, y, h = bw, n = n.grid )
  idx  <- which( dens$z == max( dens$z ), arr.ind = TRUE )[ 1L, ]
  list( x = dens$x[ idx[ 1L ] ], y = dens$y[ idx[ 2L ] ] )
}

## Octagon vertices (rectangle with cut corners) centred on (cx, cy).
.octagon.polygon <- function( cx, cy, half.x, half.y, cut.frac = 1 / 3 ) {
  cx.cut <- half.x * cut.frac
  cy.cut <- half.y * cut.frac
  data.frame(
    x = cx + c( -half.x + cx.cut,  half.x - cx.cut,  half.x,           half.x,
                half.x - cx.cut, -half.x + cx.cut, -half.x,          -half.x ),
    y = cy + c( -half.y,          -half.y,          -half.y + cy.cut,  half.y - cy.cut,
                half.y,           half.y,           half.y - cy.cut, -half.y + cy.cut )
  )
}

## Point-in-polygon test (standard PNPOLY ray-casting), vectorised over
## query points against a single fixed polygon. Dependency-free.
.point.in.polygon <- function( x, y, poly.x, poly.y ) {
  n <- length( poly.x )
  inside <- rep( FALSE, length( x ) )
  j <- n
  for ( i in seq_len( n ) ) {
    cond <- ( ( poly.y[ i ] > y ) != ( poly.y[ j ] > y ) ) &
      ( x < ( poly.x[ j ] - poly.x[ i ] ) * ( y - poly.y[ i ] ) /
          ( poly.y[ j ] - poly.y[ i ] ) + poly.x[ i ] )
    inside[ cond ] <- xor( inside[ cond ], TRUE )
    j <- i
  }
  inside
}

## Replicates get.spectra.automated()'s AF-reference resolution (external
## unstained file if available, else the lower-25%-by-peak-channel internal
## negative) and final event-selection step (top n.candidates events by
## peak-channel brightness, then the n.spectral events with lowest cosine
## similarity to that AF reference) on a singlet-/saturation-gated read of
## fcs.path.i. Used as an independent "ground truth" positive population for
## ground.truth.method = "automated" in spectra.standard.workflow.plot(), so
## panel A's highlight does not depend on the standard workflow's own
## selection. Returns a matrix (scatter + spectral columns) of the selected
## events, or NULL if too few events survive singlet gating.
.identify.true.positives.automated <- function(
    fcs.path.i, spectral.channels, scatter.channels, sat.value, singlet.quantiles,
    peak.channel, unstained.src.i, control.dir, n.candidates, n.spectral
) {
  raw <- .read.fcs.raw.singlet.flags(
    fcs.path.i, spectral.channels, scatter.channels, sat.value, singlet.quantiles
  )
  stained.mat <- raw$mat[ raw$retained, , drop = FALSE ]
  if ( nrow( stained.mat ) < 10L ) return( NULL )

  af.median <- if ( unstained.src.i$type == "file" ) {
    af.path <- file.path( control.dir, unstained.src.i$file )
    if ( !file.exists( af.path ) ) {
      NULL
    } else {
      af.raw  <- .read.fcs.raw.singlet.flags(
        af.path, spectral.channels, scatter.channels, sat.value, singlet.quantiles
      )
      af.mat  <- af.raw$mat[ af.raw$retained, , drop = FALSE ]
      af.spec <- intersect( spectral.channels, colnames( af.mat ) )
      if ( nrow( af.mat ) < 2L || length( af.spec ) == 0L ) {
        NULL
      } else {
        apply( af.mat[ , af.spec, drop = FALSE ], 2, stats::median )
      }
    }
  } else {
    NULL
  }

  if ( is.null( af.median ) ) {
    # internal negative fallback, mirroring get.spectra.automated(): the
    # lower 25% of events by peak-channel value
    peak.vals.all <- if ( peak.channel %in% colnames( stained.mat ) )
      stained.mat[ , peak.channel ] else rowMeans( stained.mat[ , spectral.channels, drop = FALSE ] )
    n.neg          <- max( 2L, floor( nrow( stained.mat ) * 0.25 ) )
    i.internal.neg <- order( peak.vals.all )[ seq_len( n.neg ) ]
    af.median      <- apply(
      stained.mat[ i.internal.neg, spectral.channels, drop = FALSE ], 2, stats::median
    )
  }

  n.cand.actual <- min( n.candidates, nrow( stained.mat ) )
  i.top   <- order( stained.mat[ , peak.channel ], decreasing = TRUE )[ seq_len( n.cand.actual ) ]
  top.mat <- stained.mat[ i.top, spectral.channels, drop = FALSE ]
  cs.vals <- .cosine.sim.rows( as.matrix( top.mat ), af.median )

  n.spec.actual <- min( n.spectral, length( i.top ) )
  i.spectral    <- i.top[ order( cs.vals )[ seq_len( n.spec.actual ) ] ]

  stained.mat[ i.spectral, , drop = FALSE ]
}

## Octagon-gated FSC-A vs SSC-A panel: linear scale, pseudocolour density
## exactly as in create.biplot() (black scattermore points + stat_density_2d
## polygon fill), with the octagon gate boundary overlaid.
.octagon.gate.panel <- function(
    x.vals, y.vals, poly, asp, x.lab, y.lab,
    gate.color = "darkgoldenrod1", density.palette = "rainbow", max.points = 5e4
) {
  df <- data.frame( x = x.vals, y = y.vals )
  if ( nrow( df ) > max.points ) {
    set.seed( asp$bird.seed )
    df <- df[ sample( seq_len( nrow( df ) ), max.points ), , drop = FALSE ]
  }

  x.limits <- c( asp$scatter.data.min.x, asp$scatter.data.max.x )
  y.limits <- c( asp$scatter.data.min.y, asp$scatter.data.max.y )
  x.breaks <- seq( asp$scatter.data.min.x, asp$scatter.data.max.x, asp$data.step )
  y.breaks <- seq( asp$scatter.data.min.y, asp$scatter.data.max.y, asp$data.step )
  x.labels <- paste0( round( x.breaks / 1e6, 1 ), "e6" )
  y.labels <- paste0( round( y.breaks / 1e6, 1 ), "e6" )

  viridis.colors <- c( "magma", "inferno", "plasma", "viridis",
                       "cividis", "rocket", "mako", "turbo" )

  p <- ggplot2::ggplot( df, ggplot2::aes( x, y ) ) +
    scattermore::geom_scattermore(
      pointsize = asp$figure.gate.point.size, color = "black", alpha = 1, na.rm = TRUE
    ) +
    ggplot2::stat_density_2d(
      ggplot2::aes( fill = ggplot2::after_stat( level ) ), geom = "polygon", na.rm = TRUE
    )

  p <- if ( density.palette %in% viridis.colors ) {
    p + ggplot2::scale_fill_viridis_c( option = density.palette )
  } else {
    p + ggplot2::scale_fill_gradientn( colors = asp$density.palette.base.color )
  }

  p +
    ggplot2::geom_polygon(
      data = poly, ggplot2::aes( x, y ), inherit.aes = FALSE,
      fill = NA, color = gate.color, linewidth = 1
    ) +
    ggplot2::scale_x_continuous(
      name = x.lab, breaks = x.breaks, labels = x.labels, limits = x.limits,
      expand = ggplot2::expansion( asp$figure.gate.scale.expand )
    ) +
    ggplot2::scale_y_continuous(
      name = y.lab, breaks = y.breaks, labels = y.labels, limits = y.limits,
      expand = ggplot2::expansion( asp$figure.gate.scale.expand )
    ) +
    .theme.biplot( asp, legend.position = "none" ) +
    ggplot2::theme( aspect.ratio = 1 )
}

## KDE-smoothed 1D histogram on the peak channel with two brackets: the
## negative-gate range (negative.quantile.min to negative.quantile quantiles,
## "Negative") and the top n.bright.events ("Positive"). The x-axis lower
## limit is supplied directly by the caller (x.min.value, typically the
## x.min.quantile quantile of the data) rather than derived from the raw
## data minimum, since flow data can have a long, sparse negative tail that
## would otherwise force a lot of unused blank space onto the plot.
.dual.selection.panel <- function(
    values, negative.range, positive.range, x.min.value, asp, x.lab,
    fill.color = "steelblue", line.color = "black",
    negative.bracket.color = "#377EB8", positive.bracket.color = "#E41A1C"
) {
  biexp   <- .biexp.transform.legacy( asp )
  x.trans <- biexp( values )

  dens    <- stats::density( x.trans, na.rm = TRUE )
  dens.df <- data.frame( x = dens$x, y = dens$y )
  y.max   <- max( dens.df$y, na.rm = TRUE )
  bar.y   <- y.max * 1.15

  make.bracket <- function( rng ) {
    bx <- biexp( rng )
    data.frame(
      x = c( bx[ 1L ], bx[ 1L ], bx[ 2L ], bx[ 2L ] ),
      y = c( bar.y * 0.95, bar.y, bar.y, bar.y * 0.95 )
    )
  }

  neg.bracket <- make.bracket( negative.range )
  pos.bracket <- make.bracket( positive.range )

  breaks      <- asp$ribbon.breaks
  limits      <- c( x.min.value, asp$expr.data.max )
  axis.labels <- sapply( breaks, function( v ) {
    if ( v == 0 ) "0" else parse( text = paste0( "10^", log10( abs( v ) ) ) )
  } )

  ggplot2::ggplot( dens.df, ggplot2::aes( x, y ) ) +
    ggplot2::geom_area( fill = fill.color, color = NA, alpha = 0.7, na.rm = TRUE ) +
    ggplot2::geom_line( color = line.color, linewidth = asp$figure.spectra.line.size, na.rm = TRUE ) +
    ggplot2::geom_path(
      data = neg.bracket, ggplot2::aes( x, y ),
      color = negative.bracket.color, linewidth = 0.9, inherit.aes = FALSE
    ) +
    ggplot2::geom_path(
      data = pos.bracket, ggplot2::aes( x, y ),
      color = positive.bracket.color, linewidth = 0.9, inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text", x = mean( neg.bracket$x ), y = bar.y * 1.12,
      label = "Negative", size = 3.2, color = negative.bracket.color
    ) +
    ggplot2::annotate(
      "text", x = mean( pos.bracket$x ), y = bar.y * 1.12,
      label = "Positive", size = 3.2, color = positive.bracket.color
    ) +
    ggplot2::scale_x_continuous(
      name = x.lab, breaks = biexp( breaks ), limits = biexp( limits ), labels = axis.labels
    ) +
    ggplot2::scale_y_continuous( name = "Density" ) +
    ggplot2::coord_cartesian( ylim = c( 0, bar.y * 1.3 ) ) +
    .theme.biplot( asp )
}

## Cowplot triple: negative fraction (left, plain black points -- matching
## the Unstained panel convention in .make.cosine.filter.biplot()), positive
## fraction (middle, coloured by cosine similarity to the negative fraction's
## median profile), and the colour-bar legend split into its own column
## (right) so both biplots keep equal width.
.make.dual.cosine.biplot <- function(
    neg.mat, pos.mat, x.dim, y.dim, pos.cs.vals, asp,
    negative.point.color = "black",
    point.size = NULL, max.points = 5e4, y.lab.suffix = "",
    gate.color = "darkgoldenrod1"
) {
  biexp   <- .biexp.transform.legacy( asp )
  y.lab   <- if ( nzchar( y.lab.suffix ) ) paste( y.dim, y.lab.suffix ) else y.dim
  pt.size <- if ( !is.null( point.size ) ) point.size else asp$figure.gate.point.size * 1.3

  # -- negative fraction (plain black points, no colour mapping)
  df.neg <- data.frame( x = neg.mat[ , x.dim ], y = neg.mat[ , y.dim ] )
  if ( nrow( df.neg ) > max.points ) {
    set.seed( asp$bird.seed )
    df.neg <- df.neg[ sample( seq_len( nrow( df.neg ) ), max.points ), , drop = FALSE ]
  }
  df.neg$x.trans <- biexp( df.neg$x )
  df.neg$y.trans <- biexp( df.neg$y )

  p.neg <- ggplot2::ggplot( df.neg, ggplot2::aes( x.trans, y.trans ) ) +
    scattermore::geom_scattermore(
      pointsize = asp$figure.gate.point.size, color = negative.point.color,
      alpha = 1, na.rm = TRUE
    ) +
    .biplot.scales( biexp, asp, x.dim, y.lab ) +
    ggplot2::labs( title = "Negative fraction" ) +
    .theme.biplot( asp ) +
    ggplot2::theme( aspect.ratio = 1 )

  # -- positive fraction (coloured by cosine similarity, continuous)
  df.pos <- data.frame( x = pos.mat[ , x.dim ], y = pos.mat[ , y.dim ], CosineSim = pos.cs.vals )
  if ( nrow( df.pos ) > max.points ) {
    set.seed( asp$bird.seed )
    idx    <- sample( seq_len( nrow( df.pos ) ), max.points )
    df.pos <- df.pos[ idx, , drop = FALSE ]
  }
  df.pos$x.trans <- biexp( df.pos$x )
  df.pos$y.trans <- biexp( df.pos$y )

  # this workflow does not subset the positive fraction any further --
  # every event shown here is carried into the population-level average --
  # so the "gate" is simply the panel's own square axis limits, boxed to
  # make that explicit
  box.limits <- biexp( c( asp$ribbon.plot.min, asp$expr.data.max ) )

  p.pos <- ggplot2::ggplot( df.pos, ggplot2::aes( x.trans, y.trans, color = CosineSim ) ) +
    scattermore::geom_scattermore( pointsize = pt.size, na.rm = TRUE ) +
    .cosine.gradient.scale( name = "Cosine similarity\nto negative" ) +
    ggplot2::annotate(
      "rect",
      xmin = box.limits[ 1L ], xmax = box.limits[ 2L ],
      ymin = box.limits[ 1L ], ymax = box.limits[ 2L ],
      fill = NA, color = gate.color, linewidth = 1
    ) +
    .biplot.scales( biexp, asp, x.dim, y.lab ) +
    ggplot2::labs( title = "Positive fraction (all events retained)" ) +
    .theme.biplot( asp, legend.position = "right" ) +
    ggplot2::theme( aspect.ratio = 1 )

  # extract the colour-bar legend into its own column so the two biplots
  # above stay the same width as each other
  legend.grob <- cowplot::get_legend( p.pos )
  p.pos       <- p.pos + ggplot2::theme( legend.position = "none" )

  cowplot::plot_grid(
    p.neg, p.pos, legend.grob, ncol = 3, rel_widths = c( 1, 1, 0.35 )
  )
}


# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' @title Plot Standard (Manual-Gating) Spectra Extraction Workflow
#'
#' @description
#' Builds a manuscript-ready, multi-panel figure illustrating a "standard"
#' single-stained-control workflow of the kind used in vendor software (e.g.
#' SpectroFlo): an octagon gate on FSC-A vs SSC-A positioned on the dominant
#' population, selection of the brightest events on the fluorophore's peak
#' channel as the positive population and the dimmest fraction as an internal
#' negative, biplots of both fractions (negative in black, positive coloured
#' by cosine similarity), and a population-
#' level background subtraction. Intended as a direct comparison against
#' [spectra.automated.steps.plot()] for the same fluorophore(s). Uses the same
#' building blocks as the rest of AutoSpectral ([create.biplot()]-style
#' pseudocolour density, [spectral.trace()]) so panel styling matches the
#' package's other figures.
#'
#' Requires `plot_spectra_automated_steps.R` to be loaded in the same package
#' namespace (several private helpers are shared).
#'
#' @param control.dir Character. Path to the directory containing the
#'   single-stained control FCS files.
#' @param control.def.file Character. Path to (or filename of) the control
#'   definition CSV -- for the plotted fluorophore itself, only the
#'   `fluorophore` / `filename` columns are used, since this workflow always
#'   derives an internal negative from the same file. If the file has a
#'   `control.type` column (as used by [define.flow.control()] /
#'   [check.control.file()]), a row with `control.type == "beads"` sharing a
#'   fluorophore name with a "cells" row is treated as a paired bead control:
#'   its `filename` is the positive bead sample, and its `universal.negative`
#'   must point at the matching unstained bead file. When present, this pair
#'   supplies the "<Fluorophore> (Beads)" trace in panel D instead of the
#'   static spectral reference library.
#' @param asp The AutoSpectral parameter list from `get.autospectral.param()`.
#' @param fluorophores Character vector of fluorophore name(s) to illustrate.
#'   Default `NULL` illustrates the first fluorophore found.
#' @param octagon.width.factor Numeric, default `3`. Controls how large the
#'   octagon gate is drawn, as a multiple of the local median-absolute-
#'   deviation spread of events around the density peak.
#' @param n.bright.events Integer, default `2000`. Number of top-expressing
#'   events on the peak channel (within the octagon gate) selected as the
#'   "positive" population.
#' @param negative.quantile Numeric in `(0, 1)`, default `0.25`. Upper
#'   quantile bound (on the peak channel, within octagon-gated events) of the
#'   internal "negative" population.
#' @param negative.quantile.min Numeric in `[0, negative.quantile)`, default
#'   `0.01`. Lower quantile bound of the internal "negative" population --
#'   the negative gate spans events between the `negative.quantile.min` and
#'   `negative.quantile` quantiles of the peak channel, excluding the extreme
#'   dim tail (typically debris) below `negative.quantile.min`.
#' @param x.min.quantile Numeric in `[0, 1)`, default `0.005`. Quantile of the
#'   peak-channel data used as the x-axis lower limit in panel B. Flow data
#'   can have a long, sparse negative tail, so this is set independently of
#'   `negative.quantile.min` to keep the histogram readable.
#' @param singlet.quantiles Numeric, default `c(0.85, 0.975)`. Quantile
#'   thresholds used only when cleaning a paired bead control (see
#'   `control.def.file`), matching [get.spectra.automated()].
#' @param gate.color Colour of the octagon gate boundary in panel A, and of
#'   the gate box drawn around the positive fraction in panel C. Default
#'   `"darkgoldenrod1"` (matching `do.gate()`'s default).
#' @param density.palette Fill palette for the pseudocolour density in panel
#'   A: one of the viridis options (`"plasma"`, `"viridis"`, etc.) or any
#'   other value to use `asp$density.palette.base.color`. Default `"rainbow"`
#'   (matching `create.biplot()`'s default).
#' @param selection.fill.color,selection.line.color Fill and line colours for
#'   the KDE histogram in panel B. Defaults `"steelblue"` / `"black"`.
#' @param negative.bracket.color,positive.bracket.color Colours for the
#'   negative- and positive-selection brackets in panel B. Defaults
#'   `"#377EB8"` (blue) / `"#E41A1C"` (red).
#' @param negative.point.color Colour for the negative-fraction events in
#'   panel C, which are plotted as plain points with no colour mapping
#'   (matching the "Unstained" convention in
#'   [spectra.automated.steps.plot()]'s panel D). Default `"black"`. Only the
#'   positive-fraction events are coloured by cosine similarity.
#' @param event.point.size Numeric or `NULL` (default). Point size for the
#'   positive-fraction events in panel C. If `NULL`, defaults to
#'   `asp$figure.gate.point.size * 1.3`.
#' @param ground.truth.method Character, one of `"automated"` (default),
#'   `"legacy"`, or `"none"`. Determines how the "true positive" events
#'   highlighted in red in panel A are identified. `"automated"` replicates
#'   [get.spectra.automated()]'s own candidate + cosine-to-AF-reference filter
#'   on this control file. `"legacy"` uses the actual AF-removal
#'   `gate.population.idx` from `remove.af()`/[clean.controls()], which
#'   requires running [define.flow.control()] and [clean.controls()] over the
#'   whole control set (see `legacy.flow.control`/`legacy.diagnostics.env` to
#'   avoid repeating this per fluorophore). `"none"` reverts to this
#'   workflow's own top-`n.highlight` cosine-similarity selection, which is
#'   not an independent ground truth and is kept only for reference.
#' @param truth.n.candidates,truth.n.spectral Integers, defaults `1000` and
#'   `200`. Only used when `ground.truth.method = "automated"`; mirror
#'   `n.candidates`/`n.spectral` in [get.spectra.automated()].
#' @param legacy.gating.system,legacy.af.remove,legacy.universal.negative,
#'   legacy.downsample,legacy.scatter.match,legacy.k.neighbors,
#'   legacy.negative.n,legacy.positive.n Only used when
#'   `ground.truth.method = "legacy"` and `legacy.flow.control`/
#'   `legacy.diagnostics.env` are not supplied; passed through to
#'   [define.flow.control()] / [clean.controls()] exactly as in
#'   [spectra.legacy.steps.plot()].
#' @param legacy.flow.control,legacy.diagnostics.env Optional, default `NULL`.
#'   Precomputed outputs of [define.flow.control()] + [clean.controls()]
#'   (the latter called with a `diagnostics.env`). Supply both together to
#'   skip re-running the legacy pipeline when illustrating multiple
#'   fluorophores with `ground.truth.method = "legacy"` -- e.g. by calling
#'   [spectra.legacy.steps.plot()] first and reusing its internal objects, or
#'   running the two functions yourself as shown in `spectra.legacy.steps.plot()`.
#' @param n.highlight Integer, default `200`. Only used when
#'   `ground.truth.method = "none"`. Number of positive-fraction events
#'   (smallest cosine similarity to the negative fraction, i.e. least
#'   AF-like) highlighted in red in panel A.
#' @param clean.positive.color Colour for the panel A ground-truth highlight.
#'   Default `"red"`.
#' @param clean.positive.point.size Numeric or `NULL` (default). Point size
#'   for the panel A highlight. If `NULL`, defaults to
#'   `asp$figure.gate.point.size * 1.5`.
#' @param cells.trace.color,beads.trace.color,af.trace.color Colours for the
#'   three traces in panel D: the population-level background-subtracted
#'   profile ("Cells"), the reference profile ("Beads"), and the negative-
#'   fraction trace. Defaults `"#D95F02"` / `"#377EB8"` / `"grey40"`.
#' @param max.points Integer. Maximum events plotted per panel (randomly
#'   downsampled beyond this for speed). Default `5e4`.
#' @param panel.width,panel.height Numeric. Width/height (inches) used per
#'   sub-panel when sizing the saved composite figure. Defaults `4` and `4`.
#' @param composite.width,composite.height Numeric or `NULL` (default).
#'   Override the overall saved figure dimensions (inches); if `NULL`, these
#'   are computed from `panel.width` / `panel.height`.
#' @param output.dir Character or `NULL` (default). Directory to save the
#'   composite figure(s). Defaults to the current working directory.
#' @param save Logical, default `TRUE`. Whether to save the composite figure
#'   for each fluorophore to `output.dir`.
#' @param file.type Character string, one of `"jpg"` (default), `"tiff"`,
#'   `"png"`, or `"pdf"`.
#' @param verbose Logical, default `TRUE`. Print progress messages.
#'
#' @return Invisibly, a named list (one entry per fluorophore), each
#'   containing:
#'   \describe{
#'     \item{`gate.panel`}{Panel A, the octagon gate on FSC-A vs SSC-A with
#'       the ground-truth positive events highlighted.}
#'     \item{`selection.panel`}{Panel B, the brightest/negative event
#'       selection histogram.}
#'     \item{`cosine.panel`}{Panel C, the negative/positive-fraction
#'       cosine-similarity biplots.}
#'     \item{`subtraction.plot`}{Panel D, the final spectral profile
#'       comparison ([spectral.trace()] of Cells / Beads / AF).}
#'     \item{`composite`}{The assembled four-panel cowplot object saved to
#'       `output.dir` when `save = TRUE`.}
#'     \item{`peak.channel`}{Character. The fluorophore's nominal peak
#'       channel, looked up from `fluorophore_database.csv` for
#'       `asp$cytometer`.}
#'     \item{`y.channel.peak`}{Character. The non-colliding peak AF channel
#'       used as the y-axis of panel C.}
#'     \item{`reference.profile`}{Named numeric vector (over
#'       `spectral.channels`) used as the "Beads" trace in panel D, or
#'       `NULL` if neither a paired bead control nor the spectral reference
#'       library had data for this fluorophore.}
#'   }
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_manual theme_bw theme geom_polygon
#' @importFrom ggplot2 stat_density_2d after_stat scale_fill_viridis_c
#' @importFrom ggplot2 scale_fill_gradientn geom_area geom_line geom_path
#' @importFrom ggplot2 annotate coord_cartesian labs ggsave expansion
#' @importFrom scattermore geom_scattermore
#' @importFrom cowplot plot_grid get_legend
#' @importFrom MASS kde2d bandwidth.nrd
#' @importFrom stats density mad median quantile setNames
#' @importFrom ragg agg_jpeg agg_tiff agg_png
#'
#' @seealso [spectra.automated.steps.plot()], [do.gate()], [create.biplot()],
#'   [spectral.trace()]
#'
#' @export

spectra.standard.workflow.plot <- function(
    control.dir,
    control.def.file,
    asp,
    fluorophores             = NULL,
    octagon.width.factor     = 3,
    n.bright.events          = 2000L,
    negative.quantile        = 0.25,
    negative.quantile.min    = 0.01,
    x.min.quantile           = 0.005,
    singlet.quantiles         = c( 0.85, 0.975 ),
    gate.color               = "darkgoldenrod1",
    density.palette          = "rainbow",
    selection.fill.color     = "steelblue",
    selection.line.color     = "black",
    negative.bracket.color   = "#377EB8",
    positive.bracket.color   = "#E41A1C",
    negative.point.color     = "black",
    event.point.size         = NULL,
    ground.truth.method       = c( "automated", "legacy", "none" ),
    truth.n.candidates         = 1000L,
    truth.n.spectral           = 200L,
    legacy.gating.system        = c( "density", "landmarks" ),
    legacy.af.remove             = TRUE,
    legacy.universal.negative    = TRUE,
    legacy.downsample             = TRUE,
    legacy.scatter.match          = TRUE,
    legacy.k.neighbors            = 3L,
    legacy.negative.n              = asp$negative.n,
    legacy.positive.n              = asp$positive.n,
    legacy.flow.control            = NULL,
    legacy.diagnostics.env         = NULL,
    n.highlight               = 200L,
    clean.positive.color      = "red",
    clean.positive.point.size = NULL,
    cells.trace.color         = "#D95F02",
    beads.trace.color         = "#377EB8",
    af.trace.color            = "grey40",
    max.points               = 5e4,
    panel.width              = 4,
    panel.height             = 4,
    composite.width          = NULL,
    composite.height         = NULL,
    output.dir                = NULL,
    save                      = TRUE,
    file.type                 = "jpg",
    verbose                   = TRUE
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

  if ( negative.quantile <= 0 || negative.quantile >= 1 )
    stop( "negative.quantile must be strictly between 0 and 1.", call. = FALSE )
  if ( negative.quantile.min < 0 || negative.quantile.min >= negative.quantile )
    stop( "negative.quantile.min must be in [0, negative.quantile).", call. = FALSE )
  if ( x.min.quantile < 0 || x.min.quantile >= 1 )
    stop( "x.min.quantile must be in [0, 1).", call. = FALSE )

  ground.truth.method  <- match.arg( ground.truth.method, c( "automated", "legacy", "none" ) )
  legacy.gating.system <- match.arg( legacy.gating.system, c( "density", "landmarks" ) )

  legacy.precomputed <- !is.null( legacy.flow.control ) || !is.null( legacy.diagnostics.env )
  if ( legacy.precomputed && ( is.null( legacy.flow.control ) || is.null( legacy.diagnostics.env ) ) )
    stop(
      "legacy.flow.control and legacy.diagnostics.env must both be supplied together, or both left NULL.",
      call. = FALSE
    )

  # -- 1. Read control table (only fluorophore/filename are used)
  ctrl.path <- if ( file.exists( control.def.file ) ) {
    control.def.file
  } else {
    file.path( control.dir, control.def.file )
  }
  ctrl.tbl <- utils::read.csv(
    ctrl.path, stringsAsFactors = FALSE, strip.white = TRUE
  )

  is.beads.row <- .is.beads.row( ctrl.tbl )

  fluor.rows.all <- which(
    !is.beads.row &
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

  # -- 2. Resolve channels + fluorophore database
  scatter.channels <- read.scatter.parameter( asp )
  if ( length( scatter.channels ) < 2L )
    stop( "At least two scatter channels (FSC/SSC) are required.", call. = FALSE )

  first.fcs.path    <- file.path( control.dir, ctrl.tbl$filename[ fluor.rows.all[ 1L ] ] )
  spectral.channels <- .derive.spectral.channels( first.fcs.path, asp )
  sat.value         <- if ( !is.null( asp$expr.data.max ) ) asp$expr.data.max else Inf

  fluor.data.path <- system.file(
    "extdata", "fluorophore_database.csv", package = "AutoSpectral"
  )
  fluorophore.database <- utils::read.csv( fluor.data.path )
  fluorophore.database[ fluorophore.database == "" ] <- NA

  db.col <- .cytometer.to.db.col( asp$cytometer )

  if ( ground.truth.method == "legacy" && !legacy.precomputed ) {
    if ( verbose )
      message( "\033[34m-- Running legacy pipeline for ground-truth positives (define.flow.control() + clean.controls()) --\033[0m" )
    legacy.diagnostics.env <- new.env( parent = emptyenv() )
    legacy.flow.control <- define.flow.control(
      control.dir = control.dir, control.def.file = control.def.file, asp = asp,
      gate = TRUE, gating.system = legacy.gating.system, parallel = FALSE, verbose = verbose
    )
    legacy.flow.control <- clean.controls(
      legacy.flow.control, asp,
      af.remove = legacy.af.remove, universal.negative = legacy.universal.negative,
      downsample = legacy.downsample, negative.n = legacy.negative.n, positive.n = legacy.positive.n,
      scatter.match = legacy.scatter.match, k.neighbors = legacy.k.neighbors,
      intermediate.figures = FALSE, main.figures = FALSE, parallel = FALSE,
      verbose = verbose, diagnostics.env = legacy.diagnostics.env
    )
  }

  fsc.a <- scatter.channels[ 1L ]
  ssc.a <- scatter.channels[ 2L ]

  composite.width.use  <- if ( !is.null( composite.width ) )  composite.width  else panel.width * 2
  row.heights           <- rep( panel.height, 4L )
  composite.height.use  <- if ( !is.null( composite.height ) ) composite.height else sum( row.heights )
  clean.positive.point.size.use <- if ( !is.null( clean.positive.point.size ) )
    clean.positive.point.size else asp$figure.gate.point.size * 1.5

  results <- list()

  # -- 3. Loop over requested fluorophores
  for ( fluor in fluorophores ) {

    if ( verbose )
      message( "\033[34m-- Building standard-workflow figure for ", fluor, " --\033[0m" )

    peak.channel <- tryCatch(
      .lookup.fluorophore.channel( fluor, asp, fluorophore.database ),
      error = function( e ) {
        warning( e$message, call. = FALSE )
        NULL
      }
    )
    if ( is.null( peak.channel ) ) next

    row.i      <- fluor.rows.all[ match( fluor, all.fluor.names ) ]
    fluor.file <- ctrl.tbl$filename[ row.i ]
    fcs.path.i <- file.path( control.dir, fluor.file )

    ref.profile.i <- .get.reference.profile(
      fluor, ctrl.tbl, control.dir, spectral.channels, scatter.channels,
      sat.value, singlet.quantiles, asp, db.col, verbose = FALSE
    )

    # an explicit external unstained control for this fluorophore
    # (universal.negative pointing at a file, exactly as in the automated
    # workflow) takes priority for the "Unstained (Autofluorescence)" trace;
    # falls back to the in-sample negative gate (panels B/C are unaffected)
    external.af.mean <- NULL
    unstained.src.i  <- .parse.unstained.source( ctrl.tbl$universal.negative[ row.i ] )

    if ( unstained.src.i$type == "file" ) {
      external.af.path <- file.path( control.dir, unstained.src.i$file )
      if ( !file.exists( external.af.path ) ) {
        warning(
          "External unstained file not found for '", fluor, "': ", unstained.src.i$file,
          "; falling back to the internal negative gate for the Autofluorescence trace.",
          call. = FALSE
        )
      } else {
        external.af.mat <- .read.fcs.clean(
          external.af.path, paste0( "Unstained (", unstained.src.i$file, ")" ),
          spectral.channels, scatter.channels, sat.value, singlet.quantiles,
          asp = asp, verbose = FALSE
        )
        af.spec.cols <- intersect( spectral.channels, colnames( external.af.mat ) )
        if ( nrow( external.af.mat ) >= 2L && length( af.spec.cols ) > 0L ) {
          external.af.mean <- stats::setNames(
            rep( NA_real_, length( spectral.channels ) ), spectral.channels
          )
          external.af.mean[ af.spec.cols ] <- colMeans( external.af.mat[ , af.spec.cols, drop = FALSE ] )
        } else {
          warning(
            "Too few clean events in the external unstained file for '", fluor,
            "'; falling back to the internal negative gate for the Autofluorescence trace.",
            call. = FALSE
          )
        }
      }
    }

    probe.cols.i <- colnames( readFCS( fcs.path.i, start.row = 1, end.row = 1 ) )
    cols.keep    <- intersect( c( scatter.channels, spectral.channels ), probe.cols.i )
    mat          <- readFCS( fcs.path.i, columns = cols.keep )

    spec.present <- intersect( spectral.channels, colnames( mat ) )
    if ( length( spec.present ) > 0 ) {
      mat <- mat[ rowSums( mat[ , spec.present, drop = FALSE ] >= sat.value ) == 0, , drop = FALSE ]
    }

    if ( !peak.channel %in% colnames( mat ) ) {
      warning(
        "Peak channel '", peak.channel, "' not found for '", fluor, "'. Skipping.",
        call. = FALSE
      )
      next
    }

    # -- Ground-truth positive identification for panel A (see
    # ground.truth.method): resolved independently of this workflow's own
    # octagon gate / brightness / cosine-similarity selection, so the panel A
    # highlight compares against (rather than reproduces) that selection.
    truth.mat <- if ( ground.truth.method == "automated" ) {
      tryCatch(
        .identify.true.positives.automated(
          fcs.path.i        = fcs.path.i,
          spectral.channels = spectral.channels,
          scatter.channels  = scatter.channels,
          sat.value         = sat.value,
          singlet.quantiles = singlet.quantiles,
          peak.channel      = peak.channel,
          unstained.src.i   = unstained.src.i,
          control.dir       = control.dir,
          n.candidates      = truth.n.candidates,
          n.spectral        = truth.n.spectral
        ),
        error = function( e ) {
          warning(
            "Automated ground-truth selection failed for '", fluor, "': ", e$message,
            call. = FALSE
          )
          NULL
        }
      )
    } else if ( ground.truth.method == "legacy" ) {
      samp.truth.i <- legacy.flow.control$sample[
        which( legacy.flow.control$fluorophore == fluor )[ 1L ]
      ]
      diag.truth.i <- if ( !is.null( samp.truth.i ) )
        legacy.diagnostics.env[[ samp.truth.i ]] else NULL

      if ( is.null( diag.truth.i ) ) {
        NULL
      } else if ( is.null( diag.truth.i$scatter.data.pos ) ) {
        warning(
          "legacy.diagnostics.env for '", fluor, "' has no scatter.data.pos -- ",
          "requires the updated remove.af() that stores scatter columns.",
          call. = FALSE
        )
        NULL
      } else {
        diag.truth.i$scatter.data.pos[ diag.truth.i$gate.population.idx, , drop = FALSE ]
      }
    } else {
      NULL
    }

    # -- A. Octagon gate (position via simplified 2D KDE peak search)
    peak.pt <- .find.density.peak.2d( mat[ , fsc.a ], mat[ , ssc.a ] )

    near.idx <- which(
      abs( mat[ , fsc.a ] - peak.pt$x ) < diff( range( mat[ , fsc.a ] ) ) * 0.15 &
        abs( mat[ , ssc.a ] - peak.pt$y ) < diff( range( mat[ , ssc.a ] ) ) * 0.15
    )
    if ( length( near.idx ) < 10L ) near.idx <- seq_len( nrow( mat ) )

    mad.x  <- stats::mad( mat[ near.idx, fsc.a ], center = peak.pt$x )
    mad.y  <- stats::mad( mat[ near.idx, ssc.a ], center = peak.pt$y )
    half.x <- octagon.width.factor * max( mad.x, 1 )
    half.y <- octagon.width.factor * max( mad.y, 1 )

    octagon.poly <- .octagon.polygon( peak.pt$x, peak.pt$y, half.x, half.y )

    gate.panel <- .octagon.gate.panel(
      mat[ , fsc.a ], mat[ , ssc.a ], octagon.poly, asp, fsc.a, ssc.a,
      gate.color = gate.color, density.palette = density.palette, max.points = max.points
    )

    # gate.row (the padded, left-aligned wrapper around gate.panel) is
    # assembled later, once panel C's cosine filter has identified the
    # "clean positive" events, so they can be highlighted here first

    in.gate    <- .point.in.polygon( mat[ , fsc.a ], mat[ , ssc.a ], octagon.poly$x, octagon.poly$y )
    gated.mat  <- mat[ in.gate, , drop = FALSE ]

    if ( nrow( gated.mat ) < 10L ) {
      warning( "Too few events inside octagon gate for '", fluor, "'. Skipping.", call. = FALSE )
      next
    }

    # -- B. Brightest / negative event selection
    peak.vals  <- gated.mat[ , peak.channel ]
    order.peak <- order( peak.vals, decreasing = TRUE )

    n.bright.actual <- min( n.bright.events, nrow( gated.mat ) )
    positive.idx    <- order.peak[ seq_len( n.bright.actual ) ]

    # negative gate: events between the negative.quantile.min and
    # negative.quantile quantiles of the peak channel (excludes the extreme
    # dim tail, typically debris, below negative.quantile.min)
    neg.bound.low  <- stats::quantile( peak.vals, negative.quantile.min )
    neg.bound.high <- stats::quantile( peak.vals, negative.quantile )
    negative.idx   <- which( peak.vals >= neg.bound.low & peak.vals <= neg.bound.high )

    if ( length( negative.idx ) < 1L ) {
      warning( "No events in the negative gate for '", fluor, "'. Skipping.", call. = FALSE )
      next
    }

    x.min.value <- stats::quantile( peak.vals, x.min.quantile )

    selection.panel <- .dual.selection.panel(
      values          = peak.vals,
      negative.range  = c( neg.bound.low, neg.bound.high ),
      positive.range  = range( peak.vals[ positive.idx ] ),
      x.min.value     = x.min.value,
      asp             = asp, x.lab = peak.channel,
      fill.color      = selection.fill.color, line.color = selection.line.color,
      negative.bracket.color = negative.bracket.color,
      positive.bracket.color = positive.bracket.color
    )

    positive.mat <- gated.mat[ positive.idx, spectral.channels, drop = FALSE ]
    negative.mat <- gated.mat[ negative.idx, spectral.channels, drop = FALSE ]

    # -- C. Cosine-similarity biplots (negative left, positive right)
    neg.median <- apply( negative.mat, 2, stats::median )
    candidate.channels <- setdiff( colnames( negative.mat ), peak.channel )
    y.channel.peak <- names( which.max( neg.median[ candidate.channels ] ) )

    pos.cs.vals <- .cosine.sim.rows( as.matrix( positive.mat ), neg.median )

    n.highlight.actual <- min( n.highlight, length( pos.cs.vals ) )
    highlight.idx <- positive.idx[ order( pos.cs.vals )[ seq_len( n.highlight.actual ) ] ]
    highlight.mat <- gated.mat[ highlight.idx, , drop = FALSE ]

    # panel A's red highlight uses the independently-identified ground-truth
    # positives (truth.mat, resolved above via ground.truth.method) rather
    # than this workflow's own brightness + cosine-similarity selection --
    # that selection is exactly what this figure is meant to critique, so
    # using it to mark "true positives" would be circular. Falls back to the
    # in-workflow selection only when ground.truth.method == "none" or the
    # ground-truth lookup failed for this fluorophore.
    highlight.source <- if ( !is.null( truth.mat ) &&
                             all( c( fsc.a, ssc.a ) %in% colnames( truth.mat ) ) ) {
      truth.mat
    } else {
      if ( ground.truth.method != "none" )
        warning(
          "No '", ground.truth.method, "' ground-truth positives available for '",
          fluor, "'; falling back to the standard workflow's own selection for the panel A highlight.",
          call. = FALSE
        )
      highlight.mat
    }

    gate.panel <- .add.highlight.layer(
      gate.panel, highlight.source[ , fsc.a ], highlight.source[ , ssc.a ],
      color = clean.positive.color, pointsize = clean.positive.point.size.use
    )

    # panel A is forced square (aspect.ratio = 1); its row cell is
    # composite.width.use wide by row.heights[1] tall, which is usually wider
    # than tall, so ggplot's own aspect-ratio padding would centre the square
    # in blank space. Instead, size the column to just fit the square and pad
    # only on the right, so the panel stays left-aligned.
    gate.square.side <- min( composite.width.use, row.heights[ 1L ] )
    gate.row <- if ( gate.square.side < composite.width.use ) {
      cowplot::plot_grid(
        gate.panel, NULL, ncol = 2,
        rel_widths = c( gate.square.side, composite.width.use - gate.square.side )
      )
    } else {
      gate.panel
    }

    cosine.panel <- .make.dual.cosine.biplot(
      neg.mat = negative.mat, pos.mat = positive.mat,
      x.dim = peak.channel, y.dim = y.channel.peak,
      pos.cs.vals = pos.cs.vals, asp = asp,
      negative.point.color = negative.point.color,
      point.size = event.point.size, max.points = max.points,
      y.lab.suffix = "(peak AF channel)", gate.color = gate.color
    )

    # -- D. Final spectral profile comparison: the population-level
    # background-subtracted signature ("Cells"), the reference profile for
    # the same fluorophore ("Beads" -- a paired bead control if provided,
    # otherwise the static spectral reference library), and the negative-
    # fraction trace as the autofluorescence estimate
    neg.mean   <- colMeans( negative.mat )
    pos.mean   <- colMeans( positive.mat )
    final.mean <- pos.mean - neg.mean

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

    af.trace <- if ( !is.null( external.af.mean ) ) external.af.mean else neg.mean
    af.trace[ !is.finite( af.trace ) ] <- 0

    subtraction.mat <- rbind(
      final.mean / max( abs( final.mean ), 1e-9 ),
      beads.trace,
      af.trace / max( abs( af.trace ), 1e-9 )
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

    # -- Assemble composite figure
    composite <- cowplot::plot_grid(
      cowplot::plot_grid( gate.row,          labels = "A" ),
      cowplot::plot_grid( selection.panel,   labels = "B" ),
      cowplot::plot_grid( cosine.panel,      labels = "C" ),
      cowplot::plot_grid( subtraction.plot,  labels = "D" ),
      ncol = 1, rel_heights = row.heights
    )

    if ( save ) {
      out.file <- file.path(
        output.dir, sprintf( "%s_standard_workflow.%s", fluor, file.type )
      )
      ggplot2::ggsave(
        out.file, plot = composite, device = plot.device,
        width = composite.width.use, height = composite.height.use, limitsize = FALSE
      )
      if ( verbose ) message( "\033[32m  Saved: ", out.file, "\033[0m" )
    }

    results[[ fluor ]] <- list(
      gate.panel        = gate.panel,
      selection.panel    = selection.panel,
      cosine.panel       = cosine.panel,
      subtraction.plot   = subtraction.plot,
      composite          = composite,
      peak.channel       = peak.channel,
      y.channel.peak     = y.channel.peak,
      reference.profile  = ref.profile.i
    )
  }

  invisible( results )
}
