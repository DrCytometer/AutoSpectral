# plot_spectra_legacy_steps.R
#
# Manuscript figure helper: illustrates the internal steps of the "legacy"
# workflow -- define.flow.control(), clean.controls(), and
# get.fluorophore.spectra() -- for one or more single-stained control
# samples, using the same biplot / gate-panel / spectral-trace building
# blocks as the rest of AutoSpectral so panels are stylistically consistent
# with the automated and standard workflow figures. Depends on private
# helpers defined in plot_spectra_automated_steps.R and
# plot_spectra_standard_workflow.R (.biexp.transform.legacy(),
# .theme.biplot(), .biplot.scales(), .cosine.gradient.scale(),
# .cosine.sim.rows(), .octagon.gate.panel(), .embed.or.placeholder(),
# .add.highlight.layer(), .get.reference.profile()) -- all three files must
# be loaded into the package together (e.g. all in R/, or all source()'d).
#
# This function calls define.flow.control() and clean.controls() itself for
# the requested control set, so it reproduces the entire legacy pipeline run
# (not just the illustrated fluorophore) and can be slow for large panels.
# clean.controls() is called with diagnostics.env set, which requires the
# diagnostics.env-aware version of remove.af()/run.af.removal()/
# clean.controls() including its `scatter.data.pos` capture (the scatter-parameter
# columns aligned to `expr.data.pos`) -- panel A's true-positive highlight reads
# from that field directly, the same way spectra.standard.workflow.plot()'s
# `ground.truth.method = "legacy"` branch already does.
#
# Panels produced per fluorophore:
#   A. FSC-A vs SSC-A pseudocolour density (as in the standard workflow's
#      panel A), with the *actual* gate boundary computed by
#      define.flow.control() overlaid (re-derived by calling
#      define.gate.landmarks()/define.gate.density() again with the same
#      control table and gating.system -- deterministic given the same
#      asp$bird.seed, and saved to a temp directory rather than
#      asp$figure.gate.dir so this doesn't duplicate the real run's output --
#      or, if `gate.list` was supplied, reusing that pre-defined gate
#      directly, matching define.flow.control()'s own priority). The
#      `n.true.positive` brightest events -- ranked by projection onto the
#      fitted RLM trend direction in panel D's (peak channel, intrusive-AF
#      channel) space, not merely the full AF-gate-excluded population --
#      are highlighted larger in red, the same way true positives are
#      marked in the other two workflow figures. Left-aligned rather than
#      centred, since the panel is forced square.
#   B. Autofluorescence removal, in the automated workflow's panel D style:
#      unstained/AF sample (left, black) and single-stained control (right,
#      all events at this stage, coloured by cosine similarity to the AF
#      median), both with the AF-exclusion gate boundary identified by
#      clean.controls() overlaid. Requires af.remove = TRUE and a paired
#      universal negative for the fluorophore; otherwise shows a
#      placeholder.
#   C. Scatter-matching of the universal negative, embedded directly from
#      the JPEG clean.controls() already saves via scatter.match.plot().
#   D. Robust linear model diagnostic: linear-scale (not biexponential)
#      biplot of the fluorophore's peak channel (x) vs. the identified
#      intrusive-AF channel (y), for the events exactly as clean.controls()
#      finalises them for this sample (post remove.af(), universal-negative
#      scatter-matching, and brightest-event selection -- see
#      run.universal.negative()/downsample.control(), and remove.af()'s own
#      internal scatter-matching when af.remove = TRUE) -- read from
#      flow.control$clean.expr, i.e. the actual population
#      get.fluorophore.spectra() fits its RLM on, which is not the same as
#      the AF-exclusion gate.population.idx panel A's highlight is ranked
#      against. The robust linear fit (as used by get.fluorophore.spectra())
#      is drawn in blue over every plotted event, and the fit's slope and
#      R^2 (squared Pearson correlation) are annotated. Left-aligned rather
#      than centred, since the panel is forced square.
#   E. Final spectral profile comparison, built the same way as the
#      automated workflow's panel F and the standard workflow's panel D: the
#      per-channel signature computed the same way get.fluorophore.spectra()
#      computes it ("<Fluorophore> (Cells)") -- fit to the same
#      flow.control$clean.expr population panel D displays, not just
#      gate.population.idx -- the reference profile for the same fluorophore
#      ("<Fluorophore> (Beads)" -- from a paired bead control in
#      control.def.file if provided, otherwise the cytometer's static
#      spectral reference library), and the matched-negative AF trace
#      ("Unstained (Autofluorescence)"). Labels, trace colours, and legend
#      are consistent with the corresponding panels in
#      spectra.automated.steps.plot() and spectra.standard.workflow.plot().

# ---------------------------------------------------------------------------
# Private helpers specific to the legacy workflow
# ---------------------------------------------------------------------------

## Resolves (and caches in cache.env, since multiple fluorophores can share a
## gate) the real gate boundary for a given gate.name. If `gate.list` was
## supplied and has an entry for `gate.name`, that pre-defined gate is reused
## directly (matching define.flow.control()'s own priority); otherwise the
## boundary is re-derived by calling the same gate-definition function
## define.flow.control() would have used, saving to a temp directory so this
## re-derivation doesn't duplicate output already written by the real
## pipeline run. Returns a data.frame(x, y), or NULL if gate definition
## failed (caller should show a placeholder).
.resolve.legacy.gate <- function(
    gate.name, control.table, control.dir, asp, gating.system, gate.list, cache.env
) {
  if ( !is.null( cache.env[[ gate.name ]] ) ) return( cache.env[[ gate.name ]] )

  if ( !is.null( gate.list ) && gate.name %in% names( gate.list ) ) {
    gate.boundary <- gate.list[[ gate.name ]]
    result <- data.frame( x = gate.boundary$x, y = gate.boundary$y )
    cache.env[[ gate.name ]] <- result
    return( result )
  }

  is.orphan <- grepl( "density_orphan", gate.name )

  gate.boundary <- tryCatch(
    if ( is.orphan || gating.system == "density" ) {
      define.gate.density(
        control.file = NULL, control.dir = control.dir, asp = asp,
        gate.name = gate.name, output.dir = tempdir(), verbose = FALSE,
        control.table = control.table
      )
    } else {
      define.gate.landmarks(
        control.file = NULL, control.dir = control.dir, asp = asp,
        gate.name = gate.name, output.dir = tempdir(), verbose = FALSE,
        control.table = control.table, check = FALSE
      )
    },
    error = function( e ) {
      warning( "Gate definition failed for '", gate.name, "': ", e$message, call. = FALSE )
      NULL
    }
  )

  result <- if ( is.null( gate.boundary ) ) NULL else
    data.frame( x = gate.boundary$x, y = gate.boundary$y )
  cache.env[[ gate.name ]] <- result
  result
}

## Cowplot triple: unstained/AF sample (left, plain black points) and the
## single-stained control (right, coloured by cosine similarity to the AF
## median, matching the automated workflow's panel D), with the AF-exclusion
## gate boundary identified by clean.controls() overlaid on both panels.
.make.af.exclusion.biplot <- function(
    unstained.mat, stained.mat, x.dim, y.dim,
    cs.vals, af.boundary, asp,
    unstained.point.color = "black", max.points = 5e4,
    cosine.point.size = NULL, af.gate.color = "black"
) {
  biexp <- .biexp.transform.legacy( asp )
  stained.pointsize <- if ( !is.null( cosine.point.size ) )
    cosine.point.size else asp$figure.gate.point.size * 1.3

  boundary.df <- NULL
  if ( !is.null( af.boundary ) && !is.null( af.boundary$x ) && length( af.boundary$x ) >= 3L ) {
    boundary.df <- data.frame(
      x = biexp( c( af.boundary$x, af.boundary$x[ 1L ] ) ),
      y = biexp( c( af.boundary$y, af.boundary$y[ 1L ] ) )
    )
  }

  # -- unstained/AF panel (plain black points)
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
    )
  if ( !is.null( boundary.df ) )
    p.unstained <- p.unstained +
      ggplot2::geom_path(
        data = boundary.df, mapping = ggplot2::aes( x, y ), inherit.aes = FALSE,
        color = af.gate.color, linewidth = asp$figure.gate.line.size
      )
  p.unstained <- p.unstained +
    .biplot.scales( biexp, asp, x.dim, y.dim ) +
    ggplot2::labs( title = "Unstained / AF" ) +
    .theme.biplot( asp ) +
    ggplot2::theme( aspect.ratio = 1 )

  # -- single-stained panel (coloured by cosine similarity, continuous)
  df.s <- data.frame(
    x = stained.mat[ , x.dim ], y = stained.mat[ , y.dim ], CosineSim = cs.vals
  )
  if ( nrow( df.s ) > max.points ) {
    set.seed( asp$bird.seed )
    df.s <- df.s[ sample( seq_len( nrow( df.s ) ), max.points ), , drop = FALSE ]
  }
  df.s$x.trans <- biexp( df.s$x )
  df.s$y.trans <- biexp( df.s$y )

  p.stained <- ggplot2::ggplot( df.s, ggplot2::aes( x.trans, y.trans, color = CosineSim ) ) +
    scattermore::geom_scattermore( pointsize = stained.pointsize, na.rm = TRUE ) +
    .cosine.gradient.scale()
  if ( !is.null( boundary.df ) )
    p.stained <- p.stained +
      ggplot2::geom_path(
        data = boundary.df, mapping = ggplot2::aes( x, y ), inherit.aes = FALSE,
        color = af.gate.color, linewidth = asp$figure.gate.line.size
      )
  p.stained <- p.stained +
    .biplot.scales( biexp, asp, x.dim, y.dim ) +
    ggplot2::labs( title = "Single-stained control (all events)" ) +
    .theme.biplot( asp, legend.position = "right" ) +
    ggplot2::theme( aspect.ratio = 1 )

  # extract the colour-bar legend into its own column so the two biplots
  # above stay the same width as each other
  legend.grob <- cowplot::get_legend( p.stained )
  p.stained   <- p.stained + ggplot2::theme( legend.position = "none" )

  cowplot::plot_grid(
    p.unstained, p.stained, legend.grob, ncol = 3, rel_widths = c( 1, 1, 0.35 )
  )
}

## Linear-scale (not biexponential, per the manuscript's request) biplot of
## the fluorophore's peak channel (x) vs. the identified intrusive-AF peak
## channel (y) for the single-stained control, with the robust linear fit
## used for spectral extraction overlaid, the "clean" (AF-gate-excluded)
## events highlighted, and the fit's slope and R^2 annotated.
.make.rlm.biplot <- function(
    x.vals, y.vals, x.lab, y.lab, clean.idx, rlm.coef, r.squared, asp,
    clean.positive.color = "red", clean.positive.point.size = 3,
    base.point.color = "grey50", base.point.size = 1.5,
    rlm.line.color = "blue"
) {
  df <- data.frame( x = x.vals, y = y.vals )
  clean.flag <- rep( FALSE, nrow( df ) )
  clean.flag[ clean.idx ] <- TRUE

  ggplot2::ggplot() +
    scattermore::geom_scattermore(
      data = df[ !clean.flag, , drop = FALSE ], mapping = ggplot2::aes( x, y ),
      pointsize = base.point.size, color = base.point.color, alpha = 0.5, na.rm = TRUE
    ) +
    scattermore::geom_scattermore(
      data = df[ clean.flag, , drop = FALSE ], mapping = ggplot2::aes( x, y ),
      pointsize = clean.positive.point.size, color = clean.positive.color,
      alpha = 1, na.rm = TRUE
    ) +
    ggplot2::geom_abline(
      intercept = rlm.coef[ 1L ], slope = rlm.coef[ 2L ],
      color = rlm.line.color, linewidth = 1
    ) +
    ggplot2::annotate(
      "text", x = -Inf, y = Inf,
      label = sprintf( "slope = %.3g\nR^2 = %.3f", rlm.coef[ 2L ], r.squared ),
      hjust = -0.15, vjust = 1.4, size = 3.5, color = "black"
    ) +
    ggplot2::scale_x_continuous( name = x.lab ) +
    ggplot2::scale_y_continuous( name = y.lab ) +
    ggplot2::labs( title = "Robust linear model: intrusive AF vs. fluorophore peak" ) +
    .theme.biplot( asp ) +
    ggplot2::theme( aspect.ratio = 1 )
}

## Full per-channel robust-linear-model spectral signature for one
## fluorophore's "clean" (AF-gate-excluded) events, replicating
## get.fluorophore.spectra()'s per-row computation (fit.robust.linear.model()
## of every spectral channel against the peak channel, then normalised so the
## peak channel reads 1) without re-running the whole-panel extraction. Used
## to build panel E's "Cells" trace.
.legacy.rlm.spectrum <- function( clean.mat, peak.channel, spectral.channels ) {
  coef <- stats::setNames( rep( 0, length( spectral.channels ) ), spectral.channels )
  coef[ peak.channel ] <- 1.0

  peak.expr <- clean.mat[ , peak.channel ]

  for ( channel in setdiff( spectral.channels, peak.channel ) ) {
    fit <- fit.robust.linear.model(
      peak.expr, clean.mat[ , channel ], peak.channel, channel
    )
    coef[ channel ] <- fit[ 2L ]
  }

  coef / max( coef )
}


# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' @title Plot Legacy Spectra Extraction Pipeline Steps
#'
#' @description
#' Builds a manuscript-ready, multi-panel figure illustrating the internal
#' steps of the legacy workflow ([define.flow.control()], [clean.controls()],
#' [get.fluorophore.spectra()]) for one or more single-stained control
#' samples: the automated scatter gate, autofluorescence-exclusion gating,
#' scatter-matched universal-negative selection, and the robust linear model
#' fit used to extract the fluorophore's spectral signature. Uses the same
#' building blocks as the rest of AutoSpectral ([create.biplot()]-style
#' biexponential biplots, [scatter.match.plot()]) so panel styling matches
#' the package's other figures.
#'
#' This function runs the full legacy pipeline itself -- it calls
#' [define.flow.control()] and [clean.controls()] internally for the entire
#' control set in `control.def.file`, not just the illustrated
#' fluorophore(s) -- so it can be slow for large panels. `parallel = FALSE`
#' is used internally throughout, since diagnostic capture is unreliable
#' under forked parallel processing.
#'
#' Requires `plot_spectra_automated_steps.R` and
#' `plot_spectra_standard_workflow.R` to be loaded in the same package
#' namespace (several private helpers are shared), and requires the
#' `diagnostics.env`-aware versions of `remove.af()` / `run.af.removal()` /
#' `clean.controls()`.
#'
#' @param control.dir Character. Path to the directory containing the
#'   single-stained control FCS files.
#' @param control.def.file Character. Path to (or filename of) the control
#'   definition CSV, in the full legacy format required by
#'   [define.flow.control()] (including `control.type`, `gate.name`,
#'   `gate.define`, etc. -- see [check.control.file()]).
#' @param asp The AutoSpectral parameter list from `get.autospectral.param()`.
#' @param fluorophores Character vector of fluorophore name(s) to illustrate.
#'   Default `NULL` illustrates the first cell-based fluorophore with a
#'   paired universal negative.
#' @param gating.system Character, one of `"density"` (default) or
#'   `"landmarks"`, matching [define.flow.control()]'s argument of the same
#'   name.
#' @param gate.list Optional named list of gates. To use this, pre-define the
#' gates using `define.gate.landmarks()` and/or `define.gate.density()`, ensure
#' that the names of the gates correspond to the names in the `control.def.file`,
#' and ensure that the `gate.name` column has been filled in for the
#' `control.def.file`. Default `NULL` will revert to creating new gates.
#' Passed through to [define.flow.control()] for the real pipeline run, and
#' also reused directly for panel A's re-derived gate boundary (rather than
#' recomputing it), so the figure shows the same gate the real run used.
#' @param af.remove Logical, default `TRUE`. Passed to [clean.controls()].
#'   Panels B and D require this to be `TRUE` and require the illustrated
#'   fluorophore to have a paired universal negative; otherwise those panels
#'   show a placeholder.
#' @param universal.negative,downsample,scatter.match,k.neighbors,negative.n,positive.n
#' Passed through to [clean.controls()]. See that function's documentation.
#' @param singlet.quantiles Numeric, default `c(0.85, 0.975)`. Quantile
#' thresholds for the two-stage FSC/SSC singlet discrimination used only when
#' cleaning a paired bead control (see `control.def.file`), matching
#' [get.spectra.automated()].
#' @param color.palette Optional character string defining the viridis color
#'   palette to be used for the fluorophore traces. Use `rainbow`
#'   to be similar to FlowJo or SpectroFlo. Other options are the viridis color
#'   options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`, `mako`
#'   and `turbo`.
#' @param gate.color Colour of the panel A gate boundary. Default
#'   `"darkgoldenrod1"` (matching `do.gate()`'s default).
#' @param density.palette Fill palette for the panel A pseudocolour density.
#'   Default `"rainbow"`.
#' @param unstained.point.color Colour for the unstained/AF events in panel
#'   B. Default `"black"`.
#' @param cosine.point.size Numeric or `NULL` (default). Point size for the
#'   single-stained control events in panel B. If `NULL`, defaults to
#'   `asp$figure.gate.point.size * 1.3`.
#' @param af.gate.color Colour of the AF-exclusion gate boundary drawn on
#'   both panel B biplots. Default `"black"`.
#' @param clean.positive.color Colour for the highlighted "true positive"
#'   events in panels A and D. Default `"red"`.
#' @param clean.positive.point.size Numeric or `NULL` (default). Point size
#'   for the panels A/D highlight. If `NULL`, defaults to
#'   `asp$figure.gate.point.size * 1.5`.
#' @param n.true.positive Integer, default `50L`. Number of "true positive"
#'   events highlighted in red in panels A and D: the brightest
#'   `n.true.positive` events among the AF-gate-excluded population, ranked
#'   by projection onto the fitted RLM trend direction in (peak channel,
#'   intrusive-AF channel) space, rather than every AF-gate-excluded event
#'   (which is simply "not AF", not "positively stained").
#' @param rlm.line.color Colour of the robust-linear-model fit line in panel
#' D. Default `"blue"`.
#' @param cells.trace.color,beads.trace.color,af.trace.color Colours for the
#' three traces in panel E: the RLM-based per-channel signature ("Cells"),
#' the reference profile ("Beads"), and the matched-negative AF trace.
#' Defaults `"#D95F02"` / `"#377EB8"` / `"grey40"`, matching
#' [spectra.automated.steps.plot()]'s panel F and
#' [spectra.standard.workflow.plot()]'s panel D.
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
#' @param verbose Logical, default `TRUE`. Print progress messages (also
#'   controls verbosity of the internal `define.flow.control()` /
#'   `clean.controls()` calls).
#' @param allow.duplicate.controls Logical, default `TRUE`. Set `TRUE` to
#'   permit multiple single-stained controls for the same fluorophore
#'   (diagnostic/QC use only). Each is tracked internally under a unique
#'   `sample` identifier. The resulting spectral reference library still needs
#'   to be reduced to one row per fluorophore before unmixing --
#'   see `check.spectra.duplicates()`.
#'
#' @return Invisibly, a named list (one entry per fluorophore), each
#'   containing:
#'   \describe{
#'     \item{`gate.panel`}{Panel A, the automated scatter gate (or a
#'       placeholder if gate definition failed), with the `n.true.positive`
#'       brightest-along-the-RLM-trend events highlighted larger in red when
#'       AF-removal diagnostics were available.}
#'     \item{`af.panel`}{Panel B, the AF-exclusion cosine-similarity biplot
#'       (or a placeholder if `af.remove = FALSE` or no paired universal
#'       negative was available for this fluorophore).}
#'     \item{`scatter.match.panel`}{Panel C, the [clean.controls()]
#'       kNN scatter-match figure embedded from its saved JPEG.}
#'     \item{`rlm.panel`}{Panel D, the robust-linear-model diagnostic, fit
#'       to and displaying `flow.control$clean.expr` for this sample (the
#'       events as clean.controls() actually finalises them, not just
#'       `gate.population.idx`), or a placeholder alongside `af.panel` when
#'       AF-removal diagnostics were unavailable or too few clean.controls()
#'       events remained.}
#'     \item{`subtraction.plot`}{Panel E, the final spectral profile
#'       comparison ([spectral.trace()] of Cells / Beads / AF), fit to the
#'       same `flow.control$clean.expr` population as `rlm.panel`, or a
#'       placeholder when AF-removal diagnostics were unavailable, too few
#'       clean.controls() events remained, or RLM extraction failed.}
#'     \item{`composite`}{The assembled five-panel cowplot object saved to
#'       `output.dir` when `save = TRUE`.}
#'     \item{`gate.name`}{Character. The `gate.name` resolved for this
#'       fluorophore's sample, or `NA` if none was assigned.}
#'     \item{`af.peak.channel`}{Character. The intrusive-AF channel used as
#'       panels B/D's y-axis, or `NA_character_` if AF-removal diagnostics
#'       were unavailable.}
#'     \item{`fluor.peak`}{Character. The fluorophore's peak channel used as
#'       panels B/D's x-axis, or `NA_character_` if AF-removal diagnostics
#'       were unavailable.}
#'     \item{`reference.profile`}{Named numeric vector (over the panel-wide
#'       spectral channels) used as the "Beads" trace in panel E, or `NULL`
#'       if neither a paired bead control nor the spectral reference library
#'       had data for this fluorophore.}
#'   }
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 geom_path geom_abline annotate labs theme theme_bw
#' @importFrom ggplot2 scale_color_manual
#' @importFrom scattermore geom_scattermore
#' @importFrom cowplot plot_grid get_legend
#' @importFrom stats cor median setNames
#'
#' @seealso [define.flow.control()], [clean.controls()],
#'   [get.fluorophore.spectra()], [spectra.automated.steps.plot()],
#'   [spectra.standard.workflow.plot()]
#'
#' @export

spectra.legacy.steps.plot <- function(
    control.dir,
    control.def.file,
    asp,
    fluorophores               = NULL,
    gating.system              = c( "density", "landmarks" ),
    gate.list                  = NULL,
    af.remove                  = TRUE,
    universal.negative         = TRUE,
    downsample                 = TRUE,
    scatter.match              = TRUE,
    k.neighbors                = 3L,
    negative.n                 = asp$negative.n,
    positive.n                 = asp$positive.n,
    singlet.quantiles           = c( 0.85, 0.975 ),
    color.palette              = NULL,
    gate.color                 = "darkgoldenrod1",
    density.palette             = "rainbow",
    unstained.point.color       = "black",
    cosine.point.size           = NULL,
    af.gate.color                = "black",
    clean.positive.color         = "red",
    clean.positive.point.size    = NULL,
    n.true.positive               = 50L,
    rlm.line.color                = "blue",
    cells.trace.color            = "#D95F02",
    beads.trace.color            = "#377EB8",
    af.trace.color                = "grey40",
    max.points                  = 5e4,
    panel.width                 = 4,
    panel.height                = 4,
    composite.width             = NULL,
    composite.height            = NULL,
    output.dir                  = NULL,
    save                        = TRUE,
    file.type                   = "jpg",
    verbose                     = TRUE,
    allow.duplicate.controls    = TRUE
) {

  # -- 0. Validate inputs
  if ( !dir.exists( control.dir ) )
    stop( "control.dir does not exist: ", control.dir, call. = FALSE )

  gating.system <- match.arg( gating.system )

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

  clean.positive.point.size.use <- if ( !is.null( clean.positive.point.size ) )
    clean.positive.point.size else asp$figure.gate.point.size * 1.5

  # -- 1. Run the full legacy pipeline (self-contained: this function calls
  # define.flow.control() and clean.controls() itself, for the whole control
  # set, not just the illustrated fluorophore(s))
  if ( verbose ) message( "\033[34m-- Running define.flow.control() --\033[0m" )
  flow.control <- define.flow.control(
    control.dir = control.dir,
    control.def.file = control.def.file,
    asp = asp,
    gate = TRUE,
    gating.system = gating.system,
    gate.list = gate.list,
    parallel = FALSE,
    verbose = verbose,
    color.palette = color.palette,
    allow.duplicate.controls = allow.duplicate.controls
  )

  diagnostics.env <- new.env( parent = emptyenv() )

  if ( verbose ) message( "\033[34m-- Running clean.controls() --\033[0m" )
  flow.control <- clean.controls(
    flow.control, asp,
    af.remove = af.remove, universal.negative = universal.negative,
    downsample = downsample, negative.n = negative.n, positive.n = positive.n,
    scatter.match = scatter.match, k.neighbors = k.neighbors,
    intermediate.figures = FALSE, main.figures = TRUE, parallel = FALSE,
    verbose = verbose, diagnostics.env = diagnostics.env
  )

  # -- 2. Resolve requested fluorophore(s) (cell-based, non-AF, non-negative,
  # with a paired universal negative, matching define.flow.control()'s own
  # negative-sample pattern)
  neg.pattern <- "(^AF$)|(\\bAF\\b)|(\\bnegative\\b)|(\\bunstained\\b)"
  fluor.rows.all <- which(
    flow.control$control.type == "cells" &
      !grepl( neg.pattern, flow.control$fluorophore, ignore.case = TRUE )
  )

  if ( length( fluor.rows.all ) == 0 )
    stop( "No cell-based fluorophore rows found in ", control.def.file, call. = FALSE )

  all.fluor.names <- flow.control$fluorophore[ fluor.rows.all ]

  if ( is.null( fluorophores ) ) {
    fluorophores <- all.fluor.names[ 1L ]
    if ( verbose )
      message( "No `fluorophores` specified; illustrating: ", fluorophores )
  }

  missing.fluor <- setdiff( fluorophores, all.fluor.names )
  if ( length( missing.fluor ) > 0 )
    stop(
      "Fluorophore(s) not found among cell-based controls: ",
      paste( missing.fluor, collapse = ", " ), call. = FALSE
    )

  # -- 3. Rebuild the assign.gates()-processed control table, purely to look
  # up each fluorophore's gate.name for panel A (define.flow.control() does
  # not return this internal object)
  ctrl.tbl.raw <- utils::read.csv(
    if ( file.exists( control.def.file ) ) control.def.file
    else file.path( control.dir, control.def.file ),
    stringsAsFactors = FALSE, strip.white = TRUE
  )
  ctrl.tbl.raw[] <- lapply( ctrl.tbl.raw, function( x ) {
    if ( is.character( x ) ) { x <- trimws( x ); x[ x == "" ] <- NA; x } else x
  } )
  ctrl.tbl.raw$sample <- .build.control.sample.names(
    ctrl.tbl.raw$fluorophore, ctrl.tbl.raw$control.type, ctrl.tbl.raw$marker
  )
  ctrl.tbl.gated <- assign.gates( ctrl.tbl.raw, gating.system, gate = TRUE, verbose = FALSE )

  scatter.channels <- read.scatter.parameter( asp )
  fsc.a <- scatter.channels[ 1L ]
  ssc.a <- scatter.channels[ 2L ]
  sat.value <- if ( !is.null( asp$expr.data.max ) ) asp$expr.data.max else Inf
  spectral.channels <- flow.control$spectral.channel
  db.col <- .cytometer.to.db.col( asp$cytometer )

  gate.cache <- new.env( parent = emptyenv() )

  composite.width.use  <- if ( !is.null( composite.width ) )  composite.width  else panel.width * 2
  row.heights           <- rep( panel.height, 5L )
  composite.height.use  <- if ( !is.null( composite.height ) ) composite.height else sum( row.heights )

  results <- list()

  # -- 4. Loop over requested fluorophores
  for ( fluor in fluorophores ) {

    if ( verbose )
      message( "\033[34m-- Building legacy pipeline-step figure for ", fluor, " --\033[0m" )

    samp.i      <- flow.control$sample[
      which( flow.control$fluorophore == fluor & flow.control$control.type == "cells" )[ 1L ]
    ]
    fluor.file  <- flow.control$filename[ samp.i ]
    gate.name.i <- ctrl.tbl.gated$gate.name[ ctrl.tbl.gated$sample == samp.i ][ 1L ]

    ref.profile.i <- .get.reference.profile(
      fluor, ctrl.tbl.raw, control.dir, spectral.channels, scatter.channels,
      sat.value, singlet.quantiles, asp, db.col, verbose = FALSE
    )

    # -- A. Automated scatter gate (re-derived; see .resolve.legacy.gate())
    fluor.file.path <- file.path( control.dir, fluor.file )
    probe.cols.i    <- colnames( readFCS( fluor.file.path, start.row = 1, end.row = 1 ) )
    cols.keep.i     <- intersect( c( scatter.channels, flow.control$spectral.channel ), probe.cols.i )
    raw.mat.i       <- readFCS( fluor.file.path, columns = cols.keep.i )
    spec.present.i  <- intersect( flow.control$spectral.channel, colnames( raw.mat.i ) )
    if ( length( spec.present.i ) > 0 )
      raw.mat.i <- raw.mat.i[
        rowSums( raw.mat.i[ , spec.present.i, drop = FALSE ] >= sat.value ) == 0, , drop = FALSE
      ]

    gate.poly.i <- if ( is.na( gate.name.i ) ) NULL else
      .resolve.legacy.gate(
        gate.name.i, ctrl.tbl.gated, control.dir, asp, gating.system, gate.list, gate.cache
      )

    gate.panel <- if ( is.null( gate.poly.i ) ) {
      .embed.or.placeholder( NULL, label = "Automated scatter gate (gate definition failed)" )
    } else {
      .octagon.gate.panel(
        raw.mat.i[ , fsc.a ], raw.mat.i[ , ssc.a ], gate.poly.i, asp, fsc.a, ssc.a,
        gate.color = gate.color, density.palette = density.palette, max.points = max.points
      )
    }

    # -- B & D. Autofluorescence-removal diagnostics (requires af.remove =
    # TRUE and a paired universal negative for this sample)
    diag.i <- diagnostics.env[[ samp.i ]]

    if ( is.null( diag.i ) ) {
      af.panel <- .embed.or.placeholder(
        NULL, label = paste(
          "AF-exclusion gate: diagnostics not available for", fluor,
          "(af.remove may be FALSE, or this sample has no paired universal negative)"
        )
      )
      rlm.panel <- .embed.or.placeholder(
        NULL, label = "Robust linear model (requires AF-removal diagnostics)"
      )
      subtraction.plot <- .embed.or.placeholder(
        NULL, label = "Final spectral profile (requires AF-removal diagnostics)"
      )
    } else {
      neg.median.i <- apply( diag.i$expr.data.neg, 2, stats::median )
      cs.vals.i    <- .cosine.sim.rows( diag.i$expr.data.pos, neg.median.i )

      af.panel <- .make.af.exclusion.biplot(
        unstained.mat = diag.i$expr.data.neg, stained.mat = diag.i$expr.data.pos,
        x.dim = diag.i$af.peak.channel, y.dim = diag.i$fluor.peak,
        cs.vals = cs.vals.i, af.boundary = diag.i$af.boundaries$upper, asp = asp,
        unstained.point.color = unstained.point.color, max.points = max.points,
        cosine.point.size = cosine.point.size, af.gate.color = af.gate.color
      )

      # x.vals.i/y.vals.i and the fit below feed *only* panel A's
      # true-positive ranking (n.true.positive) -- panel D's own biplot uses
      # a separate, later computation on flow.control$clean.expr, since
      # gate.population.idx predates clean.controls()'s later
      # scatter-matching/brightest-selection step
      x.vals.i <- diag.i$expr.data.pos[ , diag.i$fluor.peak ]
      y.vals.i <- diag.i$expr.data.pos[ , diag.i$af.peak.channel ]

      if ( length( diag.i$gate.population.idx ) >= 2L ) {
        x.fit.i <- x.vals.i[ diag.i$gate.population.idx ]
        y.fit.i <- y.vals.i[ diag.i$gate.population.idx ]
      } else {
        warning(
          "Too few AF-gate-excluded events for '", fluor, "'; ranking panel ",
          "A's true positives against the full population instead.", call. = FALSE
        )
        x.fit.i <- x.vals.i
        y.fit.i <- y.vals.i
      }

      rlm.coef.i <- fit.robust.linear.model(
        x.fit.i, y.fit.i, diag.i$fluor.peak, diag.i$af.peak.channel
      )

      # "true positive" events (red in panel A) are the brightest
      # n.true.positive events among the AF-gate-excluded population, ranked
      # by their projection onto the fitted RLM trend direction (1, slope)
      # in (peak channel, intrusive-AF channel) space -- i.e. position along
      # the fitted line, not just raw peak-channel brightness. This is a
      # much smaller, tighter population than "AF-gate-excluded", which is
      # simply "not AF" rather than "positively stained".
      trend.proj.i         <- x.fit.i + rlm.coef.i[ 2L ] * y.fit.i
      n.true.positive.i    <- min( n.true.positive, length( trend.proj.i ) )
      top.trend.i          <- order( trend.proj.i, decreasing = TRUE )[ seq_len( n.true.positive.i ) ]
      true.positive.idx.i  <- diag.i$gate.population.idx[ top.trend.i ]

      # true-positive highlight for panel A: plotted at the real FSC/SSC
      # coordinates from remove.af()'s `scatter.data.pos` capture
      # (independent of gate.panel's own raw.mat.i read -- both are just
      # measured detector values from the same file, so no row
      # correspondence between the two is needed). Mirrors
      # spectra.standard.workflow.plot()'s `ground.truth.method = "legacy"`
      # branch.
      if ( is.null( diag.i$scatter.data.pos ) ) {
        warning(
          "diagnostics.env for '", fluor, "' has no scatter.data.pos -- ",
          "requires the updated remove.af() that stores scatter columns.",
          call. = FALSE
        )
      } else if ( all( c( fsc.a, ssc.a ) %in% colnames( diag.i$scatter.data.pos ) ) ) {
        true.pos.scatter.i <- diag.i$scatter.data.pos[
          true.positive.idx.i, , drop = FALSE
        ]
        gate.panel <- .add.highlight.layer(
          gate.panel, true.pos.scatter.i[ , fsc.a ], true.pos.scatter.i[ , ssc.a ],
          color = clean.positive.color, pointsize = clean.positive.point.size.use
        )
      }

      # panel D: the events exactly as clean.controls() finalises them for
      # this sample (post remove.af(), universal-negative scatter-matching,
      # and brightest-event selection) -- flow.control$clean.expr is the
      # actual population get.fluorophore.spectra() fits its RLM on
      final.clean.mat.i <- flow.control$clean.expr[
        flow.control$clean.event.sample == samp.i, , drop = FALSE
      ]

      if ( nrow( final.clean.mat.i ) < 2L ) {
        rlm.panel <- .embed.or.placeholder(
          NULL, label = paste(
            "Robust linear model: too few clean.controls() events for", fluor
          )
        )
      } else {
        x.clean.i <- final.clean.mat.i[ , diag.i$fluor.peak ]
        y.clean.i <- final.clean.mat.i[ , diag.i$af.peak.channel ]

        rlm.coef.panel.i <- fit.robust.linear.model(
          x.clean.i, y.clean.i, diag.i$fluor.peak, diag.i$af.peak.channel
        )
        r.squared.panel.i <- stats::cor( x.clean.i, y.clean.i ) ^ 2

        rlm.panel <- .make.rlm.biplot(
          x.vals = x.clean.i, y.vals = y.clean.i,
          x.lab = diag.i$fluor.peak, y.lab = diag.i$af.peak.channel,
          clean.idx = seq_along( x.clean.i ),
          rlm.coef = rlm.coef.panel.i, r.squared = r.squared.panel.i, asp = asp,
          clean.positive.color = clean.positive.color,
          clean.positive.point.size = clean.positive.point.size.use,
          rlm.line.color = rlm.line.color
        )
      }

      # -- E. Final spectral profile comparison: the RLM-based per-channel
      # signature computed the same way get.fluorophore.spectra() computes
      # it (Cells), the reference profile for the same fluorophore (Beads --
      # a paired bead control if provided, otherwise the static spectral
      # reference library), and the matched-negative AF trace. Reuses
      # final.clean.mat.i (flow.control$clean.expr for this sample) computed
      # above for panel D, since it's the same population
      # get.fluorophore.spectra() actually fits on.
      clean.mat.i <- final.clean.mat.i

      if ( nrow( clean.mat.i ) < 3L ) {
        subtraction.plot <- .embed.or.placeholder(
          NULL, label = paste(
            "Final spectral profile: too few AF-gate-excluded events for", fluor
          )
        )
      } else {
        cells.spectrum.i <- tryCatch(
          .legacy.rlm.spectrum( clean.mat.i, diag.i$fluor.peak, spectral.channels ),
          error = function( e ) {
            warning(
              "Robust linear model spectrum extraction failed for '", fluor,
              "': ", e$message, call. = FALSE
            )
            NULL
          }
        )

        if ( is.null( cells.spectrum.i ) ) {
          subtraction.plot <- .embed.or.placeholder(
            NULL, label = paste( "Final spectral profile: extraction failed for", fluor )
          )
        } else {
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

          af.mean.i <- colMeans( diag.i$expr.data.neg )[ spectral.channels ]

          subtraction.mat <- rbind(
            cells.spectrum.i,
            beads.trace,
            af.mean.i / max( abs( af.mean.i ), 1e-9 )
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
    }

    # -- C. Scatter-matching (embed the JPEG clean.controls() already saved)
    scatter.match.file <- file.path(
      asp$figure.scatter.dir.base, paste( samp.i, asp$scatter.match.plot.filename, sep = "_" )
    )
    scatter.match.panel <- .embed.or.placeholder(
      scatter.match.file, label = "kNN scatter-matched universal negative (clean.controls())"
    )

    # panels A and D are forced square (aspect.ratio = 1); their row cell is
    # composite.width.use wide by their own row.heights entry tall, which is
    # usually wider than tall, so ggplot's own aspect-ratio padding would
    # centre the square in blank space. Instead, size the column to just fit
    # the square and pad only on the right, so the panel stays left-aligned
    # (matching spectra.standard.workflow.plot()'s panel A).
    gate.square.side <- min( composite.width.use, row.heights[ 1L ] )
    gate.row <- if ( gate.square.side < composite.width.use ) {
      cowplot::plot_grid(
        gate.panel, NULL, ncol = 2,
        rel_widths = c( gate.square.side, composite.width.use - gate.square.side )
      )
    } else {
      gate.panel
    }

    rlm.square.side <- min( composite.width.use, row.heights[ 4L ] )
    rlm.row <- if ( rlm.square.side < composite.width.use ) {
      cowplot::plot_grid(
        rlm.panel, NULL, ncol = 2,
        rel_widths = c( rlm.square.side, composite.width.use - rlm.square.side )
      )
    } else {
      rlm.panel
    }

    # -- Assemble composite figure
    composite <- cowplot::plot_grid(
      cowplot::plot_grid( gate.row,            labels = "A" ),
      cowplot::plot_grid( af.panel,             labels = "B" ),
      cowplot::plot_grid( scatter.match.panel,  labels = "C" ),
      cowplot::plot_grid( rlm.row,              labels = "D" ),
      cowplot::plot_grid( subtraction.plot,     labels = "E" ),
      ncol = 1, rel_heights = row.heights
    )

    if ( save ) {
      out.file <- file.path(
        output.dir, sprintf( "%s_legacy_pipeline_steps.%s", fluor, file.type )
      )
      ggplot2::ggsave(
        out.file, plot = composite, device = plot.device,
        width = composite.width.use, height = composite.height.use, limitsize = FALSE
      )
      if ( verbose ) message( "\033[32m  Saved: ", out.file, "\033[0m" )
    }

    results[[ fluor ]] <- list(
      gate.panel           = gate.panel,
      af.panel              = af.panel,
      scatter.match.panel   = scatter.match.panel,
      rlm.panel             = rlm.panel,
      subtraction.plot      = subtraction.plot,
      composite             = composite,
      gate.name             = gate.name.i,
      af.peak.channel        = if ( !is.null( diag.i ) ) diag.i$af.peak.channel else NA_character_,
      fluor.peak             = if ( !is.null( diag.i ) ) diag.i$fluor.peak else NA_character_,
      reference.profile      = ref.profile.i
    )
  }

  invisible( results )
}
