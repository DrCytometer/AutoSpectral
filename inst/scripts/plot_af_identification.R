# plot_af_identification.R
#
# Manuscript figure helper: illustrates the principles of single-cell
# autofluorescence (AF) identification and extraction (Supplementary Figure 1
# of the AutoSpectral manuscript), using the same building blocks as the rest
# of the package (get.af.spectra(), create.biplot(), spectral.trace(),
# unmix.autospectral(), unmix.autospectral.joint()) so the panels are
# stylistically consistent with the package's other manuscript figures.
# Self-contained -- unlike plot_spectra_standard_workflow.R, it does not
# depend on private helpers from plot_spectra_automated_steps.R.
#
# Panels produced:
#   A. Autofluorescence spectral trace from get.af.spectra(), run on the
#      unstained sample (spectral.trace()).
#   B. Raw (no-AF-extraction) biplot of the unstained sample, coloured by
#      per-cell AF Index assignment (discrete "Set2"-style palette).
#   C. Three create.biplot()-style panels (AF abundance vs. panel.c.y.dim) of
#      the unstained sample unmixed via unmix.autospectral() using residual
#      minimization, fluorophore-signal minimization, and
#      unmix.autospectral.joint() ("residuals" / "fluorophores" / "joint").
#   D. Cosine-similarity distributions: the "random assignment" baseline
#      (pairwise similarity among the AF spectra themselves) alongside the
#      per-cell similarity between each cell's raw spectrum and its assigned
#      AF spectrum for each of the three assignment methods in C, each
#      Kolmogorov-Smirnov-tested against the random baseline.
#   E. AF prediction error: raw detector-space reconstruction error (summed
#      absolute residual) using a single population-median AF signature vs.
#      per-cell AF extraction, with and without a synthetic fluorophore
#      overlay added to the unstained data. Predates sim.flow.data() -- see
#      the note in the function body.
#   F. Fluorophore recovery error: |true - unmixed| for the known synthetic
#      fluorophore signal, single-median AF vs. per-cell AF.
#   G. Stability of the per-cell AF assignment under a synthetic fluorophore
#      overlay (cosine similarity of the assigned AF spectrum, real vs.
#      synthetic; the "no synthetic" condition is definitionally 1 for every
#      cell and is shown as a reference line rather than a density).
#   H. Side-scatter (or other named scatter channel) distribution by AF
#      Index.
#
# Panel "L" of the original figure (a cell-type-by-AF-Index breakdown
# produced in FlowJo) is out of scope for this function.

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

## Formats a p-value for on-plot annotation, consistent across all panels.
.fmt.p <- function( p ) {
  if ( !is.finite( p ) ) return( "p = NA" )
  if ( p < 1e-4 ) "p < 0.0001" else sprintf( "p = %.4f", p )
}

## Corner-of-panel p-value annotation for a two-sample KS test of `x` against
## `reference`.
.ks.p.annotation <- function( x, reference ) {
  p <- stats::ks.test( x, reference )$p.value
  ggplot2::annotate(
    "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
    label = .fmt.p( p ), size = 3.2
  )
}

## Per-row cosine similarity between `mat` and the row of `ref.mat` indexed
## by the corresponding element of `idx` -- e.g. mat = raw cell spectra,
## ref.mat = af.spectra, idx = the per-cell "AF Index" column. A vectorised
## replacement for a per-cell sapply() loop; also usable with mat itself
## being an indexed subset of ref.mat (see panel G).
.cosine.sim.paired <- function( mat, ref.mat, idx ) {
  ref   <- ref.mat[ idx, , drop = FALSE ]
  dot   <- rowSums( mat * ref )
  mat.n <- sqrt( rowSums( mat^2 ) )
  ref.n <- sqrt( rowSums( ref^2 ) )
  dot / ( mat.n * ref.n + 1e-9 )
}

## Single 1D density panel for a cosine-similarity distribution, with an
## optional KS-test p-value (vs `reference`) annotated in the corner. Falls
## back to a reference vline when `x` is (near-)constant, since
## ggplot2::geom_density() requires a strictly positive bandwidth.
.cosine.density.panel <- function(
    x, title, reference = NULL, fill.color = "lightblue"
) {
  x <- x[ is.finite( x ) ]
  degenerate <- length( x ) < 2 || stats::sd( x ) < 1e-8
  
  p <- if ( degenerate ) {
    ggplot2::ggplot( data.frame( x = x ), ggplot2::aes( x ) ) +
      ggplot2::geom_vline(
        xintercept = if ( length( x ) > 0 ) mean( x ) else NA_real_,
        linetype = "dashed", color = "grey30", linewidth = 0.6
      )
  } else {
    ggplot2::ggplot( data.frame( x = x ), ggplot2::aes( x ) ) +
      ggplot2::geom_density( fill = fill.color, alpha = 0.5, na.rm = TRUE )
  }
  
  p <- p +
    ggplot2::labs( x = "Cosine similarity", y = "Density", title = title ) +
    ggplot2::xlim( c( 0, 1 ) ) +
    ggplot2::theme_minimal( base_size = 11 ) +
    ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ) )
  
  if ( !degenerate && !is.null( reference ) )
    p <- p + .ks.p.annotation( x, reference )
  
  p
}

## Violin + jittered points for an "Error by Group" comparison, matching the
## "each cell plotted as a point with a violin overlay" convention from the
## original figure. `p.labels`, if supplied, is a data.frame(x, y, label).
.error.violin.panel <- function(
    df, y.lab, title, p.labels = NULL, point.size = 0.3, point.alpha = 0.2
) {
  p <- ggplot2::ggplot( df, ggplot2::aes( Group, Error, fill = Group ) ) +
    ggplot2::geom_violin( alpha = 0.5, na.rm = TRUE, scale = "width" ) +
    ggplot2::geom_jitter(
      width = 0.15, size = point.size, alpha = point.alpha, na.rm = TRUE
    ) +
    ggplot2::labs( x = NULL, y = y.lab, title = title ) +
    ggplot2::theme_minimal( base_size = 11 ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 )
    )
  
  if ( !is.null( p.labels ) )
    p <- p + ggplot2::annotate(
      "text", x = p.labels$x, y = p.labels$y, label = p.labels$label, size = 3.2
    )
  
  p
}

## Discrete-colour biplot matching create.biplot()'s axis/transform
## conventions, but colouring points by a categorical AF Index rather than a
## continuous pseudocolour density. The 8-colour ColorBrewer "Set2" palette
## is hardcoded (matching the RdYlBu hardcoding in get_spectra_automated.R,
## so RColorBrewer doesn't need to become a package dependency) and extended
## by interpolation when more than 8 AF spectra are present, since
## get.af.spectra() commonly returns more than 8 rows.
.af.index.biplot <- function(
    plot.data, x.dim, y.dim, af.index, asp,
    x.lab = NULL, y.lab = NULL,
    x.min = -5000, x.max = asp$expr.data.max,
    y.min = -5000, y.max = asp$expr.data.max,
    x.width.basis = -1000, y.width.basis = -1000,
    max.points = 5e4, color.palette = "Set2",
    point.size = NULL
) {
  
  if ( is.null( x.lab ) ) x.lab <- x.dim
  if ( is.null( y.lab ) ) y.lab <- y.dim
  
  df <- data.frame(
    x = plot.data[ , x.dim ],
    y = plot.data[ , y.dim ],
    z = factor( af.index, levels = sort( unique( af.index ) ) )
  )
  
  if ( nrow( df ) > max.points ) {
    set.seed( asp$bird.seed )
    df <- df[ sample( seq_len( nrow( df ) ), max.points ), , drop = FALSE ]
  }
  
  x.breaks <- asp$ribbon.breaks[ asp$ribbon.breaks < x.max ]
  y.breaks <- asp$ribbon.breaks[ asp$ribbon.breaks < y.max ]
  x.axis.labels <- sapply( x.breaks, function( v ) {
    if ( v == 0 ) "0" else parse( text = paste0( "10^", log10( abs( v ) ) ) )
  } )
  y.axis.labels <- sapply( y.breaks, function( v ) {
    if ( v == 0 ) "0" else parse( text = paste0( "10^", log10( abs( v ) ) ) )
  } )
  
  biexp.trans.x <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = x.max, pos = log10( x.max ) - 1,
    neg = asp$default.transformation.param$neg,
    widthBasis = x.width.basis, inverse = FALSE
  )
  biexp.trans.y <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = y.max, pos = log10( y.max ) - 1,
    neg = asp$default.transformation.param$neg,
    widthBasis = y.width.basis, inverse = FALSE
  )
  
  df$x.trans <- biexp.trans.x( df$x )
  df$y.trans <- biexp.trans.y( df$y )
  
  set2.base <- c( "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3" )
  n.levels  <- nlevels( df$z )
  pal       <- if ( n.levels <= 8 ) set2.base[ seq_len( n.levels ) ] else
    grDevices::colorRampPalette( set2.base )( n.levels )
  
  pt.size <- if ( !is.null( point.size ) ) point.size else asp$figure.gate.point.size
  
  ggplot2::ggplot( df, ggplot2::aes( x.trans, y.trans, color = z ) ) +
    scattermore::geom_scattermore( pointsize = pt.size, alpha = 1, na.rm = TRUE ) +
    ggplot2::scale_color_manual( values = pal ) +
    ggplot2::scale_x_continuous(
      name = x.lab, breaks = biexp.trans.x( x.breaks ),
      limits = biexp.trans.x( c( x.min, x.max ) ), labels = x.axis.labels
    ) +
    ggplot2::scale_y_continuous(
      name = y.lab, breaks = biexp.trans.y( y.breaks ),
      limits = biexp.trans.y( c( y.min, y.max ) ), labels = y.axis.labels
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      legend.position = "none",
      axis.ticks = ggplot2::element_line( linewidth = asp$figure.panel.line.size ),
      axis.text = ggplot2::element_text( size = asp$figure.axis.text.size ),
      axis.title = ggplot2::element_text( size = asp$figure.axis.title.size ),
      panel.border = ggplot2::element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## Reverts a per-cell AutoSpectral unmixing result (as returned by
## unmix.autospectral() / unmix.autospectral.joint()) back to raw detector
## space: fluorophore abundances projected through `spectra`, plus each
## cell's own assigned AF spectrum (af.spectra["AF Index",]) scaled by its
## "AF" abundance.
.predict.raw.percell <- function( unmixed, spectra, af.spectra ) {
  fluor.names <- rownames( spectra )
  pred.fluor  <- unmixed[ , fluor.names, drop = FALSE ] %*% spectra
  af.idx      <- unmixed[ , "AF Index" ]
  pred.af     <- sweep( af.spectra[ af.idx, , drop = FALSE ], 1, unmixed[ , "AF" ], "*" )
  pred.fluor + pred.af
}


# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' @title Plot Autofluorescence Identification Workflow
#'
#' @description
#' Builds a manuscript-ready, multi-panel figure illustrating how AutoSpectral
#' identifies and extracts per-cell autofluorescence (AF) from an unstained
#' sample: the AF spectral library itself, raw and unmixed biplots of the
#' AF Index assignment, cosine-similarity validation of that assignment
#' against a random-assignment baseline, and AF/fluorophore prediction-error
#' comparisons between a single population-level AF signature and per-cell AF
#' extraction, with and without a synthetic fluorophore overlay. Uses the same
#' building blocks as the rest of AutoSpectral ([get.af.spectra()],
#' [create.biplot()], [spectral.trace()], [unmix.autospectral()],
#' [unmix.autospectral.joint()]) so panel styling matches the package's other
#' figures. Self-contained; does not require any other manuscript-figure
#' helper file to be loaded.
#'
#' @param unstained.sample Character path to the unstained sample FCS file.
#'   Used as a label only when `unstained.exprs` is supplied directly.
#' @param spectra Spectral signatures of fluorophores (the acquisition panel),
#'   normalized between 0 and 1, fluorophores in rows and detectors in
#'   columns. Must not contain an `"AF"` row.
#' @param asp The AutoSpectral parameter list from `get.autospectral.param()`.
#' @param unstained.exprs Optional pre-loaded raw expression matrix (cells x
#'   detectors), bypassing `readFCS()`. Must contain every column of `spectra`
#'   plus, unless `unstained.scatter` is supplied, `ssc.channel`. Default
#'   `NULL` reads from `unstained.sample`.
#' @param unstained.scatter Optional numeric vector, parallel to
#'   `unstained.exprs`, of scatter values for panel H. Only used when
#'   `unstained.exprs` is supplied directly and does not itself contain
#'   `ssc.channel`. Default `NULL`.
#' @param n.cells Integer, default `10000L`. Number of events sampled from the
#'   unstained sample for the whole figure.
#' @param sample.method Character, one of `"random"` (default) or `"first"`.
#'   How `n.cells` events are chosen.
#' @param af.spectra Optional precomputed AF spectral library (as returned by
#'   `get.af.spectra()`). Default `NULL` computes it via `get.af.spectra()`
#'   on the same `n.cells` subsample used for the rest of the figure, which
#'   can be slow -- supply a precomputed library to skip re-running it.
#' @param get.af.spectra.args Named list of additional arguments forwarded to
#'   `get.af.spectra()` (e.g. `refine`, `deduplicate`, `som.dim`) when
#'   `af.spectra` is `NULL`. `figures` and `save` default to `FALSE` here
#'   (override via this list if the CSV/plot side-effects are wanted).
#' @param panel.b.x.dim Character, default `"UV9-A"`. x-axis channel for
#'   panel B.
#' @param panel.b.y.dim Character or `NULL` (default). y-axis channel for
#'   panel B; `NULL` uses `asp$af.channel`.
#' @param panel.b.color.palette Character, default `"Set2"`. Currently only
#'   `"Set2"` is implemented (see `.af.index.biplot()`).
#' @param panel.c.x.dim,panel.c.y.dim Character, defaults `"AF"` and
#'   `"BUV615"`. Axes for the three panel-C biplots; must each be either
#'   `"AF"` or a row name of `spectra`.
#' @param ssc.channel Character or `NULL` (default). Scatter channel plotted
#'   in panel H. `NULL` uses the second entry of `read.scatter.parameter(asp)`
#'   (typically `"SSC-A"`); override to e.g. `"SSC-B-A"` for cytometers with
#'   multiple side-scatter detectors.
#' @param synthetic.fluorophores Character vector of fluorophore names (must
#'   be row names of `spectra`) used to build the synthetic-overlay data for
#'   panels E-G. Default `NULL` samples `n.synthetic.fluors` at random.
#' @param n.synthetic.fluors Integer, default `4L`. Only used when
#'   `synthetic.fluorophores` is `NULL`.
#' @param synthetic.meanlog,synthetic.sdlog Numeric, defaults `8` and `0.5`.
#'   Log-normal parameters for the per-cell synthetic fluorophore intensity
#'   (`stats::rlnorm()`).
#' @param max.points Integer, default `5e4`. Point-count downsample threshold
#'   for the biplot panels (B, C).
#' @param parallel Logical, default `TRUE`. Passed to `unmix.autospectral()` /
#'   `unmix.autospectral.joint()`.
#' @param threads Numeric or `NULL` (default). Passed to `unmix.autospectral()`
#'   / `unmix.autospectral.joint()`; `NULL` uses `asp$worker.process.n`.
#' @param panel.width,panel.height Numeric, defaults `4` and `4`. Per-panel
#'   sizing unit used to compute the composite figure dimensions.
#' @param composite.width,composite.height Numeric or `NULL` (default).
#'   Override the overall saved figure dimensions (inches); if `NULL`, these
#'   are computed from `panel.width` / `panel.height`.
#' @param output.dir Character or `NULL` (default). Directory to save the
#'   composite figure. Defaults to the current working directory.
#' @param save Logical, default `TRUE`. Whether to save the composite figure.
#' @param file.type Character string, one of `"jpg"` (default), `"tiff"`,
#'   `"png"`, or `"pdf"` -- use `"pdf"` for the manuscript.
#' @param verbose Logical, default `TRUE`. Print progress messages.
#' @param seed Integer or `NULL` (default). Random seed used for event
#'   subsampling and synthetic-data generation; `NULL` uses `asp$bird.seed`.
#'
#' @return Invisibly, a named list:
#'   \describe{
#'     \item{`af.spectra`}{The AF spectral library used throughout (as
#'       supplied, or as computed by `get.af.spectra()`).}
#'     \item{`unmixed`}{List of the three per-cell AutoSpectral unmixings of
#'       the unstained sample: `residuals`, `fluorophores`, `joint`.}
#'     \item{`synthetic`}{List describing the synthetic-overlay data used in
#'       panels E-G: `data`, `groups`, `intensities`, `unmixed.median`,
#'       `unmixed.percell`.}
#'     \item{`panels`}{List of the individual ggplot objects, named `a`
#'       through `h`.}
#'     \item{`composite`}{The assembled eight-row cowplot object saved to
#'       `output.dir` when `save = TRUE`.}
#'     \item{`panel.e.data`,`panel.f.data`}{The per-cell error data frames
#'       underlying panels E and F.}
#'   }
#'
#' @importFrom ggplot2 ggplot aes geom_density geom_vline geom_violin
#' @importFrom ggplot2 geom_jitter labs xlim theme_minimal theme element_text
#' @importFrom ggplot2 annotate scale_color_manual scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous theme_bw margin element_line
#' @importFrom ggplot2 element_rect element_blank ggsave
#' @importFrom scattermore geom_scattermore
#' @importFrom cowplot plot_grid
#' @importFrom stats ks.test wilcox.test rlnorm median sd
#' @importFrom ragg agg_jpeg agg_tiff agg_png
#'
#' @seealso [get.af.spectra()], [unmix.autospectral()],
#'   [unmix.autospectral.joint()], [create.biplot()], [spectral.trace()],
#'   [sim.flow.data()]
#'
#' @export

af.identification.plot <- function(
    unstained.sample,
    spectra,
    asp,
    unstained.exprs         = NULL,
    unstained.scatter       = NULL,
    n.cells                  = 10000L,
    sample.method             = c( "random", "first" ),
    af.spectra                = NULL,
    get.af.spectra.args       = list(),
    panel.b.x.dim              = "UV9-A",
    panel.b.y.dim              = NULL,
    panel.b.color.palette      = "Set2",
    panel.c.x.dim              = "AF",
    panel.c.y.dim              = "BUV615",
    ssc.channel                = NULL,
    synthetic.fluorophores     = NULL,
    n.synthetic.fluors          = 4L,
    synthetic.meanlog            = 8,
    synthetic.sdlog               = 0.5,
    max.points                = 5e4,
    parallel                  = TRUE,
    threads                   = NULL,
    panel.width               = 4,
    panel.height              = 4,
    composite.width           = NULL,
    composite.height          = NULL,
    output.dir                 = NULL,
    save                       = TRUE,
    file.type                  = "jpg",
    verbose                    = TRUE,
    seed                       = NULL
) {
  
  # -- 0. Validate inputs
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
  
  sample.method <- match.arg( sample.method, c( "random", "first" ) )
  if ( is.null( seed ) ) seed <- asp$bird.seed
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]
  
  if ( is.null( panel.b.y.dim ) ) panel.b.y.dim <- asp$af.channel
  if ( is.null( panel.b.y.dim ) || !nzchar( panel.b.y.dim ) )
    stop( "panel.b.y.dim was not supplied and asp$af.channel is not set.", call. = FALSE )
  
  valid.dims <- c( rownames( spectra ), "AF" )
  if ( !panel.c.x.dim %in% valid.dims || !panel.c.y.dim %in% valid.dims )
    stop(
      "panel.c.x.dim/panel.c.y.dim must each be \"AF\" or a row name of ",
      "`spectra`. Got: ", panel.c.x.dim, ", ", panel.c.y.dim, call. = FALSE
    )
  
  if ( is.null( ssc.channel ) ) {
    scatter.params <- read.scatter.parameter( asp )
    ssc.channel <- if ( length( scatter.params ) >= 2 ) scatter.params[ 2 ] else scatter.params[ 1 ]
  }
  
  spectral.channels <- colnames( spectra )
  if ( !all( c( panel.b.x.dim, panel.b.y.dim ) %in% spectral.channels ) )
    stop(
      "panel.b.x.dim/panel.b.y.dim must both be spectral channels present in ",
      "`spectra`. Got: ", panel.b.x.dim, ", ", panel.b.y.dim, call. = FALSE
    )
  
  # -- 1. Load and subsample the unstained sample
  if ( verbose ) message( "\033[34m-- Loading unstained sample --\033[0m" )
  
  if ( is.null( unstained.exprs ) ) {
    read.cols      <- unique( c( spectral.channels, ssc.channel ) )
    unstained.full <- readFCS( unstained.sample, columns = read.cols )
    n.available    <- nrow( unstained.full )
    if ( n.cells > n.available ) {
      warning(
        "n.cells (", n.cells, ") exceeds the number of events available (",
        n.available, "); using all available events.", call. = FALSE
      )
      n.cells <- n.available
    }
    idx <- if ( sample.method == "random" ) {
      set.seed( seed ); sample( seq_len( n.available ), n.cells )
    } else seq_len( n.cells )
    
    unstained.exprs <- unstained.full[ idx, spectral.channels, drop = FALSE ]
    unstained.ssc    <- unstained.full[ idx, ssc.channel ]
    
  } else {
    n.available <- nrow( unstained.exprs )
    if ( n.cells > n.available ) {
      warning(
        "n.cells (", n.cells, ") exceeds nrow(unstained.exprs) (",
        n.available, "); using all available events.", call. = FALSE
      )
      n.cells <- n.available
    }
    idx <- if ( sample.method == "random" ) {
      set.seed( seed ); sample( seq_len( n.available ), n.cells )
    } else seq_len( n.cells )
    
    if ( is.null( unstained.scatter ) ) {
      if ( !ssc.channel %in% colnames( unstained.exprs ) )
        stop(
          "`ssc.channel` ('", ssc.channel, "') not found in `unstained.exprs` ",
          "and `unstained.scatter` was not supplied.", call. = FALSE
        )
      unstained.ssc <- unstained.exprs[ idx, ssc.channel ]
    } else {
      unstained.ssc <- unstained.scatter[ idx ]
    }
    unstained.exprs <- unstained.exprs[ idx, spectral.channels, drop = FALSE ]
  }
  
  # -- 2. AF spectral library (panel A)
  if ( is.null( af.spectra ) ) {
    if ( verbose ) message( "\033[34m-- Computing AF spectra (get.af.spectra()) --\033[0m" )
    af.spectra.default.args <- list(
      unstained.sample = unstained.sample, asp = asp, spectra = spectra,
      unstained.exprs = unstained.exprs, figures = FALSE, save = FALSE,
      verbose = verbose
    )
    af.spectra <- do.call(
      get.af.spectra, utils::modifyList( af.spectra.default.args, get.af.spectra.args )
    )
  }
  
  panel.a <- suppressMessages(
    spectral.trace(
      af.spectra, asp, title = "Lung AF spectra", save = FALSE,
      split.lasers = FALSE, show.legend = FALSE,
      figure.spectra.line.size = asp$figure.spectra.line.size,
      figure.spectra.point.size = asp$figure.spectra.point.size
    )
  )
  
  # -- 3. Per-cell AF Index assignment: residuals / fluorophores / joint
  if ( verbose ) message( "\033[34m-- Assigning AF Index (residuals / fluorophores / joint) --\033[0m" )
  
  unmixed.residuals <- unmix.autospectral(
    unstained.exprs, spectra, af.spectra, asp,
    use.dist0 = FALSE, verbose = FALSE, parallel = parallel, threads = threads
  )
  unmixed.fluorophores <- unmix.autospectral(
    unstained.exprs, spectra, af.spectra, asp,
    use.dist0 = TRUE, verbose = FALSE, parallel = parallel, threads = threads
  )
  unmixed.joint <- unmix.autospectral.joint(
    unstained.exprs, spectra, af.spectra, asp,
    verbose = FALSE, parallel = parallel, threads = threads
  )
  
  # -- 4. Panel B: raw biplot coloured by AF Index (fluorophore-minimization
  # assignment, matching the original figure)
  panel.b <- .af.index.biplot(
    plot.data = unstained.exprs, x.dim = panel.b.x.dim, y.dim = panel.b.y.dim,
    af.index = unmixed.fluorophores[ , "AF Index" ], asp = asp,
    max.points = max.points, color.palette = panel.b.color.palette
  ) + ggplot2::labs( title = "Unstained lung: AF Index assignment" )
  
  # -- 5. Panel C: assignment-method biplots
  panel.c1 <- create.biplot(
    unmixed.residuals, panel.c.x.dim, panel.c.y.dim, asp,
    max.points = max.points, color.palette = "turbo", save = FALSE
  ) + ggplot2::labs( title = "Assignment: residuals" )
  
  panel.c2 <- create.biplot(
    unmixed.fluorophores, panel.c.x.dim, panel.c.y.dim, asp,
    max.points = max.points, color.palette = "turbo", save = FALSE
  ) + ggplot2::labs( title = "Assignment: fluorophores" )
  
  panel.c3 <- create.biplot(
    unmixed.joint, panel.c.x.dim, panel.c.y.dim, asp,
    max.points = max.points, color.palette = "turbo", save = FALSE
  ) + ggplot2::labs( title = "Assignment: joint" )
  
  panel.c <- cowplot::plot_grid( panel.c1, panel.c2, panel.c3, ncol = 3 )
  
  # -- 6. Panel D: cosine-similarity distributions, each KS-tested against
  # the "random assignment" (intra-AF) baseline
  random.distr    <- cosine.similarity( af.spectra )
  random.distr    <- random.distr[ lower.tri( random.distr ) ]
  
  cs.residuals    <- .cosine.sim.paired( unstained.exprs, af.spectra, unmixed.residuals[ , "AF Index" ] )
  cs.fluorophores <- .cosine.sim.paired( unstained.exprs, af.spectra, unmixed.fluorophores[ , "AF Index" ] )
  cs.joint        <- .cosine.sim.paired( unstained.exprs, af.spectra, unmixed.joint[ , "AF Index" ] )
  
  panel.d1 <- .cosine.density.panel( random.distr, "Intra-AF (random)" )
  panel.d2 <- .cosine.density.panel( cs.residuals, "Assignment: residuals", reference = random.distr )
  panel.d3 <- .cosine.density.panel( cs.fluorophores, "Assignment: fluorophores", reference = random.distr )
  panel.d4 <- .cosine.density.panel( cs.joint, "Assignment: joint", reference = random.distr )
  
  panel.d <- cowplot::plot_grid( panel.d1, panel.d2, panel.d3, panel.d4, ncol = 4 )
  
  # -- 7. Synthetic fluorophore overlay (used by panels E, F, G)
  # NOTE: this manual overlay predates sim.flow.data(). For new work
  # generating fully synthetic stained data (with proper shot/spillover/
  # detector noise), prefer sim.flow.data(spectra, asp, af.spectra =
  # af.spectra, ...) instead. Kept as-is here for figure continuity.
  if ( verbose ) message( "\033[34m-- Building synthetic fluorophore overlay --\033[0m" )
  
  if ( is.null( synthetic.fluorophores ) ) {
    set.seed( seed )
    synthetic.fluorophores <- sample( rownames( spectra ), n.synthetic.fluors )
  }
  missing.synth <- setdiff( synthetic.fluorophores, rownames( spectra ) )
  if ( length( missing.synth ) > 0 )
    stop(
      "synthetic.fluorophores not found in `spectra`: ",
      paste( missing.synth, collapse = ", " ), call. = FALSE
    )
  
  n.synth     <- length( synthetic.fluorophores )
  cells.per.f <- floor( n.cells / n.synth )
  synth.groups <- rep( synthetic.fluorophores, each = cells.per.f )
  if ( length( synth.groups ) < n.cells )
    synth.groups <- c(
      synth.groups,
      rep( synthetic.fluorophores[ n.synth ], n.cells - length( synth.groups ) )
    )
  
  set.seed( seed )
  synth.intensities <- stats::rlnorm( n.cells, meanlog = synthetic.meanlog, sdlog = synthetic.sdlog )
  
  synthetic.signal <- matrix(
    0, nrow = n.cells, ncol = ncol( unstained.exprs ), dimnames = dimnames( unstained.exprs )
  )
  for ( f in synthetic.fluorophores ) {
    f.idx <- which( synth.groups == f )
    synthetic.signal[ f.idx, ] <- outer( synth.intensities[ f.idx ], spectra[ f, ] )
  }
  synthetic.data <- unstained.exprs + synthetic.signal
  
  unmixed.percell.synth <- unmix.autospectral(
    synthetic.data, spectra, af.spectra, asp,
    use.dist0 = TRUE, verbose = FALSE, parallel = parallel, threads = threads
  )
  
  # -- 8. Panel E: AF prediction error (single-median AF vs. per-cell AF,
  # with and without the synthetic overlay)
  if ( verbose ) message( "\033[34m-- Computing AF prediction error --\033[0m" )
  
  af.median.row <- apply( unstained.exprs, 2, stats::median )
  af.median.row <- af.median.row / max( af.median.row, 1e-9 )
  combined.spectra.median <- rbind( spectra, AF = af.median.row )
  
  unmixed.median.real  <- unmix.ols( unstained.exprs, combined.spectra.median )
  unmixed.median.synth <- unmix.ols( synthetic.data,  combined.spectra.median )
  
  predicted.median.real   <- unmixed.median.real  %*% combined.spectra.median
  predicted.median.synth  <- unmixed.median.synth %*% combined.spectra.median
  predicted.percell.real  <- .predict.raw.percell( unmixed.fluorophores,  spectra, af.spectra )
  predicted.percell.synth <- .predict.raw.percell( unmixed.percell.synth, spectra, af.spectra )
  
  error.median.real   <- rowSums( abs( unstained.exprs - predicted.median.real ) )
  error.median.synth  <- rowSums( abs( synthetic.data  - predicted.median.synth ) )
  error.percell.real  <- rowSums( abs( unstained.exprs - predicted.percell.real ) )
  error.percell.synth <- rowSums( abs( synthetic.data  - predicted.percell.synth ) )
  
  panel.e.levels <- c(
    "Single median AF\n(no synthetic)", "Per-cell AF\n(no synthetic)",
    "Single median AF\n(+ synthetic)",  "Per-cell AF\n(+ synthetic)"
  )
  panel.e.df <- data.frame(
    Error = c( error.median.real, error.percell.real, error.median.synth, error.percell.synth ),
    Group = factor( rep( panel.e.levels, each = n.cells ), levels = panel.e.levels )
  )
  
  # paired (same cells on both sides of each comparison), matching the
  # violin-of-per-cell-points presentation
  p.e.real  <- stats::wilcox.test( error.median.real,  error.percell.real,  paired = TRUE )$p.value
  p.e.synth <- stats::wilcox.test( error.median.synth, error.percell.synth, paired = TRUE )$p.value
  
  panel.e <- .error.violin.panel(
    panel.e.df, y.lab = "Prediction error (\u03a3|raw - predicted|)",
    title = "AF prediction error",
    p.labels = data.frame(
      x = c( 1.5, 3.5 ), y = max( panel.e.df$Error, na.rm = TRUE ) * 1.05,
      label = c( .fmt.p( p.e.real ), .fmt.p( p.e.synth ) )
    )
  )
  
  # -- 9. Panel F: fluorophore recovery error (single-median AF vs. per-cell
  # AF) for the known synthetic fluorophore signal
  if ( verbose ) message( "\033[34m-- Computing fluorophore recovery error --\033[0m" )
  
  col.idx.median  <- match( synth.groups, colnames( unmixed.median.synth ) )
  col.idx.percell <- match( synth.groups, colnames( unmixed.percell.synth ) )
  
  recovered.median  <- unmixed.median.synth[  cbind( seq_len( n.cells ), col.idx.median ) ]
  recovered.percell <- unmixed.percell.synth[ cbind( seq_len( n.cells ), col.idx.percell ) ]
  
  error.median.recovery  <- abs( recovered.median  - synth.intensities )
  error.percell.recovery <- abs( recovered.percell - synth.intensities )
  
  panel.f.levels <- c( "Single median AF", "Per-cell AF" )
  panel.f.df <- data.frame(
    Error = c( error.median.recovery, error.percell.recovery ),
    Group = factor( rep( panel.f.levels, each = n.cells ), levels = panel.f.levels )
  )
  
  p.f <- stats::ks.test( error.median.recovery, error.percell.recovery )$p.value
  
  panel.f <- .error.violin.panel(
    panel.f.df, y.lab = "Recovery error |true - unmixed|",
    title = "Fluorophore recovery error: single AF vs. per-cell AF",
    p.labels = data.frame(
      x = 1.5, y = max( panel.f.df$Error, na.rm = TRUE ) * 1.05, label = .fmt.p( p.f )
    )
  )
  
  # -- 10. Panel G: AF assignment stability under the synthetic overlay.
  # The "without synthetic" condition compares each cell's real AF
  # assignment to itself and is therefore definitionally 1 for every cell --
  # shown as a reference line rather than a (zero-bandwidth) density.
  cs.with.synthetic <- .cosine.sim.paired(
    af.spectra[ unmixed.fluorophores[ , "AF Index" ], , drop = FALSE ],
    af.spectra, unmixed.percell.synth[ , "AF Index" ]
  )
  
  panel.g <- .cosine.density.panel(
    cs.with.synthetic, title = "AF assignment: with vs. without synthetic fluorophore"
  ) +
    ggplot2::geom_vline( xintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.6 ) +
    ggplot2::annotate(
      "text", x = 0.97, y = Inf, label = "Without synthetic (= 1)",
      angle = 90, vjust = 1.2, hjust = 1.1, size = 3
    )
  
  # -- 11. Panel H: side scatter by AF Index
  panel.h.df <- data.frame(
    SSC = unstained.ssc,
    AFIndex = factor(
      unmixed.fluorophores[ , "AF Index" ],
      levels = sort( unique( unmixed.fluorophores[ , "AF Index" ] ) )
    )
  )
  panel.h <- ggplot2::ggplot( panel.h.df, ggplot2::aes( SSC, color = AFIndex, fill = AFIndex ) ) +
    ggplot2::geom_density( alpha = 0.3, na.rm = TRUE ) +
    ggplot2::labs( title = "AF Index vs. side scatter", x = ssc.channel, y = "Density" ) +
    ggplot2::theme_minimal( base_size = 11 ) +
    ggplot2::theme( legend.position = "none" )
  
  # -- 12. Assemble composite
  composite.width.use  <- if ( !is.null( composite.width ) ) composite.width else panel.width * 4
  row.height.mult <- c( a = 0.7, b = 1, c = 1, d = 0.9, e = 1, f = 1, g = 0.9, h = 0.9 )
  composite.height.use <- if ( !is.null( composite.height ) ) composite.height else
    panel.height * sum( row.height.mult )
  
  composite <- cowplot::plot_grid(
    cowplot::plot_grid( panel.a, labels = "A" ),
    cowplot::plot_grid( panel.b, labels = "B" ),
    cowplot::plot_grid( panel.c, labels = "C" ),
    cowplot::plot_grid( panel.d, labels = "D" ),
    cowplot::plot_grid( panel.e, labels = "E" ),
    cowplot::plot_grid( panel.f, labels = "F" ),
    cowplot::plot_grid( panel.g, labels = "G" ),
    cowplot::plot_grid( panel.h, labels = "H" ),
    ncol = 1, rel_heights = row.height.mult
  )
  
  # -- 13. Save
  if ( save ) {
    out.file <- file.path( output.dir, sprintf( "AF_identification_workflow.%s", file.type ) )
    ggplot2::ggsave(
      out.file, plot = composite, device = plot.device,
      width = composite.width.use, height = composite.height.use, limitsize = FALSE
    )
    if ( verbose ) message( "\033[32m  Saved: ", out.file, "\033[0m" )
  }
  
  invisible( list(
    af.spectra = af.spectra,
    unmixed    = list( residuals = unmixed.residuals, fluorophores = unmixed.fluorophores, joint = unmixed.joint ),
    synthetic  = list(
      data = synthetic.data, groups = synth.groups, intensities = synth.intensities,
      unmixed.median = unmixed.median.synth, unmixed.percell = unmixed.percell.synth
    ),
    panels = list(
      a = panel.a, b = panel.b, c = panel.c, d = panel.d,
      e = panel.e, f = panel.f, g = panel.g, h = panel.h
    ),
    composite    = composite,
    panel.e.data = panel.e.df,
    panel.f.data = panel.f.df
  ) )
}