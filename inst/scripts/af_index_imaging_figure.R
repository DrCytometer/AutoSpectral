# af_index_imaging_figure.R

#' @title Get Imaging Parameter Columns
#'
#' @description
#' Identifies imaging-specific parameter columns (e.g. from a BD FACSDiscover
#' A8/S8 file) by name pattern, distinct from spectral fluorescence channels
#' and the two conventional scatter parameters used for gating.
#'
#' @param column.names Character vector of column names to search, e.g.
#'   `colnames(raw.data)`.
#' @param imaging.pattern Character vector of substrings identifying imaging
#'   parameter columns. Default covers the BD FACSDiscover imaging channel
#'   families: `"Radial"`, `"Correlation"`, `"Intensity"`, `"Eccentricity"`,
#'   `"Diffusivity"`, `"Center"`, `"Moment"`, `"Size"`, `"LightLoss"`,
#'   `"Imaging"` and `"Img"`. Matching is by substring, so e.g. `"LightLoss"`
#'   matches both `LightLoss (Imaging)-A` and `LightLoss (Violet)-A` --
#'   narrow this vector if only the imaging-specific variant is wanted.
#' @param exclude Optional character vector of column names to remove from
#'   the match regardless of pattern (e.g. unmixed fluorophore names, `"AF"`,
#'   `"AF Index"`). Default `NULL`.
#'
#' @return Character vector of matching column names, in their original order.
#'
#' @export

get.imaging.parameter.columns <- function(
    column.names,
    imaging.pattern = c(
      "Radial", "Correlation", "Intensity", "Eccentricity", "Diffusivity",
      "Center", "Moment", "Size", "LightLoss", "Imaging", "Img", "Delta CoM"
    ),
    exclude = NULL
) {

  pattern <- paste( imaging.pattern, collapse = "|" )
  matched <- grep( pattern, column.names, value = TRUE )

  if ( !is.null( exclude ) ) matched <- setdiff( matched, exclude )

  matched
}

#' @title Reorder AF Index by Spectral Similarity
#'
#' @description
#' Builds a mapping from the arbitrary, categorical `"AF Index"` values
#' produced by `unmix.autospectral.rcpp()` (1-based row indices into
#' `af.spectra`) to a similarity-ranked position, so that AF library entries
#' with similar spectral shape receive numerically close values. This turns
#' AF Index from an arbitrary label into a quantity that can be treated as
#' roughly continuous for plotting and summary statistics -- a smooth
#' gradient in a UMAP overlay then reflects real autofluorescence structure
#' rather than an artifact of library row order.
#'
#' The reference point is the mean AF spectrum across all rows of
#' `af.spectra`; every row is ranked by its cosine similarity to that mean,
#' most similar first. This is a 1D projection, so it's an approximation:
#' two AF profiles that are both very unlike the mean, but in different
#' directions (e.g. one has a UV spike, one has a blue spike), can still end
#' up with adjacent ranks despite not resembling each other. It's a cheap,
#' honest "how typical is this AF profile" axis, not a full embedding of AF
#' shape.
#'
#' @param af.spectra Matrix, AF spectra in rows, detectors in columns
#'   (normalized 0-1, as produced by `get.af.spectra()`).
#'
#' @return A named list:
#' \describe{
#'   \item{`rank`}{Integer vector, length `nrow(af.spectra)`, giving each
#'     original row's rank (1 = most similar to the mean spectrum), indexed
#'     by original row number -- `rank[i]` is the new value for old AF Index
#'     `i`. Use `af.index.to.continuous()` to apply this to an event vector.}
#'   \item{`order`}{Integer vector giving the original row order sorted by
#'     similarity (most to least).}
#'   \item{`cosine.similarity`}{Numeric vector, length `nrow(af.spectra)`,
#'     cosine similarity of each original row to the mean spectrum.}
#'   \item{`mean.spectrum`}{The reference spectrum used.}
#' }
#'
#' @export

reorder.af.index.by.similarity <- function( af.spectra ) {

  af.spectra    <- as.matrix( af.spectra )
  mean.spectrum <- colMeans( af.spectra )

  cosine.similarity <- apply( af.spectra, 1, function( r ) {
    denom <- sqrt( sum( r^2 ) ) * sqrt( sum( mean.spectrum^2 ) )
    if ( denom == 0 ) return( 0 )
    sum( r * mean.spectrum ) / denom
  } )

  order.idx <- order( cosine.similarity, decreasing = TRUE )
  rank      <- match( seq_along( cosine.similarity ), order.idx )

  list(
    rank              = rank,
    order              = order.idx,
    cosine.similarity  = cosine.similarity,
    mean.spectrum      = mean.spectrum
  )
}

#' @title Convert AF Index to a Similarity-Ranked Continuous Value
#'
#' @param af.index Integer vector, raw `"AF Index"` values (1-based row into
#'   the `af.spectra` used for `reorder.af.index.by.similarity()`). `0` or
#'   any value outside `1:length(af.reorder$rank)` -- unassigned events --
#'   is returned as `NA`.
#' @param af.reorder The list returned by `reorder.af.index.by.similarity()`.
#'
#' @return Numeric vector, same length as `af.index`, the similarity rank
#'   for each event, `NA` for unassigned events.
#'
#' @export

af.index.to.continuous <- function( af.index, af.reorder ) {

  n.af       <- length( af.reorder$rank )
  safe.index <- af.index
  safe.index[ !is.na( safe.index ) & ( safe.index < 1 | safe.index > n.af ) ] <- NA

  af.reorder$rank[ safe.index ]
}

# Internal helper: writes a ggplot or grid-graphics object (e.g. FlowSOM's
# PlotStars() output) to disk by opening the target device directly and
# printing, rather than relying on ggplot2::ggsave()'s object-type detection.

.save.figure <- function( plot, filename, width, height, file.type = "pdf", dpi = 300 ) {

  file.type <- tolower( file.type )
  if ( file.type == "jpeg" ) file.type <- "jpg"

  plot.device <- switch( file.type,
                         jpg  = ragg::agg_jpeg,
                         tiff = ragg::agg_tiff,
                         png  = ragg::agg_png,
                         pdf  = grDevices::pdf,
                         stop( "Unsupported file.type: ", file.type, call. = FALSE )
  )

  if ( file.type == "pdf" ) {
    plot.device( file = filename, width = width, height = height )
  } else {
    plot.device( filename = filename, width = width, height = height, units = "in", res = dpi )
  }

  on.exit( grDevices::dev.off(), add = TRUE )

  print( plot )

  invisible( filename )
}

# Internal helper: builds a qualitative color palette for a fixed number of
# discrete clusters. RColorBrewer's qualitative palettes give visually
# distinct hues, which viridis/turbo does not once many categories are
# packed close together in embedding space -- falls back to interpolating
# the palette if more colors are needed than the chosen palette provides.

.get.cluster.palette <- function( n, brewer.palette = "Paired" ) {

  if ( !requireNamespace( "RColorBrewer", quietly = TRUE ) )
    stop(
      "Cluster plotting with palette.type = \"brewer\" requires the RColorBrewer ",
      "package (not currently a package dependency -- install it separately: ",
      "install.packages(\"RColorBrewer\")).",
      call. = FALSE
    )

  max.colors <- RColorBrewer::brewer.pal.info[ brewer.palette, "maxcolors" ]

  if ( n <= 3 ) {
    base.colors <- RColorBrewer::brewer.pal( 3, brewer.palette )
    return( base.colors[ seq_len( n ) ] )
  }

  if ( n <= max.colors ) return( RColorBrewer::brewer.pal( n, brewer.palette ) )

  grDevices::colorRampPalette( RColorBrewer::brewer.pal( max.colors, brewer.palette ) )( n )
}

#' @title Cluster Imaging Parameters with FlowSOM and UMAP
#'
#' @description
#' Runs UMAP (`uwot::umap()`) on a downsampled subset of (post-cleaning)
#' events, then runs FlowSOM (self-organizing map plus consensus
#' metaclustering) on either that same 2D embedding or the higher-dimensional
#' imaging-parameter space, depending on `cluster.on`.
#'
#' Both the UMAP input and, when `cluster.on = "imaging"`, the FlowSOM
#' clustering input are by default z-scored (`z.score`) and PCA-reduced
#' (`use.pca`) first, so distance calculations aren't dominated by a handful
#' of correlated imaging features.
#'
#' Regardless of `cluster.on`, the FlowSOM object's underlying data always
#' carries the original z-scored imaging parameters alongside whatever it
#' actually clustered on, so `plot.flowsom.stars()` can look up an original
#' parameter's per-node MFI by name afterward.
#'
#' @importFrom stats sd prcomp
#'
#' @param imaging.data Numeric matrix, events x imaging parameters, with
#'   informative column names (e.g. from `get.imaging.parameter.columns()`).
#' @param asp The AutoSpectral parameter list. Used for `asp$bird.seed`.
#' @param cluster.on One of `"embedding"` (default) or `"imaging"`. Controls
#'   what FlowSOM actually computes clustering distance on. `"embedding"`
#'   clusters directly on the 2D UMAP coordinates -- useful when the point
#'   of clustering is purely to color/summarize visible groupings in the
#'   embedding, and sidesteps FlowSOM struggling in a higher-dimensional
#'   space that UMAP is already resolving better. `cluster.id`/
#'   `metacluster.id` then cover only the (possibly downsampled) `embed.idx`
#'   events, in the same row order as `embed.coords` -- no separate
#'   alignment needed for embedding plots. `"imaging"` clusters on the
#'   (PCA-reduced, if `use.pca`) z-scored imaging parameters directly,
#'   covering every event in `keep.idx`; use `align.to.embedding()` to line
#'   results up with `embed.coords` in this mode.
#' @param z.score Logical, default `TRUE`. Whether to z-score each imaging
#'   parameter (mean 0, sd 1) before PCA/clustering/embedding. Zero-variance
#'   columns are dropped automatically with a message if `TRUE`.
#' @param use.pca Logical, default `TRUE`. Whether to PCA-reduce the
#'   (z-scored) imaging parameters before UMAP, and before FlowSOM when
#'   `cluster.on = "imaging"`.
#' @param pca.components Integer, optional. Fixed number of principal
#'   components to use. Default `NULL` selects the smallest number of
#'   components reaching `pca.var.explained` cumulative variance.
#' @param pca.var.explained Numeric in (0, 1], default `0.9`. Cumulative
#'   variance threshold used to pick the component count when
#'   `pca.components` is `NULL`. Ignored if `use.pca = FALSE`.
#' @param som.xdim,som.ydim Integer, FlowSOM grid dimensions. Default `8`
#'   each (64 nodes).
#' @param n.metaclusters Integer, target number of FlowSOM metaclusters
#'   (`nClus`). Default `10`.
#' @param max.embed.events Integer, maximum events to embed with UMAP.
#'   Default `4e4`. When `cluster.on = "embedding"`, this is also the
#'   maximum number of events that receive a cluster assignment.
#' @param umap.n.neighbors,umap.min.dist,umap.metric Passed to
#'   `uwot::umap()` as `n_neighbors`, `min_dist`, `metric`.
#' @param threads Integer, threads for `uwot::umap()`. Default `1L`. `uwot`
#'   is only bit-for-bit reproducible at `threads = 1` even with a fixed
#'   seed.
#' @param verbose Logical, default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`fsom`}{The FlowSOM object.}
#'   \item{`imaging.params`}{Character vector, the *original* imaging
#'     parameter columns actually used (after dropping any zero-variance
#'     columns).}
#'   \item{`pca`}{The `prcomp()` object, or `NULL` if `use.pca = FALSE`.}
#'   \item{`n.pcs.used`}{Integer, components actually used, or `NA` if
#'     `use.pca = FALSE`.}
#'   \item{`cluster.on`}{Echoes the `cluster.on` argument.}
#'   \item{`cluster.id`}{SOM node per event covered -- see `cluster.idx`.}
#'   \item{`metacluster.id`}{Metacluster per event covered -- see
#'     `cluster.idx`.}
#'   \item{`embed.coords`}{Matrix, `length(embed.idx)` x 2, columns `UMAP1`,
#'     `UMAP2`.}
#'   \item{`keep.idx`}{Integer vector, indices into the original
#'     `imaging.data` rows that were finite and used at all.}
#'   \item{`embed.idx`}{Integer vector, indices into `keep.idx` that were
#'     embedded with UMAP.}
#'   \item{`cluster.idx`}{Integer vector, indices into `keep.idx` that
#'     `cluster.id`/`metacluster.id` cover -- identical to `embed.idx` when
#'     `cluster.on = "embedding"`, `seq_along(keep.idx)` (every event) when
#'     `cluster.on = "imaging"`. Use `cluster.event.index()` to convert this
#'     to absolute row indices into the original `imaging.data`, or
#'     `align.to.embedding()` to map a `cluster.idx`-aligned vector onto
#'     `embed.coords`.}
#' }
#'
#' @export

cluster.imaging.parameters <- function(
    imaging.data,
    asp,
    cluster.on         = c( "embedding", "imaging" ),
    z.score            = TRUE,
    use.pca            = TRUE,
    pca.components      = NULL,
    pca.var.explained   = 0.9,
    som.xdim           = 8,
    som.ydim           = 8,
    n.metaclusters     = 10,
    max.embed.events   = 4e4,
    umap.n.neighbors   = 15,
    umap.min.dist      = 0.1,
    umap.metric        = "euclidean",
    threads            = 1L,
    verbose            = TRUE
) {

  if ( !requireNamespace( "FlowSOM", quietly = TRUE ) )
    stop( "cluster.imaging.parameters() requires the FlowSOM package.", call. = FALSE )
  if ( !requireNamespace( "uwot", quietly = TRUE ) )
    stop(
      "cluster.imaging.parameters() requires the uwot package (not currently ",
      "a package dependency -- install it separately: install.packages(\"uwot\")).",
      call. = FALSE
    )

  cluster.on <- match.arg( cluster.on )

  training.data <- as.matrix( imaging.data )

  finite.mask <- apply( training.data, 1, function( r ) all( is.finite( r ) ) )
  keep.idx    <- which( finite.mask )

  if ( length( keep.idx ) < nrow( training.data ) && verbose )
    message( sprintf(
      "Dropping %d of %d events with non-finite imaging parameter values",
      nrow( training.data ) - length( keep.idx ), nrow( training.data )
    ) )

  training.data <- training.data[ keep.idx, , drop = FALSE ]

  if ( z.score ) {
    col.sd   <- apply( training.data, 2, stats::sd, na.rm = TRUE )
    zero.var <- names( col.sd )[ is.na( col.sd ) | col.sd == 0 ]

    if ( length( zero.var ) > 0 ) {
      if ( verbose )
        message( "Dropping zero-variance imaging parameter(s): ", paste( zero.var, collapse = ", " ) )
      training.data <- training.data[ , setdiff( colnames( training.data ), zero.var ), drop = FALSE ]
    }

    training.data <- scale( training.data )
  }

  imaging.params <- colnames( training.data )

  pca.result   <- NULL
  n.pcs.used   <- NA_integer_
  cluster.input <- training.data     # UMAP input: PCA scores if use.pca, else z-scored data

  if ( use.pca ) {
    pca.result    <- stats::prcomp( training.data, center = FALSE, scale. = FALSE )
    var.explained <- cumsum( pca.result$sdev^2 ) / sum( pca.result$sdev^2 )

    n.pcs.used <- if ( !is.null( pca.components ) ) {
      min( pca.components, ncol( pca.result$x ) )
    } else {
      which( var.explained >= pca.var.explained )[ 1 ]
    }
    n.pcs.used <- max( n.pcs.used, 2L )

    if ( verbose )
      message( sprintf(
        "PCA: using %d of %d components (%.1f%% variance explained)",
        n.pcs.used, ncol( pca.result$x ), 100 * var.explained[ n.pcs.used ]
      ) )

    cluster.input <- pca.result$x[ , seq_len( n.pcs.used ), drop = FALSE ]
  }

  event.n <- nrow( training.data )
  set.seed( asp$bird.seed )
  embed.idx <- if ( event.n > max.embed.events )
    sort( sample( event.n, max.embed.events ) ) else seq_len( event.n )

  if ( verbose )
    message( sprintf( "Running UMAP on %d of %d events", length( embed.idx ), event.n ) )

  embed.coords <- uwot::umap(
    cluster.input[ embed.idx, , drop = FALSE ],
    n_neighbors = umap.n.neighbors,
    min_dist    = umap.min.dist,
    metric      = umap.metric,
    n_threads   = threads,
    verbose     = verbose
  )
  colnames( embed.coords ) <- c( "UMAP1", "UMAP2" )

  # ---- what FlowSOM actually clusters on ----
  if ( cluster.on == "embedding" ) {
    cluster.basis <- embed.coords
    cluster.idx   <- embed.idx
  } else {
    cluster.basis <- cluster.input
    cluster.idx   <- seq_len( event.n )
  }

  flowsom.extra <- training.data[ cluster.idx, , drop = FALSE ]

  flowsom.input <- if ( identical( cluster.basis, flowsom.extra ) ) {
    # cluster.on = "imaging", use.pca = FALSE -- clustering basis IS the
    # original imaging parameters, nothing extra to attach.
    flowsom.extra
  } else {
    dup.names <- intersect( colnames( flowsom.extra ), colnames( cluster.basis ) )
    if ( length( dup.names ) > 0 )
      stop( "Column name collision between imaging parameters and the clustering basis: ",
            paste( dup.names, collapse = ", " ), call. = FALSE )
    cbind( flowsom.extra, cluster.basis )
  }

  if ( verbose )
    message( sprintf(
      "Clustering %d events x %d dimension(s) (%s space) with a %dx%d FlowSOM grid (%d metaclusters)",
      nrow( cluster.basis ), ncol( cluster.basis ), cluster.on, som.xdim, som.ydim, n.metaclusters
    ) )

  fsom <- FlowSOM::FlowSOM(
    input      = flowsom.input,
    compensate = FALSE,
    transform  = FALSE,
    scale      = FALSE,
    colsToUse  = colnames( cluster.basis ),
    xdim       = som.xdim,
    ydim       = som.ydim,
    nClus      = n.metaclusters,
    seed       = asp$bird.seed
  )

  cluster.id     <- FlowSOM::GetClusters( fsom )
  metacluster.id <- FlowSOM::GetMetaclusters( fsom )

  list(
    fsom           = fsom,
    imaging.params = imaging.params,
    pca             = pca.result,
    n.pcs.used      = n.pcs.used,
    cluster.on      = cluster.on,
    cluster.id     = cluster.id,
    metacluster.id = metacluster.id,
    embed.coords   = embed.coords,
    keep.idx       = keep.idx,
    embed.idx      = embed.idx,
    cluster.idx    = cluster.idx
  )
}

#' @title Absolute Event Indices for a Clustering Result
#'
#' @description
#' Returns row indices into the original data (e.g. `gated.data`) passed as
#' `imaging.data` to `cluster.imaging.parameters()`, for whichever events
#' `cluster.id`/`metacluster.id` cover (all gated events when `cluster.on =
#' "imaging"`, only the embedded subset when `cluster.on = "embedding"`).
#' Use this to pull other per-event columns (AF Index, other imaging
#' parameters, etc.) into alignment with `cluster.id`.
#'
#' @param cluster.result The list returned by `cluster.imaging.parameters()`.
#'
#' @return Integer vector of row indices.
#'
#' @export

cluster.event.index <- function( cluster.result ) {
  cluster.result$keep.idx[ cluster.result$cluster.idx ]
}

#' @title Align a Cluster-Indexed Vector to Embedding Coordinates
#'
#' @description
#' `cluster.id`/`metacluster.id` and `embed.coords` may cover different
#' (overlapping, differently-ordered) subsets of events depending on
#' `cluster.on` -- when `cluster.on = "embedding"` they cover exactly the
#' same events in the same order and this is a no-op; when `cluster.on =
#' "imaging"`, `cluster.id` covers every event and `embed.coords` covers
#' only the (downsampled) `embed.idx` subset of them. This maps a vector
#' aligned to `cluster.result$cluster.idx` onto the `embed.coords` event
#' set, so the same plotting code works regardless of which mode was used.
#'
#' @param cluster.result The list returned by `cluster.imaging.parameters()`.
#' @param values A vector aligned to `cluster.result$cluster.idx` (e.g.
#'   `cluster.result$cluster.id`, `cluster.result$metacluster.id`, or any
#'   other per-event vector pulled via `cluster.event.index()`).
#'
#' @return `values`, reordered/subset to align row-for-row with
#'   `cluster.result$embed.coords`.
#'
#' @export

align.to.embedding <- function( cluster.result, values ) {

  if ( length( values ) != length( cluster.result$cluster.idx ) )
    stop( "`values` must be the same length as `cluster.result$cluster.idx`.", call. = FALSE )

  values[ match( cluster.result$embed.idx, cluster.result$cluster.idx ) ]
}

#' @title Plot Embedding with FlowSOM Cluster Overlay
#'
#' @importFrom ggplot2 ggplot aes labs theme_bw theme margin element_line
#' @importFrom ggplot2 element_text element_rect element_blank guides guide_legend
#' @importFrom ggplot2 scale_color_manual scale_color_viridis_d
#' @importFrom scattermore geom_scattermore
#'
#' @param embed.coords Matrix, events x 2 (`UMAP1`, `UMAP2`), from
#'   `cluster.imaging.parameters()$embed.coords`.
#' @param cluster.id Vector, length `nrow(embed.coords)`, cluster/metacluster
#'   assignment for the same events (already subset to `embed.idx`).
#' @param cluster.labels Optional named character vector mapping raw cluster
#'   id (as character) to a display name, e.g.
#'   `c("1" = "Granular", "2" = "Round, low-texture")`. Ids not present in
#'   `cluster.labels` are left as their raw id. Default `NULL` (no renaming).
#' @param cluster.factor Optional vector of raw cluster ids giving the full,
#'   canonical set and order of clusters. Pass the same value here and to
#'   `plot.af.index.violin()` so colors match exactly. Default NULL derives
#'   levels from `cluster.id` alone.
#' @param cluster.levels Optional vector of raw cluster ids giving the full,
#'   canonical set and order of clusters (e.g.
#'   `sort(unique(cluster.result$metacluster.id))`). Pass the same value
#'   here and to `plot.embedding.clusters()` (and `plot.flowsom.stars()`) so
#'   violin fill colors match the embedding plot's cluster colors exactly,
#'   even if this call's `cluster.id` happens not to include every cluster.
#'   Default `NULL` derives levels from `cluster.id` alone, which is only
#'   guaranteed consistent with another plot if that plot saw the same set
#'   of clusters.
#' @param asp The AutoSpectral parameter list. Used for figure theme sizing.
#' @param palette.type One of `"brewer"` (default) or `"viridis"`. Brewer's
#'   qualitative palettes are far more distinguishable than a viridis scale
#'   once clusters sit close together in embedding space.
#' @param brewer.palette RColorBrewer qualitative palette name, used when
#'   `palette.type = "brewer"`. Default `"Paired"` (handles up to 12 clusters
#'   distinctly; interpolated beyond that).
#' @param viridis.palette Viridis option, used when `palette.type =
#'   "viridis"`. Default `"turbo"`.
#' @param point.size Point size for `scattermore::geom_scattermore()`.
#'   Default `0.5`.
#' @param legend.point.size Point size shown in the legend key (scattermore
#'   points render tiny by default). Default `3`.
#' @param title Plot title.
#' @param output.dir Output directory. Created if missing. Default
#'   `"figure_af_imaging"`.
#' @param save Logical, default `TRUE`.
#' @param file.type One of `"pdf"` (default), `"jpg"`, `"tiff"`, `"png"`.
#' @param width,height Figure dimensions in inches. Default `5`, `5`.
#'
#' @return The ggplot object (invisibly saved as a side effect if `save = TRUE`).
#'
#' @export

plot.embedding.clusters <- function(
    embed.coords,
    cluster.id,
    cluster.labels    = NULL,
    cluster.factor    = NULL,
    cluster.levels    = NULL,
    asp,
    palette.type       = c( "brewer", "viridis" ),
    brewer.palette     = "Paired",
    viridis.palette    = "turbo",
    point.size         = 0.5,
    legend.point.size  = 3,
    title              = "UMAP -- FlowSOM metaclusters",
    output.dir         = "figure_af_imaging",
    save               = TRUE,
    file.type          = "pdf",
    width              = 5,
    height             = 5
) {

  if ( nrow( embed.coords ) != length( cluster.id ) )
    stop( "`embed.coords` and `cluster.id` must have the same number of events.",
          call. = FALSE )

  palette.type <- match.arg( palette.type )

  cluster.factor <- factor( cluster.id, levels = cluster.levels )

  if ( !is.null( cluster.labels ) ) {
    lvl <- levels( cluster.factor )
    levels( cluster.factor ) <- ifelse(
      lvl %in% names( cluster.labels ), cluster.labels[ lvl ], lvl
    )
  }

  plot.df <- data.frame(
    x       = embed.coords[ , "UMAP1" ],
    y       = embed.coords[ , "UMAP2" ],
    cluster = cluster.factor
  )

  embed.plot <- ggplot( plot.df, aes( x, y, color = cluster ) ) +
    scattermore::geom_scattermore( pointsize = point.size, alpha = 0.9, na.rm = TRUE )

  embed.plot <- if ( palette.type == "brewer" )
    embed.plot + scale_color_manual(
      values = .get.cluster.palette( nlevels( cluster.factor ), brewer.palette ), name = "Cluster"
    )
  else
    embed.plot + scale_color_viridis_d( option = viridis.palette, name = "Cluster" )

  embed.plot <- embed.plot +
    guides( color = guide_legend( override.aes = list( size = legend.point.size, alpha = 1 ) ) ) +
    labs( x = "UMAP 1", y = "UMAP 2", title = title ) +
    theme_bw() +
    theme(
      plot.margin = margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      axis.ticks    = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text     = element_text( size = asp$figure.axis.text.size ),
      axis.title    = element_text( size = asp$figure.axis.title.size ),
      panel.border  = element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if ( save ) {
    if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
    .save.figure(
      embed.plot, file.path( output.dir, paste0( "embedding_clusters.", file.type ) ),
      width = width, height = height, file.type = file.type
    )
  }

  embed.plot
}

#' @title Plot FlowSOM Stars / Minimal Spanning Tree
#'
#' @param fsom A FlowSOM object, from `cluster.imaging.parameters()$fsom`.
#' @param markers Optional character vector, a subset of the imaging
#'   parameters to draw as star wedges. Default `NULL` uses every parameter
#'   the SOM was trained on -- with more than ~8-10 imaging parameters this
#'   tends to be illegible; pass a shorter, curated subset (e.g. the
#'   highest-variance ones) if stars are still too cluttered.
#' @param background One of `"metacluster"` (default) or `"none"`.
#' @param view One of `"grid"` (default) or `"MST"`. `"grid"` lays nodes out
#'   on a regular lattice, which stays legible at a fixed size; `"MST"` shows
#'   the tree topology but can crowd nodes unevenly.
#' @param equal.node.size Logical, default `TRUE`. Draws every node at the
#'   same size regardless of population, rather than scaling tiny nodes down
#'   to near-invisibility.
#' @param brewer.palette RColorBrewer qualitative palette name for the
#'   metacluster background color, matched to `plot.embedding.clusters()`'s
#'   default so the same metacluster is the same color in both panels.
#'   Default `"Paired"`.
#' @param title Plot title. Default `"FlowSOM map"`.
#' @param output.dir Output directory. Created if missing. Default
#'   `"figure_af_imaging"`.
#' @param save Logical, default `TRUE`.
#' @param file.type One of `"pdf"` (default), `"jpg"`, `"tiff"`, `"png"`.
#' @param width,height Figure dimensions in inches. Default `12`, `10` --
#'   this plot needs real canvas space; shrink after checking legibility,
#'   don't start small.
#' @param ... Additional arguments passed to `FlowSOM::PlotStars()`, e.g. a
#'   label-size argument if your installed FlowSOM version exposes one --
#'   check `?FlowSOM::PlotStars` for the exact name.
#' @param cluster.levels Optional vector of raw cluster ids giving the full,
#'   canonical set and order of clusters (e.g.
#'   `sort(unique(cluster.result$metacluster.id))`). Pass the same value
#'   here and to `plot.embedding.clusters()` (and `plot.flowsom.stars()`) so
#'   violin fill colors match the embedding plot's cluster colors exactly,
#'   even if this call's `cluster.id` happens not to include every cluster.
#'   Default `NULL` derives levels from `cluster.id` alone, which is only
#'   guaranteed consistent with another plot if that plot saw the same set
#'   of clusters.
#'
#' @return The plot object returned by `FlowSOM::PlotStars()`.
#'
#' @export

plot.flowsom.stars <- function(
    fsom,
    markers          = NULL,
    background       = c( "metacluster", "none" ),
    view             = c( "grid", "MST" ),
    equal.node.size  = TRUE,
    brewer.palette   = "Paired",
    title            = "FlowSOM map",
    output.dir       = "figure_af_imaging",
    save             = TRUE,
    file.type        = "pdf",
    width            = 12,
    height           = 10,
    cluster.levels   = NULL,
    ...
) {

  if ( !requireNamespace( "FlowSOM", quietly = TRUE ) )
    stop( "plot.flowsom.stars() requires the FlowSOM package.", call. = FALSE )

  background <- match.arg( background )
  view       <- match.arg( view )

  background.values <- NULL
  background.colors <- NULL

  if ( background == "metacluster" ) {
    background.values <- fsom$metaclustering
    background.colors <- .get.cluster.palette(
      nlevels( factor( background.values, levels = cluster.levels ) ),
      brewer.palette
    )
  }

  stars.plot <- FlowSOM::PlotStars(
    fsom,
    markers          = markers,
    view             = view,
    backgroundValues = background.values,
    backgroundColors  = background.colors,
    equalNodeSize    = equal.node.size,
    title            = title,
    ...
  )

  if ( save ) {
    if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
    .save.figure(
      stars.plot, file.path( output.dir, paste0( "flowsom_stars.", file.type ) ),
      width = width, height = height, file.type = file.type
    )
  }

  stars.plot
}

#' @title Heatmap of Median Imaging Parameters (and AF Index) by Cluster
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c labs
#' @importFrom ggplot2 theme_classic theme element_text
#'
#' @param imaging.data Numeric matrix, events x imaging parameters, subset to
#'   the same events as `cluster.id` (i.e. already indexed by `keep.idx`).
#' @param cluster.id Vector, length `nrow(imaging.data)`, cluster/metacluster
#'   assignment for the same events.
#' @param imaging.params Character vector of column names in `imaging.data`
#'   to summarize, in the desired row order (top to bottom, below the AF
#'   Index row if `af.index` is supplied).
#' @param af.index Optional numeric vector, length `nrow(imaging.data)`, a
#'   (preferably similarity-ranked, see `af.index.to.continuous()`) AF Index
#'   value per event. If supplied, its per-cluster median is added as an
#'   extra row at the top of the heatmap, on the same z-scored color scale
#'   as everything else. Default `NULL` (omit).
#' @param af.index.label Row label for the AF Index row. Default
#'   `"AF Index (median)"`.
#' @param cluster.labels Optional named character vector renaming cluster ids
#'   for display, same convention as `plot.embedding.clusters()`.
#' @param z.score Logical, default `TRUE`. Whether to z-score each row's
#'   (parameter's, or AF Index's) column of cluster medians before plotting.
#' @param asp The AutoSpectral parameter list. Used for the default
#'   `axis.text.size`.
#' @param color.palette Viridis option for the fill scale. Default `"viridis"`.
#' @param title Plot title.
#' @param output.dir Output directory. Created if missing. Default
#'   `"figure_af_imaging"`.
#' @param save Logical, default `TRUE`.
#' @param file.type One of `"pdf"` (default), `"jpg"`, `"tiff"`, `"png"`.
#' @param width,height Figure dimensions in inches. Default `NULL` auto-scales.
#' @param axis.text.size Numeric, default `NULL` uses `asp$figure.axis.text.size`.
#'
#' @return The ggplot object.
#'
#' @export

plot.cluster.imaging.heatmap <- function(
    imaging.data,
    cluster.id,
    imaging.params,
    af.index        = NULL,
    af.index.label  = "AF Index (median)",
    cluster.labels  = NULL,
    z.score         = TRUE,
    asp,
    color.palette   = "viridis",
    title           = "Median imaging parameters by cluster",
    output.dir      = "figure_af_imaging",
    save            = TRUE,
    file.type       = "pdf",
    width           = NULL,
    height          = NULL,
    axis.text.size  = NULL
) {

  if ( nrow( imaging.data ) != length( cluster.id ) )
    stop( "`imaging.data` and `cluster.id` must have the same number of events.",
          call. = FALSE )
  if ( !is.null( af.index ) && length( af.index ) != length( cluster.id ) )
    stop( "`af.index` must be the same length as `cluster.id`.", call. = FALSE )

  cluster.factor <- factor( cluster.id )

  med.mat <- t( sapply( levels( cluster.factor ), function( cl ) {
    idx       <- which( cluster.factor == cl )
    param.med <- apply( imaging.data[ idx, imaging.params, drop = FALSE ], 2, stats::median, na.rm = TRUE )
    if ( !is.null( af.index ) )
      param.med <- c( param.med, stats::setNames( stats::median( af.index[ idx ], na.rm = TRUE ), af.index.label ) )
    param.med
  } ) )
  rownames( med.mat ) <- levels( cluster.factor )

  plot.order <- if ( !is.null( af.index ) ) c( af.index.label, imaging.params ) else imaging.params

  if ( z.score ) med.mat <- scale( med.mat )

  if ( !is.null( cluster.labels ) ) {
    rn <- rownames( med.mat )
    rownames( med.mat ) <- ifelse( rn %in% names( cluster.labels ), cluster.labels[ rn ], rn )
  }

  heat.df <- as.data.frame( med.mat, check.names = FALSE )
  heat.df$cluster <- rownames( heat.df )

  heat.long <- data.frame(
    cluster   = rep( heat.df$cluster, times = ncol( heat.df ) - 1 ),
    parameter = rep( colnames( heat.df )[ -ncol( heat.df ) ], each = nrow( heat.df ) ),
    value     = as.vector( as.matrix( heat.df[ , -ncol( heat.df ), drop = FALSE ] ) ),
    stringsAsFactors = FALSE
  )

  heat.long$cluster   <- factor( heat.long$cluster, levels = rownames( med.mat ) )
  heat.long$parameter <- factor( heat.long$parameter, levels = rev( plot.order ) )

  if ( is.null( axis.text.size ) ) axis.text.size <- asp$figure.axis.text.size
  if ( is.null( width ) )  width  <- max( 6, 0.6 * nlevels( heat.long$cluster ) + 3 )
  if ( is.null( height ) ) height <- max( 5, 0.28 * length( plot.order ) )

  heatmap.plot <- ggplot( heat.long, aes( cluster, parameter, fill = value ) ) +
    geom_tile() +
    scale_fill_viridis_c( option = color.palette, name = if ( z.score ) "z-score" else "median" ) +
    labs( x = NULL, y = NULL, title = title ) +
    theme_classic() +
    theme(
      axis.text.x = element_text( angle = 45, hjust = 1, size = axis.text.size ),
      axis.text.y = element_text( size = axis.text.size )
    )

  if ( save ) {
    if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
    .save.figure(
      heatmap.plot, file.path( output.dir, paste0( "cluster_imaging_heatmap.", file.type ) ),
      width = width, height = height, file.type = file.type
    )
  }

  heatmap.plot
}

#' @title Plot Embedding with AF Index Overlay
#'
#' @importFrom ggplot2 ggplot aes labs theme_bw theme margin element_line
#' @importFrom ggplot2 element_text element_rect element_blank scale_color_viridis_c
#' @importFrom scattermore geom_scattermore
#'
#' @param embed.coords Matrix, events x 2 (`UMAP1`, `UMAP2`).
#' @param af.index Numeric/integer vector, length `nrow(embed.coords)`, the
#'   `"AF Index"` column for the same events (already subset to
#'   `embed.idx`). Plotted on a continuous viridis scale for visual purposes
#'   only -- see `test.af.index.cluster.chisq()` /
#'   `test.af.index.embedding.location()` for the categorical statistical
#'   treatment.
#' @param asp The AutoSpectral parameter list. Used for figure theme sizing.
#' @param color.palette Viridis option for the continuous scale. Default
#'   `"viridis"`.
#' @param point.size Point size for `scattermore::geom_scattermore()`.
#'   Default `0.5`.
#' @param title Plot title.
#' @param output.dir Output directory. Created if missing. Default
#'   `"figure_af_imaging"`.
#' @param save Logical, default `TRUE`.
#' @param file.type One of `"pdf"` (default), `"jpg"`, `"tiff"`, `"png"`.
#' @param width,height Figure dimensions in inches. Default `5`, `5`.
#'
#' @return The ggplot object.
#'
#' @export

plot.embedding.af.index <- function(
    embed.coords,
    af.index,
    asp,
    color.palette = "viridis",
    point.size    = 0.5,
    title         = "UMAP -- AF Index",
    output.dir    = "figure_af_imaging",
    save          = TRUE,
    file.type     = "pdf",
    width         = 5,
    height        = 5
) {

  if ( nrow( embed.coords ) != length( af.index ) )
    stop( "`embed.coords` and `af.index` must have the same number of events.",
          call. = FALSE )

  plot.df <- data.frame(
    x        = embed.coords[ , "UMAP1" ],
    y        = embed.coords[ , "UMAP2" ],
    af.index = as.numeric( af.index )
  )

  embed.plot <- ggplot( plot.df, aes( x, y, color = af.index ) ) +
    scattermore::geom_scattermore( pointsize = point.size, alpha = 0.8, na.rm = TRUE ) +
    scale_color_viridis_c( option = color.palette, name = "AF Index" ) +
    labs( x = "UMAP 1", y = "UMAP 2", title = title ) +
    theme_bw() +
    theme(
      plot.margin = margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      axis.ticks    = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text     = element_text( size = asp$figure.axis.text.size ),
      axis.title    = element_text( size = asp$figure.axis.title.size ),
      panel.border  = element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if ( save ) {
    if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
    .save.figure(
      embed.plot, file.path( output.dir, paste0( "embedding_af_index.", file.type ) ),
      width = width, height = height, file.type = file.type
    )
  }

  embed.plot
}

#' @title Test Association Between AF Index and Cluster Assignment
#'
#' @description
#' Chi-squared test of independence on the cluster x AF Index contingency
#' table, using a Monte Carlo simulated p-value (robust to the sparse cells
#' expected when many AF Index levels are spread across many clusters), plus
#' Cramer's V as an effect size.
#'
#' @param cluster.id Vector, length N, cluster/metacluster assignment.
#' @param af.index Vector, length N, the `"AF Index"` column. Treated as
#'   categorical (each value is a distinct AF library entry, not an ordered
#'   quantity).
#' @param n.sim Integer, number of Monte Carlo replicates for the simulated
#'   p-value. Default `9999`.
#' @param seed Integer, random seed. Default `1`.
#'
#' @return A named list: `contingency` (the table), `statistic` (Pearson
#'   chi-squared), `p.value` (Monte Carlo simulated), `cramers.v`.
#'
#' @export

test.af.index.cluster.chisq <- function(
    cluster.id,
    af.index,
    n.sim = 9999,
    seed  = 1
) {

  if ( length( cluster.id ) != length( af.index ) )
    stop( "`cluster.id` and `af.index` must be the same length.", call. = FALSE )

  contingency <- table(
    cluster  = factor( cluster.id ),
    af.index = factor( af.index )
  )

  set.seed( seed )
  chisq.result <- stats::chisq.test( contingency, simulate.p.value = TRUE, B = n.sim )

  n <- sum( contingency )
  k <- min( dim( contingency ) )
  cramers.v <- sqrt( as.numeric( chisq.result$statistic ) / ( n * ( k - 1 ) ) )

  list(
    contingency = contingency,
    statistic   = as.numeric( chisq.result$statistic ),
    p.value     = chisq.result$p.value,
    cramers.v   = cramers.v
  )
}

#' @title Test Whether AF Index Groups Occupy Different Embedding Locations
#'
#' @description
#' Permutation MANOVA (Pillai's trace) testing whether the 2D embedding
#' location of an event depends on its (categorical) AF Index. Complemented
#' by, for each AF Index group, its Mahalanobis distance from the grand
#' centroid in embedding space, using the pooled within-group covariance as
#' the metric -- this is the effect-size table (which groups are actually
#' displaced, and by how much), while the permutation p-value is the omnibus
#' significance call.
#'
#' @param embed.coords Matrix, events x 2 (`UMAP1`, `UMAP2`).
#' @param af.index Vector, length `nrow(embed.coords)`, the `"AF Index"`
#'   column for the same events. Treated as categorical.
#' @param n.perm Integer, permutation replicates. Default `999`.
#' @param seed Integer, random seed. Default `1`.
#' @param min.group.n Integer, AF Index groups with fewer than this many
#'   events are dropped before testing. Default `10`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return A named list: `pillai.observed`, `p.value` (permutation, primary
#'   result), `n.perm`, `null.distribution`, `mahalanobis.table` (data frame,
#'   one row per retained AF Index group, sorted by distance descending),
#'   `parametric.summary` (the `stats::manova()` Pillai-trace ANOVA table for
#'   the same data, for reference only -- embedding coordinates are not
#'   remotely Gaussian, so the permutation `p.value` above is the one to
#'   report).
#'
#' @export

test.af.index.embedding.location <- function(
    embed.coords,
    af.index,
    n.perm      = 999,
    seed        = 1,
    min.group.n = 10,
    verbose     = TRUE
) {

  if ( nrow( embed.coords ) != length( af.index ) )
    stop( "`embed.coords` and `af.index` must have the same number of events.",
          call. = FALSE )

  group   <- factor( af.index )
  group.n <- table( group )
  keep.groups <- names( group.n )[ group.n >= min.group.n ]

  if ( length( keep.groups ) < 2 )
    stop( "Fewer than two AF Index groups meet `min.group.n`; cannot test association.",
          call. = FALSE )

  keep.idx <- group %in% keep.groups
  y        <- as.matrix( embed.coords[ keep.idx, c( "UMAP1", "UMAP2" ), drop = FALSE ] )
  group    <- droplevels( group[ keep.idx ] )

  if ( verbose )
    message( sprintf(
      "Testing %d AF Index groups (%d of %d events retained, min.group.n = %d)",
      nlevels( group ), sum( keep.idx ), length( keep.idx ), min.group.n
    ) )

  pillai.trace <- function( y, group ) {
    grand.mean  <- colMeans( y )
    n.g         <- as.vector( table( group ) )
    group.means <- rowsum( y, group ) / n.g

    between <- crossprod( sqrt( n.g ) * sweep( group.means, 2, grand.mean, "-" ) )
    total   <- crossprod( sweep( y, 2, grand.mean, "-" ) )

    sum( diag( between %*% solve( total ) ) )
  }

  set.seed( seed )
  pillai.obs <- pillai.trace( y, group )

  n.events    <- nrow( y )
  pillai.null <- numeric( n.perm )
  for ( i in seq_len( n.perm ) ) {
    perm.group      <- group[ sample( n.events ) ]
    pillai.null[ i ] <- pillai.trace( y, perm.group )
  }

  p.value <- ( 1 + sum( pillai.null >= pillai.obs ) ) / ( n.perm + 1 )

  grand.mean  <- colMeans( y )
  group.means <- rowsum( y, group ) / as.vector( table( group ) )
  resid       <- y - group.means[ as.character( group ), ]
  pooled.cov  <- stats::cov( resid ) * ( n.events - 1 ) / ( n.events - nlevels( group ) )

  mahal.dist <- apply( group.means, 1, function( m ) {
    d <- m - grand.mean
    sqrt( as.numeric( t( d ) %*% solve( pooled.cov ) %*% d ) )
  } )

  mahal.table <- data.frame(
    af.index              = rownames( group.means ),
    n                      = as.integer( table( group ) ),
    mahalanobis.distance  = mahal.dist,
    row.names              = NULL
  )
  mahal.table <- mahal.table[ order( -mahal.table$mahalanobis.distance ), ]

  parametric.summary <- tryCatch(
    summary( stats::manova( y ~ group ), test = "Pillai" )$stats,
    error = function( e ) NULL
  )

  list(
    pillai.observed     = pillai.obs,
    p.value              = p.value,
    n.perm                = n.perm,
    null.distribution     = pillai.null,
    mahalanobis.table     = mahal.table,
    parametric.summary    = parametric.summary
  )
}

#' @title Test AF Index Differences Across Clusters (ANOVA)
#'
#' @description
#' One-way ANOVA testing whether the (quasi-continuous, similarity-ranked)
#' AF Index differs across imaging clusters. Uses Welch's ANOVA by default,
#' which doesn't assume equal variance or equal N across clusters --
#' reasonable here since cluster sizes and AF spread are unlikely to match.
#'
#' @param af.index Numeric vector, length N -- the similarity-ranked AF
#'   Index (see `af.index.to.continuous()`).
#' @param cluster.id Vector, length N, cluster/metacluster assignment for
#'   the same events.
#' @param var.equal Logical, default `FALSE` (Welch's ANOVA). Set `TRUE` for
#'   the classic equal-variance F-test.
#'
#' @return A named list: `statistic` (F), `df1`, `df2`, `p.value`, `method`.
#'
#' @export

test.af.index.cluster.anova <- function( af.index, cluster.id, var.equal = FALSE ) {

  if ( length( af.index ) != length( cluster.id ) )
    stop( "`af.index` and `cluster.id` must be the same length.", call. = FALSE )

  keep <- !is.na( af.index ) & !is.na( cluster.id )

  test.result <- stats::oneway.test(
    af.index[ keep ] ~ factor( cluster.id[ keep ] ), var.equal = var.equal
  )

  list(
    statistic = as.numeric( test.result$statistic ),
    df1       = as.numeric( test.result$parameter[ 1 ] ),
    df2       = as.numeric( test.result$parameter[ 2 ] ),
    p.value   = test.result$p.value,
    method    = test.result$method
  )
}

#' @title Violin Plot of AF Index Distribution by Cluster
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs
#' @importFrom ggplot2 theme_bw theme margin element_line element_text
#' @importFrom ggplot2 element_rect element_blank scale_fill_manual scale_fill_viridis_d
#'
#' @param af.index Numeric vector, length N -- the similarity-ranked, quasi-
#'   continuous AF Index for each event.
#' @param cluster.id Vector, length N, cluster/metacluster assignment for
#'   the same events.
#' @param cluster.labels Optional named character vector renaming cluster
#'   ids for display, same convention as `plot.embedding.clusters()`.
#' @param cluster.levels Optional vector of raw cluster ids giving the full,
#'   canonical set and order of clusters (e.g.
#'   `sort(unique(cluster.result$metacluster.id))`). Pass the same value
#'   here and to `plot.embedding.clusters()` (and `plot.flowsom.stars()`) so
#'   violin fill colors match the embedding plot's cluster colors exactly,
#'   even if this call's `cluster.id` happens not to include every cluster.
#'   Default `NULL` derives levels from `cluster.id` alone, which is only
#'   guaranteed consistent with another plot if that plot saw the same set
#'   of clusters.
#' @param palette.type One of `"brewer"` (default) or `"viridis"`. Must
#'   match whatever `plot.embedding.clusters()` was called with for the
#'   fill colors to agree.
#' @param brewer.palette RColorBrewer qualitative palette name, used when
#'   `palette.type = "brewer"`. Default `"Paired"`.
#' @param viridis.palette Viridis option, used when `palette.type =
#'   "viridis"`. Default `"turbo"`.
#' @param show.boxplot Logical, default `TRUE`. Overlay a narrow white
#'   boxplot for contrast against the colored violins.
#' @param var.equal Passed to `test.af.index.cluster.anova()`. Default
#'   `FALSE` (Welch's ANOVA).
#' @param title Plot title.
#' @param y.lab Y-axis label. Default `"AF Index (similarity rank)"`.
#' @param output.dir Output directory. Created if missing. Default
#'   `"figure_af_imaging"`.
#' @param save Logical, default `TRUE`.
#' @param file.type One of `"pdf"` (default), `"jpg"`, `"tiff"`, `"png"`.
#' @param width,height Figure dimensions in inches. `width` default `NULL`
#'   auto-scales to the number of clusters.
#'
#' @return A named list: `plot` (the ggplot object, with the ANOVA result as
#'   its subtitle) and `anova` (the full result of
#'   `test.af.index.cluster.anova()`).
#'
#' @export

plot.af.index.violin <- function(
    af.index,
    cluster.id,
    cluster.labels  = NULL,
    cluster.levels  = NULL,
    asp,
    palette.type     = c( "brewer", "viridis" ),
    brewer.palette   = "Paired",
    viridis.palette  = "turbo",
    show.boxplot    = TRUE,
    var.equal       = FALSE,
    title           = "AF Index by cluster",
    y.lab           = "AF Index (similarity rank)",
    output.dir      = "figure_af_imaging",
    save            = TRUE,
    file.type       = "pdf",
    width           = NULL,
    height          = 5
) {

  if ( length( af.index ) != length( cluster.id ) )
    stop( "`af.index` and `cluster.id` must be the same length.", call. = FALSE )

  palette.type <- match.arg( palette.type )

  anova.result <- test.af.index.cluster.anova( af.index, cluster.id, var.equal = var.equal )

  cluster.factor <- factor( cluster.id, levels = cluster.levels )
  n.clusters     <- nlevels( cluster.factor )

  if ( !is.null( cluster.labels ) ) {
    lvl <- levels( cluster.factor )
    levels( cluster.factor ) <- ifelse(
      lvl %in% names( cluster.labels ), cluster.labels[ lvl ], lvl
    )
  }

  plot.df <- data.frame( cluster = cluster.factor, af.index = af.index )
  plot.df <- plot.df[ !is.na( plot.df$af.index ) & !is.na( plot.df$cluster ), ]
  plot.df$cluster <- droplevels( plot.df$cluster )

  if ( is.null( width ) ) width <- max( 5, 0.6 * n.clusters + 2 )

  p.value.label <- sprintf(
    "%s: F(%.0f, %.0f) = %.2f, p = %s",
    ifelse( var.equal, "ANOVA", "Welch ANOVA" ),
    anova.result$df1, anova.result$df2, anova.result$statistic,
    format.pval( anova.result$p.value, eps = 1e-4, digits = 3 )
  )

  violin.plot <- ggplot( plot.df, aes( cluster, af.index ) ) +
    geom_violin( aes( fill = cluster ), color = "black", scale = "width", na.rm = TRUE )

  violin.plot <- if ( palette.type == "brewer" )
    violin.plot + scale_fill_manual( values = .get.cluster.palette( n.clusters, brewer.palette ) )
  else
    violin.plot + scale_fill_viridis_d( option = viridis.palette )

  if ( show.boxplot )
    violin.plot <- violin.plot +
    geom_boxplot( width = 0.12, fill = "white", outlier.shape = NA, na.rm = TRUE )

  violin.plot <- violin.plot +
    labs( x = NULL, y = y.lab, title = title, subtitle = p.value.label ) +
    theme_bw() +
    theme(
      plot.margin = margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      legend.position = "none",
      axis.ticks    = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text     = element_text( size = asp$figure.axis.text.size ),
      axis.title    = element_text( size = asp$figure.axis.title.size ),
      axis.text.x   = element_text( angle = 45, hjust = 1 ),
      panel.border  = element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )

  if ( save ) {
    if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
    .save.figure(
      violin.plot, file.path( output.dir, paste0( "af_index_violin_by_cluster.", file.type ) ),
      width = width, height = height, file.type = file.type
    )
  }

  list( plot = violin.plot, anova = anova.result )
}


