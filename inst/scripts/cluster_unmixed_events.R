# cluster_unmixed_events.R

#' @title Cluster Unmixed Events into Micro-Clusters
#'
#' @description
#' Groups events into small, homogeneous micro-clusters based on their
#' currently-unmixed abundance profile, and computes a per-cluster centroid
#' abundance vector and matched raw-detector-space centroid vector. This is
#' Phase 1 of the signature-error-correction workflow: it stabilises
#' cluster-level autofluorescence/background variation and prevents a few
#' dominant populations from swamping the downstream regression, while still
#' sampling the abundance range of every fluorophore.
#'
#' @importFrom FNN knnx.index get.knn
#'
#' @param unmixed Numeric matrix (events x fluorophores), OLS/WLS unmixed
#'   abundance for the sample under investigation. Typically the output of
#'   `unmix.ols.fast()` on `raw.data` with the current `spectra`.
#' @param raw.data Numeric matrix (events x detectors), the raw detector-space
#'   data that produced `unmixed`. Column order and names must match the
#'   detectors of `spectra`.
#' @param asp The AutoSpectral parameter list from `get.autospectral.param()`.
#'   Used for `asp$bird.seed`.
#' @param method Character, `"som"` (default), `"leiden"`, or `"kmeans"`.
#'   `"som"` uses the package's own batch SOM (`get.som.codes()`, the
#'   `AutoSpectralRcpp`-accelerated engine already used by `get.af.spectra()`
#'   and `get.fluor.variants()`, not `FlowSOM`/`EmbedSOM` directly), with
#'   events assigned to their nearest node via `FNN::knnx.index()`.
#'   `"leiden"` builds a k-nearest-neighbour graph (`FNN::get.knn()`) and
#'   clusters it with `igraph::cluster_leiden()` - `igraph` is not currently
#'   a package dependency, so this method needs it installed separately for
#'   now. `"kmeans"` gives an exact cluster count and is simplest to reason
#'   about, but slower at high event counts and high k.
#' @param cluster.id Optional integer vector, length `nrow(unmixed)`, a
#'   previously-computed cluster assignment (e.g. `$cluster.id` from an
#'   earlier call). When supplied, clustering is skipped entirely and this
#'   assignment is used directly to compute centroids - intended for reuse
#'   across correction iterations in `correct.unmixing.signatures()`, so
#'   cluster identity stays fixed while only `unmixed` (and therefore
#'   `x.clust`) changes as `spectra` is corrected. Default `NULL` (cluster
#'   fresh).
#' @param n.clusters Target cluster count for `method = "kmeans"`. Default
#'   `150`. Not used by `"som"` (see `som.dim`) or `"leiden"` (see
#'   `resolution`).
#' @param som.dim Integer, SOM grid side length for `method = "som"`
#'   (`som.dim^2` nodes). Default `20`.
#' @param k.neighbors Integer, neighbours per event for the k-NN graph in
#'   `method = "leiden"`. Default `15`.
#' @param resolution Numeric, resolution parameter passed to
#'   `igraph::cluster_leiden()` for `method = "leiden"`. Higher values give
#'   more, smaller clusters. Default `1`. Unlike `n.clusters`/`som.dim`, this
#'   does not target an exact cluster count directly - see Exploration
#'   points.
#' @param unmixed.thresholds Named numeric vector of per-fluorophore positivity
#'   thresholds in unmixed space. Required for `method = "stratify"`, ignored
#'   otherwise.
#' @param n.levels Integer, number of abundance bins per fluorophore for
#'   `method = "stratify"`. Default `10`.
#' @param stratify.exclude Optional character vector of fluorophore names that
#'   are not eligible to be an event's dominant fluorophore under
#'   `method = "stratify"`. Autofluorescence belongs here: it is present in
#'   every event and would otherwise dominate every assignment. Default `NULL`.
#' @param centroid.fun Character, `"median"` (default) or `"mean"`. Function
#'   used to summarise each cluster's abundance and raw-data vectors. Median
#'   is more robust to within-cluster outliers.
#' @param min.cluster.size Integer. Clusters with fewer than this many events
#'   are dropped before the regression step, since their centroids are
#'   unstable. Default `5`.
#' @param cluster.cofactor Numeric, asinh cofactor applied to `unmixed` before
#'   clustering only. Abundances span several decades, so linear-scale
#'   distances are dominated by the brightest population; clustering on
#'   `asinh(unmixed / cluster.cofactor)` distributes nodes across the whole
#'   dynamic range. Centroids are always computed in linear units. Default
#'   `500`.
#' @param threads Integer, OpenMP threads for `method = "som"`. Default
#'   `NULL` (all available cores, matching `get.som.codes()`'s own default).
#'   Ignored by other methods.
#' @param verbose Logical, controls messaging. Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`cluster.id`}{Integer vector, length `nrow(unmixed)`, the cluster
#'     assignment for every event (`NA` for events in clusters dropped by
#'     `min.cluster.size`).}
#'   \item{`x.clust`}{Numeric matrix (clusters x fluorophores), the
#'     per-cluster centroid abundance vector. Row names are cluster IDs.}
#'   \item{`y.clust`}{Numeric matrix (clusters x detectors), the per-cluster
#'     centroid raw-data vector, matched row-for-row to `x.clust`.}
#'   \item{`cluster.size`}{Integer vector, length `nrow(x.clust)`, the number
#'     of events contributing to each cluster centroid.}
#' }
#'
#' @section Exploration points:
#' \itemize{
#'   \item SOM vs Leiden vs kmeans vs stratify: SOM grid cells are not
#'     guaranteed to be evenly populated, which is especially visible on data
#'     with strong block structure (e.g. concatenated single-stained
#'     controls, where most of the grid may be near-empty for any one
#'     fluorophore's positive population); Leiden should adapt cluster count
#'     to the actual density structure instead of a fixed grid, at the cost
#'     of not targeting an exact cluster count directly; kmeans gives the
#'     most uniform cluster sizes but assumes roughly spherical clusters in
#'     abundance space; stratify guarantees every fluorophore the same
#'     number of populated abundance levels regardless of rarity, at the
#'     cost of needing `unmixed.thresholds` up front. Worth comparing all
#'     four on the same real data.
#'   \item `resolution` (Leiden) doesn't target an exact `n.clusters` the way
#'     `som.dim`/`n.clusters` do for the other methods. A wrapper that
#'     adjusts `resolution` until `length(unique(cluster.id))` lands near a
#'     target could be added if a specific count is needed - not implemented
#'     here.
#'   \item The Leiden graph is currently unweighted (plain k-NN edges); an
#'     inverse-distance-weighted graph is a common refinement worth trying.
#'   \item median vs mean centroids: median is the codebase default
#'     elsewhere, but for skewed clusters, mean might track a systematic
#'     error pattern more faithfully. Worth comparing both.
#' }
#'
#' @export

cluster.unmixed.events <- function(
    unmixed,
    raw.data,
    asp,
    method             = c( "som", "stratify", "leiden", "kmeans" ),
    cluster.id         = NULL,
    n.clusters         = 150,
    som.dim            = 20,
    k.neighbors        = 15,
    resolution         = 1,
    unmixed.thresholds = NULL,
    n.levels           = 10L,
    stratify.exclude   = NULL,
    centroid.fun     = c( "median", "mean" ),
    min.cluster.size = 5,
    cluster.cofactor = 500,
    threads          = NULL,
    verbose          = TRUE
) {
  
  method       <- match.arg( method )
  centroid.fun <- match.arg( centroid.fun )

  if ( !is.matrix( unmixed ) )  unmixed  <- as.matrix( unmixed )
  if ( !is.matrix( raw.data ) ) raw.data <- as.matrix( raw.data )

  if ( nrow( unmixed ) != nrow( raw.data ) )
    stop( "`unmixed` and `raw.data` must have the same number of rows (events).",
          call. = FALSE )

  event.n <- nrow( unmixed )

  # ---------------------------------------------------------------------------
  # Cluster assignment
  # ---------------------------------------------------------------------------

  if ( is.null( cluster.id ) ) {
    
    set.seed( asp$bird.seed )
    
    # Abundances span several decades, so linear-scale distances are dominated
    # entirely by the brightest population and dim events collapse into a
    # handful of nodes. Cluster on an asinh-stabilised copy; centroids below
    # are still taken in linear units, which is what the residual identity
    # requires.
    cluster.space <- asinh( unmixed / cluster.cofactor )
    
    if ( method == "som" ) {
      
      if ( verbose )
        message( sprintf( "Clustering %d events with a %dx%d SOM (%d nodes)",
                          event.n, som.dim, som.dim, som.dim^2 ) )
      
      map <- get.som.codes(
        data    = cluster.space,
        som.dim = som.dim,
        seed    = asp$bird.seed,
        threads = if ( is.null( threads ) ) 0L else as.integer( threads )
      )
      
      cluster.id <- as.integer( FNN::knnx.index( data = map$codes, query = cluster.space, k = 1 ) )
      
    } else if ( method == "stratify" ) {
      
      if ( is.null( unmixed.thresholds ) )
        stop( "method = \"stratify\" requires `unmixed.thresholds`.", call. = FALSE )
      
      strat.names <- setdiff( colnames( unmixed ), stratify.exclude )
      
      if ( !all( strat.names %in% names( unmixed.thresholds ) ) )
        stop( "`unmixed.thresholds` must cover every stratified column of `unmixed`.",
              call. = FALSE )
      
      # Assign each event to the fluorophore it is most strongly positive for,
      # measured in units of that fluorophore's own positivity threshold, then
      # split each group into n.levels quantile bins of abundance. Density-driven
      # clustering spends its nodes where the events are, which in pooled
      # single-stained controls is the shared negative population; this instead
      # guarantees every fluorophore the same number of populated abundance
      # levels however rare its positives are. Because dominance is relative,
      # an event positive for a spectrally collinear pair is still assigned to
      # whichever of the two it is more strongly positive for, rather than being
      # discarded as ambiguous.
      threshold.vec <- pmax( unmixed.thresholds[ strat.names ], .Machine$double.eps )
      score         <- sweep( unmixed[ , strat.names, drop = FALSE ], 2, threshold.vec, "/" )
      
      dominant  <- max.col( score, ties.method = "first" )
      dom.score <- score[ cbind( seq_len( event.n ), dominant ) ]
      dominant[ dom.score <= 1 ] <- 0L
      
      cluster.id <- integer( event.n )
      next.id    <- 0L
      
      for ( f in c( 0L, seq_along( strat.names ) ) ) {
        
        idx <- which( dominant == f )
        if ( length( idx ) == 0 ) next
        
        bin.value <- if ( f == 0L ) dom.score[ idx ] else unmixed[ idx, strat.names[ f ] ]
        brk       <- unique( stats::quantile(
          bin.value, probs = seq( 0, 1, length.out = n.levels + 1 ) ) )
        
        bin <- if ( length( brk ) < 2 ) rep( 1L, length( idx ) ) else
          as.integer( cut( bin.value, breaks = brk, include.lowest = TRUE ) )
        
        cluster.id[ idx ] <- next.id + bin
        next.id <- next.id + max( bin )
      }
      
      if ( verbose )
        message( sprintf(
          "Stratified %d events into %d bins across %d fluorophore group(s) plus background",
          event.n, next.id, length( unique( dominant[ dominant > 0 ] ) ) ) )
      
    } else if ( method == "leiden" ) {
      
      if ( !requireNamespace( "igraph", quietly = TRUE ) )
        stop( "method = \"leiden\" requires the igraph package (not currently ",
              "a package dependency -- install it to use this method).",
              call. = FALSE )
      
      if ( verbose )
        message( sprintf( "Clustering %d events with Leiden (k=%d neighbours, resolution=%.2f)",
                          event.n, k.neighbors, resolution ) )
      
      knn       <- FNN::get.knn( cluster.space, k = k.neighbors )
      edge.from <- rep( seq_len( event.n ), each = k.neighbors )
      edge.to   <- as.integer( t( knn$nn.index ) )
      graph     <- igraph::simplify(
        igraph::graph_from_edgelist( cbind( edge.from, edge.to ), directed = FALSE )
      )
      
      leiden.out <- igraph::cluster_leiden(
        graph, objective_function = "modularity", resolution = resolution, n_iterations = 2
      )
      cluster.id <- as.integer( igraph::membership( leiden.out ) )
      
    } else {
      
      if ( verbose )
        message( sprintf( "Clustering %d events into %d micro-clusters (kmeans)",
                          event.n, n.clusters ) )
      
      km <- stats::kmeans( cluster.space, centers = n.clusters, iter.max = 50, nstart = 1 )
      cluster.id <- km$cluster
      
    }
    
  } else {
    
    if ( length( cluster.id ) != event.n )
      stop( "`cluster.id` must have one entry per row of `unmixed`.", call. = FALSE )
    if ( verbose ) message( "Reusing supplied cluster.id (no reclustering)." )
    
  }

  # ---------------------------------------------------------------------------
  # Per-cluster centroids
  # ---------------------------------------------------------------------------

  centroid.function <- if ( centroid.fun == "median" ) stats::median else base::mean

  cluster.table <- table( cluster.id )
  keep.clusters <- as.integer( names( cluster.table )[ cluster.table >= min.cluster.size ] )

  if ( length( keep.clusters ) < ncol( unmixed ) + 1 )
    warning(
      paste0(
        "Only ", length( keep.clusters ), " clusters survived min.cluster.size=",
        min.cluster.size, "; downstream per-cluster regression diagnostics ",
        "need more clusters than fluorophores to be well-posed. Consider ",
        "lowering n.clusters, lowering min.cluster.size, or increasing event count."
      ),
      call. = FALSE
    )

  drop.mask <- !cluster.id %in% keep.clusters
  cluster.id[ drop.mask ] <- NA_integer_

  x.clust <- t( sapply( keep.clusters, function( cl ) {
    idx <- which( cluster.id == cl )
    apply( unmixed[ idx, , drop = FALSE ], 2, centroid.function )
  } ) )

  y.clust <- t( sapply( keep.clusters, function( cl ) {
    idx <- which( cluster.id == cl )
    apply( raw.data[ idx, , drop = FALSE ], 2, centroid.function )
  } ) )

  rownames( x.clust ) <- rownames( y.clust ) <- keep.clusters
  colnames( x.clust ) <- colnames( unmixed )
  colnames( y.clust ) <- colnames( raw.data )

  cluster.size <- as.integer( cluster.table[ as.character( keep.clusters ) ] )

  if ( verbose )
    message( sprintf( "Retained %d micro-clusters (min size %d)",
                       nrow( x.clust ), min.cluster.size ) )

  list(
    cluster.id   = cluster.id,
    x.clust      = x.clust,
    y.clust      = y.clust,
    cluster.size = cluster.size
  )
}
