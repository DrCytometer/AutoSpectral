# assign_af_scatter_match.R

#' @title Assign AF Spectrum By Scatter-Matched Reference Averaging
#'
#' @description
#' An experimental AF assignment approach for benchmarking against the existing
#' covariance/residual methods (\code{assign.af.fluorophores},
#' \code{assign.af.residuals}, \code{assign.af.joint.cov}).
#'
#' For each cell in \code{test.data}, the \code{k} nearest neighbours in
#' \code{ref.data} are found by Euclidean distance in scatter space (using
#' \code{FNN::knnx.index}). Their spectral channels are averaged to produce a
#' cell-specific reference spectrum. The cosine similarity between each test
#' cell's spectral signature and its scatter-matched reference average is
#' returned alongside the best-fitting AF variant index, allowing direct
#' comparison against the existing assignment approaches.
#'
#' The rationale: cells with identical scatter profiles (size, granularity)
#' should have near-identical autofluorescence. If the scatter-matched average
#' is a better proxy for the true AF than the library-based assignment, the
#' returned cosine similarities will be systematically higher than those
#' obtained from the covariance/residual pipeline.
#'
#' @param test.data Numeric matrix. Expression data from the **test** unstained
#'   FCS file. Cells in rows, channels in columns. Must contain both
#'   \code{scatter.param} channels and the spectral detector channels. Can also
#'   be the path to an FCS file (character scalar), in which case the file is
#'   read with \code{readFCS()}.
#' @param ref.data Numeric matrix. Expression data from the **reference**
#'   unstained FCS file. Same column layout as \code{test.data}. Can also be
#'   an FCS file path.
#' @param scatter.param Character vector of length >= 1. Column names that
#'   identify the scatter channels to use for kNN matching (e.g.
#'   \code{c("FSC-A", "SSC-A")}). These are excluded from the spectral
#'   similarity calculation.
#' @param k Integer. Number of nearest reference neighbours to average.
#'   Default \code{5}. Larger values stabilise the reference estimate but may
#'   blur genuine AF heterogeneity.
#' @param af.spectra Optional numeric matrix. Spectral signatures of AF
#'   variants, normalised `[0, 1]`, with variants in rows and detectors in
#'   columns. When supplied, the function also assigns each test cell to the
#'   closest AF variant (by cosine similarity to the scatter-matched average)
#'   and returns that index. Omit to skip variant assignment.
#' @param scale.scatter Logical. Whether to z-score-standardise the scatter
#'   channels before computing kNN distances (recommended when FSC and SSC
#'   have very different dynamic ranges). Default \code{TRUE}.
#' @param verbose Logical. Whether to emit progress messages. Default
#'   \code{TRUE}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{\code{cosine.similarity}}{Numeric vector, length
#'       \code{nrow(test.data)}. Per-cell cosine similarity between the test
#'       cell's spectral signature and its scatter-matched reference average.
#'       Values close to 1 indicate a strong match.}
#'     \item{\code{ref.average}}{Numeric matrix, same dimensions as the
#'       spectral portion of \code{test.data}. The scatter-matched averaged
#'       reference spectrum for each test cell (i.e., the mean of the \code{k}
#'       nearest reference neighbours).}
#'     \item{\code{af.assignment}}{Integer vector (or \code{NULL} if
#'       \code{af.spectra} was not supplied). Row index into \code{af.spectra}
#'       of the best-fitting AF variant for each test cell, determined by
#'       cosine similarity between the scatter-matched reference average and
#'       each AF library spectrum.}
#'     \item{\code{nn.index}}{Integer matrix, \code{nrow(test.data)} x
#'       \code{k}. The row indices in \code{ref.data} of the \code{k} nearest
#'       scatter neighbours for each test cell. Useful for diagnostics.}
#'     \item{\code{summary}}{A one-row data frame with aggregate statistics
#'       (mean, median, sd of cosine similarities) for quick comparison
#'       against competing methods.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Using FCS file paths:
#' result <- assign.af.scatter.match(
#'   test.data    = "path/to/test_unstained.fcs",
#'   ref.data     = "path/to/ref_unstained.fcs",
#'   scatter.param = c( "FSC-A", "SSC-A" ),
#'   k            = 5,
#'   af.spectra   = my.af.spectra   # optional
#' )
#' hist( result$cosine.similarity, main = "Scatter-match cosine similarity" )
#'
#' # Or pass matrices directly:
#' result <- assign.af.scatter.match(
#'   test.data    = test.mat,
#'   ref.data     = ref.mat,
#'   scatter.param = c( "FSC-A", "SSC-A" ),
#'   k            = 5
#' )
#' }
#'
#' @importFrom FNN knnx.index
#' @export

assign.af.scatter.match <- function(
    test.data,
    ref.data,
    scatter.param,
    k             = 5L,
    af.spectra    = NULL,
    scale.scatter = TRUE,
    verbose       = TRUE
) {

  # ---------------------------------------------------------------------------
  # 1. Input ingestion: accept file paths or pre-loaded matrices
  # ---------------------------------------------------------------------------

  if ( is.character( test.data ) ) {
    if ( verbose ) message( "Reading test FCS: ", test.data )
    test.data <- readFCS( test.data )
  }

  if ( is.character( ref.data ) ) {
    if ( verbose ) message( "Reading reference FCS: ", ref.data )
    ref.data <- readFCS( ref.data )
  }

  test.data <- as.matrix( test.data )
  ref.data  <- as.matrix( ref.data )

  # ---------------------------------------------------------------------------
  # 2. Validate columns
  # ---------------------------------------------------------------------------

  missing.scatter.test <- setdiff( scatter.param, colnames( test.data ) )
  missing.scatter.ref  <- setdiff( scatter.param, colnames( ref.data  ) )

  if ( length( missing.scatter.test ) > 0 )
    stop(
      "scatter.param columns not found in test.data: ",
      paste( missing.scatter.test, collapse = ", " ),
      call. = FALSE
    )

  if ( length( missing.scatter.ref ) > 0 )
    stop(
      "scatter.param columns not found in ref.data: ",
      paste( missing.scatter.ref, collapse = ", " ),
      call. = FALSE
    )

  # Spectral columns = everything that is NOT a scatter parameter
  spectral.cols <- setdiff( colnames( test.data ), scatter.param )

  missing.spectral <- setdiff( spectral.cols, colnames( ref.data ) )
  if ( length( missing.spectral ) > 0 )
    stop(
      "Spectral columns present in test.data are absent from ref.data: ",
      paste( missing.spectral, collapse = ", " ),
      call. = FALSE
    )

  if ( length( spectral.cols ) == 0 )
    stop( "No spectral columns remain after removing scatter.param.", call. = FALSE )

  k <- as.integer( k )
  if ( k < 1L ) stop( "'k' must be >= 1.", call. = FALSE )
  if ( k > nrow( ref.data ) )
    stop(
      "'k' (", k, ") exceeds the number of reference cells (",
      nrow( ref.data ), ").",
      call. = FALSE
    )

  # ---------------------------------------------------------------------------
  # 3. Extract sub-matrices
  # ---------------------------------------------------------------------------

  test.scatter  <- test.data[ , scatter.param,  drop = FALSE ]
  ref.scatter   <- ref.data[  , scatter.param,  drop = FALSE ]
  test.spectral <- test.data[ , spectral.cols,  drop = FALSE ]
  ref.spectral  <- ref.data[  , spectral.cols,  drop = FALSE ]

  n.test <- nrow( test.data )

  if ( verbose ) {
    message(
      sprintf(
        "Scatter-match kNN: %d test cells  |  %d reference cells  |  k = %d  |  %d spectral channels",
        n.test, nrow( ref.data ), k, length( spectral.cols )
      )
    )
  }

  # ---------------------------------------------------------------------------
  # 4. Optional scatter standardisation
  # ---------------------------------------------------------------------------

  if ( scale.scatter ) {
    scatter.means <- colMeans( ref.scatter )
    scatter.sds   <- apply( ref.scatter, 2, stats::sd )
    scatter.sds[ scatter.sds < 1e-9 ] <- 1  # guard against zero-variance channels

    test.scatter.scaled <- sweep( sweep( test.scatter, 2, scatter.means, "-" ),
                                  2, scatter.sds, "/" )
    ref.scatter.scaled  <- sweep( sweep( ref.scatter,  2, scatter.means, "-" ),
                                  2, scatter.sds, "/" )
  } else {
    test.scatter.scaled <- test.scatter
    ref.scatter.scaled  <- ref.scatter
  }

  # ---------------------------------------------------------------------------
  # 5. kNN in scatter space
  # ---------------------------------------------------------------------------

  if ( verbose ) message( "Computing kNN index in scatter space ..." )

  nn.index <- FNN::knnx.index(
    data  = ref.scatter.scaled,
    query = test.scatter.scaled,
    k     = k,
    algorithm = "kd_tree"
  )   # n.test x k integer matrix

  # ---------------------------------------------------------------------------
  # 6. Average spectral channels across the k reference neighbours
  # ---------------------------------------------------------------------------

  if ( verbose ) message( "Averaging reference spectra over k neighbours ..." )

  # Efficient: for each test cell, sum k rows of ref.spectral then divide by k.
  # Avoid a loop by using index arithmetic on the full flat index vector.
  ref.spectral.mat <- as.matrix( ref.spectral )

  # nn.index is n.test x k; flatten to a single vector for bulk row extraction
  flat.idx  <- as.vector( nn.index )                        # length n.test * k
  flat.vals <- ref.spectral.mat[ flat.idx, , drop = FALSE ] # (n.test*k) x D

  # Reshape to (n.test, k, D) by summing in blocks of k
  # Using matrix arithmetic: group rows into blocks of k and colMeans each block
  D           <- ncol( ref.spectral.mat )
  block.sums  <- matrix( 0, nrow = n.test, ncol = D )

  # Vectorised summation across neighbours
  for ( j in seq_len( k ) ) {
    block.sums <- block.sums + ref.spectral.mat[ nn.index[ , j ], , drop = FALSE ]
  }
  ref.average <- block.sums / k
  colnames( ref.average ) <- spectral.cols

  # ---------------------------------------------------------------------------
  # 7. Per-cell cosine similarity: test spectral vs scatter-matched average
  # ---------------------------------------------------------------------------

  if ( verbose ) message( "Computing cosine similarities ..." )

  # row norms with epsilon guard
  test.norms <- sqrt( rowSums( test.spectral^2 ) ) + 1e-12
  ref.norms  <- sqrt( rowSums( ref.average^2   ) ) + 1e-12

  dot.products     <- rowSums( test.spectral * ref.average )
  cosine.sim       <- dot.products / ( test.norms * ref.norms )

  # Clamp to [-1, 1] for numerical safety
  cosine.sim <- pmax( pmin( cosine.sim, 1 ), -1 )

  # ---------------------------------------------------------------------------
  # 8. Optional: assign to closest AF library variant
  # ---------------------------------------------------------------------------

  af.assignment <- NULL

  if ( !is.null( af.spectra ) ) {

    if ( verbose ) message( "Assigning AF variants by cosine similarity to reference average ..." )

    af.spectra <- as.matrix( af.spectra )

    # Align columns: use intersection, in the order of spectral.cols
    common.cols <- intersect( spectral.cols, colnames( af.spectra ) )

    if ( length( common.cols ) == 0 )
      stop(
        "No shared columns between spectral channels and af.spectra. ",
        "Check that af.spectra columns match the detector names in test.data.",
        call. = FALSE
      )

    ref.avg.sub  <- ref.average[ , common.cols, drop = FALSE ]
    af.spec.sub  <- af.spectra[  , common.cols, drop = FALSE ]

    # Cosine similarity: ref.average (n.test x D) vs af.spectra (n.af x D)
    # Result: n.test x n.af matrix
    norm.row <- function( m ) {
      nrm <- sqrt( rowSums( m^2 ) ) + 1e-12
      m / nrm
    }

    sim.matrix    <- norm.row( ref.avg.sub ) %*% t( norm.row( af.spec.sub ) )
    af.assignment <- max.col( sim.matrix, ties.method = "first" )
  }

  # ---------------------------------------------------------------------------
  # 9. Summary statistics
  # ---------------------------------------------------------------------------

  summary.df <- data.frame(
    n.test.cells   = n.test,
    k              = k,
    mean.cosine    = mean(   cosine.sim ),
    median.cosine  = stats::median( cosine.sim ),
    sd.cosine      = stats::sd(     cosine.sim ),
    pct.above.0.9  = mean(   cosine.sim > 0.9 ) * 100,
    pct.above.0.95 = mean(   cosine.sim > 0.95 ) * 100,
    stringsAsFactors = FALSE
  )

  if ( verbose ) {
    message(
      sprintf(
        "  Mean cosine = %.4f  |  Median = %.4f  |  SD = %.4f  |  >0.9: %.1f%%  |  >0.95: %.1f%%",
        summary.df$mean.cosine,
        summary.df$median.cosine,
        summary.df$sd.cosine,
        summary.df$pct.above.0.9,
        summary.df$pct.above.0.95
      )
    )
  }

  # ---------------------------------------------------------------------------
  # 10. Return
  # ---------------------------------------------------------------------------

  return( list(
    cosine.similarity = cosine.sim,
    ref.average       = ref.average,
    af.assignment     = af.assignment,
    nn.index          = nn.index,
    summary           = summary.df
  ) )
}


# ---------------------------------------------------------------------------
# Benchmarking helper
# ---------------------------------------------------------------------------

#' @title Benchmark Scatter-Match Against Existing AF Assignment Methods
#'
#' @description
#' Convenience wrapper that runs \code{assign.af.scatter.match} alongside the
#' three existing methods (\code{assign.af.fluorophores},
#' \code{assign.af.residuals}, \code{assign.af.joint.cov}) on the same test
#' and reference unstained data. For each existing method, the assigned AF
#' variant spectrum is looked up and its cosine similarity to each test cell
#' is computed, allowing direct apples-to-apples comparison with the
#' scatter-match approach.
#'
#' @param test.data Numeric matrix or FCS file path. Test unstained data
#'   (cells x channels).
#' @param ref.data Numeric matrix or FCS file path. Reference unstained data
#'   (cells x channels).
#' @param scatter.param Character vector of scatter channel names.
#' @param spectra Numeric matrix. Fluorophore spectra (fluorophores x
#'   detectors), as used by the existing assign.af.* functions.
#' @param af.spectra Numeric matrix. AF variant spectra (variants x
#'   detectors).
#' @param k Integer. Neighbours for scatter-matching. Default \code{5}.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{scatter.match}}{Full output of
#'       \code{assign.af.scatter.match}.}
#'     \item{\code{comparison}}{Data frame with one row per method and
#'       columns: \code{method}, \code{mean.cosine}, \code{median.cosine},
#'       \code{sd.cosine}, \code{pct.above.0.9}, \code{pct.above.0.95}.}
#'     \item{\code{per.cell}}{Data frame with one row per test cell containing
#'       cosine similarities from all four methods, for cell-level analysis.}
#'   }
#'
#' @export

benchmark.af.scatter.match <- function(
    test.data,
    ref.data,
    scatter.param,
    spectra,
    af.spectra,
    k       = 5L,
    verbose = TRUE
) {

  # ---- Load data if paths ----
  if ( is.character( test.data ) ) test.data <- readFCS( test.data )
  if ( is.character( ref.data  ) ) ref.data  <- readFCS( ref.data  )
  test.data <- as.matrix( test.data )

  spectral.cols <- setdiff( colnames( test.data ), scatter.param )
  test.spectral <- test.data[ , spectral.cols, drop = FALSE ]

  # ---- Helper: cosine similarity of each test cell against its assigned AF ----
  cosine.from.assignment <- function( assignments ) {
    assigned.spectra <- af.spectra[ assignments, , drop = FALSE ]
    # align columns
    common <- intersect( spectral.cols, colnames( af.spectra ) )
    ts <- test.spectral[ , common, drop = FALSE ]
    as <- assigned.spectra[ , common, drop = FALSE ]
    t.norms <- sqrt( rowSums( ts^2 ) ) + 1e-12
    a.norms <- sqrt( rowSums( as^2 ) ) + 1e-12
    pmax( pmin( rowSums( ts * as ) / ( t.norms * a.norms ), 1 ), -1 )
  }

  # ---- Scatter-match ----
  if ( verbose ) message( "\n--- Method: scatter.match ---" )
  sm.result <- assign.af.scatter.match(
    test.data     = test.data,
    ref.data      = ref.data,
    scatter.param = scatter.param,
    k             = k,
    af.spectra    = af.spectra,
    verbose       = verbose
  )

  # ---- Existing methods ----
  if ( verbose ) message( "\n--- Method: assign.af.fluorophores ---" )
  idx.fluor  <- assign.af.fluorophores( test.spectral, spectra, af.spectra )
  cs.fluor   <- cosine.from.assignment( idx.fluor )

  if ( verbose ) message( "\n--- Method: assign.af.residuals ---" )
  idx.resid  <- assign.af.residuals( test.spectral, spectra, af.spectra )
  cs.resid   <- cosine.from.assignment( idx.resid )

  if ( verbose ) message( "\n--- Method: assign.af.joint.cov ---" )
  idx.joint  <- assign.af.joint.cov( test.spectral, spectra, af.spectra )
  cs.joint   <- cosine.from.assignment( idx.joint )

  # ---- Summarise ----
  summarise.cs <- function( cs, method.name ) {
    data.frame(
      method         = method.name,
      mean.cosine    = mean(   cs ),
      median.cosine  = stats::median( cs ),
      sd.cosine      = stats::sd(     cs ),
      pct.above.0.9  = mean(   cs > 0.9  ) * 100,
      pct.above.0.95 = mean(   cs > 0.95 ) * 100,
      stringsAsFactors = FALSE
    )
  }

  comparison <- rbind(
    summarise.cs( sm.result$cosine.similarity, "scatter.match"          ),
    summarise.cs( cs.fluor,                    "assign.af.fluorophores" ),
    summarise.cs( cs.resid,                    "assign.af.residuals"    ),
    summarise.cs( cs.joint,                    "assign.af.joint.cov"    )
  )

  if ( verbose ) {
    message( "\n===== Benchmark summary =====" )
    print( comparison, digits = 4, row.names = FALSE )
  }

  per.cell <- data.frame(
    scatter.match          = sm.result$cosine.similarity,
    assign.af.fluorophores = cs.fluor,
    assign.af.residuals    = cs.resid,
    assign.af.joint.cov    = cs.joint
  )

  return( list(
    scatter.match = sm.result,
    comparison    = comparison,
    per.cell      = per.cell
  ) )
}
