# get_af_spectra_kmeans.r

#' @title Get Autofluorescence Spectra Using K-means
#'
#' @description
#' Extracts autofluorescence spectra from an unstained sample. Intended for use
#' with `unmix.autospectral`. Uses FlowSOM (EmbedSOM) clustering for rapid
#' identification of cells with similar AF profiles.
#'
#' Optionally deduplicates the resulting spectra by cosine similarity
#' (`deduplicate = TRUE`, default) to remove near-identical profiles that cause
#' spurious over-correction of near-zero events in fully stained samples. When
#' `refine = TRUE`, a second round of targeted modulation is performed on cells
#' that remain far from zero after the first-pass correction; modulated spectra
#' are screened for redundancy against each other and against the base library
#' before being appended.
#'
#' @importFrom FlowSOM SOM
#' @importFrom parallelly availableCores
#'
#' @param unstained.sample Path and file name for an unstained sample FCS file.
#'   The sample type and processing (protocol) method should match the fully
#'   stained samples to which the AF will be applied, ideally.
#' @param asp The AutoSpectral parameter list. Prepare using
#'   `get.autospectral.param`.
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and
#'   1, with fluorophores in rows and detectors in columns.
#' @param som.dim Number of x and y dimensions for the SOM. Default is `10`.
#' @param figures Logical, whether to plot the spectral traces and heatmap for
#'   the AF signatures. Default is `TRUE`.
#' @param plot.dir Directory (folder) where the plots will be saved. Default is
#'   `NULL`, which inherits from `asp$figure.af.dir`.
#' @param table.dir Directory (folder) where the spectra csv file will be saved.
#'   Default is `NULL`, which inherits from `asp$table.af.dir`.
#' @param title Title for the output spectral plots and csv file. Default is
#'   `"Autofluorescence spectra"`.
#' @param verbose Logical, controls messaging. Default is `TRUE`.
#' @param deduplicate Logical, default `TRUE`. Whether to deduplicate AF spectra
#'   by cosine similarity after the base clustering stage and again after the
#'   refinement stage. Deduplication removes near-identical spectral profiles
#'   that cause overzealous matching of near-zero events in fully stained samples,
#'   reducing over-correction artefacts. Set to `FALSE` to restore behavior
#'   identical to versions prior to deduplication support.
#' @param duplication.threshold Numeric, default `0.99`. The cosine similarity
#'   threshold used for deduplication. A spectrum is dropped if its cosine
#'   similarity to any already-retained spectrum meets or exceeds this value.
#'   Only used when `deduplicate = TRUE`.
#' @param refine Logical, default `FALSE`. Controls whether to perform a second
#'   round of autofluorescence measurement on "problem cells": those with the
#'   highest residual fluorophore signal after the first-pass per-cell AF
#'   extraction, as defined by `problem.quantile`. When `FALSE`, behavior is
#'   identical to versions of AutoSpectral prior to 1.0.0. If you are working
#'   with samples containing complex autofluorescence, e.g. tissues or tumors,
#'   using `refine = TRUE` will improve autofluorescence extraction at the cost
#'   of an increase in unmixing time.
#' @param problem.quantile Numeric, default `0.99`. The quantile for determining
#'   which cells are "problematic" after first-pass per-cell AF extraction. Cells
#'   at or above this quantile with respect to the L2 norm of their unmixed
#'   fluorophore channels (i.e. still furthest from zero) are selected for the
#'   second-round modulation. A value of `0.99` means the top 1% of cells.
#' @param remove.contaminants Logical, default `TRUE`. A QC check is performed
#'   to exclude any autofluorescence spectrum that is nearly identical to a
#'   fluorophore signature in `spectra`. This guards against low-level
#'   contamination of the unstained sample by single-stained controls.
#' @param parallel Logical, default `TRUE`, which enables parallel processing
#'   for per-cell AF identification. Used when `refine = TRUE`.
#' @param threads Numeric, defaults to a single thread for sequential processing
#'   (`parallel = FALSE`) or all available cores if `parallel = TRUE`. Used when
#'   `refine = TRUE`.
#' @param heatmap.color.palette Optional character string defining the viridis
#'   color palette for the fluorophore heatmap. Default is `"viridis"`. Options:
#'   `"magma"`, `"inferno"`, `"plasma"`, `"viridis"`, `"cividis"`, `"rocket"`,
#'   `"mako"`, `"turbo"`.
#' @param spectral.trace.color.palette Optional character string defining the
#'   color palette for the AF traces. Default is `NULL` (default R Brewer
#'   colors). Options: same as `heatmap.color.palette`.
#' @param af.fill.color Color for the shaded region indicating the range of
#'   autofluorescence variation in the variant plot. Default is `"red"`.
#' @param af.line.color Color for the median autofluorescence line in the
#'   variant plot. Default is `"black"`.
#'
#' @return A matrix of autofluorescence spectra (spectra in rows, detectors in
#'   columns). Row 1 is the population mean of the base spectra; subsequent rows
#'   are the deduplicated base spectra and, if `refine = TRUE`, modulated spectra
#'   for problem cells.
#'
#' @export
#'
#' @references
#' Van Gassen S et al. (2015). "FlowSOM: Using self-organizing maps for
#' visualization and interpretation of cytometry data." \emph{Cytometry Part A},
#' 87(7), 636-645. \doi{10.1002/cyto.a.22625}
#' Wehrens R, Kruisselbrink J (2018). "Flexible Self-Organizing Maps in kohonen
#' 3.0." \emph{Journal of Statistical Software}, \emph{87}(7), 1-18.
#' \doi{10.18637/jss.v087.i07}

get.af.spectra.kmeans <- function(
    unstained.sample,
    asp,
    spectra,
    som.dim              = 10,
    figures              = TRUE,
    plot.dir             = NULL,
    table.dir            = NULL,
    title                = "Autofluorescence spectra",
    verbose              = TRUE,
    deduplicate          = TRUE,
    duplication.threshold = 0.99,
    refine               = TRUE,
    problem.quantile     = 0.99,
    remove.contaminants  = TRUE,
    parallel             = TRUE,
    threads              = if ( parallel ) 0 else 1,
    heatmap.color.palette         = "viridis",
    spectral.trace.color.palette  = NULL,
    af.fill.color        = "red",
    af.line.color        = "black"
) {

  # ---------------------------------------------------------------------------
  # Setup
  # ---------------------------------------------------------------------------

  if ( is.null( plot.dir ) )  plot.dir  <- asp$figure.af.dir
  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir )
  if ( is.null( table.dir ) ) table.dir <- asp$table.spectra.dir
  if ( !dir.exists( table.dir ) ) dir.create( table.dir )
  if ( is.null( title ) ) title <- asp$af.file.name

  if ( is.null( threads ) ) threads <- asp$worker.process.n
  if ( parallel & threads == 0 ) threads <- parallelly::availableCores()

  # remove any pre-existing AF row from the fluorophore spectra
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  spectral.channels <- colnames( spectra )

  # ---------------------------------------------------------------------------
  # Import and prepare unstained sample
  # ---------------------------------------------------------------------------

  unstained.ff   <- readFCS( unstained.sample, return.keywords = TRUE )
  file.name      <- unstained.ff$keywords[[ "$FIL" ]]
  unstained.exprs <- unstained.ff$data[ , spectral.channels ]

  # OLS unmix without AF - combined with raw data for richer clustering features
  unmixed.no.af <- unmix.ols.fast( unstained.exprs, spectra )
  cluster.data  <- cbind( unstained.exprs, unmixed.no.af )

  # ---------------------------------------------------------------------------
  # Stage 1 - Base AF spectra via SOM
  # ---------------------------------------------------------------------------

  if ( verbose ) message( "Creating a self-organizing map of the autofluorescence" )

  cell.n <- nrow( cluster.data )

  if ( cell.n < 100 ) {
    stop( paste( "Inadequate cell numbers provided:", cell.n ) )
  } else if ( cell.n < 500 ) {
    som.dim <- max( 2, floor( sqrt( cell.n / 3 ) ) )
  }

  set.seed( asp$gate.downsample.seed )
  km <- stats::kmeans(
    cluster.data,
    centers  = som.dim * som.dim,
    nstart   = 10,
    iter.max = 100
  )

  # L-infinity normalise
  af.spectra <- t( apply( km$centers[ , spectral.channels ], 1, function( x ) x / max( x ) ) )
  af.spectra <- as.matrix( stats::na.omit( af.spectra ) )

  # Prepend population mean
  mean.af    <- colMeans( af.spectra )
  af.spectra <- rbind( mean.af, af.spectra )
  rownames( af.spectra ) <- paste0( "AF", seq_len( nrow( af.spectra ) ) )

  # Contamination QC: remove spectra resembling fluorophores
  af.spectra <- qc.af.spectra( af.spectra, spectra, plot.dir, remove.contaminants, pass = 2 )

  # Deduplication of base spectra
  if ( deduplicate ) {
    if ( verbose ) message( "Deduplicating base AF spectra by cosine similarity" )
    n.before   <- nrow( af.spectra )
    af.spectra <- deduplicate.spectra( af.spectra, threshold = duplication.threshold )
    n.after    <- nrow( af.spectra )
    if ( verbose )
      message(
        sprintf(
          "  %d base spectra retained after deduplication (dropped %d)",
          n.after, n.before - n.after
        )
      )
  }

  # Refresh row names after any removals
  rownames( af.spectra ) <- paste0( "AF", seq_len( nrow( af.spectra ) ) )

  # ---------------------------------------------------------------------------
  # Stage 1 figures - base spectra
  # ---------------------------------------------------------------------------

  if ( figures ) {
    if ( verbose ) message( "Plotting autofluorescence spectra" )

    # flip sign of negatively vectored AF spectra for plotting only (occurs on S8, A8)
    af.spectra.plot <- t( apply( af.spectra, 1, function( x ) {
      max.x <- ifelse( max( abs( x ) ) > max( x ), min( x ), max( x ) )
      x / max.x
    } ) )

    tryCatch(
      expr = {
        spectral.trace(
          spectral.matrix      = af.spectra.plot,
          asp                  = asp,
          title                = paste( title, "Autofluorescence spectra" ),
          plot.dir             = plot.dir,
          split.lasers         = FALSE,
          color.palette        = spectral.trace.color.palette
        )
        spectral.heatmap(
          spectra              = af.spectra.plot,
          title                = title,
          plot.dir             = plot.dir,
          color.palette        = heatmap.color.palette
        )
      },
      error = function( e ) {
        message( "Error in plotting AF spectra: ", e$message )
        return( NULL )
      }
    )
  }

  # ---------------------------------------------------------------------------
  # Stage 2 - Refine: targeted modulation for problem cells
  # ---------------------------------------------------------------------------

  if ( refine ) {

    if ( verbose ) message( "Refine: identifying best-fitting AF - first pass" )

    # Per-cell AF assignment on the unstained sample using the base spectra
    if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
         "assign.af.fluor.fast" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
      af.assignments <- AutoSpectralRcpp::assign.af.fluor.fast(
        raw.data  = unstained.exprs,
        spectra   = spectra,
        af.spectra = af.spectra,
        threads   = asp$worker.process.n
      )
    } else {
      af.assignments <- assign.af.fluorophores(
        raw.data   = unstained.exprs,
        spectra    = spectra,
        af.spectra = af.spectra
      )
    }

    # Unmix each cell with its assigned AF spectrum, tracking residuals and
    # projected fluorophore signal so we can compute a detector-space error
    fluor.idx        <- 2:( nrow( spectra ) + 1 )
    af.abundance     <- rep( 0, nrow( unstained.exprs ) )
    unmixed          <- cbind( af.abundance, unmixed.no.af )
    residuals        <- matrix( 0, nrow = nrow( unstained.exprs ), ncol = ncol( spectra ) )
    proj.fluor       <- matrix( 0, nrow = nrow( unstained.exprs ), ncol = ncol( spectra ) )

    combined.spectra <- matrix( NA_real_, nrow = nrow( spectra ) + 1, ncol = ncol( spectra ) )
    combined.spectra[ fluor.idx, ] <- spectra

    for ( af in seq_len( nrow( af.spectra ) ) ) {
      combined.spectra[ 1, ] <- af.spectra[ af, ]
      cell.idx <- which( af.assignments == af )

      if ( length( cell.idx ) > 0 ) {
        unmixed[ cell.idx, ] <- unmix.ols.fast(
          unstained.exprs[ cell.idx, , drop = FALSE ],
          combined.spectra
        )
        residuals[ cell.idx, ] <-
          unstained.exprs[ cell.idx, , drop = FALSE ] -
          ( unmixed[ cell.idx, , drop = FALSE ] %*% combined.spectra )
        proj.fluor[ cell.idx, ] <-
          unmixed[ cell.idx, fluor.idx, drop = FALSE ] %*%
          combined.spectra[ fluor.idx, , drop = FALSE ]
      }
    }

    # detector-space error = fluorophore projection + raw residuals
    error <- residuals + proj.fluor

    # ---- Identify problem cells (those still furthest from zero) -------------

    if ( verbose ) message( "Refine: calculating error magnitude for problem cell selection" )

    if ( length( fluor.idx ) > 1 ) {
      error.magnitude <- sqrt( rowSums( unmixed[ , fluor.idx ]^2 ) )
    } else {
      error.magnitude <- abs( unmixed[ , fluor.idx ] )
    }

    # Step the quantile down in 5 % increments until we have >= 500 cells
    while ( TRUE ) {
      threshold      <- stats::quantile( error.magnitude, problem.quantile )
      problem.idx    <- which( error.magnitude > threshold )
      problem.cell.n <- length( problem.idx )

      if ( problem.cell.n >= 500 ) break

      problem.quantile <- problem.quantile - 0.05

      if ( problem.quantile < 0.5 ) {
        threshold      <- stats::quantile( error.magnitude, problem.quantile )
        problem.idx    <- which( error.magnitude > threshold )
        problem.cell.n <- length( problem.idx )
        break
      }
    }

    if ( verbose )
      message(
        sprintf(
          "Refine: %d problem cells selected (quantile = %.2f, threshold = %.2f)",
          problem.cell.n, problem.quantile, threshold
        )
      )

    # ---- Modulate base spectra using error clusters -------------------------

    if ( problem.cell.n > 10 ) {

      # AF abundance for the problem cells (normalisation denominator)
      af.abundance.problem <- unmixed[ problem.idx, 1 ]
      af.abundance.problem[ af.abundance.problem == 0 ] <- 1e-6

      # Spill ratios: per-channel error normalised by AF abundance,
      # giving a dimensionless signature of how the current AF estimate is wrong
      spill.ratios <- error[ problem.idx, ] / af.abundance.problem

      if ( verbose )
        message(
          paste( "Refine: clustering", problem.cell.n, "problem cells by spillover error pattern" )
        )

      som.dim.error <- max( 2, floor( sqrt( problem.cell.n / 3 ) ) )

      set.seed( asp$gate.downsample.seed )
      map.error <- FlowSOM::SOM(
        spill.ratios,
        xdim   = som.dim.error,
        ydim   = som.dim.error,
        silent = TRUE
      )

      cluster.ids <- unique( map.error$mapping[ , 1 ] )

      modulated.list <- lapply( cluster.ids, function( cl ) {
        cl.sub.idx <- which( map.error$mapping[ , 1 ] == cl )
        global.idx <- problem.idx[ cl.sub.idx ]

        # median correction pattern for this error cluster
        median.ratio <- apply(
          spill.ratios[ cl.sub.idx, , drop = FALSE ],
          2,
          stats::median
        )

        # which base AF spectra were assigned to cells in this cluster?
        contributing.af.ids <- unique( af.assignments[ global.idx ] )

        # modulate each contributing base spectrum
        new.specs <- lapply( contributing.af.ids, function( id ) {
          base.spec <- af.spectra[ id, ]
          updated   <- base.spec * ( 1 + median.ratio )
          peak      <- max( updated )
          if ( peak > 1e-12 ) updated <- updated / peak
          return( updated )
        } )

        return( do.call( rbind, new.specs ) )
      } )

      modulated.af.spectra <- do.call( rbind, modulated.list )
      modulated.af.spectra <- as.matrix( stats::na.omit( modulated.af.spectra ) )

      if ( nrow( modulated.af.spectra ) > 0 && deduplicate ) {

        # Step 1: deduplicate modulated spectra against each other
        n.mod.before         <- nrow( modulated.af.spectra )
        modulated.af.spectra <- deduplicate.spectra(
          modulated.af.spectra,
          threshold = duplication.threshold
        )

        # Step 2: drop any modulated spectrum too similar to an already-kept
        # base spectrum (cross-deduplication)
        cross.sim  <- cosine.similarity.cross( modulated.af.spectra, af.spectra )
        # cross.sim is (n_modulated x n_existing); keep rows where max sim < threshold
        novel.mask <- apply( cross.sim, 1, max ) < duplication.threshold
        modulated.af.spectra <- modulated.af.spectra[ novel.mask, , drop = FALSE ]

        n.novel <- nrow( modulated.af.spectra )
        if ( verbose )
          message(
            sprintf(
              "Refine: %d novel modulated spectra retained after deduplication (dropped %d)",
              n.novel, n.mod.before - n.novel
            )
          )
      }

      if ( nrow( modulated.af.spectra ) > 0 ) {

        af.spectra <- rbind( af.spectra, modulated.af.spectra )
        af.spectra <- as.matrix( stats::na.omit( af.spectra ) )

        # Contamination QC on the expanded set
        af.spectra <- qc.af.spectra( af.spectra, spectra, plot.dir, remove.contaminants )

        rownames( af.spectra ) <- paste0( "AF", seq_len( nrow( af.spectra ) ) )

        if ( verbose )
          message(
            sprintf(
              "Refine: %d total AF spectra after modulation and QC",
              nrow( af.spectra )
            )
          )

      } else {
        if ( verbose )
          message( "Refine: all modulated spectra were redundant with base spectra - nothing appended." )
      }

      # ---- Refine figures ---------------------------------------------------

      if ( figures ) {
        if ( verbose ) message( "Refine: identifying best-fitting AF - second pass" )

        if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
             "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
          unmixed.second <- AutoSpectralRcpp::unmix.autospectral.rcpp(
            raw.data   = unstained.exprs,
            spectra    = spectra,
            af.spectra = af.spectra,
            verbose    = FALSE,
            parallel   = TRUE,
            threads    = threads
          )
        } else {
          af.assignments.second <- assign.af.fluorophores(
            raw.data   = unstained.exprs,
            spectra    = spectra,
            af.spectra = af.spectra
          )

          unmixed.second <- unmixed

          for ( af in seq_len( nrow( af.spectra ) ) ) {
            combined.spectra[ 1, ] <- af.spectra[ af, ]
            cell.idx <- which( af.assignments.second == af )
            if ( length( cell.idx ) > 0 ) {
              unmixed.second[ cell.idx, ] <- unmix.ols.fast(
                unstained.exprs[ cell.idx, , drop = FALSE ],
                combined.spectra
              )
            }
          }
        }

        if ( ncol( unmixed.no.af ) > 1 ) {
          if ( verbose ) message( "Refine: plotting impact of AF extraction" )

          channel.sd    <- apply( unmixed.no.af, 2, stats::sd )
          worst.channels <- colnames( unmixed.no.af )[ order( channel.sd, decreasing = TRUE )[ 1:2 ] ]

          tryCatch(
            expr = {
              create.biplot(
                unmixed.no.af,
                x.dim      = worst.channels[ 1 ],
                y.dim      = worst.channels[ 2 ],
                asp        = asp,
                title      = paste( file.name, "_", title, "_No_AF_Extraction" ),
                output.dir = plot.dir
              )
              create.biplot(
                unmixed,
                x.dim      = worst.channels[ 1 ],
                y.dim      = worst.channels[ 2 ],
                asp        = asp,
                title      = paste0( file.name, "_", title, "_PerCell_AF_Extraction_First_Pass" ),
                output.dir = plot.dir
              )
              create.biplot(
                unmixed.second,
                x.dim      = worst.channels[ 1 ],
                y.dim      = worst.channels[ 2 ],
                asp        = asp,
                title      = paste0( file.name, "_", title, "_PerCell_AF_Extraction_Second_Pass" ),
                output.dir = plot.dir
              )
            },
            error = function( e ) {
              message( "Error in plotting AF extraction: ", e$message )
              return( NULL )
            }
          )
        }
      }

    } else {
      message( "Refine: insufficient problem cells found - skipping modulation." )
    }

  }   # end refine

  # ---------------------------------------------------------------------------
  # Save and final figures
  # ---------------------------------------------------------------------------

  if ( is.null( title ) )
    af.file.name <- paste0( file.name, "_", asp$af.file.name, ".csv" )
  else
    af.file.name <- paste0( file.name, "_", title, ".csv" )

  utils::write.csv( af.spectra, file = file.path( table.dir, af.file.name ) )

  if ( figures ) {
    if ( verbose ) message( "Plotting autofluorescence variation" )

    tryCatch(
      expr = {
        spectral.variant.plot(
          af.spectra,
          mean.af,
          title               = paste( title, "Autofluorescence variation" ),
          save                = TRUE,
          plot.dir            = plot.dir,
          variant.fill.color  = af.fill.color,
          median.line.color   = af.line.color
        )
      },
      error = function( e ) e
    )
  }

  return( af.spectra )

}


# ---------------------------------------------------------------------------
# Helper: greedy cosine-similarity deduplication
# ---------------------------------------------------------------------------

#' @title Deduplicate Spectra by Cosine Similarity
#'
#' @description
#' Iterates through rows of a spectral matrix in order, retaining a row only
#' if its cosine similarity to every already-retained row is strictly below
#' `threshold`. Intended for internal use by `get.af.spectra`.
#'
#' @param spectra Numeric matrix, spectra in rows and detectors in columns.
#'   Assumed to be L-infinity normalised.
#' @param threshold Numeric scalar in (0, 1]. Rows at or above this similarity
#'   to any retained row are dropped. Default `0.99`.
#'
#' @return Numeric matrix with redundant rows removed.
#'
#' @keywords internal

deduplicate.spectra <- function( spectra, threshold = 0.99 ) {

  if ( nrow( spectra ) <= 1 ) return( spectra )

  # Row-normalise to unit length for cosine similarity via dot product
  norms  <- sqrt( rowSums( spectra^2 ) )
  norms  <- ifelse( norms < 1e-12, 1, norms )
  s.norm <- spectra / norms

  kept <- integer( 0 )

  for ( i in seq_len( nrow( spectra ) ) ) {
    if ( length( kept ) == 0 ) {
      kept <- c( kept, i )
    } else {
      # cosine similarities of row i against all kept rows
      sims <- s.norm[ kept, , drop = FALSE ] %*% s.norm[ i, ]
      if ( all( sims < threshold ) )
        kept <- c( kept, i )
    }
  }

  return( spectra[ kept, , drop = FALSE ] )

}


# ---------------------------------------------------------------------------
# Helper: cross cosine similarity matrix
# ---------------------------------------------------------------------------

#' @title Cross Cosine Similarity
#'
#' @description
#' Computes pairwise cosine similarity between rows of two matrices `a` and
#' `b`. Returns an (nrow(a) x nrow(b)) matrix. Intended for internal use by
#' `get.af.spectra`.
#'
#' @param a Numeric matrix.
#' @param b Numeric matrix with the same number of columns as `a`.
#'
#' @return Numeric matrix of cosine similarities, shape (nrow(a), nrow(b)).
#'
#' @keywords internal

cosine.similarity.cross <- function( a, b ) {

  norm.row <- function( m ) {
    norms <- sqrt( rowSums( m^2 ) )
    norms <- ifelse( norms < 1e-12, 1, norms )
    m / norms
  }

  norm.row( a ) %*% t( norm.row( b ) )

}
