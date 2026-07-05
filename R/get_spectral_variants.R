# get_spectral_variants.r

#' @title Get Spectral Variations for Fluorophores
#'
#' @description
#' Cycles through all fluorophores defined in \code{control.def.file},
#' identifying variation in their spectral profiles via SOM clustering on
#' scatter-matched, per-event background-corrected data.
#'
#' For each fluorophore the autofluorescence reference is derived \strong{in
#' situ} from the paired universal-negative file (or internally from the lower
#' 25\% of events when no universal negative is supplied). The AF mean vector is
#' used to project out autofluorescence and identify the empirical peak
#' detector. All positive events are scatter-matched to unstained events and
#' their per-event background is subtracted before SOM clustering. This gives a
#' comprehensive, population-level picture of true fluorophore spectral
#' variability without requiring a pre-computed \code{af.spectra} matrix.
#'
#' The output is saved as an .rds file and per-fluorophore variant plots are
#' produced if requested.
#'
#' @importFrom lifecycle deprecate_warn
#'
#' @param control.dir Character. Path to the single-stained control FCS files.
#' @param control.def.file Character. Path to the control definition CSV.
#'   Must pass \code{check.control.file()}.
#' @param asp The AutoSpectral parameter list from \code{get.autospectral.param()}.
#' @param spectra Numeric matrix. Reference spectra; fluorophores in rows,
#'   detectors in columns.
#' @param figures Logical, default \code{TRUE}. Whether to save variant-spectrum
#'   plots.
#' @param output.dir Character or \code{NULL}. Directory for figures and the
#'   .rds output file. Defaults to \code{asp$variant.dir}.
#' @param parallel Logical, default \code{FALSE}. Enable parallel processing
#'   across fluorophores.
#' @param verbose Logical, default \code{TRUE}. Set to \code{FALSE} to suppress
#'   messages.
#' @param threads Numeric or \code{NULL}. Number of parallel workers. Defaults
#'   to \code{asp$worker.process.n}.
#' @param n.cells Integer, default \code{10000}. Maximum positive events per
#'   fluorophore used for SOM clustering. Files with more events above threshold
#'   are randomly downsampled. Passed to \code{get.fluor.variants}.
#' @param som.dim Integer, default \code{10}. Side length of the square SOM
#'   grid; up to \code{som.dim^2} candidate variants per fluorophore before
#'   cosine QC. Passed to \code{get.fluor.variants}.
#' @param k.neighbors Integer, default \code{3}. Number of scatter-space
#'   nearest neighbours from the unstained pool used to estimate per-event
#'   background. Passed to \code{get.fluor.variants}.
#' @param sim.threshold Numeric, default \code{0.99}. Minimum cosine similarity
#'   to the reference spectrum for a SOM centroid to be retained as a variant.
#'   Passed to \code{get.fluor.variants}.
#' @param variant.fill.color Color for the shaded ribbon in variant plots.
#'   Default \code{"red"}.
#' @param variant.fill.alpha Alpha for \code{variant.fill.color}. Default
#'   \code{0.7}.
#' @param median.line.color Color for the reference-spectrum line. Default
#'   \code{"black"}.
#' @param median.linewidth Width of the reference-spectrum line. Default
#'   \code{1}.
#' @param stained.sample Optional file path to a representative stained FCS
#'   file. When supplied, it is read and unmixed to obtain per-fluorophore
#'   median positive signal (MFI), which weights the optimization necessity
#'   scores by fluorophore brightness. Pass `NULL` (default) to use purely
#'   geometric scores.
#' @param optimize.necessity.threshold Numeric in `[0, 1]`, default `0.01`.
#'   Passed to `calculate.optimize.necessity()`. Fluorophores whose normalised
#'   leakage score falls below this value are flagged as not requiring per-cell
#'   spectral optimisation. The result is stored in
#'   `$optimize.recommended` in the returned list and used automatically by
#'   `unmix.autospectral.rcpp()` to skip unnecessary optimisation passes.
#' @param ... Ignored. Catches and warns on previously used deprecated
#'   arguments: \code{af.spectra}, \code{refine}, \code{problem.quantile},
#'   \code{pos.quantile}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{thresholds}}{Named numeric vector of positivity thresholds in
#'     the unmixed space, one per fluorophore.}
#'   \item{\code{variants}}{Named list of variant-spectra matrices, one per
#'     fluorophore. Each matrix has variants in rows and detectors in columns.}
#'   \item{\code{delta.list}}{Named list of delta matrices (variant minus
#'     reference spectrum), one per fluorophore.}
#'   \item{\code{delta.norms}}{Named list of Euclidean norms of the deltas,
#'     one numeric vector per fluorophore.}
#' }
#' The list is also saved as an .rds file in \code{output.dir}.
#'
#' @export

get.spectral.variants <- function(
    control.dir,
    control.def.file,
    asp,
    spectra,
    figures            = TRUE,
    output.dir         = NULL,
    parallel           = FALSE,
    verbose            = TRUE,
    threads            = NULL,
    n.cells            = 10000L,
    som.dim            = 10L,
    k.neighbors        = 3L,
    sim.threshold      = 0.985,
    variant.fill.color = "red",
    variant.fill.alpha = 0.7,
    median.line.color  = "black",
    median.linewidth   = 1,
    stained.sample               = NULL,
    optimize.necessity.threshold = 0.01,
    ...
) {

  # ---------------------------------------------------------------------------
  # Deprecated argument handling
  # ---------------------------------------------------------------------------

  dots <- list( ... )

  for ( old.arg in c( "pos.quantile" ) ) {
    if ( !is.null( dots[[ old.arg ]] ) )
      lifecycle::deprecate_warn( "0.9.0",
        paste0( "get.spectral.variants(", old.arg, ")" ),
        details = "no longer used" )
  }

  for ( old.arg in c( "refine", "problem.quantile" ) ) {
    if ( !is.null( dots[[ old.arg ]] ) )
      lifecycle::deprecate_warn( "1.6.0",
        paste0( "get.spectral.variants(", old.arg, ")" ),
        details = paste0(
          "The second-pass refinement is superseded by in-situ ",
          "scatter-matched background subtraction."
        ) )
  }

  if ( !is.null( dots$af.spectra ) )
    lifecycle::deprecate_warn( "1.6.0",
      "get.spectral.variants(af.spectra)",
      details = paste0(
        "AF is now derived in situ from the universal-negative files ",
        "listed in the control table. The af.spectra argument is ignored."
      ) )

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  # Catches the common failure mode of a positional argument shifting into
  # the wrong parameter (e.g. passing the now-deprecated `af.spectra` as the
  # 5th positional argument, which silently lands in `figures`).

  .type.err <- function( arg.name, expected, x ) {
    stop(
      paste0(
        "`", arg.name, "` must be ", expected, ", but got an object of class `",
        paste( class( x ), collapse = "/" ), "` with length ", length( x ), ".\n",
        "  If you're passing a spectra matrix or other object positionally, ",
        "check that it lines up with the current argument order for ",
        "get.spectral.variants() -- deprecated arguments like `af.spectra` ",
        "are ignored and must not be passed positionally."
      ),
      call. = FALSE
    )
  }

  if ( !is.character( control.dir ) || length( control.dir ) != 1 || is.na( control.dir ) )
    .type.err( "control.dir", "a single character path", control.dir )

  if ( !is.character( control.def.file ) || length( control.def.file ) != 1 || is.na( control.def.file ) )
    .type.err( "control.def.file", "a single character path", control.def.file )

  if ( !is.list( asp ) )
    .type.err( "asp", "a list (from get.autospectral.param())", asp )

  if ( !is.matrix( spectra ) && !is.data.frame( spectra ) )
    .type.err( "spectra", "a numeric matrix (fluorophores in rows, detectors in columns)", spectra )
  spectra <- as.matrix( spectra )
  storage.mode( spectra ) <- "double"
  if ( is.null( rownames( spectra ) ) )
    stop( "`spectra` must have rownames giving fluorophore names (including \"AF\").", call. = FALSE )

  if ( !is.logical( figures ) || length( figures ) != 1 || is.na( figures ) )
    .type.err( "figures", "a single TRUE/FALSE value", figures )

  if ( !is.null( output.dir ) && ( !is.character( output.dir ) || length( output.dir ) != 1 ) )
    .type.err( "output.dir", "NULL or a single character path", output.dir )

  if ( !is.logical( parallel ) || length( parallel ) != 1 || is.na( parallel ) )
    .type.err( "parallel", "a single TRUE/FALSE value", parallel )

  if ( !is.logical( verbose ) || length( verbose ) != 1 || is.na( verbose ) )
    .type.err( "verbose", "a single TRUE/FALSE value", verbose )

  if ( !is.null( threads ) && ( !is.numeric( threads ) || length( threads ) != 1 ) )
    .type.err( "threads", "NULL or a single number", threads )

  if ( !is.numeric( n.cells ) || length( n.cells ) != 1 || n.cells <= 0 )
    .type.err( "n.cells", "a single positive number", n.cells )

  if ( !is.numeric( som.dim ) || length( som.dim ) != 1 || som.dim <= 0 )
    .type.err( "som.dim", "a single positive number", som.dim )

  if ( !is.numeric( k.neighbors ) || length( k.neighbors ) != 1 || k.neighbors <= 0 )
    .type.err( "k.neighbors", "a single positive number", k.neighbors )

  if ( !is.numeric( sim.threshold ) || length( sim.threshold ) != 1 ||
       sim.threshold < 0 || sim.threshold > 1 )
    .type.err( "sim.threshold", "a single number in [0, 1]", sim.threshold )

  if ( !is.null( stained.sample ) && ( !is.character( stained.sample ) || length( stained.sample ) != 1 ) )
    .type.err( "stained.sample", "NULL or a single character path", stained.sample )

  if ( !is.numeric( optimize.necessity.threshold ) || length( optimize.necessity.threshold ) != 1 ||
       optimize.necessity.threshold < 0 || optimize.necessity.threshold > 1 )
    .type.err( "optimize.necessity.threshold", "a single number in [0, 1]", optimize.necessity.threshold )

  # ---------------------------------------------------------------------------
  # Setup
  # ---------------------------------------------------------------------------

  if ( is.null( output.dir ) ) output.dir <- asp$variant.dir
  if ( !dir.exists( output.dir ) ) dir.create( output.dir )

  if ( som.dim > 20 ) {
    n.cells <- min( 5000, n.cells )
    warning(
      paste(
        "Argument `som.dim` has been set to", som.dim, "which will produce",
        som.dim^2, "spectral variants per fluorophore.", "\n",
        "This requires proprotionally more cells in `n.cells` as input,",
        "and may trigger failure.",
        "`n.cells` has been automatically adjusted to a minimum of 5000."
      ),
      call. = FALSE
    )
  }

  fluorophores     <- rownames( spectra )[ rownames( spectra ) != "AF" ]
  spectra <- spectra[ fluorophores, , drop = FALSE ]
  spectral.channel <- colnames( spectra )

  # ---------------------------------------------------------------------------
  # Read and validate control file
  # ---------------------------------------------------------------------------

  if ( !file.exists( control.def.file ) )
    stop( paste( "Unable to locate control.def.file:", control.def.file ),
          call. = FALSE )

  if ( verbose ) message( "\033[32mChecking control file for errors \033[0m" )
  check.control.file( control.dir, control.def.file, asp, strict = TRUE )

  control.table <- utils::read.csv(
    control.def.file, stringsAsFactors = FALSE, strip.white = TRUE
  )
  control.table[] <- lapply( control.table, function( x ) {
    if ( is.character( x ) ) { x <- trimws( x ); x[ x == "" ] <- NA; x } else x
  } )

  # scatter channels (needed for KNN matching in get.fluor.variants)
  scatter.channel          <- read.scatter.parameter( asp )
  spectral.channel         <- colnames( spectra )

  if ( grepl( "Discover", asp$cytometer ) )
    spectral.channel <- spectral.channel[ grep( asp$spectral.channel, spectral.channel ) ]

  # per-sample metadata
  table.fluors       <- control.table$fluorophore
  table.fluors       <- table.fluors[ !is.na( table.fluors ) ]
  universal.negative <- control.table$universal.negative
  universal.negative[ is.na( universal.negative ) ] <- "FALSE"
  names( universal.negative ) <- table.fluors
  flow.channel       <- control.table$channel
  names( flow.channel ) <- table.fluors
  flow.file.name     <- control.table$filename
  names( flow.file.name ) <- table.fluors
  control.type <- control.table$control.type
  names( control.type ) <- table.fluors

  if ( !( "AF" %in% table.fluors ) )
    stop(
      "Unable to locate `AF` control in control file. An unstained cell control is required.",
      call. = FALSE
    )

  # ensure spectra columns match channel order
  spectra.cols <- colnames( spectra )
  if ( !identical( spectral.channel, spectra.cols ) ) {
    if ( all( spectra.cols %in% spectral.channel ) &&
         length( spectra.cols ) == length( spectral.channel ) ) {
      spectra <- spectra[ , spectral.channel ]
      message( "Columns of spectra reordered to match data" )
    } else {
      stop( "Column names in spectra and data do not match.", call. = FALSE )
    }
  }

  # reconcile fluorophores
  if ( !all( table.fluors %in% fluorophores ) ) {
    fluor.to.match  <- table.fluors[ !grepl( "Negative|^AF$", table.fluors ) ]
    matching.fluors <- fluor.to.match %in% fluorophores
    if ( !any( matching.fluors ) )
      stop( "No matching fluorophores between `spectra` and the control file.",
            call. = FALSE )
    if ( !all( matching.fluors ) )
      warning( "Some fluorophores in the control file are absent from `spectra`.",
               call. = FALSE )
    table.fluors <- fluor.to.match[ matching.fluors ]
  }

  # ---------------------------------------------------------------------------
  # Positivity thresholds from the unstained file
  # ---------------------------------------------------------------------------

  if ( verbose )
    message( paste0( "\033[32m", "Measuring background in unstained samples", "\033[0m" ) )

  unstained <- readFCS( file.path( control.dir, flow.file.name[ "AF" ] ) )

  if ( nrow( unstained ) > asp$gate.downsample.n.cells ) {
    set.seed( asp$bird.seed )
    unstained.idx <- sample( nrow( unstained ), asp$gate.downsample.n.cells )
    unstained     <- unstained[ unstained.idx, spectral.channel, drop = FALSE ]
  } else {
    unstained <- unstained[ , spectral.channel, drop = FALSE ]
  }

  raw.thresholds <- apply( unstained, 2, function( col )
    stats::quantile( col, 0.995 ) )

  # get AF spectra in place
  af.spectra <- get.af.spectra(
    file.path( control.dir, flow.file.name[ "AF" ] ),
    asp,
    spectra,
    som.dim = 10,
    figures = FALSE,
    save = FALSE,
    refine = FALSE
  )

  # derive per-file AF PCs for all unique unstained cell files used as negatives
  cell.fluors <- names( control.type )[ control.type == "cells" ]
  univ.neg.files <- unique( universal.negative[
    names( universal.negative ) %in% cell.fluors &
      universal.negative != "FALSE" &
      !is.na( universal.negative ) &
      grepl( "\\.fcs$", universal.negative, ignore.case = TRUE )
  ] )

  af.pcs.list <- lapply( univ.neg.files, function( fn ) {
    dat <- readFCS( file.path( control.dir, fn ) )
    if ( nrow( dat ) > asp$gate.downsample.n.cells ) {
      set.seed( asp$bird.seed )
      dat <- dat[ sample( nrow( dat ), asp$gate.downsample.n.cells ),
                  spectral.channel, drop = FALSE ]
    } else {
      dat <- dat[ , spectral.channel, drop = FALSE ]
    }
    sv <- svd( dat, nu = 0, nv = 4 )
    t( sv$v )
  } )
  names( af.pcs.list ) <- univ.neg.files

  # find the likely positivity thresholds for determining what needs refinement
  unstained.unmixed <- unmix.autospectral(
    unstained,
    spectra,
    af.spectra,
    verbose = FALSE
  )
  unmixed.thresholds <- apply(
    unstained.unmixed[ , fluorophores, drop = FALSE ], 2, function( col )
      stats::quantile( col, 0.995 )
  )

  if ( is.null( names( table.fluors ) ) ) names( table.fluors ) <- table.fluors

  # ---------------------------------------------------------------------------
  # Parallel setup
  # ---------------------------------------------------------------------------

  if ( is.null( threads ) ) threads <- asp$worker.process.n

  args.list <- list(
    file.name          = flow.file.name,
    control.dir        = control.dir,
    asp                = asp,
    spectra            = spectra,
    figures            = figures,
    output.dir         = output.dir,
    verbose            = verbose,
    spectral.channel   = spectral.channel,
    scatter.channel    = scatter.channel,
    universal.negative = universal.negative,
    control.type       = control.type,
    raw.thresholds     = raw.thresholds,
    unmixed.thresholds = unmixed.thresholds,
    flow.channel       = flow.channel,
    af.pcs             = af.pcs.list,
    n.cells            = n.cells,
    som.dim            = som.dim,
    k.neighbors        = k.neighbors,
    sim.threshold      = sim.threshold,
    variant.fill.color = variant.fill.color,
    variant.fill.alpha = variant.fill.alpha,
    median.line.color  = median.line.color,
    median.linewidth   = median.linewidth
  )

  if ( parallel ) {
    internal.functions <- c(
      "get.fluor.variants",
      "cosine.similarity",
      "spectral.variant.plot.dens",
      "unmix.ols",
      ".cosine.sim.rows"
    )
    exports <- c( "args.list", "table.fluors", internal.functions )

    result <- create.parallel.lapply(
      asp,
      exports,
      parallel   = parallel,
      threads    = threads,
      export.env = environment()
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # ---------------------------------------------------------------------------
  # Main loop
  # ---------------------------------------------------------------------------

  if ( verbose )
    message( paste0( "\033[34m", "Identifying spectral variation", "\033[0m" ) )

  # initialise with base spectra as safe fallback
  spectral.variants <- lapply( table.fluors, function( fl )
    spectra[ fl, , drop = FALSE ] )
  names( spectral.variants ) <- table.fluors

  updated.variants <- tryCatch(
    expr = {
      lapply.function( table.fluors, function( f ) {
        tryCatch(
          expr = {
            if ( is.na( args.list$flow.channel[ f ] ) )
              stop( paste( "No flow channel mapped for", f ) )
            do.call( get.fluor.variants, c( list( f ), args.list ) )
          },
          error = function( e ) {
            list( is.error = TRUE, msg = conditionMessage( e ) )
          }
        )
      } )
    },
    finally = {
      if ( !is.null( result$cleanup ) ) result$cleanup()
    }
  )

  names( updated.variants ) <- table.fluors

  for ( f in table.fluors ) {
    res <- updated.variants[[ f ]]
    if ( is.list( res ) && isTRUE( res$is.error ) ) {
      warning( paste( "Variant calculation failed for:", f, "| Error:", res$msg ) )
    } else if ( !is.null( res ) ) {
      spectral.variants[[ f ]] <- res
    }
  }

  # ---------------------------------------------------------------------------
  # Deltas
  # ---------------------------------------------------------------------------
  # calculate deltas for each fluorophore's variants
  delta.list <- lapply( names( spectral.variants ), function( fl ) {
    spectral.variants[[ fl ]] - matrix(
      spectra[ fl, ],
      nrow = nrow( spectral.variants[[ fl ]] ),
      ncol = ncol( spectra ),
      byrow = TRUE
    )
  } )
  names( delta.list ) <- names( spectral.variants )

  delta.norms <- lapply( delta.list, function( d ) sqrt( rowSums( d^2 ) ) )
  names( delta.norms ) <- names( spectral.variants )

  ### calculate optimization necessity scores ###

  # spectra matrix without AF row for scoring
  spectra.no.af <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  # optionally derive MFI weights from a representative stained sample
  mu.weights <- NULL

  if ( !is.null( stained.sample ) ) {

    if ( !file.exists( stained.sample ) ) {
      warning(
        paste( "stained.sample file not found:", stained.sample,
               "- proceeding with geometric scores only." ),
        call. = FALSE
      )
    } else {

      if ( verbose )
        message( paste0(
          "\033[34m",
          "Computing per-fluorophore MFI weights from stained sample",
          "\033[0m"
        ) )

      stained.raw <- readFCS( stained.sample )

      if ( nrow( stained.raw ) > 5000 ) {
        set.seed( asp$bird.seed )
        stained.raw <- stained.raw[
          sample( nrow( stained.raw ), 5000 ),
          spectral.channel,
          drop = FALSE
        ]
      } else {
        stained.raw <- stained.raw[ , spectral.channel, drop = FALSE ]
      }
      # TBD: add autospectralRcpp option here
      stained.unmixed <- unmix.autospectral(
        stained.raw,
        spectra,
        af.spectra,
        verbose = FALSE
      )

      fluor.cols <- rownames( spectra.no.af )

      mu.weights <- apply(
        stained.unmixed[ , fluor.cols, drop = FALSE ],
        2,
        function( x ) {
          pos <- x[ x > 0 ]
          if ( length( pos ) == 0 ) 0 else stats::median( pos )
        }
      )
    }
  }

  necessity <- calculate.optimize.necessity(
    spectra    = spectra.no.af,
    delta.list = delta.list,
    mu         = mu.weights,
    threshold  = optimize.necessity.threshold,
    verbose    = verbose
  )

  if ( verbose )
    message( paste0( "\033[34m", "Spectral variation computed!", "\033[0m" ) )

  variants <- list(
    thresholds  = unmixed.thresholds,
    variants    = spectral.variants,
    delta.list  = delta.list,
    delta.norms = delta.norms,
    optimize.scores      = necessity$scores.norm,
    optimize.recommended = necessity$optimize.recommended
  )

  saveRDS( variants, file = file.path( output.dir, asp$variant.filename ) )

  return( variants )
}
