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
#' A cell-based unstained sample is still required at this stage, since it
#' anchors the positivity thresholds and per-node AF library. It is normally
#' read from the control file's \code{"AF"} row; if the control file has no
#' \code{"AF"} row (e.g. a bead-only or negative-free control setup), supply
#' one directly via \code{unstained.sample} instead.
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
#' @param parallel Logical, default \code{TRUE}. Enable parallel processing
#'   for SOM clustering (requires AutoSpectralRcpp).
#' @param verbose Logical, default \code{TRUE}. Set to \code{FALSE} to suppress
#'   messages.
#' @param threads Numeric or \code{NULL}. Number of parallel workers. Defaults
#'   to \code{asp$worker.process.n}.
#' @param n.cells Integer, default \code{10000}. Maximum positive events per
#'   fluorophore used for SOM clustering. Files with more events above threshold
#'   are randomly downsampled. Passed to \code{get.fluor.variants}.
#' @param som.dim Integer, default \code{5}. Side length of the square SOM
#'   grid; up to \code{som.dim^2} candidate variants per fluorophore before
#'   cosine QC. Passed to \code{get.fluor.variants}.
#' @param k.neighbors Integer, default \code{3}. Number of scatter-space
#'   nearest neighbours from the unstained pool used to estimate per-event
#'   background. Passed to \code{get.fluor.variants}.
#' @param sim.threshold Numeric, default \code{0.99}. Minimum cosine similarity
#'   to the reference spectrum for a SOM centroid to be retained as a variant.
#'   Passed to \code{get.fluor.variants}.
#' @param sim.threshold.floor Numeric, default \code{0.90}. Lower bound for
#'   adaptive relaxation of \code{sim.threshold} when the initial cutoff
#'   retains fewer than 20 events. Relaxation is logged via \code{warning()}
#'   and the threshold actually used is returned as the
#'   \code{"cosine.threshold.used"} attribute.
#' @param af.collinear.threshold Numeric, default \code{0.95}. Minimum
#'   cosine similarity between \code{fluor}'s reference spectrum and any of
#'   its paired unstained file's AF principal directions (\code{af.pcs}) at
#'   or above which the AF-component projection step is skipped, since a
#'   joint OLS fit against near-collinear AF and fluorophore directions can
#'   push real fluorophore signal into the AF term. Recorded as the
#'   \code{"af.collinear"} attribute.
#' @param noise.floor.tail.fraction Numeric in (0, 1), default \code{0.20}.
#'   Fraction of each detector's raw values (lowest end) used to estimate the
#'   per-control noise floor. Passed to \code{get.fluor.variants}.
#' @param variant.fill.color Color for the shaded ribbon in variant plots.
#'   Default \code{"red"}.
#' @param variant.fill.alpha Alpha for \code{variant.fill.color}. Default
#'   \code{0.7}.
#' @param median.line.color Color for the reference-spectrum line. Default
#'   \code{"black"}.
#' @param median.linewidth Width of the reference-spectrum line. Default
#'   \code{1}.
#' @param use.unmixed Logical, default \code{TRUE}. Whether AF extraction
#'   (\code{get.af.spectra()}) and fluorophore variant assessment
#'   (\code{get.fluor.variants()}) may use full-spectra OLS unmixing as part
#'   of their SOM clustering input, positivity selection, and Spillover
#'   Spreading Matrix construction. Set to \code{FALSE} when \code{spectra}
#'   contains several similar or collinear fluorophores (e.g. a bead-cell
#'   comparison panel), where a full-spectra unmix is itself unstable or
#'   unsolvable and would corrupt rather than inform those steps. When
#'   \code{FALSE}, clustering falls back to raw detector space only, the
#'   unstained-sample positivity thresholds used internally by
#'   \code{get.fluor.variants()} are not computed, and the returned
#'   \code{spillover.spread} is always \code{NULL}.
#' @param unstained.sample Optional file path to a cell-based unstained FCS
#'   file, used as the autofluorescence reference when the control file has
#'   no \code{"AF"} row. Required in that case; ignored (with a message) if
#'   the control file does have an \code{"AF"} row, since the in-situ
#'   unstained paired with the single-stained controls is used instead.
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
#'   \item{\code{noise.floor}}{Named numeric vector, per-detector electronic
#'     noise floor in signal units (SD), pooled by minimum across controls.
#'     Matches the units of `noise.floor` elsewhere in the package
#'     (`unmix.fcs()`, the C++ pipeline). Square it before passing to
#'     `estimate.noise.model(read.var.floor = ...)`, which expects a
#'     variance.}
#'   \item{\code{spillover.spread}}{Matrix (source fluorophore x target
#'     channel), the Spillover Spreading Matrix: increase in unmixed
#'     variance a source fluorophore's positive population contributes to
#'     each other channel, per unit of its own on-channel signal. Diagonal
#'     entries are `NA`. `NULL` if no control supplied enough positive
#'     events. Saved as a heatmap when `figures = TRUE`.}
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
    parallel           = TRUE,
    verbose            = TRUE,
    threads            = NULL,
    n.cells            = 10000L,
    som.dim            = 5L,
    k.neighbors        = 3L,
    sim.threshold      = 0.985,
    sim.threshold.floor    = 0.90,
    af.collinear.threshold = 0.95,
    noise.floor.tail.fraction = 0.20,
    variant.fill.color = "red",
    variant.fill.alpha = 0.7,
    median.line.color  = "black",
    median.linewidth   = 1,
    use.unmixed                  = TRUE,
    unstained.sample             = NULL,
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
  # captured before any coercion below, since as.matrix()/storage.mode<- are
  # not guaranteed to preserve the "fluorophore" attribute
  spectra.fluorophore <- attr( spectra, "fluorophore" )
  spectra <- as.matrix( spectra )
  storage.mode( spectra ) <- "double"
  attr( spectra, "fluorophore" ) <- spectra.fluorophore
  if ( is.null( rownames( spectra ) ) )
    stop( "`spectra` must have rownames giving fluorophore names (including \"AF\").", call. = FALSE )

  # `spectra` must contain exactly one row per fluorophore. Checked here,
  # independently of the control-file check below, since `spectra` is a
  # separate argument that may have been built earlier, cached, or loaded
  # from a CSV, and so is not guaranteed to correspond to `control.def.file`.
  # Unlike check.spectra.duplicates() elsewhere, the fallback path here (used
  # when `spectra` carries no "fluorophore" attribute) is a hard stop, not a
  # warning: unmixing against two near-duplicate reference rows lets a solve
  # push arbitrarily large, opposite-signed signal into each of them, and the
  # resulting per-fluorophore "variance" this function computes would be
  # numerical noise rather than biology -- there's no safe way to proceed.
  if ( !is.null( spectra.fluorophore ) ) {

    check.spectra.duplicates( spectra )

  } else {

    # No identity attribute (e.g. `spectra` was loaded via read.spectra()).
    # A `sample`-disambiguated rowname (e.g. "PE (cells)", "PE (cells) (CD4)")
    # is unique by construction and so cannot reveal a duplicate directly;
    # recover true identity instead by matching the leading fluorophore name
    # in each rowname against the fluorophore database.
    non.af.rows <- rownames( spectra )[ rownames( spectra ) != "AF" ]
    fluor.db.path <- system.file(
      "extdata", "fluorophore_database.csv", package = "AutoSpectral"
    )
    fluor.db <- utils::read.csv( fluor.db.path, stringsAsFactors = FALSE )

    check.fluor <- if ( verbose ) {
      match.fluorophores( non.af.rows, fluor.db )
    } else {
      suppressMessages( match.fluorophores( non.af.rows, fluor.db ) )
    }
    names( check.fluor ) <- NULL

    unmatched <- non.af.rows[ check.fluor == "No match" ]
    if ( length( unmatched ) > 0 )
      warning(
        paste0(
          "Could not resolve a fluorophore identity for the following ",
          "`spectra` row(s) against the fluorophore database, so they ",
          "could not be checked for duplicate controls: ",
          paste( unmatched, collapse = ", " ), "."
        ),
        call. = FALSE
      )

    dup.fluor <- unique(
      check.fluor[ duplicated( check.fluor ) & check.fluor != "No match" ]
    )

    if ( length( dup.fluor ) > 0 )
      stop(
        paste0(
          "get.spectral.variants() does not support multiple controls per ",
          "fluorophore. `spectra` appears to contain more than one row for: ",
          paste( dup.fluor, collapse = ", " ), " (rows: ",
          paste( non.af.rows[ check.fluor %in% dup.fluor ], collapse = ", " ),
          "). Reduce `spectra` to a single row per fluorophore before ",
          "calling get.spectral.variants()."
        ),
        call. = FALSE
      )
  }

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

  if ( !is.numeric( sim.threshold.floor ) || length( sim.threshold.floor ) != 1 ||
       sim.threshold.floor < 0 || sim.threshold.floor > sim.threshold )
    .type.err( "sim.threshold.floor",
               "a single number in [0, sim.threshold]", sim.threshold.floor )

  if ( !is.numeric( af.collinear.threshold ) || length( af.collinear.threshold ) != 1 ||
       af.collinear.threshold < 0 || af.collinear.threshold > 1 )
    .type.err( "af.collinear.threshold", "a single number in [0, 1]", af.collinear.threshold )

  if ( !is.logical( use.unmixed ) || length( use.unmixed ) != 1 || is.na( use.unmixed ) )
    .type.err( "use.unmixed", "a single TRUE/FALSE value", use.unmixed )

  if ( !is.null( unstained.sample ) && ( !is.character( unstained.sample ) || length( unstained.sample ) != 1 ) )
    .type.err( "unstained.sample", "NULL or a single character path", unstained.sample )

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
  # get.spectral.variants() does not support multiple controls per
  # fluorophore (see the `table.fluors` duplicate check below for why), so
  # this is always strict regardless of what the rest of the pipeline allows
  check.control.file(
    control.dir, control.def.file, asp, strict = TRUE,
    allow.duplicate.controls = FALSE
  )

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

  # per-sample metadata. get.spectral.variants() unmixes against `spectra`
  # and computes per-fluorophore variance from the result, so two controls
  # for the same fluorophore is not just a labeling ambiguity here: unmixing
  # against two near-identical reference rows lets a solve push arbitrarily
  # large, opposite-signed signal into each of them, and the resulting
  # "variance" is numerical noise, not biology. `sample`-level disambiguation
  # (used elsewhere in the package to permit multiple controls) is
  # deliberately not applied here.
  table.fluors <- control.table$fluorophore
  table.fluors <- table.fluors[ !is.na( table.fluors ) ]

  if ( anyDuplicated( table.fluors ) != 0 )
    stop(
      paste0(
        "get.spectral.variants() does not support multiple controls per ",
        "fluorophore. Duplicated in `control.def.file`: ",
        paste( unique( table.fluors[ duplicated( table.fluors ) ] ), collapse = ", " ), "."
      ),
      call. = FALSE
    )

  universal.negative <- control.table$universal.negative
  universal.negative[ is.na( universal.negative ) ] <- "FALSE"
  names( universal.negative ) <- table.fluors
  flow.channel       <- control.table$channel
  names( flow.channel ) <- table.fluors
  flow.file.name     <- control.table$filename
  names( flow.file.name ) <- table.fluors
  control.type <- control.table$control.type
  names( control.type ) <- table.fluors

  has.af.row <- "AF" %in% table.fluors

  if ( !has.af.row && is.null( unstained.sample ) )
    stop(
      "An unstained cell control is required for get.spectral.variants(): ",
      "either include an `AF` row in the control file, or supply the ",
      "`unstained.sample` argument.",
      call. = FALSE
    )

  if ( has.af.row && !is.null( unstained.sample ) && verbose )
    message(
      "\033[33mBoth an `AF` row in the control file and `unstained.sample` ",
      "were supplied; using the in-situ `AF` control from the control ",
      "file.\033[0m"
    )

  unstained.file <- if ( has.af.row )
    file.path( control.dir, flow.file.name[ "AF" ] )
  else
    unstained.sample

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

  unstained <- readFCS( unstained.file, columns = spectral.channel )

  if ( nrow( unstained ) > asp$gate.downsample.n.cells ) {
    set.seed( asp$bird.seed )
    unstained.idx <- sample( nrow( unstained ), asp$gate.downsample.n.cells )
    unstained     <- unstained[ unstained.idx, , drop = FALSE ]
  }

  raw.thresholds <- apply( unstained, 2, function( col )
    stats::quantile( col, 0.995 ) )

  # get AF spectra in place
  af.spectra <- get.af.spectra(
    unstained.file,
    asp,
    spectra,
    som.dim = 10,
    figures = FALSE,
    save = FALSE,
    use.unmixed = use.unmixed,
    refine = FALSE,
    parallel = parallel,
    threads = threads
  )

  # derive per-file AF PCs for all unique unstained cell files used as negatives
  cell.fluors <- names( control.type )[ control.type == "cells" ]
  univ.neg.files <- unique( universal.negative[
    names( universal.negative ) %in% cell.fluors &
      universal.negative != "FALSE" &
      !is.na( universal.negative ) &
      grepl( "\\.fcs$", universal.negative, ignore.case = TRUE )
  ] )

  # read each unique universal-negative file once
  neg.cache <- lapply( univ.neg.files, function( fn ) {
    dat <- readFCS(
      file.path( control.dir, fn ),
      columns = union( spectral.channel, scatter.channel )
    )
    list(
      spectral = dat[ , spectral.channel, drop = FALSE ],
      scatter  = dat[ , scatter.channel,  drop = FALSE ]
    )
  } )
  names( neg.cache ) <- univ.neg.files

  af.pcs.list <- lapply( neg.cache, function( nc ) {
    dat <- nc$spectral
    if ( nrow( dat ) > asp$gate.downsample.n.cells ) {
      set.seed( asp$bird.seed )
      dat <- dat[ sample( nrow( dat ), asp$gate.downsample.n.cells ), , drop = FALSE ]
    }
    sv <- svd( dat, nu = 0, nv = 4 )
    t( sv$v )
  } )
  names( af.pcs.list ) <- univ.neg.files

  # Parallel setup
  threads <- if ( isTRUE( parallel ) ) {
    if ( is.null( threads ) ) 0L else as.integer( threads )
  } else {
    1L
  }

  # find the likely positivity thresholds for determining what needs
  # refinement. Skipped when `use.unmixed = FALSE`: a full-pipeline unmix
  # against several similar or collinear fluorophores is exactly the
  # operation `use.unmixed = FALSE` is meant to avoid, and the resulting
  # thresholds are unused downstream in that mode -- get.fluor.variants()
  # falls back to raw-threshold selection instead.
  if ( use.unmixed ) {

    unstained.unmixed <- if (
      requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
      "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
    ) {
      AutoSpectralRcpp::unmix.autospectral.rcpp(
        raw.data = unstained,
        spectra = spectra,
        af.spectra = af.spectra,
        verbose = FALSE,
        parallel = TRUE,
        threads = threads
      )
    } else {
      unmix.autospectral(
        raw.data = unstained,
        spectra = spectra,
        af.spectra = af.spectra,
        asp = asp,
        verbose = FALSE,
        parallel = TRUE,
        threads = threads
      )
    }

    unmixed.thresholds <- apply(
      unstained.unmixed[ , fluorophores, drop = FALSE ], 2, function( col )
        stats::quantile( col, 0.995 )
    )

  } else {

    unstained.unmixed  <- NULL
    unmixed.thresholds <- stats::setNames(
      rep( NA_real_, length( fluorophores ) ), fluorophores )

  }

  # set up main loop call
  if ( is.null( names( table.fluors ) ) ) names( table.fluors ) <- table.fluors

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
    neg.cache          = neg.cache,
    use.unmixed        = use.unmixed,
    n.cells            = n.cells,
    som.dim            = som.dim,
    k.neighbors        = k.neighbors,
    sim.threshold      = sim.threshold,
    sim.threshold.floor    = sim.threshold.floor,
    af.collinear.threshold = af.collinear.threshold,
    noise.floor.tail.fraction = noise.floor.tail.fraction,
    variant.fill.color = variant.fill.color,
    variant.fill.alpha = variant.fill.alpha,
    median.line.color  = median.line.color,
    median.linewidth   = median.linewidth,
    parallel           = parallel,
    threads            = threads
  )

  # ---------------------------------------------------------------------------
  # Main loop
  # ---------------------------------------------------------------------------

  if ( verbose )
    message( paste0( "\033[34m", "Identifying spectral variation", "\033[0m" ) )

  # initialise with base spectra as safe fallback
  spectral.variants <- lapply( table.fluors, function( fl )
    spectra[ fl, , drop = FALSE ] )
  names( spectral.variants ) <- table.fluors

  updated.variants <- lapply( table.fluors, function( f ) {
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

  # ---------------------------------------------------------------------------
  # Noise floor, pooled across controls
  # ---------------------------------------------------------------------------
  # Each control supplies an upper bound per detector, in signal units (SD),
  # matching every other `noise.floor` in this package. Every detector is far
  # from SOME control's emission, so the minimum across controls is the
  # tightest available bound. With few controls the minimum is noisy; with
  # many, consider a low quantile instead of the strict minimum.

  noise.floor <- NULL
  floor.list  <- lapply( spectral.variants, function( v ) attr( v, "noise.floor" ) )
  floor.list  <- floor.list[ !vapply( floor.list, is.null, logical( 1 ) ) ]

  if ( length( floor.list ) > 0 ) {

    floor.mat   <- do.call( rbind, floor.list )
    noise.floor <- apply( floor.mat, 2, function( v ) {
      v <- v[ is.finite( v ) & v > 0 ]
      if ( length( v ) == 0 ) NA_real_ else min( v )
    } )

    if ( anyNA( noise.floor ) )
      warning( "Noise floor could not be estimated at ",
               sum( is.na( noise.floor ) ), " detector(s).", call. = FALSE )

    if ( verbose )
      message( sprintf(
        "Noise floor from %d control(s): median SD %.1f (range %.1f - %.1f)",
        nrow( floor.mat ),
        stats::median( noise.floor, na.rm = TRUE ),
        min( noise.floor, na.rm = TRUE ),
        max( noise.floor, na.rm = TRUE ) ) )
  }

  # ---------------------------------------------------------------------------
  # Spillover Spreading Matrix
  # ---------------------------------------------------------------------------
  # For source fluorophore a and target channel b, the increase in unmixed
  # variance that a's positive population contributes to channel b, per unit
  # of a's own on-channel signal:
  #
  #   SS(a, b) = ( mad(pos_a(b))^2 - mad(neg(b))^2 ) / ( median(pos_a(a)) - median(neg(a)) )
  #
  # Computed in unmixed (compensated) units, so it reflects what actually
  # limits resolution after spectral unmixing rather than raw detector
  # crosstalk. Diagonal entries (a source fluorophore against its own
  # channel) are set to NA -- that is the width of the positive peak itself,
  # not spillover, and left in it dominates the colour scale.

  spillover.spread <- NULL
  spread.list <- lapply( spectral.variants, function( v ) attr( v, "spillover.spread" ) )
  mfi.list    <- lapply( spectral.variants, function( v ) attr( v, "on.channel.mfi" ) )
  have.spread <- !vapply( spread.list, is.null, logical( 1 ) )

  if ( any( have.spread ) ) {

    spread.fluors <- names( spectral.variants )[ have.spread ]
    neg.mad       <- apply( unstained.unmixed[ , fluorophores, drop = FALSE ], 2, stats::mad )
    neg.median    <- apply( unstained.unmixed[ , fluorophores, drop = FALSE ], 2, stats::median )

    spillover.spread <- matrix(
      NA_real_, nrow = length( spread.fluors ), ncol = length( fluorophores ),
      dimnames = list( spread.fluors, fluorophores )
    )

    for ( a in spread.fluors ) {
      pos.mad <- spread.list[[ a ]][ fluorophores ]
      denom   <- mfi.list[[ a ]] - neg.median[ a ]
      if ( !is.finite( denom ) || denom <= 0 ) next
      spillover.spread[ a, ] <- ( pos.mad^2 - neg.mad^2 ) / denom
      if ( a %in% colnames( spillover.spread ) ) spillover.spread[ a, a ] <- NA_real_
    }

    if ( verbose )
      message( sprintf(
        "Spillover spread matrix computed for %d of %d fluorophore(s)",
        length( spread.fluors ), length( fluorophores ) ) )

    if ( figures ) {
      tryCatch(
        expr = {
          create.heatmap(
            matrix        = spillover.spread,
            title         = "spillover_spread",
            legend.label  = "Spread (var / on-channel MFI)",
            plot.dir      = output.dir,
            color.palette = "viridis"
          )
        },
        error = function( e ) {
          warning( "Spillover spread heatmap failed: ", conditionMessage( e ),
                   call. = FALSE )
        }
      )
    }
  }

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

      stained.raw <- readFCS( stained.sample, columns = spectral.channel )

      if ( nrow( stained.raw ) > 5000 ) {
        set.seed( asp$bird.seed )
        stained.raw <- stained.raw[ sample( nrow( stained.raw ), 5000 ), , drop = FALSE ]
      }

      stained.unmixed <- if (
        requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
        "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
      ) {
        AutoSpectralRcpp::unmix.autospectral.rcpp(
          raw.data = stained.raw,
          spectra = spectra,
          af.spectra = af.spectra,
          verbose = FALSE,
          parallel = TRUE,
          threads = threads
        )
      } else {
        unmix.autospectral(
          raw.data = stained.raw,
          spectra = spectra,
          af.spectra = af.spectra,
          asp = asp,
          verbose = FALSE,
          parallel = TRUE,
          threads = threads
        )
      }

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
    noise.floor = noise.floor,
    spillover.spread  = spillover.spread,
    optimize.scores      = necessity$scores.norm,
    optimize.recommended = necessity$optimize.recommended
  )

  saveRDS( variants, file = file.path( output.dir, asp$variant.filename ) )

  return( variants )
}
