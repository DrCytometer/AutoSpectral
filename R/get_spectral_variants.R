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
#' Uses the Spillover Spreading Matrix  built from \code{get.fluor.variants()}'s
#' Residual Model regression rather than an empirical MAD ratio
#' against a separately-unmixed unstained baseline. There is no
#' \code{spread.denom.min.mad}/\code{"snr"} concept here: each control's own
#' regression already separates its baseline (intercept) from its
#' abundance-scaled spread (slope), estimated from its own dim-to-bright
#' event range rather than compared against a different sample's estimate.
#' A row is instead trusted once its regression used at least
#' \code{spread.min.events} events; below that, the hotspot-matrix fallback
#' (unchanged from \code{get.spectral.variants()}) fills the row when enough
#' trusted rows exist to calibrate against.
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
#'   fluorophore used for SOM clustering. Passed to \code{get.fluor.variants}.
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
#'   retains fewer than 20 events.
#' @param af.collinear.threshold Numeric, default \code{0.95}. Minimum
#'   cosine similarity between a fluorophore's reference spectrum and any of
#'   its paired unstained file's AF principal directions at or above which
#'   the AF-component projection step (and the low-rank fit feeding the
#'   Residual Model) is skipped.
#' @param noise.floor.tail.fraction Numeric in (0, 1), default \code{0.20}.
#'   Fraction of each detector's raw values (lowest end) used to estimate the
#'   per-control noise floor. Passed to \code{get.fluor.variants}.
#' @param noise.spillover.floor Numeric in (0, 1), default \code{0.002}. A
#'   detector is excluded from a control's contribution to the pooled
#'   noise-model regression when that control's own reference spectrum
#'   falls below this fraction of its own peak there, in addition to the
#'   upper exclusion already applied for detectors it dominates (see
#'   \code{get.fluor.variants(noise.mask.threshold)}).
#' @param spread.min.events Integer, default \code{50}. Minimum number of
#'   events behind a source fluorophore's Residual Model regression
#'   (\code{"spillover.spread.n"}) before its Spillover Spreading Matrix row
#'   is trusted. A 2-parameter regression needs more support than the old
#'   MAD point estimate did; rows below this are left blank unless filled by
#'   the hotspot-matrix fallback.
#' @param spread.hotspot.fallback Logical, default \code{TRUE}. When a
#'   source fluorophore's Spillover Spreading Matrix row fails the
#'   \code{spread.min.events} check (a weak or under-titrated control),
#'   fill that row from \code{calculate.hotspot.matrix(spectra)} instead of
#'   leaving it blank. The hotspot matrix is a purely geometric measure of
#'   pairwise spread susceptibility from the reference spectra alone;
#'   filling uses a single calibration constant (the median ratio of
#'   trusted \code{spillover.spread} entries to their hotspot-matrix
#'   counterparts), so it requires at least \code{spread.hotspot.min.pairs}
#'   trusted entries to calibrate against. Filled rows are tagged
#'   \code{"hotspot"} in the returned matrix's \code{"source"} attribute.
#' @param spread.hotspot.min.pairs Integer, default \code{20}. Minimum
#'   number of trusted (source, target) entries required to calibrate the
#'   hotspot-matrix fallback. Below this, weak controls are left blank as
#'   before and a message explains why.
#' @param huber.k Numeric, default \code{1.345}. Huber tuning constant
#'   passed to \code{get.fluor.variants()}'s spillover-spread
#'   regression.
#' @param huber.max.iter Integer, default \code{100L}. Maximum IRLS
#'   iterations passed to the same regression.
#' @param variant.fill.color Color for the shaded ribbon in variant plots.
#'   Default \code{"red"}.
#' @param variant.fill.alpha Alpha for \code{variant.fill.color}. Default
#'   \code{0.7}.
#' @param median.line.color Color for the reference-spectrum line. Default
#'   \code{"black"}.
#' @param median.linewidth Width of the reference-spectrum line. Default
#'   \code{1}.
#' @param use.unmixed Logical, default \code{TRUE}. Whether AF extraction and
#'   fluorophore variant assessment may use full-spectra OLS unmixing as
#'   part of their SOM clustering input, positivity selection, and Spillover
#'   Spreading Matrix construction. Set to \code{FALSE} when \code{spectra}
#'   contains several similar or collinear fluorophores (e.g. a bead-cell
#'   comparison panel). When \code{FALSE}, the returned \code{spillover.spread}
#'   is always \code{NULL}.
#' @param unstained.sample Optional file path to a cell-based unstained FCS
#'   file, used as the autofluorescence reference when the control file has
#'   no \code{"AF"} row.
#' @param stained.sample Optional file path to a representative stained FCS
#'   file, weighting the optimization necessity scores by fluorophore
#'   brightness. Pass `NULL` (default) to use purely geometric scores.
#' @param optimize.necessity.threshold Numeric in `[0, 1]`, default `0.01`.
#'   Passed to `calculate.optimize.necessity()`.
#' @param ... Ignored. Catches and warns on previously used deprecated
#'   arguments: \code{af.spectra}, \code{refine}, \code{problem.quantile},
#'   \code{pos.quantile}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{thresholds}}{Named numeric vector of positivity thresholds in
#'     the unmixed space, one per fluorophore.}
#'   \item{\code{neg.thresholds}}{Named numeric vector, the 0.5th percentile of
#'     each fluorophore's unstained unmixed distribution.}
#'   \item{\code{variants}}{Named list of variant-spectra matrices, one per
#'     fluorophore.}
#'   \item{\code{delta.list}}{Named list of delta matrices (variant minus
#'     reference spectrum), one per fluorophore.}
#'   \item{\code{delta.norms}}{Named list of Euclidean norms of the deltas.}
#'   \item{\code{noise.floor}}{Named numeric vector, per-detector electronic
#'     noise floor in signal units (SD), pooled by minimum across controls.}
#'   \item{\code{spillover.spread}}{Matrix (source fluorophore x target
#'     channel), the Residual Model Spillover Spreading Matrix: the
#'     Huber-robust slope of each target's squared residual-projection
#'     against the source's own recovered abundance -- added unmixed
#'     variance per unit of the source's on-channel abundance. Diagonal
#'     entries are `NA`. Rows below `spread.min.events` are left `NA` across
#'     the row unless filled by the hotspot-matrix fallback. Carries
#'     `"n.events"` and `"source"` attributes (named integer/character
#'     vectors, one per source fluorophore): event count behind the
#'     regression, and whether the row came from the regression
#'     (`"residual"`) or the hotspot-matrix fallback (`"hotspot"`). `NULL`
#'     if no control supplied enough positive events. Saved as a heatmap
#'     when `figures = TRUE`.}
#'   \item{\code{spillover.spread.intercept}}{Matrix, same shape as
#'     \code{spillover.spread}: each source's regression intercept per
#'     target channel -- that source's own estimate of the target's
#'     baseline unmixed variance at zero abundance. Not filled by the
#'     hotspot fallback (the hotspot matrix has no calibrated intercept
#'     term); `NA` wherever \code{spillover.spread} came from that
#'     fallback or was left blank.}
#' }
#' The list is also saved as an .rds file in \code{output.dir}.
#'
#' @references
#' Cai X et al. (2026). Residual Model for unmixed spread prediction in
#' spectral flow cytometry. \emph{bioRxiv} 2026.01.27.701929.
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
    noise.spillover.floor = 0.002,
    spread.min.events        = 50L,
    spread.hotspot.fallback  = TRUE,
    spread.hotspot.min.pairs = 20,
    huber.k         = 1.345,
    huber.max.iter  = 100L,
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

  # `spectra` must contain exactly one row per fluorophore.
  if ( !is.null( spectra.fluorophore ) ) {

    check.spectra.duplicates( spectra )

  } else {

    non.af.rows <- rownames( spectra )[ rownames( spectra ) != "AF" ]
    fluor.db.path <- system.file(
      "extdata", "fluorophore_database.csv", package = "AutoSpectral"
    )
    fluor.db <- utils::read.csv( fluor.db.path, stringsAsFactors = FALSE )

    check.fluor <- match.fluorophores( non.af.rows, fluor.db, verbose = FALSE )
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

  if ( verbose ) message( "\033[32mChecking control file and fluorophore labels for errors \033[0m" )
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
  fluor.to.match <- table.fluors[ !grepl( "Negative|^AF$", table.fluors ) ]
  if ( !all( fluor.to.match %in% fluorophores ) ) {
    matching.fluors <- fluor.to.match %in% fluorophores
    if ( !any( matching.fluors ) )
      stop( "No matching fluorophores between `spectra` and the control file.",
            call. = FALSE )
    if ( !all( matching.fluors ) )
      warning(
        sprintf(
          "Some fluorophores in the control file are absent from `spectra`: %s.",
          paste( fluor.to.match[ !matching.fluors ], collapse = ", " )
        ),
        call. = FALSE
      )
    table.fluors <- fluor.to.match[ matching.fluors ]
  } else {
    table.fluors <- fluor.to.match
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
  # refinement. Skipped when `use.unmixed = FALSE`.
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

    neg.thresholds <- apply(
      unstained.unmixed[ , fluorophores, drop = FALSE ], 2, function( col )
        stats::quantile( col, 0.005 )
    )

  } else {

    unstained.unmixed  <- NULL
    unmixed.thresholds <- stats::setNames(
      rep( NA_real_, length( fluorophores ) ), fluorophores )
    neg.thresholds <- stats::setNames(
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
    huber.k            = huber.k,
    huber.max.iter     = huber.max.iter,
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
  # Noise model (read.var, kappa), pooled across single-stained controls
  # ---------------------------------------------------------------------------
  if ( verbose )
    message( paste0( "\033[34m", "Modelling detector noise", "\033[0m" ) )

  noise.model <- NULL
  events.list <- lapply( spectral.variants, function( v ) attr( v, "noise.events" ) )
  mask.list   <- lapply( spectral.variants, function( v ) attr( v, "noise.mask" ) )
  have.noise  <- !vapply( events.list, is.null, logical( 1 ) )

  if ( sum( have.noise ) >= 2L ) {

    noise.fluors <- names( spectral.variants )[ have.noise ]

    fit.list <- lapply( noise.fluors, function( fl ) {
      ev  <- events.list[[ fl ]]
      ref <- spectra[ fl, , drop = FALSE ]
      co  <- unmix.ols( ev, ref )
      list( fitted = co %*% ref, resid = ev - ( co %*% ref ) )
    } )
    names( fit.list ) <- noise.fluors

    pool.fitted  <- do.call( rbind, lapply( fit.list, function( x ) x$fitted ) )
    pool.resid   <- do.call( rbind, lapply( fit.list, function( x ) x$resid ) )
    pool.file.id <- unlist( lapply( noise.fluors, function( fl )
      rep( fl, nrow( events.list[[ fl ]] ) ) ) )

    row.start <- 1L
    for ( fl in noise.fluors ) {
      n.fl     <- nrow( events.list[[ fl ]] )
      rows     <- row.start:( row.start + n.fl - 1L )
      low.mask <- spectra[ fl, ] < noise.spillover.floor * max( spectra[ fl, ] )
      excl     <- mask.list[[ fl ]] | low.mask
      pool.resid[ rows, excl ] <- NA_real_
      row.start <- row.start + n.fl
    }

    noise.model <- tryCatch(
      suppressWarnings( .fit.noise.regression(
        y.hat          = pool.fitted,
        resid          = pool.resid,
        det.names      = spectral.channel,
        read.var.floor = if ( length( floor.list ) > 0 ) noise.floor^2 else NULL,
        unstained.data = unstained,
        file.id        = pool.file.id,
        verbose        = FALSE
      ) ),
      error = function( e ) {
        warning( "Pooled single-stained-control noise model failed: ",
                 conditionMessage( e ), call. = FALSE )
        NULL
      }
    )

    if ( verbose && !is.null( noise.model ) )
      message( sprintf(
        "Noise model from %d control(s): median read SD %.1f, median counts.per.unit %.3g",
        length( noise.fluors ),
        stats::median( sqrt( noise.model$read.var ) ),
        stats::median( noise.model$counts.per.unit ) ) )

  } else if ( verbose ) {
    message( "Fewer than 2 controls carried noise-model events; ",
             "skipping pooled noise-model estimation." )
  }

  # ---------------------------------------------------------------------------
  # Spillover Spreading Matrix (Residual Model)
  # ---------------------------------------------------------------------------
  # Each source fluorophore's row comes directly from get.fluor.variants()'s
  # regression: "spillover.spread" is already the slope -- SS(source, .) in
  # variance-per-unit-abundance units -- with "spillover.spread.intercept" the
  # matching baseline. No separate unstained-population baseline is combined
  # in here; that step, and the estimator mismatch it introduced, no longer
  # exists in this version.

  spillover.spread           <- NULL
  spillover.spread.intercept <- NULL

  slope.list     <- lapply( spectral.variants, function( v ) attr( v, "spillover.spread" ) )
  intercept.list <- lapply( spectral.variants, function( v ) attr( v, "spillover.spread.intercept" ) )
  n.list         <- lapply( spectral.variants, function( v ) attr( v, "spillover.spread.n" ) )
  source.list    <- lapply( spectral.variants, function( v ) attr( v, "spillover.spread.source" ) )
  have.spread    <- !vapply( slope.list, is.null, logical( 1 ) )

  if ( any( have.spread ) ) {

    spread.fluors <- names( spectral.variants )[ have.spread ]

    spillover.spread <- do.call( rbind, lapply( spread.fluors, function( a )
      slope.list[[ a ]][ fluorophores ] ) )
    dimnames( spillover.spread ) <- list( spread.fluors, fluorophores )

    spillover.spread.intercept <- do.call( rbind, lapply( spread.fluors, function( a )
      intercept.list[[ a ]][ fluorophores ] ) )
    dimnames( spillover.spread.intercept ) <- list( spread.fluors, fluorophores )

    for ( a in spread.fluors ) {
      if ( a %in% colnames( spillover.spread ) ) {
        spillover.spread[ a, a ]           <- NA_real_
        spillover.spread.intercept[ a, a ] <- NA_real_
      }
    }

    spread.n      <- stats::setNames( rep( NA_integer_, length( spread.fluors ) ), spread.fluors )
    spread.source <- stats::setNames( rep( NA_character_, length( spread.fluors ) ), spread.fluors )
    for ( a in spread.fluors ) {
      spread.n[ a ]      <- if ( is.null( n.list[[ a ]] ) ) NA_integer_ else as.integer( n.list[[ a ]] )
      spread.source[ a ] <- if ( is.null( source.list[[ a ]] ) ) NA_character_ else source.list[[ a ]]
    }

    attr( spillover.spread, "n.events" ) <- spread.n
    attr( spillover.spread, "source" )   <- spread.source

    trusted <- spread.fluors[ !is.na( spread.n ) & spread.n >= spread.min.events ]
    weak    <- setdiff( spread.fluors, trusted )

    if ( verbose )
      message( sprintf(
        "Spillover spread (Residual Model) computed for %d of %d fluorophore(s); %d below spread.min.events = %d",
        length( trusted ), length( fluorophores ), length( weak ), spread.min.events ) )

    # -------------------------------------------------------------------
    # Hotspot-matrix fallback for weak controls
    # -------------------------------------------------------------------
    # Unchanged in spirit from get.spectral.variants(): calibrate a single
    # constant from the panel's own trusted rows against
    # calculate.hotspot.matrix(spectra), then use it to fill rows below
    # spread.min.events. No intercept fallback is attempted -- the hotspot
    # matrix has no calibrated baseline term, so
    # spillover.spread.intercept is left NA for any row filled this way.

    if ( spread.hotspot.fallback && length( trusted ) > 0 && length( weak ) > 0 ) {

      hotspot <- calculate.hotspot.matrix( spectra )
      hotspot <- hotspot[ rownames( spillover.spread ), colnames( spillover.spread ), drop = FALSE ]

      cal.ratio <- spillover.spread[ trusted, , drop = FALSE ] / hotspot[ trusted, , drop = FALSE ]
      cal.ratio <- cal.ratio[ is.finite( cal.ratio ) & hotspot[ trusted, , drop = FALSE ] > 1e-6 ]

      if ( length( cal.ratio ) >= spread.hotspot.min.pairs ) {

        k <- stats::median( cal.ratio )

        for ( a in weak ) {
          filled <- k * hotspot[ a, ]
          filled[ !is.finite( filled ) ] <- NA_real_
          spillover.spread[ a, ] <- filled
          if ( a %in% colnames( spillover.spread ) ) spillover.spread[ a, a ] <- NA_real_
          spread.source[ a ] <- "hotspot"
        }

        attr( spillover.spread, "source" ) <- spread.source

        if ( verbose )
          message( sprintf(
            "Filled %d weak control(s) from the hotspot matrix (calibration k = %.3g from %d trusted pair(s)): %s",
            length( weak ), k, length( cal.ratio ), paste( weak, collapse = ", " ) ) )

      } else if ( verbose ) {
        message( "Too few trusted pairs to calibrate a hotspot-matrix fallback; ",
                 "weak control(s) left blank: ", paste( weak, collapse = ", " ) )
      }
    }

    # A source's own on-channel signal cannot suppress variance elsewhere on
    # average; the true slope and intercept are both >= 0. Small negative
    # entries are Huber-fit noise around zero, not real values -- clip them
    # rather than carry them into downstream thresholds and weighting.
    spillover.spread[ !is.na( spillover.spread ) & spillover.spread < 0 ] <- 0
    spillover.spread.intercept[
      !is.na( spillover.spread.intercept ) & spillover.spread.intercept < 0 ] <- 0
    # use the unit-normalized spillover spread for visualization
    ssm <- l2.normalize.spectra( spillover.spread )
    # save the data as CSV
    utils::write.csv(
      ssm,
      file = file.path(
        asp$figure.similarity.heatmap.dir,
        paste0( "Normalized_", asp$spillover.spread.file.name, ".csv" )
      )
    )
    utils::write.csv(
      spillover.spread,
      file = file.path(
        asp$figure.similarity.heatmap.dir,
        paste0( asp$spillover.spread.file.name, ".csv" )
      )
    )

    if ( figures ) {

      tryCatch(
        expr = {
          spectral.heatmap(
            spectra       = ssm,
            title         = asp$spillover.spread.file.name,
            plot.dir      = asp$figure.similarity.heatmap.dir,
            legend.label  = "Normalized Spillover Spread",
            color.palette = "magma"
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
    thresholds     = unmixed.thresholds,
    neg.thresholds = neg.thresholds,
    variants    = spectral.variants,
    delta.list  = delta.list,
    delta.norms = delta.norms,
    noise.floor = noise.floor,
    noise.model = noise.model,
    spillover.spread           = spillover.spread,
    spillover.spread.intercept = spillover.spread.intercept,
    optimize.scores      = necessity$scores.norm,
    optimize.recommended = necessity$optimize.recommended
  )

  saveRDS( variants, file = file.path( output.dir, asp$variant.filename ) )

  return( variants )
}
