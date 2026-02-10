# get_spectral_variants.r

#' @title Get Spectral Variations for Fluorophores
#'
#' @description
#' This function cycles through all the fluorophores defined in `control.def.file`,
#' identifying variations in spectral profiles. It does this by performing SOM
#' clustering on the positive events in the cleaned control data. The output is
#' saved as an .rds file, and figures summarizing the variation are saved, if
#' desired. Note that the .rds file contains all the needed information for
#' downstream processing (per-cell unmixing), so you can just load that using
#' the `readRDS` function) rather than re-running this process.
#'
#' @importFrom lifecycle deprecate_warn
#' @importFrom flowCore read.FCS exprs
#'
#' @param control.dir File path to the single-stained control FCS files.
#' @param control.def.file CSV file defining the single-color control file names,
#' fluorophores they represent, marker names, peak channels, and gating requirements.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param spectra A matrix containing the spectral data. Fluorophores in rows,
#' detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param n.cells Numeric, default `2000`. Number of cells to use for defining
#' the variation in spectra. Up to `n.cells` cells will be selected as positive
#' events in the peak channel for each fluorophore, above the 99.5th percentile
#' level in the unstained sample.
#' @param som.dim Numeric, default `10`. Number of x and y dimensions to use in
#' the SOM for clustering the spectral variation. The number of spectra returned
#' for each fluorophore will increase with the quadratic of `som.dim`, so for 10,
#' you will get up to 100 variants. Somewhere between 4 and 7 appears to be
#' sufficient, but with the pruning of variants implemented in
#' `unmix.autospectral()` in v1.0.0, this is less important.
#' @param figures Logical, controls whether the variation in spectra for each
#' fluorophore is plotted in `output.dir`. Default is `TRUE`.
#' @param output.dir File path to whether the figures and .rds data file will be
#' saved. Default is `NULL`, in which case `asp$variant.dir` will be used.
#' @param parallel Logical, default is `FALSE`, in which case sequential processing
#' will be used. The new parallel processing should always be faster.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#' @param threads Numeric, default is `NULL`, in which case `asp$worker.process.n`
#' will be used. `asp$worker.process.n` is set by default to be one less than the
#' available cores on the machine. Multi-threading is only used if `parallel` is
#' `TRUE`.
#' @param refine Logical, default is `TRUE`. Controls whether to perform a second
#' round of variation measurement on "problem cells", which are those with the
#' highest spillover, as defined by `problem.quantile`.
#' @param problem.quantile Numeric, default `0.95`. The quantile for determining
#' which cells will be considered "problematic" after unmixing with per-cell AF
#' extraction. Cells in the `problem.quantile` or above with respect to total
#' signal in the fluorophore (non-AF) channels after per-cell AF extraction will
#' be used to determine additional autofluorescence spectra, using a second round
#' of clustering and modulation of the previously selected autofluroescence
#' spectra. A value of `0.95` means the top 5% of cells, those farthest from zero,
#' will be selected for further investigation.
#' @param ... Ignored. Previously used for deprecated arguments such as
#' `pos.quantile` and `sim.threshold`, which are now fixed internally and no
#' longer user-settable.
#'
#' @return A vector with the indexes of events inside the initial gate.
#'
#' @export

get.spectral.variants <- function(
    control.dir,
    control.def.file,
    asp,
    spectra,
    af.spectra,
    n.cells = 2000,
    som.dim = 10,
    figures = TRUE,
    output.dir = NULL,
    parallel = FALSE,
    verbose = TRUE,
    threads = NULL,
    refine = TRUE,
    problem.quantile = 0.95,
    ...
) {

  dots <- list( ... )

  if ( !is.null( dots$sim.threshold ) ) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "get.spectral.variants(sim.threshold)",
      details = "no longer used"
    )
  }
  if ( !is.null( dots$pos.quantile ) ) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "get.spectral.variants(pos.quantile)",
      details = "no longer used"
    )
  }

  # check that af.spectra is a matrix and has rows
  if ( is.null( af.spectra ) ) {
    stop(
      "Multiple AF in `af.spectra` must be provided in a matrix. See `?get.af.spectra`.",
      call. = FALSE
    )
  }
  if ( !is.matrix( af.spectra ) || nrow( af.spectra ) < 2 ) {
    stop(
      "Multiple AF in `af.spectra` must be provided in a matrix. See `?get.af.spectra`.",
      call. = FALSE
    )
  }

  if ( is.null( output.dir ) )
    output.dir <- asp$variant.dir
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

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

  fluorophores <- rownames( spectra )[ rownames( spectra ) != "AF" ]
  spectral.channel <- colnames( spectra )

  # read control file
  if ( !file.exists( control.def.file ) ) {
    stop(
      paste( "Unable to locate control.def.file:", control.def.file ),
      call. = FALSE
    )
  }

  if ( verbose ) message( "\033[32mChecking control file for errors \033[0m" )
  check.control.file( control.dir, control.def.file, asp, strict = TRUE )

  control.table <- utils::read.csv(
    control.def.file,
    stringsAsFactors = FALSE,
    strip.white = TRUE
  )

  control.table[] <- lapply( control.table, function( x ) {
    if ( is.character( x ) ) {
      x <- trimws( x )
      x[ x == "" ] <- NA
      x
    } else x
  } )

  # set channels to be used
  spectral.channel <- colnames( spectra )
  flow.scatter.parameter <- read.scatter.parameter( asp )
  flow.scatter.and.channel.spectral <- c(
    asp$default.time.parameter,
    flow.scatter.parameter,
    spectral.channel
  )

  if ( grepl( "Discover", asp$cytometer ) ) {
    spectral.channel <- spectral.channel[
      grep( asp$spectral.channel, spectral.channel ) ]
  }

  # define list of samples
  flow.sample <- control.table$sample
  table.fluors <- control.table$fluorophore
  table.fluors <- table.fluors[ !is.na( table.fluors ) ]
  universal.negative <- control.table$universal.negative
  universal.negative[ is.na( universal.negative ) ] <- FALSE
  names( universal.negative ) <- table.fluors
  flow.channel <- control.table$channel
  names( flow.channel ) <- table.fluors
  flow.file.name <- control.table$filename
  names( flow.file.name ) <- table.fluors
  control.type <- control.table$control.type
  names( control.type ) <- table.fluors

  # stop if "AF" sample is not present, fluorophore mismatch
  if ( !( "AF" %in% table.fluors ) ) {
    stop(
      "Unable to locate `AF` control in control file. An unstained cell control is required.",
      call. = FALSE
    )
  }

  # check for data/spectra column matching
  spectra.cols <- colnames( spectra )

  if ( !identical( spectral.channel, spectra.cols ) ) {

    # ensure both actually have the same columns before reordering
    if ( all( spectra.cols %in% spectral.channel ) &&
        length( spectra.cols ) == length( spectral.channel ) ) {
      # reorder raw.data to match the order of spectra
      spectra <- spectra[ , spectral.channel ]
      message( "Columns of spectra reordered to match data" )
    } else {
      stop(
        "Column names in spectra and raw.data do not match perfectly;
        cannot reorder by name alone."
      )
    }
  }

  if ( ! all( table.fluors %in% fluorophores ) ) {
    # check for 'Negative', 'AF', check for match again
    fluor.to.match <- table.fluors[ !grepl( "Negative", table.fluors ) ]
    fluor.to.match <- fluor.to.match[ !fluor.to.match == "AF" ]
    matching.fluors <- fluor.to.match %in% fluorophores

    if ( !all( matching.fluors ) ) {
      warning( "The fluorophores in your control file don't match those in your spectra.",
               call. = FALSE )
      if ( !any( matching.fluors ) )
        stop( "No matching fluorophores between provided `spectra` and the control file.",
              call. = FALSE )
    }
    table.fluors <- fluor.to.match[ matching.fluors ]
  }

  # get thresholds for positivity
  if ( verbose ) message( paste0( "\033[32m", "Calculating positivity thresholds", "\033[0m" ) )
  unstained <- suppressWarnings(
    flowCore::read.FCS(
      file.path( control.dir, flow.file.name[ "AF" ] ),
      transformation = NULL,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )
  )

  # read exprs for spectral channels only
  if ( nrow( unstained ) > asp$gate.downsample.n.cells ) {
    set.seed( asp$gate.downsample.seed )
    unstained.idx <- sample( nrow( unstained ), asp$gate.downsample.n.beads )
    unstained <- flowCore::exprs( unstained )[ unstained.idx, spectral.channel ]
  } else {
    unstained <- flowCore::exprs( unstained )[ , spectral.channel ]
  }

  raw.thresholds <- apply( unstained, 2, function( col ) stats::quantile( col, 0.995 ) )

  unstained.unmixed <- unmix.autospectral(
    unstained,
    spectra,
    af.spectra,
    verbose = FALSE
  )
  unmixed.thresholds <- apply(
    unstained.unmixed[ , fluorophores ], 2, function( col )
      stats::quantile( col, 0.995 )
  )

  if ( is.null( names( table.fluors ) ) ) names( table.fluors ) <- table.fluors

  # main loop
  if ( parallel & is.null( threads ) ) threads <- asp$worker.process.n

  # construct list of arguments
  args.list <- list(
    file.name = flow.file.name,
    control.dir = control.dir,
    asp = asp,
    spectra = spectra,
    af.spectra = af.spectra,
    n.cells = n.cells,
    som.dim = som.dim,
    figures = figures,
    output.dir = output.dir,
    verbose = verbose,
    spectral.channel = spectral.channel,
    universal.negative = universal.negative,
    control.type = control.type,
    raw.thresholds = raw.thresholds,
    unmixed.thresholds = unmixed.thresholds,
    flow.channel = flow.channel,
    refine = refine,
    problem.quantile = problem.quantile
  )

  # Set up parallel processing
  if ( parallel ) {
    internal.functions <- c(
      "get.fluor.variants",
      "cosine.similarity",
      "spectral.variant.plot",
      "unmix.ols",
      "assign.af.residuals"
    )
    exports <- c( "args.list", "table.fluors", internal.functions )

    result <- create.parallel.lapply(
      asp,
      exports,
      parallel = parallel,
      threads = threads,
      export.env = environment()
    )
    lapply.function <- result$lapply

  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  if ( verbose ) message( paste0( "\033[34m", "Identifying spectral variation", "\033[0m" ) )

  # initialize with base spectra for safety
  spectral.variants <- lapply( table.fluors, function( fl ) {
    spectra[ fl, , drop = FALSE ]
  } )
  names( spectral.variants ) <- table.fluors

  # update with real variants where possible
  tryCatch(
    expr = {
      updated.variants <- lapply.function( table.fluors, function( f ) {
        tryCatch(
          expr = {
            do.call(
              get.fluor.variants,
              c( list( f ), args.list )
            )
          },
          error = function( e ) {
            message( "Error for fluorophore ", f, ": ", conditionMessage( e ) )
            NULL
          }
        )
      } )
      names( updated.variants ) <- table.fluors

      # replace base spectra only if updated.variants[[f]] is not NULL
      for ( f in table.fluors ) {
        if ( !is.null( updated.variants[[ f ]] ) ) {
          spectral.variants[[ f ]] <- updated.variants[[ f ]]
        } else {
          warning(
            paste( "Calculation failed for:", f, "- using base spectrum as fallback." )
          )
        }
      }
    },
    finally = {
      if ( !is.null( result$cleanup ) )
        result$cleanup()
    }
  )

  ### calculate deltas, delta.norms

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

  # precompute delta norms
  delta.norms <- lapply( delta.list, function( d ) {
    sqrt( rowSums( d^2 ) )
  } )
  names( delta.norms ) <- names( spectral.variants )

  if ( verbose )
    message( paste0( "\033[34m", "Spectral variation computed!", "\033[0m" ) )

  variants <- list(
    thresholds = unmixed.thresholds,
    variants = spectral.variants,
    delta.list = delta.list,
    delta.norms = delta.norms
  )

  saveRDS( variants, file = file.path( output.dir, asp$variant.filename ) )

  return( variants )

}
