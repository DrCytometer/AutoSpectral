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
#' @param som.dim Numeric, default `7`. Number of x and y dimensions to use in
#' the SOM for clustering the spectral variation. The number of spectra returned
#' for each fluorophore will increase with the quadratic of `som.dim`, so for 7,
#' you will get up to 49 variants. Increasing the SOM dimensions further does
#' not help. Somewhere between 4 and 7 appears to be sufficient, but with the
#' pruning of variants implemented in `unmix.autospectral()` in v1.0.0, this is
#' less important.
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
#' @param ... Ignored. Previously used for deprecated arguments such as
#' `pos.quantile` and `sim.threshold`, which are now fixed internally and no
#' longer user-settable.
#'
#' @return A vector with the indexes of events inside the initial gate.
#'
#' @export

get.spectral.variants <- function(
    control.dir, control.def.file,
    asp, spectra, af.spectra,
    n.cells = 2000,
    som.dim = 7,
    figures = TRUE,
    output.dir = NULL,
    parallel = FALSE,
    verbose = TRUE,
    threads = NULL,
    ...
) {

  dots <- list( ... )

  if ( !is.null( dots$sim.threshold ) ) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "get.spectral.variants(sim.threshold)",
      "no longer used"
    )
  }
  if ( !is.null( dots$pos.quantile ) ) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "get.spectral.variants(pos.quantile)",
      "no longer used"
    )
  }

  if ( is.null( af.spectra ) )
    stop( "Multiple AF spectra must be provided." )
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided." )

  if ( is.null( output.dir ) )
    output.dir <- asp$variant.dir
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  if ( som.dim > 12 )
    warning( paste(
      "Argument `som.dim` has been set to", som.dim, "which will produce",
      som.dim^2, "spectral variants per fluorophore.", "\n",
      "More spectral variants means slower unmixing and will also require",
      "proprotionally more cells in `n.cells` as input."
    ) )

  fluorophores <- rownames( spectra )[ rownames( spectra ) != "AF" ]
  spectral.channel <- colnames( spectra )

  # read control file
  if ( !file.exists( control.def.file ) )
    stop( paste( "Unable to locate control.def.file:", control.def.file ) )

  if ( verbose ) message( "\033[34mChecking control file for errors \033[0m" )
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
  flow.set.resolution <- asp$expr.data.max
  flow.set.channel.table <- read.channel( control.dir, control.def.file, asp )
  flow.set.channel <- flow.set.channel.table[[ 1 ]]
  non.spectral.channel <- asp$non.spectral.channel
  non.spectral.channel <- paste0( non.spectral.channel, collapse = "|" )

  flow.spectral.channel <- flow.set.channel[ !grepl( non.spectral.channel,
                                                     flow.set.channel ) ]
  flow.scatter.parameter <- read.scatter.parameter( asp )
  flow.scatter.and.channel.spectral <- c( asp$default.time.parameter,
                                          flow.scatter.parameter,
                                          flow.spectral.channel )

  if ( grepl( "Discover", asp$cytometer ) )
    flow.spectral.channel <- flow.spectral.channel[ grep( asp$spectral.channel,
                                                          flow.spectral.channel ) ]

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
  if ( !( "AF" %in% table.fluors ) )
    stop( "Unable to locate `AF` control in control file. An unstained cell control is required" )

  if ( ! all( table.fluors %in% fluorophores ) ) {
    # check for 'Negative', 'AF', check for match again
    fluor.to.match <- table.fluors[ !grepl( "Negative", table.fluors ) ]
    fluor.to.match <- fluor.to.match[ !fluor.to.match == "AF" ]
    matching.fluors <- fluor.to.match %in% fluorophores

    if ( !all( matching.fluors ) ) {
      warning( "The fluorophores in your control file don't match those in your spectra." )
      if ( !any( matching.fluors ) )
        stop()
    }
    table.fluors <- fluor.to.match[ matching.fluors ]
  }

  # get thresholds for positivity
  if ( verbose ) message( paste0( "\033[33m", "Calculating positivity thresholds", "\033[0m" ) )
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
    spectral.channel = flow.spectral.channel,
    universal.negative = universal.negative,
    control.type = control.type,
    raw.thresholds = raw.thresholds,
    unmixed.thresholds = unmixed.thresholds,
    flow.channel = flow.channel
  )

  # Set up parallel processing
  if ( parallel ) {
    internal.functions <- c( "get.fluor.variants", "cosine.similarity",
                             "spectral.variant.plot", "unmix.ols" )
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

  if ( verbose ) message( paste0( "\033[33m", "Getting spectral variants", "\033[0m" ) )

  spectral.variants <- tryCatch(
    expr = {
      lapply.function( table.fluors, function( f ) {
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
    },
    finally = {
      if ( !is.null( result$cleanup ) )
        result$cleanup()
    }
  )

  # drop any that are NULL due to error
  failed <- sapply( spectral.variants, is.null )
  failed.variants <- names( spectral.variants )[ failed ]
  spectral.variants <- spectral.variants[ !failed ]

  if ( any( failed ) ) {
    warning(
      paste0(
        "Calculation of spectral variation failed for ",
        paste( failed.variants, collapse = "\n" )
      )
    )
  }

  if ( verbose ) message( paste0( "\033[33m", "Spectral variation computed!", "\033[0m" ) )

  variants <- list(
    thresholds = unmixed.thresholds,
    variants = spectral.variants
  )

  saveRDS( variants, file = file.path( output.dir, asp$variant.filename ) )

  return( variants )

}
