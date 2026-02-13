# clean_controls.r

#' @title Clean Controls
#'
#' @description
#' A multi-part function to clean single-color controls in order to extract
#' fluorophore signatures. Any part can be run independently:
#'
#' - **Stage 1**: Autofluorescence noise removal using PCA unmixing on matching
#' unstained (cells only).
#' - **Stage 2**: Brightest event selection from positive, universal negative
#' from matching negative, and downsampling to speed up RLM spectra optimization.
#'
#' @importFrom lifecycle deprecate_warn
#'
#' @param flow.control A list prepared using `define.flow.control`, containing
#' the data and essential information about the cytometer and data structure.
#' @param asp The AutoSpectral parameter list, prepared using
#' `get.autospectral.param`.
#' @param time.clean Logical, default is `FALSE`. Whether to run PeacoQC to
#' remove time-based inconsistencies in the controls.
#' @param trim Logical, default is `FALSE`. Whether to remove extreme events
#' (positive and negative) from controls. Deprecated.
#' @param trim.factor Numeric. Default is `asp$rlm.trim.factor`. Required if
#' `trim = TRUE`. Deprecated.
#' @param af.remove Logical, default is `TRUE`. Whether to remove intrusive
#' autofluorescence contamination from cell controls using PCA-based
#' identification and gating. Requires universal negatives to be defined in the
#' control file and in `flow.control`.
#' @param universal.negative Logical, default is `TRUE`. Whether to use a
#' universal negative sample as the negative for spectral extraction. Requires
#' universal negatives to be defined in the control file and in `flow.control`.
#' @param downsample Logical, default is `TRUE`. Whether to reduce cell and bead
#' control events to speed up processing.
#' @param negative.n Integer. Number of events to include in the downsampled
#' negative population. Default is `asp$negative.n`.
#' @param positive.n Integer. Number of events to include in the downsampled
#' positive population. Default is `asp$positive.n`.
#' @param scatter.match Logical, default is `TRUE`. Whether to select negative
#' events based on scatter profiles matching the positive events. Defines a
#' region of FSC and SSC based on the distribution of selected positive events.
#' @param scrub Logical, if `TRUE` allows for re-cleaning of already cleaned
#' data, provided there are clean data in `flow.control`. Deprecated.
#' @param intermediate.figures Logical, if `TRUE` returns additional figures to
#' show the inner workings of the cleaning, including definition of low-AF cell
#' gates on the PCA-unmixed unstained and spectral ribbon plots of the AF
#' exclusion from the unstained. Default is `FALSE` to speed up processing.
#' @param main.figures Logical, if `TRUE` creates the main figures to show the
#' impact of intrusive autofluorescent event removal and scatter-matching for
#' the negatives.
#' @param parallel Logical, default is `FALSE`, in which case parallel processing
#' will not be used. Parallel processing will likely be faster when many small
#' files are read in. If the data is larger, parallel processing may not
#' accelerate the process much.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#' @param threads Numeric, number of threads to use for parallel processing.
#' Default is `NULL` which will revert to `asp$worker.process.n` if
#' `parallel=TRUE`.
#'
#' @return
#' Returns a modified `flow.control` with the original data intact. New, cleaned
#' data and corresponding factors are stored in new slots.
#'
#' @export

clean.controls <- function(
    flow.control,
    asp,
    time.clean = FALSE,
    trim = FALSE,
    trim.factor = NULL,
    af.remove = TRUE,
    universal.negative = TRUE,
    downsample = TRUE,
    negative.n = asp$negative.n,
    positive.n = asp$positive.n,
    scatter.match = TRUE,
    scrub = FALSE,
    intermediate.figures = FALSE,
    main.figures = TRUE,
    parallel = FALSE,
    verbose = TRUE,
    threads = NULL
  ) {

  if ( intermediate.figures & !main.figures ) main.figures <- TRUE

  if ( parallel & is.null( threads ) )
    threads <- asp$worker.process.n

  if ( parallel )
    threads <- max( floor( threads / 2 ), 1 )

  if ( parallel )
    warning( "Parallelization has not been fully optimized for `clean.controls`." )

  # warn user if deprecated arguments are called
  # these will be phased out later
  if ( time.clean ) {
    lifecycle::deprecate_warn(
      "0.9.1",
      "clean.controls(time.clean)",
      details = "will be deprecated going forward"
    )
  }
  if ( trim ) {
    lifecycle::deprecate_warn(
      "0.9.1",
      "clean.controls(trim)",
      details = "will be deprecated going forward"
    )
  }
  if ( !is.null( trim.factor ) ) {
    lifecycle::deprecate_warn(
      "0.9.1",
      "clean.controls(trim.factor)",
      details = "will be deprecated going forward"
    )
  }
  if ( scrub ) {
    lifecycle::deprecate_warn(
      "0.9.1",
      "clean.controls(scrub)",
      details = "will be deprecated going forward"
    )
  }

  flow.sample <- flow.control$sample
  flow.sample.n <- length( flow.sample )
  flow.control.type <- flow.control$control.type
  clean.event.type <- flow.control$event.type

  if ( scrub & !is.null( flow.control$clean.expr ) ) {

    clean.expr <- flow.control$clean.expr
    clean.event.sample <- flow.control$clean.event.sample
    if ( ! is.null( flow.control$clean.universal.negative ) ) {
      flow.negative <- flow.control$clean.universal.negative
      clean.universal.negative <- flow.control$clean.universal.negative
    } else {
      flow.negative <- flow.control$universal.negative
      clean.universal.negative <- NULL
    }

    asp$min.cell.warning.n <- asp$min.cell.warning.n / 2

  } else {

    clean.expr <- flow.control$expr.data
    clean.event.sample <- flow.control$event.sample
    flow.negative <- flow.control$universal.negative
    clean.universal.negative <- NULL
  }

  # split expression data by sample into a list
  clean.expr <- lapply( flow.sample, function( fs ) {
    clean.expr[ clean.event.sample == fs, ]
  } )

  names( clean.expr ) <- flow.sample
  spectral.channel <- flow.control$spectral.channel

  all.channels <- flow.control$scatter.and.channel.original

  ### Stage 1: Use PeacoQC to clean up flow  -----------------
  # clean based on irregularities in flow by TIME parameter
  # this can be slow and may not show much effect

  if ( time.clean ) {
    if ( !dir.exists( asp$figure.peacoqc.dir ) & main.figures )
      dir.create( asp$figure.peacoqc.dir )

    clean.expr <- run.peacoQC(
      expr.data = clean.expr,
      spectral.channel = spectral.channel,
      all.channels = all.channels,
      asp = asp,
      figures = main.figures,
      parallel = parallel,
      threads = threads,
      verbose = verbose
    )
  }

  ### Stage 2: Trimming -----------------
  # clean by trimming extreme events
  # to be used in case of aggregates
  # not recommended for regular use due to loss of brightest events

  if ( trim ) {

    if ( is.null( trim.factor ) ) {
      trim.factor <- asp$rlm.trim.factor
    }

    # trim fluorophore controls only
    trim.sample <- flow.control$fluorophore[
      !grepl( "AF|negative", flow.control$fluorophore, ignore.case = TRUE ) ]

    trim.peak.channels <- flow.control$channel[
      flow.control$fluorophore %in% trim.sample ]

    trim.sample.data <- clean.expr[ trim.sample ]

    trimmed.expr <- run.trim.events(
      trim.sample.data, trim.sample,
      trim.peak.channels, trim.factor, asp
    )

    rm( trim.sample.data )

    # merge in trimmed data
    clean.expr[ names( trimmed.expr ) ] <- trimmed.expr

  }

  ### Stage 3: Remove Autofluorescence intrusions -----------------
  # clean by removing autofluorescence contamination
  # recommended for controls from tissues (e.g., mouse splenocytes)
  # not needed for PBMCs

  if ( af.remove ) {
    if ( main.figures ) {
      if ( !dir.exists( asp$figure.clean.control.dir ) )
        dir.create( asp$figure.clean.control.dir )
      if ( !dir.exists( asp$figure.spectral.ribbon.dir ) )
        dir.create( asp$figure.spectral.ribbon.dir )
      if ( !dir.exists( asp$figure.scatter.dir.base ) )
        dir.create( asp$figure.scatter.dir.base )
    }

    # identify universal negative cell samples
    univ.neg <- unique( flow.negative )
    univ.neg <- univ.neg[ !is.null( univ.neg ) ]
    univ.neg <- univ.neg[ !is.na( univ.neg ) ]
    univ.neg <- univ.neg[ univ.neg != FALSE ]

    univ.neg.sample <- flow.control$sample[
      flow.control.type == "cells" & flow.sample %in% univ.neg ]
    # create new negative slot
    clean.universal.negative <- flow.negative
    flow.control$clean.universal.negative <- clean.universal.negative

    if ( ( length( univ.neg ) > 0 ) ) {

      # select cell-based single-stained samples to be used
      # must be cells and must have a corresponding universal negative
      af.removal.sample <- flow.sample[
        flow.control.type == "cells" & flow.negative %in% univ.neg ]

      # don't remove AF from negatives (to allow scatter matching)
      af.removal.sample <- af.removal.sample[
        !( af.removal.sample  %in% univ.neg.sample ) ]

      # check that we have samples to work with
      if ( length( af.removal.sample ) > 0 ) {
        af.remove.peak.channels <- flow.control$channel[
          flow.control$fluorophore %in% af.removal.sample ]

        # remove identified AF from single-color controls
        af.removed.expr <- run.af.removal(
          clean.expr = clean.expr,
          af.removal.sample = af.removal.sample,
          spectral.channel = spectral.channel,
          peak.channel = af.remove.peak.channels,
          universal.negative = flow.negative,
          asp = asp,
          scatter.param = flow.control$scatter.parameter,
          negative.n = negative.n,
          positive.n = positive.n,
          scatter.match = scatter.match,
          intermediate.figures = intermediate.figures,
          main.figures = main.figures,
          parallel = parallel,
          threads = threads,
          verbose = verbose
        )

        # store cleaned data
        cleaned.samples <- names( af.removed.expr )
        clean.expr[ cleaned.samples ] <- af.removed.expr[ cleaned.samples ]

      } else {
        warning(
          "No cell-based universal negative samples could be identified for `af.remove`."
          )
        message(
          "\033[31m",
          paste(
            "No cell-based universal negative samples could be identified.",
            "To perform autofluorescence removal, you must specify a universal negative in the fcs_control_file,",
            "and you must have single-stained cell controls.",
            "If you only have bead-based controls, set `af.remove` to FALSE and try again.",
            "Skipping autofluorescence removal.",
            sep = "\n"
          ),
          "\033[0m"
        )
      }
    } else {
      warning(
      "No cell-based universal negative samples could be identified.
      To perform autofluorescence removal, you must have single-stained cells,
      and specify a universal negative in the fcs_control_file."
      )
      message(
        paste0(
          "\033[31m",
          "No cell-based universal negative samples could be identified.",
          "\n",
          "Skipping autofluorescence removal.",
          "\033[0m"
        )
      )
    }
  }

  ### Stage 4: Universal Negatives and Downsampling -----------------
  # clean by selecting universal negative (AF-removed if af.remove performed)
  # also selects brightest positive events
  # and downsamples to speed up subsequent calculations
  # recommended whenever you have an appropriate universal negative

  if ( universal.negative ) {
    if ( main.figures ) {
      if ( !dir.exists( asp$figure.scatter.dir.base ) )
        dir.create( asp$figure.scatter.dir.base )
      if ( !dir.exists( asp$figure.spectral.ribbon.dir ) )
        dir.create( asp$figure.spectral.ribbon.dir )
    }

    # use AF-removed universal negative samples if available
    if ( !is.null( clean.universal.negative ) )
      flow.negative <- clean.universal.negative
    else
      flow.negative

    # select fluorophore samples to be used
    univ.sample <- flow.control$fluorophore[
      !grepl(
        "negative",
        flow.control$fluorophore,
        ignore.case = TRUE
      ) & flow.negative != FALSE ]

    # exclude AF
    univ.sample <- univ.sample[ univ.sample != "AF" ]

    # if AF removal has been done with scatter-matching, run only on beads
    if ( af.remove ) {
      univ.sample.type <- flow.control.type[
        names( flow.control.type ) %in% univ.sample ]
      univ.sample <- univ.sample[ univ.sample.type == "beads" ]
    }

    # check for remaining samples
    if ( length( univ.sample ) != 0 ) {

      univ.peak.channels <- flow.control$channel[
        flow.control$fluorophore %in% univ.sample ]

      univ.neg.expr <- run.universal.negative(
        clean.expr = clean.expr,
        univ.sample = univ.sample,
        universal.negatives = flow.negative,
        scatter.param = flow.control$scatter.parameter,
        peak.channels = univ.peak.channels,
        downsample = downsample,
        negative.n = negative.n,
        positive.n = positive.n,
        spectral.channel = spectral.channel,
        asp = asp,
        control.type = flow.control.type,
        scatter.match = scatter.match,
        intermediate.figures = intermediate.figures,
        main.figures = main.figures,
        verbose = verbose
      )

      # merge in cleaned data
      clean.expr[ names( univ.neg.expr ) ] <- univ.neg.expr
    }

    # simply downsample AF to speed up calculations
    af.expr <- clean.expr[[ "AF" ]]

    if ( "AF" %in% flow.sample & nrow( af.expr ) > negative.n ) {
      set.seed( asp$gate.downsample.seed )
      downsample.idx <- sample( 1:nrow( af.expr ), negative.n )
      clean.expr[[ "AF" ]] <- af.expr[ downsample.idx, ]
    }
  }

  # downsample only
  if ( !universal.negative & downsample ) {
    if ( main.figures ) {
      if ( !dir.exists( asp$figure.scatter.dir.base ) )
        dir.create( asp$figure.scatter.dir.base )
      if ( !dir.exists( asp$figure.spectral.ribbon.dir ) )
        dir.create( asp$figure.spectral.ribbon.dir )
    }

    # select fluorophore samples to be used
    downsample.sample <- flow.control$fluorophore[
      !grepl( "AF|negative", flow.control$fluorophore, ignore.case = TRUE ) ]

    downsample.peak.channels <- flow.control$channel[
      flow.control$fluorophore %in% downsample.sample ]

    downsample.expr <- run.downsample(
      clean.expr.data = clean.expr,
      downsample.sample = downsample.sample,
      peak.channels = downsample.peak.channels,
      negative.n = negative.n,
      positive.n = positive.n,
      verbose = verbose
    )

    # merge in cleaned data
    clean.expr[ names( downsample.expr ) ] <- downsample.expr
  }

  # merge data and re-establish corresponding factors
  names( clean.expr ) <- flow.sample

  # get maximum number of events per sample to adjust event numbering
  flow.sample.event.number.max <- 0

  for ( fs.idx in 1 : flow.sample.n )
  {
    flow.sample.event.number <- nrow( clean.expr[[ fs.idx ]]  )

    rownames( clean.expr[[ fs.idx ]] ) <- paste(
      flow.sample[ fs.idx ],
      seq_len( flow.sample.event.number ),
      sep = "_"
    )

    if ( flow.sample.event.number < asp$min.cell.warning.n ) {
      warning(
        paste0(
          "\033[31m",
          "Fewer than ",
          asp$min.cell.warning.n,
          " gated events in ",
          names( clean.expr )[ fs.idx ],
          "\033[0m",
          "\n"
        )
      )
    }

    if ( flow.sample.event.number > flow.sample.event.number.max )
      flow.sample.event.number.max <- flow.sample.event.number
  }

  flow.event.number.width <-
    floor( log10( flow.sample.event.number.max ) ) + 1
  flow.event.regexp <- sprintf( "\\.[0-9]{%d}$", flow.event.number.width )

  # set rownames
  for ( fs.idx in 1 : flow.sample.n )
  {
    flow.sample.event.number <- nrow( clean.expr[[ fs.idx ]]  )
    flow.the.sample <- flow.sample[ fs.idx ]

    flow.the.event <- sprintf(
      "%s.%0*d", flow.the.sample,
      flow.event.number.width, 1 : flow.sample.event.number
    )
    rownames( clean.expr[[ fs.idx ]] ) <- flow.the.event
  }

  clean.expr <- do.call( rbind, clean.expr )

  # set events
  flow.event <- rownames( clean.expr )
  flow.event.n <- length( flow.event )

  flow.event.sample <- sub( flow.event.regexp, "", flow.event )
  flow.event.sample <- factor( flow.event.sample, levels = flow.sample )

  event.type.factor <- flow.sample
  names( event.type.factor ) <- flow.control$control.type
  flow.event.type <- factor(
    flow.event.sample,
    levels = event.type.factor,
    labels = names( event.type.factor )
    )

  # store in flow.control
  flow.control$clean.expr <- clean.expr
  flow.control$clean.event.sample <- flow.event.sample
  flow.control$clean.event.type <- flow.event.type

  return( flow.control )
}
