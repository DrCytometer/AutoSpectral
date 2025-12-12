# check_control_file.r

#' @title Check Control File
#' @description
#' Attempts to find potential failure points and input errors in the control file
#' `control.def.file` used to define the single-stained control setup for
#' AutoSpectral.
#'
#' @importFrom utils read.csv
#' @importFrom flowCore read.FCSheader
#' @importFrom dplyr filter
#'
#' @param control.dir File path to the single-stained control FCS files.
#' @param control.def.file CSV file defining the single-color control file names,
#' fluorophores they represent, marker names, peak channels, and gating requirements.
#' @param asp The AutoSpectral parameter list defined
#' using `get.autospectral.param`.
#' @param strict Logical. Controls whether the function triggers a break or
#' continues and outputs a list of errors. Default is `FALSE`.
#'
#' @return A named list of errors and warnings intended to help the user fix
#' problems with the `control.def.file`.
#'
#' @export

check.control.file <- function( control.dir, control.def.file, asp,
                                strict = FALSE ) {

  errors <- list()
  warning.list <- list()

  # check that control file exists
  if ( !file.exists( control.def.file ) )
    stop( paste( "Unable to locate control.def.file:", control.def.file ) )

  # check that control.dir exists
  if ( !dir.exists( control.dir ) )
    stop( paste( "Unable to locate control.dir:", control.dir ) )

  # check that fcs files exist
  control.files <- list.files( path = control.dir, pattern = ".fcs",
                               ignore.case = TRUE )

  if ( length( control.files ) < 1 )
    stop( paste( "Unable to find FCS files in the control.dir:", control.dir ) )

  # read control file
  control.table <- read.csv( control.def.file, na.strings = "",
                             stringsAsFactors = FALSE )

  # check that columns are all present
  expected.colnames <- c( "filename", "fluorophore", "marker", "channel",
                         "control.type", "universal.negative", "large.gate" )
  actual.colnames <- colnames( control.table )

  colnames.check <- expected.colnames %in% actual.colnames

  if ( !all( colnames.check ) ) {
    missing.columns <- expected.colnames[ !colnames.check ]
    errors$missing.columns <- missing.columns
    warning( paste( "\033[31mRequired column:", missing.columns, "was not found in the
                 control.def.file.\033[0m", collapse = "\n" ) )

    if ( !strict )
      return( errors )
  }

  # check for missing data in filename column
  missing.filenames <- dplyr::filter( control.table, filename == "" )

  if ( nrow( missing.filenames ) > 0 ) {
    warning( "Some rows are missing information for filename." )
    message( "\033[32mSome rows are missing information for filename:\033[0m" )
    message( paste( missing.filenames, collapse = "\n" ) )
    warning.list$missing.filenames <- missing.filenames
    control.table <- dplyr::filter( control.table, filename != "" )
  }

  # check that it has rows of data > 1
  if ( nrow( control.table ) < 2 ) {
    warning( "Control table only has one row. Please add more information." )
    message( "\033[31mControl table only has one row. Please add more information.\033[0m" )
    errors$low.sample.number <- "Control table only has one row. Please add more information."

    if ( !strict )
      return( errors )
  }

  if ( nrow( control.table ) < 3 ) {
    message( "\033[32mControl table only has two rows. Please check that it is complete.\033[0m" )
    warning.list$low.sample.number <- paste( "Only", nrow( control.table ), "samples provided" )
  }

  # check for duplicated filenames
  duplicated.file.check <- duplicated( control.table$filename )

  if ( all( duplicated.file.check ) ) {
    duplicate.filenames <- control.table$filename[ duplicated.file.check ]
    warning( "One or more filenames are duplicated." )
    message( "\033[31mOne or more filenames are duplicated:\033[0m" )
    message( paste( duplicate.filenames, collapse = "\n" ) )
    errors$duplicate.filenames <- duplicate.filenames
  }

  # check that file names contain .fcs
  filename.check <- grepl( ".fcs", control.table$filename, ignore.case = TRUE )

  if ( !all( filename.check ) ) {
    non.fcs.filenames <- control.table$filename[ !filename.check ]
    warning( "One or more filenames are not .fcs files." )
    message( "\033[31mOne or more filenames are not .fcs files:\033[0m" )
    message( paste( non.fcs.filenames, collapse = "\n" ) )
    errors$non.fcs.filenames <- non.fcs.filenames

    if ( !strict )
      return( errors )
  }

  # check that fluorophore has been filled in
  missing.fluorophore.1 <- is.na( control.table$fluorophore )
  missing.fluorophore.2 <- control.table$fluorophore[ control.table$fluorophore == "" ]

  if ( any( missing.fluorophore.1 ) || length( missing.fluorophore.2 ) > 0 ) {
    missing.fluorophore <- control.table$filename[ missing.fluorophore.1 ]
    warning( "Fluorophores have not been filled in for some controls." )
    message( "\033[31mFluorophores have not been filled in for the following controls:\033[0m" )
    message( paste( missing.fluorophore, collapse = "\n" ) )
    errors$missing.fluorophore <- missing.fluorophore
  }

  # check for duplicated fluorophore names
  duplicated.fluorophores <- duplicated( control.table$fluorophore )
  if ( any( duplicated.fluorophores ) ) {
    dup.fluor <- control.table$fluorophore[ duplicated.fluorophores ]
    warning( "Duplicate entries exist for fluorophores." )
    message( "\033[31mThe following fluorophores are duplicated:\033[0m" )
    message( paste( dup.fluor, collapse = "\n" ) )
    message( "\033[31mTo run multiple controls per fluorophore, name them uniquely, e.g.,
             FITC_1, FITC_2, etc.\033[0m" )
    errors$duplicated.fluorophore <- dup.fluor
  }

  # check that an unstained cell sample has been provided and is named "AF"
  if ( ! any( control.table$fluorophore == "AF" ) ) {
    message( "\033[32mNo Autofluorescence `AF` control sample has been provided.\033[0m" )
    errors$missing.af <- "No `AF` sample listed"
  } else {
    # check that the sample has been marked as `cells`
    af.control <- control.table$filename[ control.table$fluorophore == "AF" ]
    if ( length( af.control ) > 1 ) {
      message( "\033[32mDuplicate `AF` control samples have been provided. Provide unique names.\033[0m" )
      errors$duplicate.af <- af.control
    }
    if ( ! control.table$control.type[ control.table$fluorophore == "AF" ] == "cells" ) {
      message( "\033[32mThe `AF` control must be `cells` in the `control.type` column.\033[0m" )
      errors$wrong.af.type <- af.control
    }
  }

  # check that marker has been filled in
  missing.marker.1 <- is.na( control.table$marker )
  missing.marker.2 <- control.table$fluorophore[ control.table$marker == "" ]

  acceptable.missing <- c( "AF", ".*negative.*" )

  if ( any( missing.marker.1 ) || length( missing.marker.2 ) > 0 ) {
    missing.marker <- control.table$filename[ missing.marker.1 ]
    missing.marker.fluor <- control.table$fluorophore[ missing.marker.1 ]

    #problem.missing.marker <- !sapply( missing.marker.fluor, function( x )
    #  any( x == acceptable.missing ) )
    problem.missing.marker <- !sapply( missing.marker.fluor, function( x )
      any( grepl( paste( acceptable.missing, collapse = "|" ), x, ignore.case = TRUE ) ) )
    problem.missing.marker <- names( problem.missing.marker )[ problem.missing.marker ]
    problem.missing.marker <- control.table$filename[ control.table$fluorophore %in% problem.missing.marker ]

    if ( length( problem.missing.marker ) > 0 ) {
      message( "\033[32mMarkers have not been filled in for the following controls:\033[0m" )
      message( paste( problem.missing.marker, collapse = "\n" ) )
      message( paste( "\033[32mNote that markers are not required for unstained controls,",
                      "negative controls or the autofluorescence control. Markers can be",
                      "left out entirely, but your FCS files will not be labeled.\033[0m",
                      "\n", sep = "\n" ) )
      warning.list$missing.marker <- missing.marker
    }
  }

  # check that peak channels have been filled in
  missing.channel.1 <- is.na( control.table$channel )
  missing.channel.2 <- control.table$fluorophore[ control.table$channel == "" ]

  if ( any( missing.channel.1 ) || length( missing.channel.2 ) > 0 ) {
    missing.channel <- control.table$filename[ missing.channel.1 ]
    missing.channel.fluor <- control.table$fluorophore[ missing.channel.1 ]

    problem.missing.channel <- !sapply( missing.channel.fluor, function( x )
      any( grepl( paste( acceptable.missing, collapse = "|" ), x, ignore.case = TRUE ) ) )
    problem.missing.channel <- names( problem.missing.channel )[ problem.missing.channel ]
    problem.missing.channel <- control.table$filename[ control.table$fluorophore %in% problem.missing.channel ]

    if ( length( problem.missing.channel ) > 0 ) {
      warning( "Channels have not been filled in for some controls." )
      message( "\033[31mChannels have not been filled in for the following controls:\033[0m" )
      message( paste( problem.missing.channel, collapse = "\n" ) )

      errors$missing.channel <- problem.missing.channel
    }
  }

  # check that control type has been filled in
  missing.type.1 <- is.na( control.table$control.type )
  missing.type.2 <- control.table$fluorophore[ control.table$control.type == "" ]

  if ( any( missing.type.1 ) || length( missing.type.2 ) > 0 ) {
    missing.type <- control.table$filename[ missing.type.1 ]
    warning( "Control type has not been filled in for the some controls" )
    message( "\033[31mControl type has not been filled in for the following controls:\033[0m" )
    message( paste( missing.type, collapse = "\n" ) )
    message( "\033[31mOptions are `beads` or `cells`.\033[0m" )
    errors$missing.type <- missing.type
  }

  # check that it matches "beads" or "cells"
  wrong.type <- !grepl( "beads|cells", control.table$control.type )

  if ( any( wrong.type ) ) {
    wrong.type <- control.table$filename[ wrong.type ]
    warning( "Control type has been filled in incorrectly for some controls" )
    message( "\033[31mControl type has been filled in incorrectly for the following controls:\033[0m" )
    message( paste( wrong.type, collapse = "\n" ) )
    message( "\033[31mControl type must be `beads` or `cells`.\033[0m" )
    errors$wrong.type <- wrong.type
  }

  # check that "No Match" has been replaced with real data
  no.match <- grepl( "No match", control.table$fluorophore, ignore.case = TRUE )

  if ( any( no.match ) ) {
    no.match <- control.table$filename[ no.match ]
    warning( "Your control table still lists `No match` for some files" )
    message( "\033[31mYour control table still lists `No match` for these files:\033[0m" )
    message( paste( no.match, collapse = "\n" ) )
    message( paste( "`\033[31mNo match` means AutoSpectral couldn't find a match for the",
             "fluorophore from the FCS file's name. Fill it in manually.\033[0m",
             sep = "\n" ) )
    errors$no.match <- no.match
  }

  # check that universal negative have been provided
  missing.univ.neg.1 <- is.na( control.table$universal.negative )
  missing.univ.neg.2 <- control.table$universal.negative[ control.table$universal.negative == "" ]

  if ( any( missing.univ.neg.1 ) || length( missing.univ.neg.2 ) > 0 ) {
    missing.univ.neg <- control.table$filename[ missing.univ.neg.1 ]
    message( "\033[32mUniversal negatives have not been provided for these files:\033[0m" )
    message( paste( missing.univ.neg, collapse = "\n" ) )
    message( "\033[32mUniversal negatives are not required, but are recommended.\033[0m" )
    warning.list$lacking.univ.neg <- missing.univ.neg
  }

  if ( any( !missing.univ.neg.1 ) ) {

    # check that universal negative names match a file in the file list
    univ.neg <- unique( control.table$universal.negative )
    univ.neg.samples <- sapply( control.table$filename, function( x ) any( x == univ.neg ) )
    univ.neg.present <- univ.neg %in% names( univ.neg.samples[ univ.neg.samples == TRUE ] )

    if ( !all( univ.neg.present ) ) {
      missing.univ.neg <- univ.neg[ !univ.neg.present ]
      warning( "One or more files listed as `universal.negative` were not
               included as samples under `filename`" )
      message( paste( "\033[31mThe following are listed as universal negatives but are",
                      "not listed as samples under filename:\033[0m", sep = "\n" ) )
      message( paste( missing.univ.neg, collapse = "\n" ) )
      message( paste( "\033[31mAll samples used as universal negatives must be included",
               "as a separate line in the control file. Fluorophore should be listed",
               "as `Negative`.\033[0m", sep = "\n" ) )
      errors$missing.univ.neg <- missing.univ.neg

      missing.univ.neg.grep <- paste( missing.univ.neg, collapse = "|" )
      matching.samples <- control.table$filename[ grepl( missing.univ.neg.grep,
                                                         control.table$universal.negative ) ]
      message( paste( "\033[32mIn case this is a typo, here are the samples listing those as",
               "universal negatives:\033[0m", sep = "\n" ) )
      message( paste( matching.samples, collapse = "\n" ) )
      errors$matching.samples <- matching.samples
    }

    # check that universal negative file is named Negative or AF
    negative.names <- "AF|Negative*"
    universal.negative.fluor <- control.table$fluorophore[ control.table$filename %in% univ.neg ]
    neg.name.match <- grepl( negative.names, universal.negative.fluor )

    if ( any( !neg.name.match ) ) {
      non.matching.neg.names <- univ.neg[ !neg.name.match ]
      warning( "Universal negative samples should have either `AF` or `Negative` as the fluorophore.
               Use `AF` for the unstained cell control and `Negative` for unstained beads. If you have
               multiple unstained samples, use `AF` for the first (or most releveant) cell-based unstained,
               then `Negative`, `Negative1`, `Negative2`, etc." )
      message( paste( "\033[31mUniversal negative samples should have either `AF` or `Negativeas the fluorophore.",
               "Use `AF` for the unstained cell control and `Negative` for unstained beads.",
               "If you have multiple unstained samples, use `AF` for the first (or most releveant)
               cell-based unstained, then `Negative`, `Negative1`, `Negative2`, etc.",
               "These samples are listed as universal negatives but have a different fluorophore indicated:\033[0m",
               sep = "\n" ) )
      message( paste( non.matching.neg.names, collapse = "\n" ) )
      errors$non.matching.neg.names <- non.matching.neg.names
    }

    # check that universal negative type matches sample type (beads vs cells)
    samples.with.negatives <- control.table$filename[ !is.na( control.table$universal.negative ) ]
    univ.neg.type.match <- sapply( samples.with.negatives, function( samp ) {
      sample.type <- control.table$control.type[ control.table$filename == samp ]
      sample.neg <- control.table$universal.negative[ control.table$filename == samp ]
      neg.sample.type <- control.table$control.type[ control.table$filename == sample.neg ]
      sample.type == neg.sample.type
    } )

    if ( !all( univ.neg.type.match ) ) {
      mismatched.negatives <- names( univ.neg.type.match )[ !univ.neg.type.match ]
      warning( paste( "Mismatch detected for sample type between the sample and",
                      "the negative for these samples.", sep = "\n" ) )
      message( paste( "\033[31mMismatch detected for sample type between the sample and",
                      "the negative for these samples\033[0m.", sep = "\n" ) )
      message( paste( mismatched.negatives, collapse = "\n" ) )
      message( paste( "\033[31mPlease ensure that `beads` samples are matched with a bead negative",
               "and `cells` samples with a cell negative.\033[0m", "\n", sep = "\n" ) )
      errors$mismatched.negatives <- mismatched.negatives
    }
  }

  # check for is.viability and large.gate presence
  missing.large.gate <- is.na( control.table$large.gate )
  missing.viability <- is.na( control.table$is.viability )

  if ( all( missing.large.gate ) ) {
    message( paste( "\033[32mThe `large.gate` option has not been set for any sample. Check",
             "that this is appropriate. Larger gates are appropriate if the cells",
             "bearing the marker are bigger than the main population, e.g., monocytes.\033[0m",
             "\n", sep = "\n" ) )
    warning.list$missing.large.gate <- "No large.gate set"
  }

  if ( all( missing.viability ) ) {
    message( paste( "\033[32mThe `is.viability` option has not been set for any sample.",
                    "Check that this is appropriate. If you have a viability dye in",
                    "your panel, set `is.viability` to `TRUE` for that sample.\033[0m",
                    "\n", sep = "\n" )  )
    warning.list$missing.viability <- "No viability marker set"
  }

  # check that is.viability and large.gate are either absent, TRUE or FALSE
  viability.input <- control.table$is.viability
  viability.input.valid.idx <- !is.na( viability.input )
  large.gate.input <- control.table$large.gate
  acceptable.input <- c( TRUE, FALSE )

  # Logical validity: TRUE/FALSE/NA allowed, anything else triggers error
  valid.logical <- is.na( large.gate.input ) | large.gate.input %in% acceptable.input
  non.standard.idx <- !valid.logical

  # check for non-standard entries
  if ( any( non.standard.idx ) ) {
    warning( "Only TRUE and FALSE are acceptable options for `large.gate`." )
    message( paste(
      "\033[31mOnly `TRUE` and `FALSE` are acceptable options for `large.gate`.",
      "Enter TRUE, FALSE or leave the field blank (equivalent to FALSE).\033[0m",
      sep = "\n"
    ) )
    message( "Problematic entries:" )

    # Return filenames of problematic entries
    large.gate.error <- control.table$filename[ non.standard.idx ]
    message( paste( large.gate.error, collapse = "\n" ) )
    errors$large.gate.error <- large.gate.error
  }

  # Logical validity: only TRUE, FALSE, or NA allowed
  valid.logical <- is.na( viability.input ) | ( viability.input %in% acceptable.input )
  non.standard.idx <- !valid.logical

  if ( any( non.standard.idx ) ) {
    warning( "Only TRUE and FALSE are acceptable options for `is.viability`." )
    message( paste(
      "\033[31mOnly `TRUE` and `FALSE` are acceptable options for `is.viability`.",
      "Enter TRUE, FALSE or leave the field blank (equivalent to FALSE).\033[0m",
      sep = "\n"
    ) )
    message( "Problematic entries:" )
    viability.error <- control.table$filename[ non.standard.idx ]
    message( paste( viability.error, collapse = "\n" ) )
    errors$viability.error <- viability.error
  }

  # check for misallocation of `large.gate` or `is.viability` to bead controls
  # For entries marked TRUE, ensure they correspond to cell-based controls
  true.large.gate.idx <- which( !is.na( large.gate.input ) & large.gate.input == TRUE )

  if ( length( true.large.gate.idx ) > 0 ) {
    # TRUE is allowed only when control.type == "cells"
    invalid.true.idx <- true.large.gate.idx[
      control.table$control.type[ true.large.gate.idx ] != "cells" ]

    if ( length( invalid.true.idx ) > 0) {
      warning( "Only cell-based controls should be marked as TRUE for `large.gate`." )
      message( paste(
        "\033[31mOnly cell-based controls should be marked as `TRUE` for `large.gate`.",
        "For bead-based controls, leave `large.gate` blank or set it to FALSE.\033[0m",
        sep = "\n"
      ) )
      message( "Problematic entries:" )
      large.gate.error2 <- control.table$filename[ invalid.true.idx ]
      message( paste( large.gate.error2, collapse = "\n" ) )
      errors$large.gate.error2 <- large.gate.error2
    }
  }

  # For entries that are TRUE, ensure they correspond to cell-based controls
  true.viability.idx <- which( !is.na( viability.input ) & viability.input == TRUE )

  if ( length( true.viability.idx ) > 0) {
    invalid.true.idx <- true.viability.idx[
      control.table$control.type[ true.viability.idx ] != "cells" ]

    if ( length( invalid.true.idx ) > 0 ) {
      warning( "Only cell-based controls should be marked as TRUE for `is.viability`." )
      message( paste(
        "\033[31mOnly cell-based controls should be marked as `TRUE` for `is.viability`.",
        "For bead-based viability controls (ArC beads), leave `is.viability` blank or FALSE.\033[0m",
        sep = "\n"
      ) )
      message( "Problematic entries:")
      viability.gate.error <- control.table$filename[ invalid.true.idx ]
      message( paste( viability.gate.error, collapse = "\n" ) )
      errors$viability.gate.error <- viability.gate.error
    }
  }

  # check FCS files for consistency
  header.check <- lapply( control.table$filename, function( fcs ) {

    fcs.header <- suppressWarnings( flowCore::read.FCSheader( fcs, path = control.dir ) )
    header.text <- fcs.header[[ 1 ]]
    n.events <- as.integer( header.text[ "$TOT" ] )
    param.n <- as.integer( header.text[ "$PAR" ] )
    param.names <- sapply( 1:param.n, function( i )
      header.text[[ paste0( "$P", i, "N" ) ]] )

    cyt <- header.text[ "$CYT" ]

    list( parameters = param.names,
          cytometer = cyt,
          events = n.events )
  } )

  # check for low event counts
  n.events <- sapply( header.check, function( x ) x$events )
  low.event.idx <- which( n.events < 5000 )

  if ( length( low.event.idx ) > 0 ) {
    low.count.files <- control.table$filename[ low.event.idx ]
    message( paste( "\033[31mLow event count (<5000) detected in one or more files.",
                    "Concerning files:\033[0m", sep = "\n" ) )
    message( paste( low.count.files, collapse = "\n" ) )
    warning.list$low.event.count <- low.count.files
  }

  very.low.event.idx <- which( n.events < 1000 )

  if ( length( very.low.event.idx ) > 0 ) {
    low.count.files <- control.table$filename[ very.low.event.idx ]
    warning( "Some FCS files contain fewer than 1000 events:" )
    message( paste( "\033[31mLow event count (<1000) detected in one or more files.",
                    "Problematic files:\033[0m", sep = "\n" ) )
    message( paste( low.count.files, collapse = "\n" ) )
    errors$very.low.event.count <- low.count.files
  }

  param.names.list <- lapply( header.check, `[[`, "parameters" )
  cytometer.list <- sapply( header.check, `[[`, "cytometer" )

  if ( asp$cytometer == "FACSDiscover S8" | asp$cytometer == "FACSDiscover A8" ) {
    remove.non.spectral <- function( params, patterns ) {
      params[ !sapply( params, function( p )
        any( grepl( paste( patterns, collapse = "|" ), p ) ) ) ]
    }

    non.spectral <- asp$non.spectral.channel[ 4:length( asp$non.spectral.channel ) ]

    param.names.list <- lapply( param.names.list, remove.non.spectral,
                                        patterns = non.spectral )
  }

  matching.parameters <- all( sapply( param.names.list[ -1 ], function( x )
    identical( x, param.names.list[[ 1 ]] ) ) )

  if ( !matching.parameters ) {
    warning( paste( "Inconsistencies found in parameter names across FCS files.",
                    "Please inspect or this will likely fail.", sep = "\n" ) )
    message( paste( "\033[31mInconsistencies found in parameter names across FCS files.",
                    "Please inspect or this will likely fail.\033[0m", sep = "\n" ) )
    errors$parameter.mismatch <- param.names.list[
        sapply( param.names.list, function( x ) !identical( x, param.names.list[[ 1 ]] ) )
      ]
  }

  matching.cytometers <- length( unique( cytometer.list ) ) == 1

  if ( !matching.cytometers ) {
    warning( "More than one cytometer name has been detected among the FCS files." )
    message( "\033[31mMore than one cytometer name has been detected among the FCS files.\033[0m" )
    message( "\033[31mSomething is very wrong.\033[0m" )
    errors$matching.cytometers <- unique( cytometer.list )

  } else {
    cytometer.match <- unique( cytometer.list ) == asp$cytometer

    if ( !cytometer.match ) {
      cytometer.match <- grepl( asp$cytometer, unique( cytometer.list ), ignore.case = TRUE )
      if ( !cytometer.match ) {
        warning( paste( "The name of the cytomter in your FCS files does not match the",
                        "name in the `asp` parameter list. Please ensure you have selected",
                        "the appropriate option when calling `get.autospectral.param`.",
                        sep = "\n" ) )
        message( paste( "\033[32mThe name of the cytomter in your FCS files does not match the",
                        "name in the `asp` parameter list. Please ensure you have selected",
                        "the appropriate option when calling `get.autospectral.param`.\033[0m",
                        sep = "\n" ) )
        warning.list$cytometer.asp.match <- list( AutoSpectral = asp$cytometer,
                                                  FCS = unique( cytometer.list ) )
      }
    }
  }

  # check for matches with peak channel names
  parameter.set <- unique( param.names.list )[[ 1 ]]
  peak.channels <- control.table$channel[ !is.na( control.table$channel ) ]
  parameter.matches <- peak.channels %in% parameter.set

  if ( !all( parameter.matches ) ) {
    mismatched.peak.channels <- peak.channels[ !parameter.matches ]
    warning( paste( "The following channels in your control table do not appear",
                    "in the parameter names of the FCS files.", sep = "\n" ) )
    message( paste( "\033[31mThe following channels in your control table do not appear",
                    "in the parameter names of the FCS files:\033[0m", sep = "\n" ) )
    message( paste( mismatched.peak.channels, collapse = "\n" ) )
    message( "\033[31mThis will fail.\033[0m" )
    errors$mismatched.peak.channels <- mismatched.peak.channels
  }

  if ( length( errors ) > 0 ) {
    if ( strict ) {
      stop( paste( errors, collapse = "\n" ) )
    } else {
      warning( "Critical errors found in control file" )
    }

    if ( length( warning.list ) == 0 )
      warning.list <- "See Errors"
  } else {
    message( "\033[34m No critical errors found in control file.\033[0m" )
    errors <- "No Errors Found"
    if ( length( warning.list ) == 0 )
      warning.list <- "No Issues Detected"
  }

  error.checks <- list( Errors = errors, Warnings = warning.list )

  return( error.checks )

}
