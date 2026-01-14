# validate_control_file.r

#' @title Validate Control File
#' @description
#' Checks whether the control file used to define the single-stained control
#' setup for AutoSpectral conforms to expectations and follows rules required for
#' successful running of AutoSpectral.
#'
#' @importFrom utils read.csv
#' @importFrom flowCore read.FCSheader
#'
#' @param control.dir File path to the single-stained control FCS files.
#' @param control.def.file CSV file defining the single-color control file names,
#' fluorophores they represent, marker names, peak channels, and gating requirements.
#' @param asp The AutoSpectral parameter list.
#'
#' @return A dataframe of errors and warnings intended to help the user fix
#' problems with the `control.def.file`.

validate.control.file <- function( control.dir, control.def.file, asp ) {

  issues <- list()

  ## ---------- filesystem ----------
  if ( !file.exists( control.def.file ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "missing_control_file",
                  message = paste( "Unable to locate control.def.file:", control.def.file ) )
    return( do.call( rbind, issues ) )
  }

  if ( !dir.exists( control.dir ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "missing_control_dir",
                  message = paste( "Unable to locate control.dir:", control.dir ) )
    return( do.call( rbind, issues ) )
  }

  fcs.files <- list.files( control.dir, pattern = "\\.fcs$", ignore.case = TRUE )
  if ( length( fcs.files ) < 1 ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue("error", "no_fcs_files",
                 message = paste( "No FCS files found in:", control.dir ) )
    return( do.call( rbind, issues ) )
  }

  ## ---------- load table ----------
  ct <- read.csv(
    control.def.file,
    stringsAsFactors = FALSE,
    strip.white = TRUE
  )

  # trim white space, convert blanks to NAs
  ct[] <- lapply( ct, function( x ) {
    if ( is.character( x ) ) {
      x <- trimws( x )
      x[ x == "" ] <- NA
      x
    } else x
  } )

  required.cols <- c(
    "filename", "fluorophore", "marker", "channel",
    "control.type", "universal.negative",
    "large.gate", "is.viability"
  )

  missing.cols <- setdiff( required.cols, colnames( ct ) )
  if ( length( missing.cols ) > 0) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "missing_columns",
                  column = paste( missing.cols, collapse = ", " ),
                  message = "Required columns missing from control file" )
    return( do.call( rbind, issues ) )
  }

  ## ---------- row count ----------
  if ( nrow( ct ) < 2 ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue("error", "too_few_rows",
                 message = "Control table contains fewer than two rows" )
  } else if ( nrow( ct ) < 3 ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "warning", "low_row_count",
                  message = paste( "Only", nrow(ct), "rows present in control table" ) )
  }

  ## ---------- filename ----------
  if ( any( is.na( ct$filename ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "warning", "missing_filename",
                  filename = ct$filename[ is.na( ct$filename ) ],
                  column = "filename",
                  message = "Filename missing" )
  }

  if ( any( duplicated( ct$filename ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "duplicate_filename",
                  filename = ct$filename[ duplicated( ct$filename ) ],
                  column = "filename",
                  message = "Duplicate filenames detected" )
  }

  bad.ext <- !grepl( "\\.fcs$", ct$filename, ignore.case = TRUE )
  if ( any( bad.ext ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "non_fcs_filename",
                  filename = ct$filename[ bad.ext ],
                  column = "filename",
                  message = "Filename does not end in .fcs" )
  }

  ## ---------- fluorophore ----------
  if ( any( is.na( ct$fluorophore ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "missing_fluorophore",
                  filename = ct$filename[ is.na( ct$fluorophore ) ],
                  column = "fluorophore",
                  message = "Fluorophore missing" )
  }

  if ( any( duplicated( ct$fluorophore ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "duplicate_fluorophore",
                  column = "fluorophore",
                  message = "Duplicate fluorophore names detected" )
  }

  no.match <- grepl( "No match", ct$fluorophore, ignore.case = TRUE )
  if ( any( no.match ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "no_match_present",
                  filename = ct$filename[ no.match ],
                  column = "fluorophore",
                  message = "`No match` still present in fluorophore column" )
  }

  ## ---------- AF ----------
  af.idx <- ct$fluorophore[ !is.na( ct$fluorophore ) ] == "AF"
  if ( !any( af.idx ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "missing_af",
                  column = "filename",
                  message = "No AF (autofluorescence) control provided" )
  } else {
    if ( sum( af.idx ) > 1 ) {
      issues[[ length( issues ) + 1 ]] <-
        .new_issue( "error", "duplicate_af",
                    filename = ct$filename[ af.idx ],
                    message = "Multiple AF controls detected" )
    }
    if ( any( ct$control.type[ af.idx ] != "cells" ) ) {
      issues[[ length( issues ) + 1 ]] <-
        .new_issue( "error", "af_wrong_type",
                    filename = ct$filename[ af.idx ],
                    column = "control.type",
                    message = "AF control must be of type 'cells'" )
    }
  }

  ## ---------- marker ----------
  acceptable.missing <- c( "AF", "negative" )
  missing.marker <- is.na( ct$marker )
  if ( any( missing.marker ) ) {
    bad <- !grepl( paste( acceptable.missing, collapse = "|" ),
                   ct$fluorophore[ missing.marker ],
                   ignore.case = TRUE )
    if ( any( bad ) ) {
      issues[[ length( issues ) + 1 ]] <-
        .new_issue( "warning", "missing_marker",
                    filename = ct$filename[ missing.marker ][ bad] ,
                    column = "marker",
                    message = "Marker missing where required" )
    }
  }

  ## ---------- channel ----------
  missing.channel <- is.na( ct$channel )
  if ( any( missing.channel ) ) {
    bad <- !grepl( paste( acceptable.missing, collapse = "|" ),
                   ct$fluorophore[missing.channel ],
                   ignore.case = TRUE )
    if ( any( bad ) ) {
      issues[[ length( issues ) + 1 ]] <-
        .new_issue( "error", "missing_channel",
                    filename = ct$filename[ missing.channel ][ bad ],
                    column = "channel",
                    message = "Channel missing where required" )
    }
  }

  ## ---------- control.type ----------
  if ( any( is.na( ct$control.type ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "missing_control_type",
                  filename = ct$filename[ is.na( ct$control.type ) ],
                  column = "control.type",
                  message = "control.type missing" )
  }

  bad.type <- !ct$control.type %in% c( "cells", "beads" )
  if ( any( bad.type, na.rm = TRUE ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "invalid_control_type",
                  filename = ct$filename[ bad.type ],
                  column = "control.type",
                  message = "control.type must be 'cells' or 'beads'" )
  }

  ## ---------- universal negative ----------
  missing.un <- is.na( ct$universal.negative )
  if ( all( missing.un ) ) {
    .new_issue( "warning", "no_universal_negative",
                column = "universal.negative",
                message = "No universal negative specified" )
  }

  un <- na.omit( unique( ct$universal.negative ) )
  absent <- setdiff( un, ct$filename )
  if ( length( absent ) > 0) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "universal_negative_missing",
                  filename = absent,
                  column = "filename",
                  message = "Universal negative not listed as a sample" )
  }

  neg.fluor <- ct$fluorophore[ match( un, ct$filename ) ]
  bad.neg.name <- !grepl( "^AF$|^Negative", neg.fluor )
  if ( any( bad.neg.name ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "invalid_negative_fluorophore",
                  filename = un[ bad.neg.name ],
                  column = "fluorophore",
                  message = "Universal negative fluorophore must be AF or Negative" )
  }

  sample.with.neg <- ct$filename[ !is.na( ct$universal.negative ) ]
  # can only run check if filenames are not duplicated
  if ( !any( duplicated( ct$filename ) ) ) {
    mismatch <- vapply( sample.with.neg, function( s ) {
      ct$control.type[ ct$filename == s ] !=
        ct$control.type[ ct$filename == ct$universal.negative[ ct$filename == s ] ]
    }, logical( 1 ) )

    if ( any( mismatch[ !is.na( mismatch ) ] ) ) {
      issues[[ length( issues ) + 1 ]] <-
        .new_issue( "error", "negative_type_mismatch",
                    filename = sample.with.neg[ mismatch ],
                    column = "universal.negative",
                    message = "Sample and universal negative control.type mismatch" )
    }
  }

  ## ---------- logical fields ----------
  lg <- as.logical( ct$large.gate )
  iv <- as.logical( ct$is.viability )

  if ( all( is.na( lg ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "warning", "no_large_gate_set",
                  message = "No large.gate set for any sample" )
  }

  if ( all( is.na( iv ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "warning", "no_viability_set",
                  column = "is.viability",
                  message = "No is.viability set for any sample" )
  }

  if ( any( !is.na( lg ) & ct$control.type != "cells" & lg ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "large_gate_non_cell",
                  filename = ct$filename[ !is.na( lg ) & lg & ct$control.type != "cells" ],
                  column = "large.gate",
                  message = "large.gate TRUE only allowed for cell controls" )
  }

  if ( any( !is.na( iv ) & ct$control.type != "cells" & iv ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "viability_non_cell",
                  filename = ct$filename[ !is.na( iv ) & iv & ct$control.type != "cells" ],
                  column = "is.viability",
                  message = "is.viability TRUE only allowed for cell controls" )
  }

  ## ---------- FCS headers ----------
  headers <- lapply( ct$filename, function( f ) {
    tryCatch( {
      h <- suppressWarnings(
        flowCore::read.FCSheader( f, path = control.dir,
                                  emptyValue = FALSE )[[ 1 ]]
      )
      n.par <- as.integer( h[ "$PAR" ] )
      list(
        ok = TRUE,
        file = f,
        parameters = vapply( seq_len( n.par ),
                             function( i ) h[[ paste0( "$P", i, "N" ) ]],
                             character( 1 ) ),
        cytometer = h[ "$CYT" ],
        events = as.integer( h[ "$TOT" ] ),
        error = NULL
      )
    }, error = function( e ) {
      list( ok = FALSE, file = f, error = conditionMessage( e ) )
    } )
  } )

  failed <- vapply( headers, `[[`, logical( 1 ), "ok" ) == FALSE
  if ( any( failed ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "fcs_header_read_failed",
                  filename = vapply( headers[ failed ], `[[`, character( 1 ), "file" ),
                  message = vapply( headers[ failed ], `[[`, character( 1 ), "error" ) )
  }

  valid <- headers[ !failed ]
  if ( length( valid ) == 0 ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "no_valid_fcs_headers",
                  message = "No valid FCS headers could be read" )
    return( do.call( rbind, issues ) )
  }

  ## ---------- low events ----------
  ev <- vapply( valid, `[[`, integer( 1 ), "events" )

  if ( any( ev < 5000 ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "warning", "low_event_count",
                  filename = ct$filename[ ev < 5000 ],
                  message = "Event count < 5000" )
  }

  if (any(ev < 1000)) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "very_low_event_count",
                  filename = ct$filename[ ev < 1000 ],
                  message = "Event count < 1000" )
  }

  ## ---------- parameter consistency ----------
  params <- lapply( valid, `[[`, "parameters")

  if ( asp$cytometer %in% c( "FACSDiscover S8", "FACSDiscover A8" ) ) {
    nonspec <- asp$non.spectral.channel[ 4:length( asp$non.spectral.channel ) ]
    params <- lapply( params, function( p )
      p[!grepl(paste(nonspec, collapse = "|"), p ) ] )
  }

  if ( !all( vapply( params[ -1 ], identical, logical( 1 ), params[[ 1 ]] ) ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "parameter_mismatch",
                  message = "Parameter names inconsistent across FCS files" )
  }

  ## ---------- cytometer ----------
  cyts <- vapply( valid, `[[`, character( 1 ), "cytometer" )

  if ( length( unique( cyts ) ) > 1 ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "multiple_cytometers",
                  message = "Multiple cytometers detected" )
  } else if ( !grepl( asp$cytometer, unique( cyts ), ignore.case = TRUE ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "warning", "cytometer_mismatch",
                  message = "Cytometer name does not match asp$cytometer" )
  }

  ## ---------- channel vs parameters ----------
  peak.channels <- ct$channel[ !is.na( ct$channel ) ]
  if ( !all( peak.channels %in% params[[ 1 ]] ) ) {
    issues[[ length( issues ) + 1 ]] <-
      .new_issue( "error", "channel_parameter_mismatch",
                  column = "channel",
                  message = "One or more channels not found in FCS parameters" )
  }

  if ( length( issues ) == 0) {
    return( data.frame() )
  }

  return( do.call( rbind, issues ) )
}

