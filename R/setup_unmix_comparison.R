# setup_unmix_comparison.R

## Internal helper. Resolves `fluorophore.list` to a plain character vector,
## whether the caller supplied one directly or pointed to an existing
## AutoSpectral control-file CSV (as produced by `create.control.file()`).
## `AF`, `Negative`, and `No match` rows are dropped automatically when
## reading from a control file, since these are not stained fluorophores.
##
## @keywords internal
.resolve.fluorophore.list <- function( fluorophore.list ) {
  
  is.control.file <- is.character( fluorophore.list ) &&
    length( fluorophore.list ) == 1 &&
    grepl( "\\.csv$", fluorophore.list, ignore.case = TRUE )
  
  if ( !is.control.file ) {
    return( unique( fluorophore.list ) )
  }
  
  if ( !file.exists( fluorophore.list ) ) {
    stop(
      paste0( "Control-file CSV not found: '", fluorophore.list, "'" ),
      call. = FALSE
    )
  }
  
  ct <- utils::read.csv(
    fluorophore.list, stringsAsFactors = FALSE, strip.white = TRUE
  )
  
  if ( !"fluorophore" %in% colnames( ct ) ) {
    stop(
      paste0(
        "'", fluorophore.list, "' does not contain a 'fluorophore' column; ",
        "is this an AutoSpectral control file?"
      ),
      call. = FALSE
    )
  }
  
  excluded  <- c( "AF", "Negative", "No match" )
  is.usable <- !is.na( ct$fluorophore ) &
    ct$fluorophore != "" &
    !ct$fluorophore %in% excluded
  
  unique( ct$fluorophore[ is.usable ] )
}


## Internal helper. Builds the combined setup table for a single folder:
## one row per fluorophore in `fluorophore.list`, giving the matched
## single-stained control filename and matched unmixed detector channel,
## plus one row for the unstained control. Any ambiguity is recorded in the
## `flag` column rather than raising an error, so a single problematic
## fluorophore doesn't block the rest of the folder from being scanned.
##
## @keywords internal
.setup.one.folder <- function(
    folder.path,
    label,
    fluorophore.list,
    panel.database,
    unstained.regex,
    pattern,
    verbose
) {
  
  fcs.files <- list.files( folder.path, pattern = pattern, ignore.case = TRUE )
  
  if ( length( fcs.files ) == 0 ) {
    stop(
      paste0(
        "No FCS files found in folder '", label, "' (", folder.path, ")."
      ),
      call. = FALSE
    )
  }
  
  # --- identify the unstained file -----------------------------------------
  
  is.unstained <- grepl( unstained.regex, fcs.files, ignore.case = TRUE )
  
  if ( sum( is.unstained ) == 0 ) {
    stop(
      paste0(
        "No file matching `unstained.regex` ('", unstained.regex,
        "') found in folder '", label, "'."
      ),
      call. = FALSE
    )
  }
  
  if ( sum( is.unstained ) > 1 ) {
    unstained.file <- sort( fcs.files[ is.unstained ] )[ 1 ]
    warning(
      paste0(
        "\033[31mFolder '", label, "': multiple candidate unstained files ",
        "found; using '", unstained.file, "'. Edit the setup CSV if this is ",
        "wrong.\033[0m"
      ),
      call. = FALSE
    )
  } else {
    unstained.file <- fcs.files[ is.unstained ]
  }
  
  single.stain.files <- fcs.files[ !fcs.files %in% unstained.file ]
  
  # --- match filenames to fluorophores -------------------------------------
  
  fluorophore.matches <- suppressMessages(
    match.fluorophores( single.stain.files, panel.database )
  )
  
  # --- read headers, check $PnN consistency against the unstained file ----
  
  headers <- lapply(
    fcs.files,
    function( f ) readFCSheader( file.path( folder.path, f ) )[[ 1 ]]
  )
  names( headers ) <- fcs.files
  
  param.sets <- lapply( headers, function( h ) {
    unname( unlist( h[ grep( "^\\$P\\d+N$", names( h ) ) ] ) )
  } )
  
  ref.params <- param.sets[[ unstained.file ]]
  
  is.consistent <- vapply(
    param.sets, function( p ) setequal( p, ref.params ), logical( 1 )
  )
  inconsistent.files <- names( param.sets )[ !is.consistent ]
  
  if ( length( inconsistent.files ) > 0 ) {
    warning(
      paste0(
        "\033[31mFolder '", label, "': the following file(s) have detector ",
        "parameters inconsistent with the unstained control and may not be ",
        "comparable: ", paste( inconsistent.files, collapse = ", " ),
        "\033[0m"
      ),
      call. = FALSE
    )
  }
  
  # --- match reference channel names to fluorophores (the channel map) ----
  
  channel.matches <- suppressMessages(
    match.fluorophores( unlist( ref.params, use.names = FALSE ), panel.database )
  )
  
  matched.idx <- which( channel.matches %in% fluorophore.list )
  channel.map <- stats::setNames(
    ref.params[ matched.idx ], channel.matches[ matched.idx ]
  )
  
  dup.fluor <- unique( names( channel.map )[ duplicated( names( channel.map ) ) ] )
  
  # --- assemble one row per fluorophore ------------------------------------
  
  rows <- lapply( fluorophore.list, function( fl ) {
    
    file.idx <- which( fluorophore.matches == fl )
    channel  <- if ( fl %in% names( channel.map ) ) channel.map[[ fl ]] else NA_character_
    
    flag <- character( 0 )
    if ( length( file.idx ) == 0 ) flag <- c( flag, "No file match" )
    if ( length( file.idx ) > 1 )  flag <- c( flag, "Multiple file matches" )
    if ( is.na( channel ) )        flag <- c( flag, "No channel match" )
    if ( fl %in% dup.fluor )       flag <- c( flag, "Multiple channel matches" )
    
    matched.filename <- if ( length( file.idx ) >= 1 ) {
      paste( single.stain.files[ file.idx ], collapse = "; " )
    } else {
      "No match"
    }
    
    data.frame(
      fluorophore = fl,
      channel     = if ( is.na( channel ) ) "No match" else channel,
      filename    = matched.filename,
      flag        = if ( length( flag ) == 0 ) "OK" else paste( flag, collapse = "; " ),
      stringsAsFactors = FALSE
    )
  } )
  
  setup.table <- do.call( rbind, rows )
  
  unstained.row <- data.frame(
    fluorophore = "Unstained",
    channel     = NA_character_,
    filename    = unstained.file,
    flag        = "OK",
    stringsAsFactors = FALSE
  )
  
  setup.table           <- rbind( unstained.row, setup.table )
  rownames( setup.table ) <- NULL
  
  setup.table
}


#' @title Set Up an Unmixing Comparison
#'
#' @description
#' Scans one or more folders of already-unmixed FCS files (one folder per
#' operator or unmixing method to be compared) and builds a per-folder table
#' mapping every fluorophore in a supplied panel to its matched single-stained
#' control filename and matched unmixed detector channel, plus the unstained
#' control shared as the negative reference. One CSV is written per folder for
#' manual review before calling `compare.unmix.folders()`.
#'
#' Filenames are matched to fluorophore identity with `match.fluorophores()`.
#' Detector channels are matched to fluorophore identity the same way, applied
#' to the `$PnN` parameter names read from each FCS header and restricted to
#' the supplied panel, so that unrelated substring matches elsewhere in the
#' fluorophore database cannot be picked up. Any row with `flag != "OK"` (a
#' missing or ambiguous file/channel match, or parameter inconsistency within
#' the folder) must be corrected by hand in the output CSV before proceeding.
#'
#' @param folders Named character vector of directory paths, one per unmixed
#' result set to compare (e.g. one per operator, plus one for AutoSpectral).
#' Names are used as labels throughout downstream analysis and plots. If
#' `folders` is unnamed (or partially named), directory basenames are used
#' instead and a warning is issued.
#' @param fluorophore.list Either a character vector of fluorophore names
#' (matching the `fluorophore` column of `fluorophore_database.csv`) defining
#' the panel, or a character string giving the path to an existing
#' AutoSpectral control-file CSV (as produced by `create.control.file()`),
#' from which the `fluorophore` column is used (excluding `AF`, `Negative`,
#' and `No match` rows).
#' @param unstained.regex Character. Case-insensitive regular expression used
#' to identify the unstained control file in each folder. Default
#' `"unstained"`.
#' @param fluorophore.database Data frame of fluorophore names/synonyms passed
#' to `match.fluorophores()`. Default `NULL` loads the bundled
#' `fluorophore_database.csv`.
#' @param pattern Character. Regular expression used to find FCS files in each
#' folder. Default `"\\.fcs$"`.
#' @param output.dir Character. Directory in which the per-folder setup CSVs
#' are written. Created if absent. Default `"./unmix_comparison_setup"`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Invisibly, a named list (one element per folder, using the names of
#' `folders`) of the setup data frames written to `output.dir`. Each data
#' frame has one row per fluorophore in `fluorophore.list` plus one row for
#' the unstained control, with columns `fluorophore`, `channel`, `filename`,
#' and `flag`.
#'
#' @export

setup.unmix.comparison <- function(
    folders,
    fluorophore.list,
    unstained.regex      = "unstained",
    fluorophore.database  = NULL,
    pattern               = "\\.fcs$",
    output.dir            = "./unmix_comparison_setup",
    verbose               = TRUE
) {
  
  # --- resolve folder labels ------------------------------------------------
  
  if ( is.null( names( folders ) ) || any( names( folders ) == "" ) ) {
    warning(
      "`folders` is unnamed (or partially named); using directory basenames as labels.",
      call. = FALSE
    )
    names( folders ) <- basename( folders )
  }
  
  if ( any( duplicated( names( folders ) ) ) ) {
    stop( "`folders` names (labels) must be unique.", call. = FALSE )
  }
  
  missing.dirs <- folders[ !dir.exists( folders ) ]
  if ( length( missing.dirs ) > 0 ) {
    stop(
      paste0(
        "The following folders do not exist: ",
        paste( missing.dirs, collapse = ", " )
      ),
      call. = FALSE
    )
  }
  
  # --- resolve fluorophore panel --------------------------------------------
  
  fluorophore.list <- .resolve.fluorophore.list( fluorophore.list )
  
  if ( is.null( fluorophore.database ) ) {
    fluor.data.path <- system.file(
      "extdata", "fluorophore_database.csv", package = "AutoSpectral"
    )
    fluorophore.database <- utils::read.csv( fluor.data.path )
    fluorophore.database[ fluorophore.database == "" ] <- NA
  }
  
  # `fluorophore.list` entries are matched against the full database (not yet
  # restricted) using the same fuzzy matching used for filenames, so minor
  # formatting drift (whitespace, case, or a synonym instead of the primary
  # name) is resolved instead of silently producing an empty `panel.database`
  # and a channel map that can never match anything.
  
  canonical.fluor <- suppressMessages(
    match.fluorophores( fluorophore.list, fluorophore.database )
  )
  
  unmatched.panel <- fluorophore.list[ canonical.fluor == "No match" ]
  if ( length( unmatched.panel ) > 0 ) {
    stop(
      paste0(
        "The following `fluorophore.list` entries could not be matched to ",
        "any fluorophore (or synonym) in `fluorophore.database`: ",
        paste( unmatched.panel, collapse = ", " ),
        ". Check spelling/formatting against fluorophore_database.csv, or ",
        "pass a corrected `fluorophore.database`."
      ),
      call. = FALSE
    )
  }
  
  fluorophore.list <- unique( canonical.fluor )
  
  panel.database <- fluorophore.database[
    fluorophore.database$fluorophore %in% fluorophore.list, , drop = FALSE
  ]
  
  if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )
  
  # --- scan each folder ------------------------------------------------------
  
  setup.tables <- list()
  
  for ( label in names( folders ) ) {
    
    if ( verbose ) {
      message( sprintf(
        "\033[34mScanning folder: %s (%s)\033[0m", label, folders[[ label ]]
      ) )
    }
    
    setup.tables[[ label ]] <- .setup.one.folder(
      folder.path      = folders[[ label ]],
      label            = label,
      fluorophore.list = fluorophore.list,
      panel.database   = panel.database,
      unstained.regex  = unstained.regex,
      pattern          = pattern,
      verbose          = verbose
    )
    
    out.file <- file.path(
      output.dir, paste0( label, "_unmix_comparison_setup.csv" )
    )
    utils::write.csv( setup.tables[[ label ]], out.file, row.names = FALSE )
    
    if ( verbose ) message( sprintf( "  Wrote: %s", out.file ) )
  }
  
  n.flagged <- sum( vapply(
    setup.tables, function( tb ) sum( tb$flag != "OK" ), integer( 1 )
  ) )
  
  if ( n.flagged > 0 ) {
    warning(
      paste0(
        "\033[31m", n.flagged, " row(s) across all folders need manual ",
        "review (flag != 'OK'). Edit the CSVs in '", output.dir,
        "' before running compare.unmix.folders().\033[0m"
      ),
      call. = FALSE
    )
  } else if ( verbose ) {
    message( "\033[32mAll rows matched cleanly (flag == 'OK') in every folder.\033[0m" )
  }
  
  return( invisible( setup.tables ) )
}