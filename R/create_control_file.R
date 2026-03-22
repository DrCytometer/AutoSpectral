# create_control_file.r

#' @title Create Control File
#'
#' @description
#' A helper function to draft a description of your single stained control
#' files such that AutoSpectral can understand and process them correctly.
#' Given a set of single stained control fcs files, `create.control.file` will
#' produce a csv file listing the matching peak detector channels for your
#' fluorophores (if known). If your files contain bead or cell tags in the filename,
#' it will assign your controls as cells or beads. You will need to fill in any
#' "No Match" results manually. You will need to set universal negatives manually.
#' You will need to add marker names manually.
#'
#' @param control.dir file path to the single stained control fcs files
#' @param asp The AutoSpectral parameter list. Generate using
#' `get.autospectral.param`
#' @param fill.gate.name Logical, default is `TRUE`. Will attempt to automatically
#' assign gate names for the `gate.name` column if `TRUE`.
#' @param filename Character string defining the output filename. Default is
#' "fcs_control_file", to which .csv will be appended.
#'
#' @return No returns. Outputs a csv file called fcs_control_file.csv
#' @export

create.control.file <- function(
    control.dir,
    asp,
    fill.gate.name = TRUE,
    filename = "fcs_control_file"
  ) {

  # check for existing control file and generate a new name if it exists
  control.file.name <- paste0( filename, ".csv" )
  file.count <- 1

  while ( file.exists( control.file.name ) ) {
    control.file.name <- paste0( filename, "_", file.count, ".csv")
    file.count <- file.count + 1
  }

  # find the files
  control.files <- list.files( control.dir, pattern = ".fcs", ignore.case = TRUE )

  if ( is.null( control.files ) || length( control.files ) <= 1 ) {
    stop( "Single-stained control files not found. Check directory.", call. = FALSE )
  }

  control.colnames <- c(
    "filename", "fluorophore", "marker", "channel", "control.type",
    "universal.negative", "large.gate", "is.viability", "gate.name", "gate.define"
  )

  control.table <- data.frame(
    matrix(
      ncol = length( control.colnames ),
      nrow = length( control.files )
    )
  )
  colnames( control.table ) <- control.colnames

  # define filenames
  control.table$filename <- control.files

  # find corresponding fluorophores if possible
  fluor.data.path <- system.file(
    "extdata", "fluorophore_database.csv", package = "AutoSpectral"
  )
  fluorophore.database <- utils::read.csv( fluor.data.path )
  fluorophore.database[ fluorophore.database == "" ] <- NA
  fluorophore.matches <- match.fluorophores( control.files, fluorophore.database )

  control.table$fluorophore <- fluorophore.matches[ control.table$filename ]

  # find markers if possible
  marker.data.path <- system.file(
    "extdata", "marker_database.csv", package = "AutoSpectral"
  )
  marker.database <- utils::read.csv( marker.data.path )
  marker.database[ marker.database == "" ] <- NA
  marker.matches <- match.markers( control.files, marker.database )
  control.table$marker <- marker.matches[ control.table$filename ]

  # set corresponding peak detectors based on cytometer
  if ( asp$cytometer == "Aurora" ) {
    if ( asp$cytometer.version == "NL" ) {
      detectors <- stats::setNames(
        fluorophore.database$channel.NL, fluorophore.database$fluorophore )
    } else {
      detectors <- stats::setNames(
        fluorophore.database$channel.Aurora, fluorophore.database$fluorophore )
    }
  } else if ( asp$cytometer == "ID7000" ) {
    detectors <- stats::setNames(
      fluorophore.database$channel.ID7000, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "FACSDiscover A8" ) {
    detectors <- stats::setNames(
      fluorophore.database$channel.s8, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "FACSDiscover S8" ) {
    detectors <- stats::setNames(
      fluorophore.database$channel.s8, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Opteon" ) {
    detectors <- stats::setNames(
      fluorophore.database$channel.opteon, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Mosaic" ) {
    detectors <- stats::setNames(
      fluorophore.database$channel.mosaic, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Xenith" ) {
    detectors <- stats::setNames(
      fluorophore.database$channel.xenith, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Symphony" ) {
    detectors <- stats::setNames(
      fluorophore.database$channel.A5SE, fluorophore.database$fluorophore )
  } else {
    stop( "Unsupported cytometer" )
  }

  detector.idx <- match( control.table$fluorophore, names( detectors ) )

  control.table$channel <- detectors[ detector.idx ]

  control.table$control.type <- sapply(
    control.table$filename, function( filename ) {
    if ( grepl( "cells", filename, ignore.case = TRUE ) ){
      type <- "cells"
    } else if ( grepl( "beads", filename, ignore.case = TRUE ) ){
      type <- "beads"
    } else {
      type <- ""
    }
    type
  } )

  # reorder list by wavelength column in database file
  control.table.merged <- merge(
    control.table, fluorophore.database, by = "fluorophore", all.x = TRUE )

  laser.order <- c( "DeepUV", "UV", "Violet", "Blue", "YellowGreen", "Red", "IR" )

  control.table.merged$excitation.laser <- factor(
    control.table.merged$excitation.laser,
    levels = laser.order
    )

  control.table.merged <- control.table.merged[
    order( control.table.merged$excitation.laser,
           control.table.merged$nominal.wavelength ), ]

  # desired.col <- c( control.colnames, "is.viability" )
  control.table <- control.table.merged[ , control.colnames ]

  # automatically fill gate names if requested
  if ( !"gate.name" %in% colnames( control.table ) ) control.table$gate.name <- ""
  if ( !"gate.define" %in% colnames( control.table ) ) control.table$gate.define <- TRUE
  if ( fill.gate.name ) {
    control.table <- assign.gates(
      control.table,
      ggating.system = "density",
      gate = TRUE
    )
  }

  # fill AF for unstained cells, Negative for unstained beads
  control.table$fluorophore[ grepl( "Unstained", control.table$filename ) ] <-
    ifelse(
      control.table$control.type[
        grepl( "Unstained", control.table$filename ) ] == "cells",
           "AF",
           "Negative" )

  # fill Negative for Negative
  control.table$fluorophore[
    grepl( "Negative", control.table$filename ) ] <- "Negative"

  # replace any NAs
  control.table[ is.na( control.table ) ] <- ""
  control.table[ control.table == "NA" ] <- ""


  utils::write.csv(
    control.table,
    file = control.file.name,
    row.names = FALSE
  )

  # check for "No Match"
  no.match <- grepl( "No match", control.table$fluorophore )

  if ( any( no.match ) ) {
    warning(
      paste0(
        "\033[31m",
        "One or more fluorophores could not be matched.",
        "\n",
        "Edit any 'No match' values before proceeding.",
        "\033[0m"
      )
    )
  }

  # check for duplicate fluorophores
  duplicate.fluorophores <- anyDuplicated( control.table$fluorophore[ !no.match ] )

  if ( duplicate.fluorophores != 0 ) {
    warning(
      paste0(
        "\033[31m",
        "Duplicated fluorophore names appear in the control file.",
        "\n",
        "Inspect and remove any extra single color control files ",
        "or edit the control file to be accurate.",
        "\n",
        "Only one control may be used per fluorophore.",
        "\033[0m"
      )
    )
  }

}
