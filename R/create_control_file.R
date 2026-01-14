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
#' @importFrom stats setNames
#'
#' @param control.dir file path to the single stained control fcs files
#' @param asp The AutoSpectral parameter list. Generate using
#' `get.autospectral.param`
#'
#' @return No returns. Outputs a csv file called fcs_control_file.csv
#' @export

create.control.file <- function( control.dir, asp ){

  # check for existing control file and generate a new name if it exists
  control.file.name.base <- "fcs_control_file"
  control.file.name <- paste0( control.file.name.base, ".csv" )
  file.count <- 1

  while ( file.exists( control.file.name ) ) {
    control.file.name <- paste0( control.file.name.base, "_", file.count, ".csv")
    file.count <- file.count + 1
  }

  # find the files
  control.files <- list.files( control.dir, pattern = ".fcs", ignore.case = TRUE )

  if ( is.null( control.files ) || length( control.files ) <= 1 ) {
    stop( "Single-stained control files not found. Check directory.", call. = FALSE )
  }

  control.colnames <- c( "filename", "fluorophore", "marker", "channel",
                         "control.type", "universal.negative", "large.gate" )

  control.def.file <- data.frame( matrix( ncol = length( control.colnames ),
                                                nrow = length( control.files ) ) )

  colnames( control.def.file ) <- control.colnames

  # define filenames
  control.def.file$filename <- control.files

  # find corresponding fluorophores if possible
  fluor.data.path <- system.file(
    "extdata", "fluorophore_database.csv", package = "AutoSpectral"
    )
  fluorophore.database <- read.csv( fluor.data.path )
  fluorophore.database[ fluorophore.database == "" ] <- NA
  fluorophore.matches <- match.fluorophores( control.files, fluorophore.database )

  control.def.file$fluorophore <- fluorophore.matches[ control.def.file$filename ]

  # find markers if possible
  marker.data.path <- system.file(
    "extdata", "marker_database.csv", package = "AutoSpectral"
    )
  marker.database <- read.csv( marker.data.path )
  marker.database[ marker.database == "" ] <- NA
  marker.matches <- match.markers( control.files, marker.database )
  control.def.file$marker <- marker.matches[ control.def.file$filename ]

  # set corresponding peak detectors based on cytometer
  if ( asp$cytometer == "Aurora" ) {
    if ( asp$cytometer.version == "NL" ) {
      detectors <- setNames(
        fluorophore.database$channel.NL, fluorophore.database$fluorophore )
    } else {
      detectors <- setNames(
        fluorophore.database$channel.Aurora, fluorophore.database$fluorophore )
    }
  } else if ( asp$cytometer == "ID7000" ) {
    detectors <- setNames(
      fluorophore.database$channel.ID7000, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "FACSDiscover A8" ) {
    detectors <- setNames(
      fluorophore.database$channel.s8, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "FACSDiscover S8" ) {
    detectors <- setNames(
      fluorophore.database$channel.s8, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Opteon" ) {
    detectors <- setNames(
      fluorophore.database$channel.opteon, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Mosaic" ) {
    detectors <- setNames(
      fluorophore.database$channel.mosaic, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Xenith" ) {
    detectors <- setNames(
      fluorophore.database$channel.xenith, fluorophore.database$fluorophore )
  } else if ( asp$cytometer == "Symphony" ) {
    detectors <- setNames(
      fluorophore.database$channel.A5SE, fluorophore.database$fluorophore )
  } else {
    stop( "Unsupported cytometer" )
  }

  detector.idx <- match( control.def.file$fluorophore, names( detectors ) )

  control.def.file$channel <- detectors[ detector.idx ]

  control.def.file$control.type <- sapply(
    control.def.file$filename, function( filename ) {
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
  control.def.file.merged <- merge(
    control.def.file, fluorophore.database, by = "fluorophore", all.x = TRUE )

  laser.order <- c( "DeepUV", "UV", "Violet", "Blue", "YellowGreen", "Red", "IR" )

  control.def.file.merged$excitation.laser <- factor(
    control.def.file.merged$excitation.laser,
    levels = laser.order
    )

  control.def.file.merged <- control.def.file.merged[
    order( control.def.file.merged$excitation.laser,
           control.def.file.merged$nominal.wavelength ), ]

  desired.col <- c( control.colnames, "is.viability" )

  control.def.file <- control.def.file.merged[, desired.col ]

  # fill AF for unstained cells, Negative for unstained beads
  control.def.file$fluorophore[ grepl( "Unstained", control.def.file$filename ) ] <-
    ifelse(
      control.def.file$control.type[ grepl(
        "Unstained", control.def.file$filename ) ] == "cells",
           "AF",
           "Negative" )

  # fill Negative for Negative
  control.def.file$fluorophore[ grepl(
    "Negative", control.def.file$filename ) ] <- "Negative"

  # replace any NAs
  control.def.file[ is.na( control.def.file ) ] <- ""
  control.def.file[ control.def.file == "NA" ] <- ""

  write.csv( control.def.file, file = control.file.name, row.names = FALSE )

  # check for duplicate fluorophores
  duplicate.fluorophores <- anyDuplicated( control.def.file$fluorophore )

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
