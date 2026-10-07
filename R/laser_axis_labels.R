# laser_axis_labels.R
#
# Shared helpers for colouring detector axis labels by excitation laser.
#
# The colours are drawn as a geom_text() layer just below the panel (one colour
# per label, given as a layer parameter), because vectorised colour input to
# element_text() is not supported by ggplot2. The plot's own axis text is kept
# but made transparent, so it reserves exactly the space the original labels
# used. Plots using this need `clip = "off"` on their coordinate system.

utils::globalVariables( c( "label.x", "label.text" ) )

# Default laser colours (names match the *_laser columns of
# cytometer_database.csv).
.laser.colors <- c(
  DeepUV      = "#4B2A7B",
  UV          = "#8E44AD",
  Violet      = "#3F51B5",
  Blue        = "#1E88E5",
  YellowGreen = "#7CB342",
  Red         = "#D32F2F",
  IR          = "#8D6E63"
)

# Column of cytometer_database.csv holding the detector names for `asp`
# (same mapping as create.control.file()). NULL if not recognised.
.cytometer.db.column <- function( asp ) {
  if ( is.null( asp ) || is.null( asp$cytometer ) ) return( NULL )
  switch(
    asp$cytometer,
    "Aurora"          = if ( identical( asp$cytometer.version, "NL" ) ) "NorthernLights" else "Aurora",
    "ID7000"          = "ID7000",
    "FACSDiscover A8" = "Discover",
    "FACSDiscover S8" = "Discover",
    "FACSDiscover"    = "Discover",
    "Opteon"          = "Opteon",
    "Mosaic"          = "Mosaic",
    "Xenith"          = "Xenith",
    "Symphony"        = "A5SE",
    "CytoStellar"     = "CytoStellar",
    NULL
  )
}

# Laser of each detector from one database column; NA where not found.
# Handles ID7000 alternative channel spellings and CytoStellar base names.
.lasers.from.column <- function( detectors, db.col, database ) {
  db.detectors <- database[[ db.col ]]
  db.lasers    <- database[[ paste0( db.col, "_laser" ) ]]
  keep <- !is.na( db.detectors ) & db.detectors != ""
  db.detectors <- db.detectors[ keep ]
  db.lasers    <- db.lasers[ keep ]

  query <- detectors
  if ( db.col == "ID7000" )
    db.detectors <- .align.id7000.channels( db.detectors, detectors, database = database )
  if ( db.col == "CytoStellar" )
    query <- sub( "-[AH]$", "", query )

  db.lasers[ match( query, db.detectors ) ]
}

#' @title Excitation Laser of Each Detector
#'
#' @description
#' Looks up the excitation laser of each detector in
#' `cytometer_database.csv`. With `asp`, the database column of that
#' cytometer is used; without it, the column whose detectors cover all of
#' `detectors` is used (if several columns match, they must agree on every
#' laser). Every detector must map; otherwise `NULL` is returned. A message
#' is given only if some, but not all, labels match detectors: an axis where
#' nothing matches (e.g. fluorophore names on a spillover-spread heatmap) is
#' treated as not being a detector axis and returns `NULL` silently.
#'
#' @param detectors Character vector of detector (channel) names.
#' @param asp Optional AutoSpectral parameter list.
#' @param verbose Logical, report why lasers could not be assigned.
#' @return Named character vector of laser names (names = `detectors`), or
#'   `NULL` if the detectors do not all map to the instrument's detectors.
#' @noRd
.detector.lasers <- function( detectors, asp = NULL, verbose = TRUE ) {

  detectors <- as.character( detectors )
  if ( length( detectors ) == 0 ) return( NULL )

  db.path <- system.file( "extdata", "cytometer_database.csv", package = "AutoSpectral" )
  if ( db.path == "" ) {
    if ( verbose ) message( "cytometer_database.csv not found; detector labels will be black." )
    return( NULL )
  }
  database <- utils::read.csv( db.path, stringsAsFactors = FALSE, check.names = FALSE )

  candidate.cols <- grep( "_(laser|alt|description)$", colnames( database ),
                          value = TRUE, invert = TRUE )
  candidate.cols <- candidate.cols[ paste0( candidate.cols, "_laser" ) %in% colnames( database ) ]

  asp.col <- .cytometer.db.column( asp )

  if ( !is.null( asp.col ) ) {
    lasers <- .lasers.from.column( detectors, asp.col, database )
    if ( anyNA( lasers ) ) {
      if ( verbose && !all( is.na( lasers ) ) )
        message( "Detector labels will be black: ", sum( is.na( lasers ) ), " of ",
                 length( lasers ), " detector(s) are not ", asp.col, " detectors (",
                 paste( utils::head( detectors[ is.na( lasers ) ], 5 ), collapse = ", " ),
                 if ( sum( is.na( lasers ) ) > 5 ) ", ..." else "", ")." )
      return( NULL )
    }
    return( stats::setNames( lasers, detectors ) )
  }

  # no (recognised) cytometer: find the column(s) covering every detector
  full <- lapply( candidate.cols, function( cc ) .lasers.from.column( detectors, cc, database ) )
  names( full ) <- candidate.cols
  any.match <- any( vapply( full, function( l ) any( !is.na( l ) ), logical( 1 ) ) )
  full <- full[ vapply( full, function( l ) !anyNA( l ), logical( 1 ) ) ]
  # CytoStellar is matched on base names (suffix stripped), so e.g. Aurora
  # "UV1-A" would also match its "UV1"; use it only if nothing else matches
  if ( length( full ) > 1 ) full <- full[ names( full ) != "CytoStellar" ]

  if ( length( full ) == 0 ) {
    if ( verbose && any.match )
      message( "Detector labels will be black: the detectors do not all match any ",
               "cytometer in cytometer_database.csv." )
    return( NULL )
  }

  agree <- all( vapply( full, identical, logical( 1 ), full[[ 1 ]] ) )
  if ( !agree ) {
    if ( verbose )
      message( "Detector labels will be black: the detectors match several cytometers (",
               paste( names( full ), collapse = ", " ), ") with different laser ",
               "assignments. Supply `asp` to resolve." )
    return( NULL )
  }

  stats::setNames( full[[ 1 ]], detectors )
}

#' @title Laser Colours for Detector Labels
#'
#' @param detectors Character vector of detector names, in label order.
#' @param asp Optional AutoSpectral parameter list.
#' @param laser.colors Named colour vector (names = laser names).
#' @param fallback Colour used when lasers cannot be assigned. Default "black".
#' @param verbose Logical.
#' @return list( colors = one colour per detector, mapped = TRUE if laser
#'   colours were assigned, FALSE if the fallback was used ).
#' @noRd
.laser.label.colors <- function(
    detectors, asp = NULL, laser.colors = .laser.colors, fallback = "black", verbose = TRUE
) {
  lasers <- .detector.lasers( detectors, asp = asp, verbose = verbose )
  if ( is.null( lasers ) )
    return( list( colors = rep( fallback, length( detectors ) ), mapped = FALSE ) )

  colors <- unname( laser.colors[ lasers ] )
  colors[ is.na( colors ) ] <- fallback
  list( colors = colors, mapped = TRUE )
}

#' @title Laser-Coloured Detector Axis Labels
#'
#' @description
#' Returns ggplot components that colour the x-axis detector labels by
#' excitation laser: a geom_text() layer at the bottom edge of the panel plus
#' a theme that keeps the original axis text transparent (so it reserves the
#' same space). If lasers cannot be assigned (or `color.by.laser = FALSE`),
#' only a theme setting the axis text to black is returned. The calling plot
#' must use `clip = "off"` in its coordinate system when `mapped` is TRUE.
#'
#' @param x Label positions on the x axis (factor/character values for a
#'   discrete axis, numbers for a continuous one).
#' @param labels Detector names (same order as `x`).
#' @param asp Optional AutoSpectral parameter list.
#' @param color.by.laser Logical. Default TRUE.
#' @param angle,hjust,vjust Text angle and justification, as in the plot's
#'   axis.text.x.
#' @param size Label size in points. Default 8.8 (theme_minimal() /
#'   theme_classic() axis text at base size 11).
#' @param facet Optional one-row data.frame of facet variable values; labels
#'   are then drawn in that facet only (e.g. the bottom panel).
#' @param laser.colors Named colour vector.
#' @param verbose Logical.
#' @return list( components = list of ggplot components, mapped = logical ).
#' @noRd
.laser.axis.labels <- function(
    x, labels, asp = NULL, color.by.laser = TRUE,
    angle = 45, hjust = 1, vjust = 1, size = 8.8,
    facet = NULL, laser.colors = .laser.colors, verbose = TRUE
) {

  black.axis <- list( ggplot2::theme( axis.text.x = ggplot2::element_text( colour = "black" ) ) )
  if ( !isTRUE( color.by.laser ) ) return( list( components = black.axis, mapped = FALSE ) )

  cols <- .laser.label.colors( labels, asp = asp, laser.colors = laser.colors, verbose = verbose )
  if ( !cols$mapped ) return( list( components = black.axis, mapped = FALSE ) )

  # trailing space gives a small gap between the panel edge and the text
  label.df <- data.frame( label.x = x, label.text = paste0( labels, " " ), stringsAsFactors = FALSE )
  if ( !is.null( facet ) )
    label.df <- cbind( label.df, facet[ rep( 1, nrow( label.df ) ), , drop = FALSE ] )

  components <- list(
    ggplot2::geom_text(
      data = label.df,
      mapping = ggplot2::aes( x = label.x, y = -Inf, label = label.text ),
      inherit.aes = FALSE,
      colour = cols$colors,
      angle = angle, hjust = hjust, vjust = vjust,
      size = size / ggplot2::.pt
    ),
    ggplot2::theme( axis.text.x = ggplot2::element_text( colour = "transparent" ) )
  )

  list( components = components, mapped = TRUE )
}
