# convert_flowjo_workspace.r

#' @title Convert FlowJo Workspace
#'
#' @description
#' Rewrites the parameter names in a FlowJo workspace (`.wsp`) or template
#' (`.wspt`) so that its gates, plots, statistics and table columns match the
#' parameter names in an AutoSpectral-unmixed FCS file. Instrument software
#' and AutoSpectral name unmixed parameters differently (e.g.
#' "StarBright Violet 440-A" versus "SBV440-A"), so a gating template built on
#' instrument-unmixed files does not apply to AutoSpectral output without this
#' conversion.
#'
#' Each parameter referenced in the workspace is matched to a parameter in
#' `fcs.file` in this order:
#' \enumerate{
#'   \item an entry in `name.map`;
#'   \item an identical parameter name (e.g. scatter and time);
#'   \item the same fluorophore, resolved through the synonyms in
#'   `fluorophore_database.csv`, ignoring a trailing "-A" and differences in
#'   case, spacing and punctuation.
#' }
#' Parameters that cannot be matched, or match more than one parameter, are
#' left unchanged and reported in a warning; supply them via `name.map`.
#'
#' FlowJo prefixes compensated parameters (by default with "Comp-"). A
#' compensated reference keeps the prefix only if the new parameter is listed
#' in the `$SPILLOVER` keyword of `fcs.file`, since FlowJo creates compensated
#' parameters only for those. Table Editor columns ("name :: marker") are
#' given the new name and the marker (`$PnS`) from `fcs.file`.
#'
#' Samples in a `.wsp` refer to the original FCS files, so they are removed,
#' together with their acquisition-defined compensation matrices. Group-level
#' gating is retained; add the AutoSpectral-unmixed files to the groups in
#' FlowJo to apply it.
#'
#' Gate coordinates are not altered. They remain appropriate only if the
#' AutoSpectral-unmixed data are on the same scale as the original data.
#'
#' @param workspace.file File path to the FlowJo workspace (`.wsp`) or
#' template (`.wspt`) to convert.
#' @param fcs.file File path to an AutoSpectral-unmixed FCS file (an unmixed
#' sample or control) providing the target parameter names.
#' @param output.file File path for the converted workspace. Default `NULL`
#' writes `<workspace name>_AutoSpectral.wsp` next to `workspace.file`.
#' @param name.map Optional named character vector of manual matches, with
#' workspace parameter names as names and `fcs.file` parameter names as values,
#' e.g. `c("AF-A" = "AF Abundance")`. Use names without the compensation
#' prefix. Takes precedence over automatic matching.
#' @param comp.prefix The prefix FlowJo uses for compensated parameters.
#' Default `"Comp-"`.
#' @param verbose Logical. Whether to print progress messages. Default `TRUE`.
#'
#' @return Invisibly, a data frame with one row per parameter referenced in
#' the workspace: `workspace.name`, `fcs.name` (`NA` if unmatched), `marker`
#' (the `$PnS` in `fcs.file`) and `match` (how it was matched).
#'
#' @export

convert.flowjo.workspace <- function( workspace.file, fcs.file,
                                      output.file = NULL, name.map = NULL,
                                      comp.prefix = "Comp-", verbose = TRUE ) {

  if ( !file.exists( workspace.file ) )
    stop( "Workspace file not found: ", workspace.file, call. = FALSE )
  if ( !file.exists( fcs.file ) )
    stop( "FCS file not found: ", fcs.file, call. = FALSE )

  if ( is.null( output.file ) ) {
    output.file <- file.path(
      dirname( workspace.file ),
      paste0( sub( "\\.wspt?$", "", basename( workspace.file ),
                   ignore.case = TRUE ), "_AutoSpectral.wsp" )
    )
  }

  # target parameters; FCS keywords are case-insensitive
  keys <- readFCSheader( fcs.file )[[ 1 ]]
  names( keys ) <- toupper( names( keys ) )
  n.par <- as.integer( keys[[ "$PAR" ]] )
  fcs.names <- unname( keys[ paste0( "$P", seq_len( n.par ), "N" ) ] )
  fcs.markers <- unname( keys[ paste0( "$P", seq_len( n.par ), "S" ) ] )
  fcs.markers[ is.na( fcs.markers ) ] <- ""
  names( fcs.markers ) <- fcs.names

  # parameters FlowJo will compensate, from the spillover keyword
  spill.key <- intersect( c( "$SPILLOVER", "SPILL", "$SPILL", "SPILLOVER" ),
                          names( keys ) )
  spill.params <- character( 0 )
  if ( length( spill.key ) > 0 ) {
    spill <- trimws( strsplit( keys[[ spill.key[ 1 ] ]], ",", fixed = TRUE )[[ 1 ]] )
    n.spill <- as.integer( spill[ 1 ] )
    if ( !is.na( n.spill ) && n.spill > 0 && length( spill ) == 1 + n.spill + n.spill^2 ) {
      spill.matrix <- matrix( suppressWarnings( as.numeric( spill[ -seq_len( 1 + n.spill ) ] ) ),
                              nrow = n.spill, byrow = TRUE )
      # FlowJo ignores an identity spillover matrix and creates no compensated
      # parameters from it
      off.diagonal <- spill.matrix[ row( spill.matrix ) != col( spill.matrix ) ]
      if ( any( off.diagonal != 0, na.rm = TRUE ) )
        spill.params <- spill[ 1 + seq_len( n.spill ) ]
    }
  }

  doc <- xml2::read_xml( workspace.file )
  ns <- xml2::xml_ns( doc )

  # remove samples, which point to the original FCS files
  sample.nodes <- xml2::xml_find_all( doc, "//*[local-name()='SampleList']/*" )
  if ( length( sample.nodes ) > 0 ) {
    xml2::xml_remove( sample.nodes )
    xml2::xml_remove( xml2::xml_find_all( doc, "//*[local-name()='SampleRefs']/*" ) )
    xml2::xml_remove( xml2::xml_find_all(
      doc, "//*[local-name()='Matrices']/*[@matrixType='acquisitionDefined']"
    ) )
    if ( verbose )
      message( sprintf( "Removed %d sample(s) referring to the original FCS files.",
                        length( sample.nodes ) ) )
  }

  # attributes holding parameter references
  ref.nodes <- xml2::xml_find_all( doc, paste(
    "//*[local-name()='Axis' or local-name()='fcs-dimension' or local-name()='parameter' or local-name()='Statistic']",
    "//*[@heatMapStatParameter or @lastYAxisName or @parameterName or @userProvidedCompInfix or @*[local-name()='parameter']]",
    sep = " | "
  ) )
  is.ref <- function( node, a ) {
    el <- xml2::xml_name( node )
    attr.local <- sub( "^.*:", "", names( a ) )
    ( attr.local == "name" & el %in% c( "Axis", "fcs-dimension", "parameter" ) ) |
      ( attr.local == "id" & el == "Statistic" ) |
      attr.local %in% c( "heatMapStatParameter", "lastYAxisName", "parameterName",
                         "userProvidedCompInfix", "parameter" )
  }

  # the workspace parameter, without the compensation prefix or a marker
  ref.base <- function( value ) {
    value <- sub( " :: .*$", "", value )
    ifelse( startsWith( value, comp.prefix ),
            substring( value, nchar( comp.prefix ) + 1 ), value )
  }

  ref.values <- unlist( lapply( ref.nodes, function( node ) {
    a <- xml2::xml_attrs( node, ns = ns )
    a[ is.ref( node, a ) ]
  } ), use.names = FALSE )
  ws.names <- sort( unique( ref.base( ref.values[ ref.values != "" ] ) ) )

  # fluorophore synonym lookup: normalized name -> normalized fluorophore
  fluor.db <- utils::read.csv(
    system.file( "extdata", "fluorophore_database.csv", package = "AutoSpectral" ),
    colClasses = "character"
  )
  normalize <- function( x ) tolower( gsub( "[^A-Za-z0-9]", "", x ) )
  fluor.cols <- c( "fluorophore",
                   grep( "^synonym[0-9]+$", colnames( fluor.db ), value = TRUE ) )
  synonym <- unlist( lapply( fluor.cols, function( col ) {
    keep <- !is.na( fluor.db[[ col ]] ) & fluor.db[[ col ]] != ""
    stats::setNames( normalize( fluor.db$fluorophore[ keep ] ),
                     normalize( fluor.db[[ col ]][ keep ] ) )
  } ) )
  synonym <- synonym[ !duplicated( names( synonym ) ) ]
  canonical <- function( x ) {
    k <- normalize( sub( "-A$", "", x ) )
    ifelse( k %in% names( synonym ), synonym[ k ], k )
  }

  fcs.key <- canonical( fcs.names )
  mapping <- data.frame( workspace.name = ws.names, fcs.name = NA_character_,
                         match = "unmatched", stringsAsFactors = FALSE )

  for ( i in seq_along( ws.names ) ) {
    w <- ws.names[ i ]
    if ( !is.null( name.map ) && w %in% names( name.map ) ) {
      if ( !name.map[[ w ]] %in% fcs.names )
        stop( sprintf( "name.map value '%s' is not a parameter in the FCS file.",
                       name.map[[ w ]] ), call. = FALSE )
      mapping[ i, c( "fcs.name", "match" ) ] <- c( name.map[[ w ]], "name.map" )
    } else if ( w %in% fcs.names ) {
      mapping[ i, c( "fcs.name", "match" ) ] <- c( w, "identical" )
    } else {
      hit <- fcs.names[ fcs.key == canonical( w ) ]
      if ( length( hit ) == 1 ) {
        mapping[ i, c( "fcs.name", "match" ) ] <- c( hit, "fluorophore" )
      } else if ( length( hit ) > 1 ) {
        mapping$match[ i ] <- "ambiguous"
      }
    }
  }
  mapping$marker <- ifelse( is.na( mapping$fcs.name ), NA_character_,
                            fcs.markers[ mapping$fcs.name ] )
  mapping <- mapping[ , c( "workspace.name", "fcs.name", "marker", "match" ) ]
  new.name <- stats::setNames( mapping$fcs.name, mapping$workspace.name )

  # translate one attribute value, preserving the compensation prefix and marker
  convert.value <- function( value ) {
    has.marker <- grepl( " :: ", value, fixed = TRUE )
    ref <- sub( " :: .*$", "", value )
    is.comp <- startsWith( ref, comp.prefix )
    base <- if ( is.comp ) substring( ref, nchar( comp.prefix ) + 1 ) else ref
    target <- new.name[ base ]
    if ( is.na( target ) ) return( value )
    out <- if ( is.comp && target %in% spill.params ) paste0( comp.prefix, target ) else target
    if ( has.marker && fcs.markers[[ target ]] != "" )
      out <- paste0( out, " :: ", fcs.markers[[ target ]] )
    unname( out )
  }

  n.changed <- 0L
  for ( node in ref.nodes ) {
    a <- xml2::xml_attrs( node, ns = ns )
    edit <- is.ref( node, a )
    for ( key in names( a )[ edit ] ) {
      new.value <- convert.value( a[[ key ]] )
      if ( !identical( new.value, a[[ key ]] ) ) {
        xml2::xml_set_attr( node, key, new.value, ns = ns )
        n.changed <- n.changed + 1L
      }
    }
    # Table Editor columns carry the marker separately
    if ( "stain" %in% names( a ) && "parameterName" %in% names( a ) ) {
      target <- new.name[ ref.base( a[[ "parameterName" ]] ) ]
      if ( !is.na( target ) )
        xml2::xml_set_attr( node, "stain", fcs.markers[[ target ]] )
    }
  }

  unmatched <- mapping$workspace.name[ is.na( mapping$fcs.name ) ]
  if ( length( unmatched ) > 0 ) {
    warning(
      sprintf( paste( "%d workspace parameter(s) not matched to %s and left",
                      "unchanged: %s. Supply matches via `name.map`." ),
               length( unmatched ), basename( fcs.file ),
               paste( unmatched, collapse = ", " ) ),
      call. = FALSE
    )
  }

  xml2::write_xml( doc, output.file )

  if ( verbose ) {
    message( sprintf( "Matched %d of %d workspace parameters; updated %d references.",
                      sum( !is.na( mapping$fcs.name ) ), nrow( mapping ), n.changed ) )
    message( "Writing: ", output.file )
  }

  invisible( mapping )
}
