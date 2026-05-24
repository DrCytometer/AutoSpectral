# get_cytometer_param.R
#
# Data-driven cytometer parameter resolution.
# Reads the $CYT (and optionally $CYTSN / CREATOR) FCS keyword from a single
# representative file and returns a list of instrument-specific parameters
# used to identify fluorescence detector columns, scatter columns, the primary
# AF channel, the saturation ceiling, and the canonical detector ordering.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

## Extract $PnN parameter names from the keyword vector, in channel order.
.extract_pnn <- function( kw ) {
  par.count <- as.integer( kw[ "$PAR" ] )
  if ( is.na( par.count ) || par.count < 1L )
    stop( "$PAR keyword missing or invalid." )
  pnn.keys <- paste0( "$P", seq_len( par.count ), "N" )
  vals <- kw[ pnn.keys ]
  vals[ is.na( vals ) ] <- ""
  unname( vals )
}

## Null-coalescing helper
`.cyt_or` <- function( a, b ) {
  if ( !is.null( a ) && !is.na( a ) && nchar( trimws( a ) ) > 0 ) a else b
}


# ---------------------------------------------------------------------------
# Per-cytometer parameter table
# ---------------------------------------------------------------------------
# Recognition logic uses $CYT first (where values are confirmed from real data),
# then CREATOR as a fallback.
#
# Confirmed $CYT values (from real instrument files):
#   Cytek Aurora / NL : "Aurora"             (distinguished later by UV detectors)
#   Sony ID7000       : "ID7000"
#   BD FACSDiscover   : "FACSDiscover A8"     (A8; S8 expected similar)
#   BD FACSymphony    : "FACSymphony Analyzer"
#   NovoCyte Opteon   : "NovoCyte Opteon"
#   CytoFLEX Mosaic   : "CytoFLEX mosaic 88 & CytoFLEX LX"
#   Attune Xenith     : "Attune Xenith"

.CYTOMETER_PARAMS <- list(

  Aurora = list(
    cyt.label        = "Cytek Aurora",
    cyt.kw.pattern   = "^Aurora$",
    creator.pattern  = "SpectroFlo",
    scatter.param    = c( "FSC-A", "SSC-A" ),
    non.spectral.pat = c( "FSC", "SSC", "Time", "-H$", "-W$" ),
    sat.value        = 4194304L,
    af.channel       = "V7-A",
    db.col           = "Aurora"
  ),

  NorthernLights = list(
    cyt.label        = "Cytek Northern Lights",
    cyt.kw.pattern   = "^Aurora$",          # same $CYT as Aurora; use UV presence
    creator.pattern  = "SpectroFlo",
    scatter.param    = c( "FSC-A", "SSC-A" ),
    non.spectral.pat = c( "FSC", "SSC", "Time", "-H$", "-W$" ),
    sat.value        = 4194304L,
    af.channel       = "V7-A",
    db.col           = "NorthernLights"
  ),

  ID7000 = list(
    cyt.label        = "Sony ID7000",
    cyt.kw.pattern   = "ID7000",
    creator.pattern  = NULL,
    scatter.param    = c( "FSC-A", "SSC-A" ),
    non.spectral.pat = c( "FSC", "SSC", "TIME", "-H$", "-W$" ),
    sat.value        = 1048576L,
    af.channel       = "405CH7-A",
    db.col           = "ID7000"
  ),

  FACSDiscover = list(
    cyt.label        = "BD FACSDiscover (S8 / A8)",
    cyt.kw.pattern   = "FACSDiscover",
    creator.pattern  = "FACSDiva|FACSSuite",
    scatter.param    = c( "FSC-A", "SSC (Violet)-A" ),
    non.spectral.pat = c(
      "FSC", "SSC", "Time", "LightLoss", "Delta", "Plate",
      "Radial", "Correlation", "Intensity", "Eccentricity",
      "Diffusivity", "Center", "Moment", "Size", "Saturated",
      "Sorted", "Row", "Column", "Img", "Protocol", "EventLabel",
      "Region", "Gate", "Index", "Phase", "Event", "Drop",
      "Spectral", "Waveform", "Merged", "Flow", "Packet", "Reserved",
      "-H$", "-W$", "-T$"
    ),
    spectral.pat     = "\\([0-9]+\\)-A$",
    sat.value        = 24140237L,
    af.channel       = "V6 (515)-A",
    db.col           = "Discover"
  ),

  Opteon = list(
    cyt.label        = "Agilent NovoCyte Opteon",
    cyt.kw.pattern   = "NovoCyte|Opteon",
    creator.pattern  = "NovoExpress|Opteon",
    scatter.param    = c( "FSC-A", "VSSC-A" ),
    non.spectral.pat = c( "FSC", "SSC", "VSSC", "Time", "-H$", "-W$", "Width" ),
    sat.value        = 16777216L,
    af.channel       = "UV508-A",
    db.col           = "Opteon"
  ),

  Mosaic = list(
    cyt.label        = "Beckman Coulter CytoFLEX Mosaic",
    cyt.kw.pattern   = "CytoFLEX",
    creator.pattern  = "CytExpert",
    scatter.param    = c( "FSC-A", "BSSC-A" ),
    non.spectral.pat = c( "FSC", "SSC", "BSSC", "Time", "-H$" ),
    sat.value        = 16777215L,
    af.channel       = "V8-A",
    db.col           = "Mosaic"
  ),

  Xenith = list(
    cyt.label        = "ThermoFisher Attune Xenith",
    cyt.kw.pattern   = "Attune|Xenith",
    creator.pattern  = "Attune|VitesseSQ",
    scatter.param    = c( "FSC51-A", "SSC52-A" ),
    non.spectral.pat = c(
      "FSC", "SSC", "Time", "Event", "Gate", "Sort",
      "Comp", "-H$", "-W$"
    ),
    sat.value        = 100000L,
    af.channel       = "FL13-A",
    db.col           = "Xenith"
  ),

  Symphony = list(
    cyt.label        = "BD FACSymphony A5 SE",
    cyt.kw.pattern   = "FACSymphony",
    creator.pattern  = "FACSDiva|FACSSuite",
    scatter.param    = c( "FSC-A", "SSC-A" ),
    non.spectral.pat = c( "FSC", "SSC", "Time", "-H$", "-W$" ),
    sat.value        = 262144L,
    af.channel       = "UV515-A",
    db.col           = "A5SE"
  )
)


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

#' @title Resolve Cytometer Parameters from an FCS File
#'
#' @description
#' Reads the FCS keyword metadata from a single representative file and
#' returns a list of instrument-specific parameters needed to correctly identify
#' fluorescence detector columns, scatter columns, the primary autofluorescence
#' channel, the saturation ceiling, and the canonical detector ordering.
#'
#' Recognition is fully data-driven: no cytometer name needs to be supplied by
#' the user. The function inspects `$CYT`, `CREATOR`, and the set of `$PnN`
#' parameter names encoded in the FCS TEXT segment.
#'
#' Confirmed `$CYT` values (from real instrument FCS files):
#' * Cytek Aurora / Northern Lights: `"Aurora"` (distinguished by UV detectors)
#' * Sony ID7000: `"ID7000"`
#' * BD FACSDiscover A8: `"FACSDiscover A8"`
#' * BD FACSymphony A5 SE: `"FACSymphony Analyzer"`
#' * NovoCyte Opteon: `"NovoCyte Opteon"`
#' * CytoFLEX Mosaic: `"CytoFLEX mosaic 88 & CytoFLEX LX"`
#' * Attune Xenith: `"Attune Xenith"`
#'
#' @param fcs.path Character scalar. Path to a single raw `.fcs` file
#'   representative of the experiment.
#' @param db.path Character scalar or `NULL`. Path to `cytometer_database.csv`.
#'   When `NULL` the function looks in the installed AutoSpectral package;
#'   if absent, detector ordering follows FCS acquisition order.
#' @param verbose Logical, default `TRUE`. Print a message confirming which
#'   cytometer was identified.
#'
#' @return A named list:
#' \describe{
#'   \item{`cyt.label`}{Human-readable cytometer name.}
#'   \item{`scatter.param`}{Character vector of FSC/SSC column names.}
#'   \item{`cols.detector`}{Fluorescence detector column names in
#'     excitation/emission order.}
#'   \item{`af.channel`}{Primary autofluorescence detector name.}
#'   \item{`sat.value`}{Numeric raw saturation ceiling.}
#'   \item{`db.col`}{Column name used in `cytometer_database.csv`.}
#'   \item{`keywords`}{Full keyword list from the FCS TEXT segment.}
#' }
#'
#' @examples
#' \dontrun{
#' p <- get.cytometer.param( "path/to/control.fcs" )
#' p$cyt.label      # "Cytek Aurora"
#' p$af.channel     # "V7-A"
#' p$cols.detector  # c("UV1-A", "UV2-A", ..., "R8-A")
#' }
#'
#' @export

get.cytometer.param <- function(
    fcs.path,
    db.path = NULL,
    verbose = TRUE
) {

  if ( !file.exists( fcs.path ) )
    stop( "File not found: ", fcs.path, call. = FALSE )

  kw           <- readFCSheader( fcs.path )[[ 1 ]]
  names( kw )  <- toupper( names( kw ) )
  all.par      <- .extract_pnn( kw )

  cyt.kw          <- `.cyt_or`( kw[ "$CYT" ],       "" )
  creator.kw      <- `.cyt_or`( kw[ "CREATOR" ],    "" )
  creator.kw2     <- `.cyt_or`( kw[ "$CREATOR" ],   "" )
  creator.combined <- paste( creator.kw, creator.kw2 )

  matched <- .match_cytometer( all.par, cyt.kw, creator.combined )

  if ( verbose )
    message( "\033[34mCytometer identified: \033[32m", matched$cyt.label, "\033[0m" )

  cols.detector <- .derive_detector_cols( all.par, matched )
  cols.detector <- .order_detectors( cols.detector, matched$db.col, db.path )

  if ( !is.null( matched$af.channel ) &&
       !matched$af.channel %in% cols.detector ) {
    warning(
      "Primary AF channel '", matched$af.channel, "' not found in this file's ",
      "detector list. AF characterisation will fall back to the channel with the ",
      "highest mean signal across the unstained control.",
      call. = FALSE
    )
    matched$af.channel <- NULL
  }

  list(
    cyt.label     = matched$cyt.label,
    scatter.param = .resolve_scatter( all.par, matched$scatter.param ),
    cols.detector = cols.detector,
    af.channel    = matched$af.channel,
    sat.value     = matched$sat.value,
    db.col        = matched$db.col,
    keywords      = kw
  )
}


# ---------------------------------------------------------------------------
# Internal matching logic
# ---------------------------------------------------------------------------

## Match FCS keywords / parameter names to a .CYTOMETER_PARAMS entry.
## Uses confirmed $CYT values as the primary discriminant; falls back to
## CREATOR for instruments where $CYT is not reliable.
.match_cytometer <- function( all.par, cyt.kw, creator.combined ) {

  ## --- Priority 1: $CYT-based matching (confirmed values) ---

  ## Sony ID7000
  if ( grepl( "ID7000", cyt.kw, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$ID7000 )

  ## BD FACSDiscover (A8/S8)
  if ( grepl( "FACSDiscover", cyt.kw, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$FACSDiscover )

  ## BD FACSymphony (A5 SE) — must come before generic FACSDiva CREATOR check
  if ( grepl( "FACSymphony", cyt.kw, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$Symphony )

  ## Agilent NovoCyte Opteon
  if ( grepl( "NovoCyte|Opteon", cyt.kw, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$Opteon )

  ## Beckman CytoFLEX Mosaic
  if ( grepl( "CytoFLEX", cyt.kw, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$Mosaic )

  ## ThermoFisher Attune Xenith
  if ( grepl( "Attune|Xenith", cyt.kw, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$Xenith )

  ## Cytek Aurora / Northern Lights — $CYT = "Aurora" for both
  ## Distinguish by presence of UV detector channels
  if ( grepl( "^Aurora$", cyt.kw, ignore.case = TRUE ) ) {
    has.uv <- any( grepl( "^UV[0-9]+-A$", all.par ) )
    return( if ( has.uv ) .CYTOMETER_PARAMS$Aurora else .CYTOMETER_PARAMS$NorthernLights )
  }

  ## --- Priority 2: CREATOR-based fallback ---

  ## BD instruments via FACSDiva / FACSSuite
  if ( grepl( "FACSDiva|FACSSuite", creator.combined, ignore.case = TRUE ) ) {
    has.discover.style <- any( grepl( "\\([0-9]+\\)-A$", all.par ) )
    return(
      if ( has.discover.style ) .CYTOMETER_PARAMS$FACSDiscover
      else                      .CYTOMETER_PARAMS$Symphony
    )
  }

  ## Cytek via SpectroFlo (Aurora/NL)
  if ( grepl( "SpectroFlo", creator.combined, ignore.case = TRUE ) ) {
    has.uv <- any( grepl( "^UV[0-9]+-A$", all.par ) )
    return( if ( has.uv ) .CYTOMETER_PARAMS$Aurora else .CYTOMETER_PARAMS$NorthernLights )
  }

  ## Opteon (NovoExpress CREATOR)
  if ( grepl( "NovoExpress", creator.combined, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$Opteon )

  ## Mosaic (CytExpert CREATOR)
  if ( grepl( "CytExpert", creator.combined, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$Mosaic )

  ## Xenith (Attune / VitesseSQ CREATOR)
  if ( grepl( "Attune|VitesseSQ|Xenith", creator.combined, ignore.case = TRUE ) )
    return( .CYTOMETER_PARAMS$Xenith )

  stop(
    "Could not identify cytometer from FCS keywords.\n",
    "  $CYT     : '", cyt.kw, "'\n",
    "  CREATOR  : '", creator.combined, "'\n",
    "Please open an issue at https://github.com/DrCytometer/AutoSpectral/issues ",
    "or supply cytometer parameters manually.",
    call. = FALSE
  )
}


## Derive fluorescence detector column names from the full $PnN list.
.derive_detector_cols <- function( all.par, matched ) {

  ## Exclusion pass: remove channels whose name contains any non-spectral pattern
  excl.pat <- paste( matched$non.spectral.pat, collapse = "|" )
  cols     <- all.par[ !grepl( excl.pat, all.par, perl = TRUE ) ]
  cols     <- cols[ nchar( cols ) > 0L ]

  ## Positive-match pass (FACSDiscover): additionally keep only spectral channels
  ## matching "XX (NNN)-A" to exclude any remaining imaging channels
  if ( !is.null( matched$spectral.pat ) ) {
    cols <- cols[ grepl( matched$spectral.pat, cols, perl = TRUE ) ]
    if ( length( cols ) == 0L )
      stop(
        "No spectral detector columns found using pattern '",
        matched$spectral.pat, "'. ",
        "Check that the FCS file is raw/unprocessed.",
        call. = FALSE
      )
    return( cols )
  }

  if ( length( cols ) == 0L )
    stop(
      "No fluorescence detector columns identified. ",
      "Check that the FCS file is raw/unprocessed and from a supported cytometer.",
      call. = FALSE
    )
  cols
}


## Order detectors using cytometer_database.csv canonical sequence.
.order_detectors <- function( cols.detector, db.col, db.path ) {

  if ( is.null( db.path ) ) {
    db.path <- tryCatch(
      system.file( "extdata", "cytometer_database.csv",
                   package = "AutoSpectral", mustWork = TRUE ),
      error = function( e ) ""
    )
  }

  if ( !nchar( db.path ) || !file.exists( db.path ) ) {
    warning(
      "cytometer_database.csv not found; detector ordering will follow ",
      "acquisition order from the FCS file.",
      call. = FALSE
    )
    return( cols.detector )
  }

  db        <- utils::read.csv( db.path, stringsAsFactors = FALSE )
  ref.order <- db[[ db.col ]]
  ref.order <- ref.order[ nchar( trimws( ref.order ) ) > 0L ]

  ordered   <- intersect( ref.order, cols.detector )
  remainder <- setdiff( cols.detector, ref.order )

  c( ordered, remainder )
}


## Verify canonical scatter channel names are present; fall back to grep if not.
.resolve_scatter <- function( all.par, canonical ) {
  present <- canonical[ canonical %in% all.par ]
  if ( length( present ) == length( canonical ) )
    return( present )

  fallback <- grep( "^FSC|^SSC|^BSSC|^VSSC", all.par, value = TRUE )
  fallback  <- grep( "-A$", fallback, value = TRUE )

  if ( length( fallback ) == 0L )
    warning(
      "No scatter channels found; scatter-matching will be unavailable.",
      call. = FALSE
    )

  fallback
}
