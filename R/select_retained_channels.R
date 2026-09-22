# select_retained_channels.r

#' @title Select Retained Channels
#'
#' @description
#' Internal helper for `unmix.fcs()`. Selects the parameters of a raw FCS file
#' that are carried through unchanged into the unmixed file, alongside the newly
#' unmixed fluorophore parameters.
#'
#' Some instruments (e.g. the BD FACSDiscover family and the FACSymphony A5 SE)
#' write "raw" FCS files that also hold the instrument's own unmixed
#' fluorophore parameters. These are never carried through: they would
#' duplicate, and could be mistaken for, the AutoSpectral output. A parameter
#' is retained only if it matches one of the named parameter families in
#' `asp$non.spectral.channel` (time, scatter, imaging and other acquisition
#' parameters). Bare suffix entries in that vector (`-H`, `-W`, `-T`) exist to
#' keep detector companions out of the spectral channel set and are not used to
#' retain anything, since fluorophore names such as `APC-H7` would otherwise
#' match them.
#'
#' Raw detector data are written only when `include.raw = TRUE`. For each
#' detector used in unmixing this is the `-A` (or `-H`) parameter that was
#' unmixed, plus its `-T` companion where the file has one. The `-H` and `-W`
#' companions of the detectors are never written.
#'
#' @param original.param Character vector of the parameter names (`$PnN`) of the
#' raw FCS file, in file order.
#' @param spectral.channel Character vector of the detector channels used for
#' unmixing (the column names of the spectra matrix).
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param`.
#' @param include.raw Logical, whether to retain the raw detector data. Default
#' is `FALSE`.
#' @param include.imaging Logical, whether to retain imaging parameters. Only
#' applies to the BD FACSDiscover family, where `FALSE` restricts the retained
#' parameters to time and scatter. Default is `TRUE`.
#' @param verbose Logical, whether to report the parameters that are not carried
#' through. Default is `TRUE`.
#'
#' @return Character vector of the parameter names to carry through, in the
#' order they are to be written.
#'
#' @keywords internal

.select.retained.channels <- function(
    original.param,
    spectral.channel,
    asp,
    include.raw = FALSE,
    include.imaging = TRUE,
    verbose = TRUE
) {

  # other acquisition suffixes recorded for the detectors used in unmixing
  # (a CytoStellar spectra matrix may be built on -H rather than -A)
  detector.base <- sub(
    "-A$|-H$", "",
    spectral.channel[ grepl( "-A$|-H$", spectral.channel ) ]
  )
  raw.companion <- intersect(
    original.param,
    as.vector( outer( detector.base, c( "-A", "-H", "-W", "-T" ), paste0 ) )
  )
  raw.channel <- union( spectral.channel, raw.companion )
  is.raw <- original.param %in% raw.channel

  # named parameter families only; bare suffix entries are excluded
  family.pattern <- asp$non.spectral.channel
  family.pattern <- family.pattern[ !grepl( "^-[A-Za-z]\\$?$", family.pattern ) ]

  if ( length( family.pattern ) > 0 ) {
    is.family <- grepl( paste0( family.pattern, collapse = "|" ), original.param )
  } else {
    is.family <- !is.raw
  }

  retained <- original.param[ is.family & !is.raw ]
  excluded <- original.param[ !is.family & !is.raw ]

  if ( verbose && length( excluded ) > 0 ) {
    shown <- utils::head( excluded, 10 )
    message(
      sprintf(
        "Original parameters not carried through to the unmixed file (%d): %s%s",
        length( excluded ),
        paste( shown, collapse = ", " ),
        if ( length( excluded ) > length( shown ) )
          sprintf( " and %d more", length( excluded ) - length( shown ) ) else ""
      )
    )
  }

  # imaging parameters are optional on the Discover; time and scatter stay
  if ( grepl( "Discover", asp$cytometer ) && !include.imaging ) {
    keep <- union( asp$time.and.scatter, asp$default.time.parameter )
    retained <- intersect( retained, keep )
  }

  # raw detector data and their -T companions, only on request
  if ( include.raw ) {
    retained <- c(
      retained,
      spectral.channel,
      raw.companion[ grepl( "-T$", raw.companion ) ]
    )
  }

  return( retained )
}
