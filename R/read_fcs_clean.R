# read_fcs_clean.r

#' @title Read And Clean An FCS File
#'
#' @description
#' Reads the scatter and spectral channels of an FCS file and removes
#' events that saturate any spectral detector, events at the scatter
#' maxima (`asp$scatter.data.max.x`/`asp$scatter.data.max.y`, when set) and,
#' optionally, doublets. Doublets are removed with a two-pass Area/Height
#' scatter-ratio gate: the FSC ratio cutoff first, then the SSC ratio cutoff
#' computed only on events that passed the FSC cutoff. This is the cleaning
#' applied to every control file by `get.spectra.automated()`, exported so
#' that other packages can pre-process related samples identically.
#'
#' @param path File name and path of the FCS file.
#' @param label Character, name used in the progress message.
#' @param spectral.channels Character vector of spectral detector names.
#' @param scatter.channels Character vector of scatter channel names (Area).
#' Matching Height channels are read, if present, for doublet removal.
#' @param sat.value Numeric, the detector saturation value. Events with any
#' spectral channel at or above this value are removed. Typically
#' `asp$expr.data.max`.
#' @param singlet.quantiles Numeric, length 2. Quantile cutoffs for the FSC
#' and SSC Area/Height ratios, e.g. `c( 0.85, 0.975 )`.
#' @param remove.doublets Logical, default `TRUE`. Skipped automatically if
#' Height channels are not present in the file.
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param`.
#' @param verbose Logical, default `TRUE`.
#'
#' @return Numeric matrix of retained events, with the scatter channels and
#' spectral channels present in the file (Height channels dropped).
#'
#' @export

read.fcs.clean <- function(
    path,
    label,
    spectral.channels,
    scatter.channels,
    sat.value,
    singlet.quantiles,
    remove.doublets = TRUE,
    asp,
    verbose = TRUE
) {
  fsc.a <- asp$default.scatter.parameter[ 1L ]
  ssc.a <- asp$default.scatter.parameter[ 2L ]
  fsc.h <- sub( "-A$", "-H", fsc.a )
  ssc.h <- sub( "-A$", "-H", ssc.a )

  height.channels <- sub( "-A$", "-H", scatter.channels )
  cols.keep       <- c( scatter.channels, height.channels, spectral.channels )

  probe.cols <- colnames( readFCS( path, start.row = 1, end.row = 1 ) )
  present    <- intersect( cols.keep, probe.cols )
  mat        <- readFCS( path, columns = present )

  # -- remove spectral-saturating events
  spec.present <- intersect( spectral.channels, colnames( mat ) )
  if ( length( spec.present ) > 0 ) {
    keep <- rowSums( mat[ , spec.present, drop = FALSE ] >= sat.value ) == 0
    mat  <- mat[ keep, , drop = FALSE ]
  }

  # -- remove scatter-saturating events--we may not actually want this
  if ( fsc.a %in% colnames( mat ) && !is.null( asp$scatter.data.max.x ) )
    mat <- mat[ mat[ , fsc.a ] < asp$scatter.data.max.x, , drop = FALSE ]
  if ( ssc.a %in% colnames( mat ) && !is.null( asp$scatter.data.max.y ) )
    mat <- mat[ mat[ , ssc.a ] < asp$scatter.data.max.y, , drop = FALSE ]

  # -- remove doublets (two-pass scatter-ratio, mirrors flowstate::select_singlets)
  if ( remove.doublets && all( c( fsc.a, fsc.h ) %in% colnames( mat ) ) ) {
    fsc.ratio <- mat[ , fsc.a ] / ( mat[ , fsc.h ] + 1e-9 )
    mat       <- mat[ fsc.ratio < stats::quantile( fsc.ratio, probs = singlet.quantiles[ 1L ] ), ,
                      drop = FALSE ]

    if ( all( c( ssc.a, ssc.h ) %in% colnames( mat ) ) ) {
      ssc.ratio <- mat[ , ssc.a ] / ( mat[ , ssc.h ] + 1e-9 )
      mat       <- mat[ ssc.ratio < stats::quantile( ssc.ratio, probs = singlet.quantiles[ 2L ] ), ,
                        drop = FALSE ]
    }
  }

  # drop height channels -- not needed downstream
  mat <- mat[ , intersect( c( scatter.channels, spectral.channels ), colnames( mat ) ),
              drop = FALSE ]

  if ( verbose )
    message( sprintf( "\033[32m  %-40s  %d events\033[0m", label, nrow( mat ) ) )
  mat
}
