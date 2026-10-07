# unmix_titration.R
#
# Reads a titration series of single-stained controls plus the matched
# unstained control, applies one scatter gate to every file, unmixes each file
# with each supplied spectral matrix and selects the positive events.
# Secondary Stain Index (SSI) values are calculated with AutoSpectral's
# calculate.ssi().

#' @param stained.files Named character vector (names = antibody amount),
#'   ordered from dimmest to brightest.
#' @param unstained.file Path of the unstained control.
#' @param spectra.list Named list of full-panel spectral matrices, e.g.
#'   list(cells = <matrix with cell-derived spectrum>, beads = <matrix with
#'   bead-derived spectrum>). All must have the same detector columns.
#' @param fluorophore Stained fluorophore (row name in every matrix).
#' @param asp AutoSpectral parameter list.
#' @param gate Scatter gate polygon (list with x, y).
#' @param unmix.fun Unmixing function (default `unmix.ols`).
#' @param positive.quantile Events are plotted/analysed only if the unmixed
#'   value in `fluorophore` is above this quantile of the unmixed unstained
#'   control (0.995 = 99.5th percentile), multiplied by `positive.multiplier`.
#' @param positive.multiplier Numeric, default 1. Multiplication factor for
#'   the quantile threshold (e.g. 2 = twice the 99.5th percentile). Note that a
#'   negative quantile value would move towards more negative values.
#' @param positive.space Which unmixing defines positivity: the name of an
#'   entry of `spectra.list` (so that the SAME events are shown under every
#'   unmixing, default the first entry), or "each" to use the threshold of the
#'   unmixing being plotted.
#' @param y.from Name of the entry of `spectra.list` under which the channel
#'   with the largest |SSI| in the brightest file is chosen as y-axis
#'   (default the last entry, i.e. the incorrect/bead-derived unmixing).
#' @return list(
#'   raw           = gated raw events per file (scatter + spectral channels),
#'   unmixed       = unmixed events, list per spectra.list entry, per file,
#'   keep          = indices of positive events, per entry and file,
#'   y.channel     = channel with largest |SSI|,
#'   ssi           = data.frame (entry, amount, n.events, SSI, Spillover),
#'   condition.number = condition number of each matrix)
unmix.titration <- function(
    stained.files,
    unstained.file,
    spectra.list,
    fluorophore,
    asp,
    gate,
    unmix.fun = unmix.ols,
    positive.quantile = 0.995,
    positive.multiplier = 1,
    positive.space = names( spectra.list )[ 1 ],
    y.from = names( spectra.list )[ length( spectra.list ) ]
) {

  spec.ch <- colnames( spectra.list[[ 1 ]] )
  scatter <- asp$default.scatter.parameter
  files <- c( stained.files, Unstained = unstained.file )

  # read, remove saturated events, apply the scatter gate
  raw <- lapply( files, function( f ) {
    d <- readFCS( f, columns = union( scatter, spec.ch ) )
    d <- d[ rowSums( d[ , spec.ch, drop = FALSE ] >= asp$expr.data.max ) == 0, , drop = FALSE ]
    apply.gate( d, gate, scatter.param = scatter, asp = asp )
  } )
  names( raw ) <- names( files )

  # unmix every file with every matrix
  unmixed <- lapply( spectra.list, function( spectra ) {
    lapply( raw, function( d ) unmix.fun( d[ , spec.ch, drop = FALSE ], spectra ) )
  } )

  # positive events
  keep.for <- function( space, entry ) {
    thr <- stats::quantile( unmixed[[ space ]][[ "Unstained" ]][ , fluorophore ], positive.quantile ) *
      positive.multiplier
    message( sprintf( "Positivity threshold (%s space) for '%s': %.1f", space, entry, thr ) )
    lapply( unmixed[[ entry ]][ names( stained.files ) ], function( u ) which( u[ , fluorophore ] > thr ) )
  }
  keep <- lapply( names( spectra.list ), function( entry ) {
    space <- if ( identical( positive.space, "each" ) ) entry else positive.space
    keep.for( space, entry )
  } )
  names( keep ) <- names( spectra.list )

  # y-axis channel: largest |SSI| in the brightest file under `y.from`
  brightest <- names( stained.files )[ length( stained.files ) ]
  ssi.bright <- calculate.ssi(
    unmixed[[ y.from ]][[ "Unstained" ]],
    unmixed[[ y.from ]][[ brightest ]][ keep[[ y.from ]][[ brightest ]], , drop = FALSE ],
    fluorophore
  )
  y.channel <- rownames( ssi.bright )[ which.max( abs( ssi.bright$SSI ) ) ]
  message( "Channel with the largest |SSI| under '", y.from, "' unmixing: ", y.channel )

  # SSI for every population at that channel
  ssi <- do.call( rbind, lapply( names( spectra.list ), function( entry ) {
    do.call( rbind, lapply( names( stained.files ), function( a ) {
      idx <- keep[[ entry ]][[ a ]]
      if ( length( idx ) < 2 )
        return( data.frame( entry = entry, amount = a, n.events = length( idx ),
                            SSI = NA_real_, Spillover = NA_real_ ) )
      res <- suppressWarnings( calculate.ssi(
        unmixed[[ entry ]][[ "Unstained" ]],
        unmixed[[ entry ]][[ a ]][ idx, , drop = FALSE ],
        fluorophore
      ) )
      data.frame( entry = entry, amount = a, n.events = length( idx ),
                  SSI = res[ y.channel, "SSI" ], Spillover = res[ y.channel, "Spillover" ] )
    } ) )
  } ) )
  rownames( ssi ) <- NULL

  list(
    raw = raw,
    unmixed = unmixed,
    keep = keep,
    y.channel = y.channel,
    ssi = ssi,
    condition.number = vapply( spectra.list, calculate.condition.number, numeric( 1 ) )
  )
}
