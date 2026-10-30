# wrong_background_biplots.R
#
# Panel e of Figure 3: one single-stained cell control unmixed in a full panel
# (library spectra + one AF spectrum) with its own spectrum extracted after
# subtracting the CORRECT background (unstained cells in the same landmark
# gate) or the WRONG background (unstained cells in another cell type's gate).
# Biplots of the stained channel against (1) AF and (2) the panel fluorophore
# most collinear with the spectral error, with create.biplot()-style
# spillover-spread boundaries (calculate.residual.spread(), correct matrix)
# and the secondary stain index (calculate.ssi()) of the positive population.

#' Build the unmixing matrix: panel library spectra + AF, stained dye replaced
build.panel.matrix <- function( library.spectra, panel.fluors, af.spectrum, stained, spectrum ) {
  ch <- colnames( library.spectra )
  m <- rbind( library.spectra[ panel.fluors, ch, drop = FALSE ], AF = af.spectrum[ ch ] )
  m[ stained, ] <- spectrum[ ch ]
  m
}

#' @param stained.raw Gated raw events of the stained control (scatter + detectors).
#' @param unstained.raw Unstained events in the same gate.
#' @param matrices list(correct = , incorrect = ) unmixing matrices.
#' @param fluorophore Stained dye (row name).
#' @param peak.channel Its peak detector.
#' @param y.channels Character vector of unmixed channels for the y axes.
#' @param marker,wrong.background Labels.
#' @param positive.quantile,positive.multiplier Positive events: unmixed value
#'   in `fluorophore` above multiplier x quantile of the unmixed unstained
#'   (correct matrix), so that the same events are shown in every plot.
#' @return list(plots = 2 x length(y.channels) named list, ssi = data.frame,
#'   spread = list, condition.number, n.positive)
wrong.background.biplots <- function(
    stained.raw, unstained.raw, matrices, fluorophore, peak.channel, y.channels,
    asp, marker, wrong.background,
    positive.quantile = 0.995, positive.multiplier = 1,
    positive.color = "#21908C", base.title.size = NULL
) {
  ch <- colnames( matrices$correct )

  unmixed <- lapply( matrices, function( m ) list(
    stained   = unmix.ols( stained.raw[ , ch, drop = FALSE ], m ),
    unstained = unmix.ols( unstained.raw[ , ch, drop = FALSE ], m )
  ) )

  thr <- stats::quantile( unmixed$correct$unstained[ , fluorophore ], positive.quantile ) *
    positive.multiplier
  keep <- which( unmixed$correct$stained[ , fluorophore ] > thr )
  if ( length( keep ) < 20 ) stop( "Only ", length( keep ), " positive events for the biplots.", call. = FALSE )

  # SSI under each matrix (calculate.ssi() selects its own top 500 positives)
  ssi <- do.call( rbind, lapply( names( matrices ), function( e ) {
    res <- suppressWarnings( calculate.ssi(
      unmixed[[ e ]]$unstained, unmixed[[ e ]]$stained[ keep, , drop = FALSE ], fluorophore ) )
    data.frame( matrix = e, channel = y.channels,
                SSI = res[ y.channels, "SSI" ], Spillover = res[ y.channels, "Spillover" ] )
  } ) )

  # spillover spread: Residual Model on the correct matrix
  spread.fit <- calculate.residual.spread(
    stained.data = stained.raw, unstained.data = unstained.raw,
    spectra = matrices$correct, fluorophore = fluorophore,
    peak.channel = peak.channel, asp = asp
  )

  # axis limits shared by the two matrices
  x.all <- unlist( lapply( unmixed, function( u ) u$unstained[ , fluorophore ] ) )
  x.min <- min( -1000, 2 * stats::quantile( x.all, 0.01 ) )

  entry.title <- c(
    correct   = sprintf( "%s bg (correct)", marker ),
    incorrect = sprintf( "%s bg (incorrect)", wrong.background )
  )

  plots <- list()
  for ( y in y.channels ) {
    y.all <- unlist( lapply( unmixed, function( u )
      c( u$unstained[ , y ], u$stained[ keep, y ] ) ) )
    y.min <- min( -1000, 2 * stats::quantile( y.all, 0.01 ) )
    sp <- list(
      ss       = spread.fit$spread[[ y ]],
      flat.pos = unname( stats::quantile( unmixed$correct$unstained[ , y ], 0.995 ) ),
      flat.neg = unname( stats::quantile( unmixed$correct$unstained[ , y ], 0.005 ) )
    )
    for ( e in names( matrices ) ) {
      s <- ssi$SSI[ ssi$matrix == e & ssi$channel == y ]
      plots[[ paste( e, y, sep = "|" ) ]] <- create.brightness.biplot(
        populations = list(
          Negative = unmixed[[ e ]]$unstained,
          Positive = unmixed[[ e ]]$stained[ keep, , drop = FALSE ]
        ),
        x.dim = fluorophore,
        y.dim = y,
        asp = asp,
        population.colors = c( Negative = "black", Positive = positive.color ),
        legend.labels = c(
          Negative = sprintf( "Unstained (%s gate)", marker ),
          Positive = sprintf( "%s+  (SSI %s)", marker, formatC( s, digits = 1, format = "f" ) )
        ),
        spread = sp,
        x.min = x.min,
        y.min = y.min,
        y.lab = if ( y == "AF" ) "AF (unmixed)" else y,
        title = paste0( fluorophore, ", ", entry.title[[ e ]] )
      )
    }
  }

  list(
    plots = plots,
    ssi = ssi,
    spread = spread.fit,
    condition.number = vapply( matrices, calculate.condition.number, numeric( 1 ) ),
    n.positive = length( keep ),
    threshold = unname( thr )
  )
}
