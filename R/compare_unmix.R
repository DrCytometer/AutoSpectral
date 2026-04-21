# compare_unmix.R

#' @title Compare Unmixing Quality Across Two Spectral References
#'
#' @description
#' A diagnostic utility for comparing unmixing performance when two different
#' spectral references are available for the same fluorophore — for example,
#' a cell-derived reference versus a bead-derived reference, or two different
#' lot measurements. Both references are used to unmix the same single-stained
#' FCS file and the resulting Secondary Stain Index (SSI) and secondary spillover
#' values are compared side-by-side.
#'
#' The function:
#' \enumerate{
#'   \item Reads the single-stained and unstained FCS files.
#'   \item Applies a FSC/SSC scatter gate (auto-detected or user-supplied).
#'   \item Unmixes both gated datasets with each spectral reference matrix using
#'     OLS (or WLS for the ID7000).
#'   \item Identifies the fluorophore channel with the largest absolute SSI
#'     across both unmixings and uses it as the y-axis for biplot visualisation.
#'   \item Annotates each biplot with the median ± rSD of the unstained
#'     distribution and the SSI / spillover value.
#'   \item Saves a side-by-side JPEG to `plot.dir` and returns a summary table.
#' }
#'
#' **Note:** This function calls `calculate.ssi()`, which must be available in
#' the package namespace.
#'
#' @param single.stained.fcs Character string. Path to the single-stained FCS
#' file for the fluorophore of interest.
#' @param unstained.fcs Character string. Path to the matched unstained (negative
#' control) FCS file.
#' @param fluorophore Character string. Name of the target fluorophore (must
#' match a row name in `spectra`, or will be appended automatically).
#' @param spectra Numeric matrix of spectral references (fluorophores × detectors,
#' values normalised 0–1) used as the starting point for constructing
#' `ref.spectra` and `test.spectra` when those are not supplied directly.
#' @param ref.spectrum Named numeric vector or single-row matrix. The reference
#' spectrum for `fluorophore` used in the first unmixing.
#' @param test.spectrum Named numeric vector or single-row matrix. The test
#' spectrum for `fluorophore` used in the second unmixing.
#' @param cytometer Character string identifying the cytometer model, passed to
#' [get.autospectral.param()]. Also determines the unmixing algorithm:
#' `"id7000"` uses WLS, all others use OLS.
#' @param ref.spectra Optional numeric matrix. A pre-built full spectral
#' reference matrix (fluorophores × detectors) for the reference condition.
#' When `NULL` (default), this is constructed from `spectra` and `ref.spectrum`.
#' @param test.spectra Optional numeric matrix. A pre-built full spectral
#' reference matrix for the test condition. When `NULL` (default), this is
#' constructed from `spectra` and `test.spectrum`.
#' @param x.channel Unused. Reserved for future use.
#' @param y.channel Unused. Reserved for future use.
#' @param x.min Numeric. Minimum x-axis value (in data units) passed to
#' [create.biplot()]. Default `-1000`.
#' @param y.min Numeric. Minimum y-axis value (in data units) passed to
#' [create.biplot()]. Default `-1000`.
#' @param x.width.basis Numeric. Width basis for the biexponential x-axis
#' transform passed to [create.biplot()]. Default `-1000`.
#' @param y.width.basis Numeric. Width basis for the biexponential y-axis
#' transform passed to [create.biplot()]. Default `-1000`.
#' @param gate Logical. Reserved for future gating control; currently the
#' scatter gate is always applied. Default `TRUE`.
#' @param gate.bound Optional list with elements `$x` and `$y` defining a
#' polygon gate boundary in FSC/SSC space. When `NULL` (default), the gate is
#' detected automatically via [do.gate()].
#' @param ref.label Character string. Display label for the reference condition
#' (e.g., `"Cells"`). Also passed as `control.type` to [do.gate()].
#' Default `"Cells"`.
#' @param test.label Character string. Display label for the test condition
#' (e.g., `"Beads"`). Default `"Beads"`.
#' @param biplot.width Numeric. Width of each individual biplot in inches.
#' The saved figure is twice this wide (two panels). Default `5`.
#' @param biplot.height Numeric. Height of the saved figure in inches.
#' Default `5`.
#' @param title Character string. Title displayed on the figure and used as
#' the JPEG filename stem. Default is constructed from `ref.label`,
#' `test.label`, and `fluorophore`.
#' @param plot.dir Character string. Directory for saving the output JPEG.
#' Created automatically if absent. Default `"./figure_compare_unmix"`.
#'
#' @return A data frame with one row per fluorophore in the reference spectra
#' matrix and four columns:
#' \describe{
#'   \item{`Fluorophore`}{Fluorophore name.}
#'   \item{`Reference.SSI`}{SSI under the reference spectrum.}
#'   \item{`Test.SSI`}{SSI under the test spectrum.}
#'   \item{`Reference.Spill`}{Fractional spillover under the reference spectrum.}
#'   \item{`Test.Spill`}{Fractional spillover under the test spectrum.}
#' }
#' The figure is saved to `plot.dir` and the combined plot is also printed to
#' the active graphics device.
#'
#' @importFrom ggplot2 geom_hline labs ggsave
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom sp point.in.polygon
#' @importFrom cowplot plot_grid
#' @importFrom ragg agg_jpeg
#'
#' @export

compare.unmix <- function(
    single.stained.fcs,
    unstained.fcs,
    fluorophore,
    spectra,
    ref.spectrum,
    test.spectrum,
    cytometer,
    ref.spectra        = NULL,
    test.spectra       = NULL,
    x.channel          = NULL,
    y.channel          = NULL,
    x.min              = -1000,
    y.min              = -1000,
    x.width.basis      = -1000,
    y.width.basis      = -1000,
    gate               = TRUE,
    gate.bound         = NULL,
    ref.label          = "Cells",
    test.label         = "Beads",
    biplot.width       = 5,
    biplot.height      = 5,
    title              = paste( ref.label, "vs", test.label, fluorophore ),
    plot.dir           = "./figure_compare_unmix"
  ) {

  # --- input validation ---
  if ( !file.exists( single.stained.fcs ) ) {
    stop(
      paste0( "File not found: '", single.stained.fcs, "'" ),
      call. = FALSE
    )
  }
  if ( !file.exists( unstained.fcs ) ) {
    stop(
      paste0( "File not found: '", unstained.fcs, "'" ),
      call. = FALSE
    )
  }

  # --- build spectral reference matrices if not supplied ---
  if ( is.null( ref.spectra ) ) {
    if ( fluorophore %in% rownames( spectra ) ) {
      ref.spectra                 <- spectra
      ref.spectra[ fluorophore, ] <- ref.spectrum
    } else {
      ref.spectra                                    <- rbind( spectra, ref.spectrum )
      rownames( ref.spectra )[ nrow( ref.spectra ) ] <- fluorophore
    }
  }

  if ( is.null( test.spectra ) ) {
    if ( fluorophore %in% rownames( spectra ) ) {
      test.spectra                 <- spectra
      test.spectra[ fluorophore, ] <- test.spectrum
    } else {
      test.spectra                                     <- rbind( spectra, test.spectrum )
      rownames( test.spectra )[ nrow( test.spectra ) ] <- fluorophore
    }
  }

  # --- condition numbers ---
  ref.cn  <- calculate.condition.number( ref.spectra )
  test.cn <- calculate.condition.number( test.spectra )
  message( "\033[34mReference spectra matrix condition number: ", ref.cn, "\033[0m" )
  message( "\033[34mTest spectra matrix condition number:      ", test.cn, "\033[0m" )

  # --- cytometer parameters ---
  asp <- get.autospectral.param( cytometer )

  # --- read FCS files ---
  ss.expr.data <- readFCS( single.stained.fcs )
  un.expr.data <- readFCS( unstained.fcs )

  # check that spectral channels are present in FCS
  missing.channels <- setdiff( colnames( spectra ), colnames( ss.expr.data ) )
  if ( length( missing.channels ) > 0 ) {
    stop(
      paste0(
        "The following spectral channels are absent from the FCS file: ",
        paste( missing.channels, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  # --- scatter gate ---
  gate.data <- ss.expr.data[ , asp$default.scatter.parameter ]

  if ( is.null( gate.bound ) ) {
    gate.bound <- do.gate(
      gate.data,
      viability.gate              = FALSE,
      large.gate                  = FALSE,
      samp                        = fluorophore,
      scatter.and.channel.label   = asp$default.scatter.parameter,
      control.type                = tolower( ref.label ),
      asp                         = asp
    )
  }

  # gate single-stained events
  in.gate.ss <- which(
    sp::point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      gate.bound$x, gate.bound$y
    ) != 0
  )
  ss.raw.data <- ss.expr.data[ in.gate.ss, colnames( spectra ), drop = FALSE ]

  # gate unstained events using the same boundary
  in.gate.un <- which(
    sp::point.in.polygon(
      un.expr.data[ , asp$default.scatter.parameter[ 1 ] ],
      un.expr.data[ , asp$default.scatter.parameter[ 2 ] ],
      gate.bound$x, gate.bound$y
    ) != 0
  )
  un.raw.data <- un.expr.data[ in.gate.un, colnames( spectra ), drop = FALSE ]

  # --- select unmixing algorithm ---
  unmix <- if ( cytometer == "id7000" ) unmix.wls else unmix.ols

  # --- unmix ---
  ss.ref.unmix  <- unmix( ss.raw.data, ref.spectra )
  ss.test.unmix <- unmix( ss.raw.data, test.spectra )
  un.ref.unmix  <- unmix( un.raw.data, ref.spectra )
  un.test.unmix <- unmix( un.raw.data, test.spectra )

  # --- SSI and spillover ---
  ref.error  <- calculate.ssi( un.ref.unmix,  ss.ref.unmix,  fluorophore )
  test.error <- calculate.ssi( un.test.unmix, ss.test.unmix, fluorophore )

  # auto-select y-channel as the one with the largest absolute SSI across both
  all.ssi    <- c( ref.error[ , "SSI" ], test.error[ , "SSI" ] )
  fluor.names <- c( rownames( ref.error ), rownames( test.error ) )
  max.idx    <- which.max( abs( all.ssi ) )
  y.channel  <- fluor.names[ max.idx ]
  message( "\033[34mChannel with the maximum SSI: ", y.channel, "\033[0m" )

  # --- axis limits ---
  x.plot.min <- min(
    stats::quantile( un.test.unmix[ , fluorophore ], 0.01 ) * 2,
    stats::quantile( un.ref.unmix[  , fluorophore ], 0.01 ) * 2,
    x.min
  )
  y.plot.min <- min(
    stats::quantile( ss.test.unmix[ , y.channel ], 0.01 ) * 2,
    stats::quantile( ss.ref.unmix[  , y.channel ], 0.01 ) * 2,
    y.min
  )

  # --- biplots ---
  ref.biplot <- create.biplot(
    ss.ref.unmix,
    x.dim         = fluorophore,
    y.dim         = y.channel,
    asp           = asp,
    x.min         = x.plot.min,
    y.min         = y.plot.min,
    x.width.basis = x.width.basis,
    y.width.basis = y.width.basis,
    save          = FALSE,
    title         = ref.label,
    width         = biplot.width,
    height        = biplot.height
  )

  test.biplot <- create.biplot(
    ss.test.unmix,
    x.dim         = fluorophore,
    y.dim         = y.channel,
    asp           = asp,
    x.min         = x.plot.min,
    y.min         = y.plot.min,
    x.width.basis = x.width.basis,
    y.width.basis = y.width.basis,
    save          = FALSE,
    title         = test.label,
    width         = biplot.width,
    height        = biplot.height
  )

  # --- biexponential y-axis transform for reference lines ---
  if ( y.width.basis < -1000 ) {
    y.excess.width.basis <- y.width.basis + 1000
    y.pos.log.delta      <- log10( abs( y.excess.width.basis ) )
    y.pos.log            <- log10( asp$expr.data.max ) - 1 - y.pos.log.delta
    y.pos.log            <- pmax( y.pos.log, 2 )
  } else {
    y.pos.log <- log10( asp$expr.data.max ) - 1
  }

  biexp.transform.y <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$expr.data.max,
    pos          = y.pos.log,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = y.width.basis,
    inverse      = FALSE
  )

  # --- annotate test biplot ---
  unstained.median.test <- stats::median( un.test.unmix[ , y.channel ] )
  unstained.rsd.test    <- stats::mad(    un.test.unmix[ , y.channel ] )
  test.ssi   <- round( test.error[ y.channel, "SSI" ],      2 )
  test.spill <- round( test.error[ y.channel, "Spillover" ], 2 ) * 100

  test.biplot <- test.biplot +
    ggplot2::geom_hline(
      yintercept = biexp.transform.y( unstained.median.test ),
      color      = "black"
    ) +
    ggplot2::geom_hline(
      yintercept = biexp.transform.y( unstained.median.test + unstained.rsd.test ),
      color      = "red"
    ) +
    ggplot2::geom_hline(
      yintercept = biexp.transform.y( unstained.median.test - unstained.rsd.test ),
      color      = "red"
    ) +
    ggplot2::labs(
      subtitle = paste( "SSI:", test.ssi, "  Spillover:", test.spill, "%" )
    )

  # --- annotate reference biplot ---
  unstained.median.ref <- stats::median( un.ref.unmix[ , y.channel ] )
  unstained.rsd.ref    <- stats::mad(    un.ref.unmix[ , y.channel ] )
  ref.ssi   <- round( ref.error[ y.channel, "SSI" ],       2 )
  ref.spill <- round( ref.error[ y.channel, "Spillover" ],  2 ) * 100

  ref.biplot <- ref.biplot +
    ggplot2::geom_hline(
      yintercept = biexp.transform.y( unstained.median.ref ),
      color      = "black"
    ) +
    ggplot2::geom_hline(
      yintercept = biexp.transform.y( unstained.median.ref + unstained.rsd.ref ),
      color      = "red"
    ) +
    ggplot2::geom_hline(
      yintercept = biexp.transform.y( unstained.median.ref - unstained.rsd.ref ),
      color      = "red"
    ) +
    ggplot2::labs(
      subtitle = paste( "SSI:", ref.ssi, "  Spillover:", ref.spill, "%" )
    )

  # --- combine and save ---
  combined.plot <- cowplot::plot_grid(
    ref.biplot,
    test.biplot,
    nrow  = 1,
    align = "v",
    axis  = "lr"
  )

  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

  ggplot2::ggsave(
    file.path( plot.dir, sprintf( "%s.jpg", title ) ),
    combined.plot,
    device    = ragg::agg_jpeg,
    width     = biplot.width * 2,
    height    = biplot.height,
    limitsize = FALSE
  )

  print( combined.plot )

  # --- result table ---
  result.table <- data.frame(
    Fluorophore    = rownames( ref.spectra ),
    Reference.SSI  = ref.error$SSI,
    Test.SSI       = test.error$SSI,
    Reference.Spill = ref.error$Spillover,
    Test.Spill     = test.error$Spillover
  )

  return( result.table )
}
