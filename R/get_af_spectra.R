# get_af_spectra.r

#' @title Get Autofluorescence Spectra
#'
#' @description
#' Extracts autofluorescence spectra from an unstained samples. Intended for use
#' with `unmix.autospectral`. Uses FlowSOM (EmbedSOM) clustering for rapid
#' identification of cells with similar AF profiles.
#'
#' @importFrom EmbedSOM SOM
#' @importFrom stats na.omit
#' @importFrom flowCore read.FCS exprs
#' @importFrom utils write.csv
#'
#' @param unstained.sample Path and file name for a unstained sample FCS file.
#' The sample type and processing (protocol) method should match the fully
#' stained samples to which the AF will be applied, ideally.
#' @param asp The AutoSpectral parameter list.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param threads Numeric. Number of threads to use for parallel processing in
#' the creation of the SOM. Default `NULL` reverts to `asp$worker.process.n`.
#' @param som.dim Number of x and y dimensions for the SOM. Default is `10`.
#' @param figures Logical, whether to plot the spectral traces and heatmap for
#' the AF signatures. Default is `TRUE`.
#' @param plot.dir Directory (folder) where the plots will be saved. Default is
#' `NULL`, which inherits from `asp$figure.af.dir`.
#' @param table.dir Directory (folder) where the spectra csv file will be saved.
#' Default is `NULL`, which inherits from `asp$table.af.dir`.
#' @param title Title for the output spectral plots and csv file. Default is
#' `Autofluorescence spectra`.
#' @param verbose Logical, controls messaging. Default is `TRUE`.
#'
#' @return A matrix of autofluorescence spectra.
#'
#' @export

get.af.spectra <- function( unstained.sample,
                            asp, spectra,
                            threads = NULL,
                            som.dim = 10,
                            figures = TRUE,
                            plot.dir = NULL,
                            table.dir = NULL,
                            title = NULL,
                            verbose = TRUE ) {

  # set defaults
  if ( is.null( plot.dir ) )
    plot.dir <- asp$figure.af.dir

  if ( is.null( threads ) )
    threads <- asp$worker.process.n

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  spectral.channels <- colnames( spectra )

  # import unstained sample
  unstained.ff <- suppressWarnings(
    flowCore::read.FCS(
      unstained.sample,
      transformation = FALSE,
      truncate_max_range = FALSE,
      emptyValue = FALSE )
  )

  unstained.exprs <- flowCore::exprs( unstained.ff )[ , spectral.channels ]

  unmixed <- unmix.ols( unstained.exprs, spectra )

  cluster.data <- cbind( unstained.exprs, unmixed )

  # get cluster of AF from unstained
  if ( verbose ) message( "Creating a self-organizing map of the autofluorescence" )

  set.seed( 42 )
  map <- EmbedSOM::SOM(
    cluster.data,
    xdim = som.dim, ydim = som.dim,
    batch = TRUE, parallel = TRUE,
    threads = threads
    )

  af.spectra <- t(
    apply( map$codes[ , spectral.channels ], 1, function( x ) {
      max.x <- ifelse( max( abs( x ) ) > max( x ), min( x ), max( x ) )
      x / max.x } )
    )
  rownames( af.spectra ) <- paste0( "AF", 1:nrow( af.spectra ) )

  af.spectra <- as.matrix( na.omit( af.spectra ) )

  # save as CSV
  if ( is.null( title ) )
    af.file.name <- paste0( asp$af.file.name, ".csv" )
  else
    af.file.name <- paste0( title, ".csv" )

  if ( is.null( table.dir ) )
    table.dir <- asp$table.spectra.dir

  write.csv( af.spectra, file = file.path( table.dir, af.file.name ) )

  if ( figures ) {
    if ( verbose ) message( "Plotting spectra" )

    if ( !dir.exists( plot.dir ) )
      dir.create( plot.dir )

    if ( is.null( title ) )
      title <- asp$af.file.name

    spectral.trace(
      spectral.matrix = af.spectra,
      asp = asp,
      title = title,
      plot.dir = plot.dir,
      split.lasers = FALSE
      )

    spectral.heatmap( af.spectra, title, plot.dir )
  }

  return( af.spectra )
}
