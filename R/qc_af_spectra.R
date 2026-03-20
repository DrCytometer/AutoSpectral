# qc_af_spectra.r

#' @title QC Autofluorescence Spectra
#'
#' @description
#' This function performs a comparison between the user's fluorophore spectra,
#' from `spectra`, and identified autofluorescence spectra to check for any
#' contamination. Contamination may occur if the unstained sample used to generate
#' the AF spectra was not actually unstained, contained fluorescent protein
#' reporter constructs, or was contaminated during acquisition, for example by
#' crossover from nearby wells on a plate. Quality control is performed using
#' cosine similarity.
#'
#' @param af.spectra Matrix or dataframe containing autofluorescence spectral
#' signatures. This should be in format AFn x detectors. Row names will be used
#' as the names.
#' @param spectra Matrix or dataframe containing spectral data. This
#' should be in format fluorophores x detectors. Row names will be used as the
#' fluorophore names. Column names will be used as the detectors (channels).
#' @param output.dir Directory where the files will be saved.
#' Default is `./figure_autofluorescence`.
#' @param remove Logical. When `TRUE`, identified highly similar AF spectra will
#' be removed, returning only AF spectra that at least 0.995 or less from any
#' fluorophore spectrum.
#' @param pass Numeric, default `1`. Counter to separate multiple passes of AF
#' extraction, such as occur when `refine=TRUE` in `get.af.spectra()`.
#'
#' @seealso
#' * [af.qc.plot()]
#' * [get.af.spectra()]
#'
#' @return None. Plots are saved to a PDF in `plot.dir`.
#'
#' @export

qc.af.spectra <- function(
    af.spectra,
    spectra,
    output.dir = "./figure_autofluorescence",
    remove = TRUE,
    pass = 1
  ) {

  fluorophore.n <- nrow( spectra )

  sim.values <- numeric( nrow( af.spectra ) )
  sim.pairs  <- character( nrow( af.spectra ) )

  # loop for each af.spectrum vs spectra, get max cosine similarity value and fluor
  for ( af in seq_len( nrow( af.spectra ) ) ) {
    combined.spectra <- rbind( af.spectra[ af, ], spectra )
    # similarity result, excluding self, for this AF
    sim.matrix <- cosine.similarity( combined.spectra )[ 1, 2:( fluorophore.n + 1 ), drop = FALSE ]

    # which is the most similar?
    most.similar <- which.max( sim.matrix )
    sim.values[ af ] <- sim.matrix[ , most.similar ]
    sim.pairs[ af ] <- colnames( sim.matrix )[ most.similar ]
  }

  # identify any above 0.995
  contaminants <- which( sim.values > 0.995 ) # test this on some data sets
  contaminant.n <- length( contaminants )

  # if any, remove
  if ( contaminant.n > 0 ) {

    if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

    message.body <- paste(
      contaminant.n, "autofluorescence spectra were identified as being extremely",
      "similar to fluorophore signatures in `spectra`. This usually indicates",
      "experimental contamination of the unstained sample with stained cells,",
      "often from single-color controls. This can occur when samples are run in",
      "plates with orbital shaking, causing physical transfer of well contents",
      "from one to another. For future experiments, position the unstained sample",
      "apart from any stained samples, and acquire multiple replicates to ensure",
      "clean data."
    )

    message.footer <- paste(
      "These 'contaminating' autofluorescence spectra have been removed and can",
      "be inspected in the QC plots in folder:", output.dir,
      "\n\nIf you believe these to be real autofluorescence profiles, re-run",
      "`get.af.spectra()`, setting `remove.contaminants = FALSE`."
    )

    # wrap the text to the width of the console and issue the warning
    warning(
      paste(
        strwrap( message.body, width = 80 ),
        collapse = "\n"
      ),
      "\n\n",
      paste(
        strwrap( message.footer, width = 80 ),
        collapse = "\n"
      ),
      call. = FALSE
    )

    # table, print and save with similar pairs
    contaminant.table <- data.frame(
      AF = rownames( af.spectra )[ contaminants ],
      Fluorophore = sim.pairs[ contaminants ],
      Similarity = sim.values[ contaminants ],
      stringsAsFactors = FALSE
    )
    print( contaminant.table )
    filename <- paste0( "Contaminant_table_AF_spectra_QC_pass_", pass, ".csv" )
    utils::write.csv( contaminant.table, file.path( output.dir, filename ) )

    # print pdf with QC table and plots
    tryCatch(
      expr = {
        af.qc.plot( af.spectra, spectra, contaminant.table,
                    plot.dir = output.dir )
      },
      error = function( e ) {
        message( "Error in plotting AF spectral QC: ", e$message )
        return( NULL )
      }
    )

    # remove problematic AF spectra
    af.spectra <- af.spectra[ which( sim.values < 0.995 ), , drop = FALSE ]
  }

  # return only dissimilar af.spectra
  return( af.spectra )
}
