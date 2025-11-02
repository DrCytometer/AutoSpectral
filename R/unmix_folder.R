# unmix_folder.r

#' @title Unmix All FCS Files in a Directory
#'
#' @description
#' This function unmixes all FCS files in a specified directory using the
#' provided spectra and method, and saves the unmixed FCS files to an output
#' directory of the user's choice.
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @param fcs.dir Directory containing FCS files to be unmixed.
#' @param spectra Matrix containing spectra information.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method A character string specifying the unmixing method to use.
#' The default is `Automatic`, which uses `AutoSpectral` for AF extraction if
#' af.spectra are provided and automatically selects `OLS` or `WLS` depending
#' on which is normal for the given cytometer in `asp$cytometer`. This means
#' that files from the ID7000, A8 and S8 will be unmixed using `WLS` while
#' others will be unmixed with `OLS`. Any option can be set manually.
#' Manual options are `OLS`, `WLS`, `AutoSpectral`, `Poisson` and `FastPoisson`.
#' Default is `OLS`. `FastPoisson` requires installation of `AutoSpectralRcpp`.
#' @param weighted Logical, whether to use ordinary or weighted least squares
#' unmixing as the base algorithm in AutoSpectral unmixing.
#' Default is `FALSE` and will use OLS.
#' @param weights Optional numeric vector of weights: one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means. Only used for `WLS`
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`. Required for `AutoSpectral` unmixing. Default is
#' `NULL` and will thus provoke failure if no spectra are provided and
#' `AutoSpectral` is selected.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants`. Default is `NULL`. Used for
#' AutoSpectral unmixing. Required for per-cell fluorophore optimization.
#' @param output.dir Directory to save the unmixed FCS files
#' (default is asp$unmixed.fcs.dir).
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`
#' @param include.raw Logical indicating whether to include raw data in the
#' written FCS file. Default is `FALSE`
#' @param include.imaging Logical indicating whether to include imaging data in
#' the written FCS file: relevant for S8 and A8. Default is `FALSE`
#' @param calculate.error Logical, whether to calculate the RMSE unmixing model
#' accuracy and include it as a keyword in the FCS file.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#' signature for each cell is determined by which unmixing brings the fluorophore
#' signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing minimizes the
#' per-cell residual (`use.dist0` = `FALSE`). Default is `TRUE`. Used for
#' AutoSpectral unmixing.
#' @param divergence.threshold Numeric. Used for `FastPoisson` only. Threshold
#' to trigger reversion towards WLS unmixing when Poisson result diverges.
#' Default is `1e4`
#' @param divergence.handling String. How to handle divergent cells from Poisson
#' IRLS. Options are `NonNeg`, in which case non-negativity will be enforced,
#' `WLS`, where values will revert to the WLS initial unmix or `Balance`,
#' where `WLS` and `NonNeg` results will be averaged. Default is `Balance`
#' @param balance.weight Numeric. Weighting to average non-convergent cells.
#' Used for `Balance` option under `divergence.handling`. Default is `0.5`
#' @param speed Selector for the precision-speed trade-off in AutoSpectral per-cell
#' fluorophore optimization. Options are the default `fast`, which selects the
#' best spectral fit per cell by updating the predicted values for each
#' fluorophore independently without repeating the unnmixing, `medium` which uses
#' a Woodbury-Sherman-Morrison rank-one updating of the unnmixing matrix for
#' better results and a moderate slow-down, or `slow`, which explicitly
#' recomputes the unmixing matrix for each variant for maximum precision. The
#' `fast` method is only one available in the `AutoSpectral` package and will be
#' slow in the pure R implementation. Installation of `AutoSpectralRcpp` is
#' strongly encouraged.
#'
#' @return None. Saves the unmixed FCS files to the specified output directory.
#'
#' @export

unmix.folder <- function( fcs.dir, spectra, asp, flow.control,
                          method = "Automatic",
                          weighted = FALSE,
                          weights = NULL,
                          af.spectra = NULL,
                          spectra.variants = NULL,
                          output.dir = NULL,
                          file.suffix = NULL,
                          include.raw = FALSE,
                          include.imaging = FALSE,
                          calculate.error = FALSE,
                          use.dist0 = TRUE,
                          divergence.threshold = 1e4,
                          divergence.handling = "Balance",
                          balance.weight = 0.5,
                          speed = "fast" ){

  if ( is.null( output.dir ) )
    output.dir <- asp$unmixed.fcs.dir

  files.to.unmix <- list.files( fcs.dir, pattern = ".fcs", full.names = TRUE )

  # set up parallel processing
  if ( asp$parallel && ( method == "OLS" || method == "WLS" ) ) {
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  lapply.function( files.to.unmix, FUN = unmix.fcs, spectra, asp, flow.control,
                   method, weighted, weights, af.spectra, spectra.variants,
                   output.dir, file.suffix, include.raw,include.imaging,
                   calculate.error, use.dist0,
                   divergence.threshold, divergence.handling, balance.weight,
                   speed )
}
