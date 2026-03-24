# define_keywords.r

#' @title Define Keywords
#'
#' @description
#' Updates FCS file keywords after unmixing to define the new parameters. Tracks
#' existing keywords from the input FCS file for metadata compatibility.
#'
#' @param fcs.keywords The input keywords obtained from the read FCS file.
#' @param final.matrix The expression data, containing both unmixed data and any
#' retained parameters such as scatter and time.
#' @param original.param The original parameter (column) names of the input FCS
#' expression data.
#' @param spectra A matrix containing the spectral data.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns.
#' @param flow.control A list containing flow cytometry control parameters.
#' @param asp The AutoSpectral parameter list.
#' @param method A character string specifying the unmixing method used.
#' @param file.name The name of the FCS file to be written.
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector).
#' @param spectral.channel Optional character vector of the channels used for
#' unmixing. Should match `weights` in length (one weight per channel).
#'
#' @return The updated keyword list for writing the FCS file.
#'
#' @export

define.keywords <- function(
    fcs.keywords,
    final.matrix,
    original.param,
    spectra,
    af.spectra,
    flow.control,
    asp,
    method,
    file.name,
    weights = NULL,
    spectral.channel = NULL
  ) {

  # Identify non-parameter keywords
  non.param.keys <- fcs.keywords[!grepl("^\\$?P\\d+", names(fcs.keywords))]
  if (asp$cytometer == "Mosaic") {
    non.param.keys <- non.param.keys[!grepl("^\\$?CH\\d+", names(non.param.keys))]
  }

  # Build parameter lookup from original file
  pN.keys <- grep("^\\$?P\\d+N$", names(fcs.keywords), value = TRUE)
  param.lookup <- lapply(pN.keys, function(k) {
    p.idx <- sub("\\$?P(\\d+)N", "\\1", k)
    matches <- grep(paste0("^\\$?P", p.idx, "(?:[A-Z]+)$"), names(fcs.keywords), value = TRUE)
    stats::setNames(fcs.keywords[matches], matches)
  })
  names(param.lookup) <- sapply(pN.keys, function(k) fcs.keywords[[k]])

  # Build new parameter keywords
  param.keywords <- list()
  n.param <- ncol(final.matrix)
  p.names <- colnames(final.matrix)
  bit.depth <- if (!is.null(asp$bit.depth)) asp$bit.depth else "32"

  for (i in seq_len(n.param)) {
    p.name <- p.names[i]
    p.prefix <- paste0("$P", i)

    # AF Index
    if (p.name == "AF Index") {
      param.keywords[[paste0(p.prefix, "N")]] <- "AF Index"
      param.keywords[[paste0(p.prefix, "S")]] <- "Autofluorescence Index"
      param.keywords[[paste0(p.prefix, "B")]] <- as.character(bit.depth)
      param.keywords[[paste0(p.prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p.prefix, "R")]] <- as.character(nrow(af.spectra))
      param.keywords[[paste0(p.prefix, "G")]] <- "1"
      param.keywords[[paste0(p.prefix, "DISPLAY")]] <- "LIN"
      param.keywords[[paste0(p.prefix, "TYPE")]] <- "Fluorescence"

      # Autofluorescence
    } else if (p.name == "AF-A") {
      param.keywords[[paste0(p.prefix, "N")]] <- p.name
      param.keywords[[paste0(p.prefix, "S")]] <- "Autofluorescence"
      param.keywords[[paste0(p.prefix, "B")]] <- as.character(bit.depth)
      param.keywords[[paste0(p.prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p.prefix, "R")]] <- as.character(asp$expr.data.max)
      param.keywords[[paste0(p.prefix, "G")]] <- "1"
      param.keywords[[paste0(p.prefix, "DISPLAY")]] <- "LOG"
      param.keywords[[paste0(p.prefix, "TYPE")]] <- "Fluorescence"

      # Existing parameters (Scatter, Time, etc.)
    } else if (p.name %in% names(param.lookup)) {
      old.entry <- param.lookup[[p.name]]
      # Rename all keywords for this parameter to the new index 'i'
      new.names <- gsub("^\\$?P\\d+", p.prefix, names(old.entry))
      for (k in seq_along(old.entry)) {
        param.keywords[[new.names[k]]] <- old.entry[[k]]
      }

      # New Unmixed Fluorophores
    } else {
      param.keywords[[paste0(p.prefix, "N")]] <- p.name
      param.keywords[[paste0(p.prefix, "B")]] <- as.character(bit.depth)
      param.keywords[[paste0(p.prefix, "E")]] <- "0,0"
      param.keywords[[paste0(p.prefix, "R")]] <- as.character(asp$expr.data.max)
      param.keywords[[paste0(p.prefix, "G")]] <- "1"
      param.keywords[[paste0(p.prefix, "DISPLAY")]] <- "LOG"
      param.keywords[[paste0(p.prefix, "TYPE")]] <- "Fluorescence"

      # Map Marker/Stain from flow.control
      clean.name <- sub("-A$", "", p.name)
      f.idx <- match(clean.name, flow.control$fluorophore)
      marker <- if (!is.na(f.idx)) as.character(flow.control$antigen[f.idx]) else p.name
      param.keywords[[paste0(p.prefix, "S")]] <- marker
    }
  }

  # Format Spectra for Keywords
  format.matrix.string <- function(m) {
    vals <- as.vector(t(m))
    formatted <- formatC(vals, digits = 8, format = "fg", flag = "#")
    paste(c(nrow(m), ncol(m), rownames(m), colnames(m), formatted), collapse = ",")
  }

  # Combine everything
  new.keywords <- utils::modifyList(non.param.keys, param.keywords)

  # Add package versions
  asp.ver <- as.character(utils::packageVersion("AutoSpectral"))
  rcpp.ver <- if (requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    as.character(utils::packageVersion("AutoSpectralRcpp"))
  } else {
    "0"
  }

  new.keywords <- utils::modifyList(new.keywords, list(
    "$FIL" = file.name,
    "$PAR" = as.character(n.param),
    "$TOT" = as.character(nrow(final.matrix)),
    "UNMIXINGMETHOD" = method,
    "$BYTEORD" = "1,2,3,4",
    "$DATATYPE" = "F",
    "SPECTRA" = format.matrix.string(spectra),
    "FLUOROCHROMES" = paste(rownames(spectra), collapse = ","),
    "AUTOSPECTRAL" = asp.ver,
    "AUTOSPECTRALRCPP" = rcpp.ver
  ))

  # add AF spectra if used
  if (!is.null(af.spectra)) {
    new.keywords[["AUTOFLUORESCENCE"]] <- format.matrix.string(af.spectra)
  }

  # add weights if used
  if (!is.null(weights) && !is.null(spectral.channel)) {
    weights.str <- paste(
      c(length(spectral.channel), spectral.channel,
        formatC(weights, digits = 8, format = "fg")),
      collapse = ","
    )
    new.keywords[["WEIGHTS"]] <- weights.str
  }

  return(new.keywords)
}
