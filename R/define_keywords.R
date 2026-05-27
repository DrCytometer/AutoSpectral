# define_keywords.r

#' @title Define Keywords
#'
#' @description
#' Updates FCS file keywords after unmixing to define the new parameters. Tracks
#' existing keywords from the input FCS file for metadata compatibility.
#'
#' Per [Data File Standard for Flow Cytometry, version FCS 3.1](https://doi.org/10.1002/cyto.a.20825):
#' \itemize{
#'   \item $LAST_MODIFIED/dd-mmm-yyyy hh:mm:ss.cc
#'   \item $LAST_MODIFIER/string/
#'   \item $ORIGINALITY/string
#' }
#'
#' are added as additional keywords to denote modification of the FCS file.
#'
#' @param fcs.keywords The input keywords obtained from the read FCS file.
#' @param final.matrix The expression data, containing both unmixed data and any
#' retained parameters such as scatter and time.
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
#' @param include.imaging A logical value indicating whether to include imaging
#' parameters in the written FCS file. Default is `FALSE` to provide smaller
#' output files.
#'
#' @return The updated keyword list for writing the FCS file.
#'
#' @export

define.keywords <- function(
    fcs.keywords,
    final.matrix,
    spectra,
    af.spectra,
    flow.control,
    asp,
    method,
    file.name,
    weights = NULL,
    spectral.channel = NULL,
    include.imaging = FALSE
) {

  # Identify non-parameter keywords
  non.param.keys <- fcs.keywords[!grepl("^\\$?P\\d+", names(fcs.keywords))]
  if (asp$cytometer == "Mosaic") {
    non.param.keys <- non.param.keys[!grepl("^\\$?CH\\d+", names(non.param.keys))]
  }
  non.param.keys <- non.param.keys[!names(non.param.keys) %in% c(
    "BDCHORUSDATARECORD"
  )]

  # Build parameter lookup from original file
  pN.keys <- grep("^\\$?P\\d+N$", names(fcs.keywords), value = TRUE)
  param.lookup <- lapply(pN.keys, function(k) {
    p.idx <- sub("\\$?P(\\d+)N", "\\1", k)
    matches <- grep(paste0("^\\$?P", p.idx, "(?:[A-Z]+)$"), names(fcs.keywords), value = TRUE)
    stats::setNames(fcs.keywords[matches], matches)
  })
  names(param.lookup) <- sapply(pN.keys, function(k) fcs.keywords[[k]])

  # Whitelist: only channels genuinely carried over from the original file may
  # match param.lookup. BD instruments (A5SE, Discover S8/A8) write pre-unmixed
  # FCS files whose $PnN values are fluorophore names (e.g. "FITC-A", "PE-A")
  # that are identical to AutoSpectral's unmixed output column names. Without
  # this filter, param.lookup would match those names and inject metadata from
  # the wrong original parameter index into the new file, scrambling channels.
  #
  # asp$time.and.scatter is an explicit whitelist defined for Discover cytometers.
  # For other cytometers it is NULL, so we fall back to matching
  # asp$non.spectral.channel against the original parameter names — those are the
  # scatter/time channels on those instruments and are safe to carry through.
  if (!is.null(asp$time.and.scatter)) {
    whitelist.channels <- asp$time.and.scatter
    if (isTRUE(include.imaging) && !is.null(asp$non.spectral.channel)) {
      # for Discover + include.imaging, also retain imaging channels that match
      # the non-spectral pattern, deduplicated against time.and.scatter
      ns.pattern <- paste0(asp$non.spectral.channel, collapse = "|")
      ns.matches <- grep(ns.pattern, names(param.lookup), value = TRUE)
      whitelist.channels <- union(whitelist.channels, ns.matches)
    }
  } else if (!is.null(asp$non.spectral.channel)) {
    # no explicit whitelist: derive it by matching non-spectral patterns against
    # the original parameter names (covers scatter and time on Aurora, A5SE, etc.)
    ns.pattern <- paste0(asp$non.spectral.channel, collapse = "|")
    whitelist.channels <- grep(ns.pattern, names(param.lookup), value = TRUE)
  } else {
    whitelist.channels <- character(0)
  }
  # restrict lookup to whitelisted channel names; everything else is treated as
  # a new unmixed fluorophore and receives freshly generated metadata.
  # Raw spectral channels are also safe to carry through when retained
  # (include.raw = TRUE): they exist in the original file under the same names
  # and their voltage/range metadata is correct and worth preserving.
  if (!is.null(spectral.channel)) {
    whitelist.channels <- union(whitelist.channels, spectral.channel)
  }
  param.lookup <- param.lookup[names(param.lookup) %in% whitelist.channels]

  # Build new parameter keywords
  param.keywords <- list()
  n.param <- ncol(final.matrix)
  p.names <- colnames(final.matrix)
  bit.depth <- if (!is.null(asp$bit.depth)) asp$bit.depth else "32"

  for (i in seq_len(n.param)) {
    p.name <- p.names[i]
    # FCS 3.1 3.2.23: commas are not allowed in $PnN values as they conflict
    # with keywords such as $SPILLOVER and $TR. Replace any with underscores.
    if (grepl(",", p.name, fixed = TRUE)) {
      warning(sprintf(
        "Parameter name '%s' contains a comma, which is invalid in $PnN (FCS 3.1 3.2.23). Replacing with '_'.",
        p.name
      ), call. = FALSE)
      p.name <- gsub(",", "_", p.name, fixed = TRUE)
    }
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

      # Whitelisted carry-through parameters (Scatter, Time, and optionally imaging)
    } else if (p.name %in% names(param.lookup)) {
      old.entry <- param.lookup[[p.name]]
      # whitelist keywords
      keep.fields <- c("N", "S", "B", "E", "R", "G", "V", "DISPLAY", "TYPE" )

      old.names <- names(old.entry)
      field.types <- sub("^\\$?P\\d+", "", old.names)

      keep.idx <- field.types %in% keep.fields

      filtered.entry <- old.entry[keep.idx]
      new.names <- paste0(p.prefix, field.types[keep.idx])

      for (k in seq_along(filtered.entry)) {
        param.keywords[[new.names[k]]] <- filtered.entry[[k]]
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

  # FCS 3.1 3.2.18: $PnE and $PnR are required for every parameter.
  # 3.2.22: when $DATATYPE/F/, all parameters must have $PnE/0,0/.
  # Carry-through parameters copied from source files may be missing either;
  # apply safe fallbacks here
  for (i in seq_len(n.param)) {
    p.prefix <- paste0("$P", i)
    e.key <- paste0(p.prefix, "E")
    r.key <- paste0(p.prefix, "R")
    if (is.null(param.keywords[[e.key]])) {
      param.keywords[[e.key]] <- "0,0"
    }
    if (is.null(param.keywords[[r.key]])) {
      param.keywords[[r.key]] <- as.character(asp$expr.data.max)
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

  # if other aspects of FCS3.2 keywords prove problematic, activate this section:
  #sanitize.value <- function(x) {
  #  x <- as.character(x)
  #  x <- gsub("\\|", "/", x)   # replace delimiter
  #  x <- gsub("[\r\n]", " ", x) # remove line breaks
  #  x
  #}
  #new.keywords <- lapply(new.keywords, sanitize.value)

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
    "AUTOSPECTRALRCPP" = rcpp.ver,
    '$LAST_MODIFIED' = toupper(format(Sys.time(), "%d-%b-%Y %H:%M:%OS2")),
    '$LAST_MODIFIER' = sprintf("AutoSpectral_%s", utils::packageVersion("AutoSpectral")),
    '$ORIGINALITY' = "DataModified"
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
