# readFCSheader.r

#' @title Read FCS Header
#'
#' @description
#' A light weight FCS file header reader based on `flowstate`, trying to mimic
#' `flowCore`.
#'
#' @param fcs.path A character string specifying the file path (directory and
#' file name) for the .fcs file to be read.
#' @param keyword Optional argument specifying which keyword(s) to return if the
#' whole set is not desired.
#'
#' @return A named list of keywords (metadata)
#'
#' @export
#'
#' @seealso \code{\link[flowCore:read.FCSheader]{read.FCSheader}}
#'
#' @references
#' Laniewski, Nathan. \emph{flowstate}.
#' \url{https://github.com/nlaniewski/flowstate}
#'
#' Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J, Jiang M,
#' Finak G (2025). \emph{flowCore: Basic structures for flow cytometry data}.
#' \doi{10.18129/B9.bioc.flowCore}
#' \url{https://bioconductor.org/packages/flowCore}

readFCSheader <- function(fcs.path, keyword = NULL) {
  con <- file(fcs.path, open = "rb")
  on.exit(close(con))

  # Read the 58-byte header
  header.raw <- readChar(con, 58)
  if (length(header.raw) == 0) return(NULL)

  # Extract FCS version
  fcs.version <- trimws(substr(header.raw, 4, 10))
  txt.st <- as.numeric(trimws(substr(header.raw, 11, 18)))
  txt.en <- as.numeric(trimws(substr(header.raw, 19, 26)))

  # Jump to and read the text segment
  seek(con, txt.st)
  txt.raw <- readBin(con, "raw", n = txt.en - txt.st + 1)
  txt <- rawToChar(txt.raw[txt.raw != as.raw(0)])

  # Parse keywords
  delim <- substr(txt, 1, 1)
  kv <- strsplit(txt, delim, fixed = TRUE)[[1]]

  # Remove the first element if empty (common with leading delimiters)
  if (length(kv) > 0 && kv[1] == "") kv <- kv[-1]
  if (length(kv) %% 2 != 0) kv <- c(kv, "")

  keys.names <- kv[seq(1, length(kv), by = 2)]
  keys.values <- kv[seq(2, length(kv), by = 2)]

  # Prepend the FCSversion to match flowCore's 3340 count
  all.names <- c("FCSversion", keys.names)
  all.values <- c(fcs.version, keys.values)

  # Format as a named character vector inside a named list
  result.vector <- stats::setNames(all.values, all.names)

  # keyword filtering
  if (!is.null(keyword)) {
    result.vector <- result.vector[all.names %in% keyword]
    if (length(result.vector) == 0) result.vector <- NULL
  }

  out <- list(result.vector)
  names(out) <- fcs.path

  return(out)
}

