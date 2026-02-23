# readFCS.r

#' @title Read FCS File
#'
#' @description
#' A light weight FCS file binary reader based on `flowstate`. Lower memory usage
#' and faster than `flowCore`.
#'
#' @param fcs.path A character string specifying the file path (directory and
#' file name) for the .fcs file to be read.
#' @param return.keywords Logical, default `FALSE`. Controls whether keywords are
#' returned as well as the expression data matrix (useful for writing a new file).
#'
#' @return If `return.keywords = TRUE`, a list containing two elements:
#' 1) the expression data in a matrix, and 2) the keywords.
#' If `return.keywords = FALSE`, only the expression data is returned as a matrix.
#'
#' @export
#'
#' @seealso \code{\link[flowCore:read.FCS]{read.FCS}}
#'
#' @references
#' Laniewski, Nathan. \emph{flowstate}.
#' \url{https://github.com/nlaniewski/flowstate}
#'
#' Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J, Jiang M,
#' Finak G (2025). \emph{flowCore: Basic structures for flow cytometry data}.
#' \doi{10.18129/B9.bioc.flowCore}
#' \url{https://bioconductor.org/packages/flowCore}


readFCS <- function(fcs.path, return.keywords = FALSE) {
  con <- file(fcs.path, open = "rb")
  on.exit(close(con))

  # Header & Text Offsets
  header <- readChar(con, 58)
  txt.st <- as.numeric(trimws(substr(header, 11, 18)))
  txt.en <- as.numeric(trimws(substr(header, 19, 26)))

  # Extract Keywords
  seek(con, txt.st)
  txt.raw <- readBin(con, "raw", n = txt.en - txt.st + 1)
  # Remove null bytes
  txt <- rawToChar(txt.raw[txt.raw != as.raw(0)])

  # Identify delimiter (first character)
  delim <- substr(txt, 1, 1)
  kv <- strsplit(txt, delim, fixed = TRUE)[[1]]

  # Remove the first element if it's empty (since txt starts with delim)
  if (length(kv) > 0 && kv[1] == "") {
    kv <- kv[-1]
  }

  # handle empty keywords (Opteon)
  if (length(kv) %% 2 != 0) {
    kv <- c(kv, "") # Append an empty string to balance the pairs
  }

  # Map keys to values
  keys_names <- kv[seq(1, length(kv), by = 2)]
  keys_values <- kv[seq(2, length(kv), by = 2)]

  keywords <- stats::setNames(as.list(keys_values), keys_names)

  # Data Offsets & Metadata
  data.st <- as.numeric(keywords[["$BEGINDATA"]])
  data.en <- as.numeric(keywords[["$ENDDATA"]])
  n.par <- as.numeric(keywords[["$PAR"]])

  # Determine Endianness
  # 1,2,3,4 is little-endian (Intel); 4,3,2,1 is big-endian (Sun/PowerPC)
  byte_ord <- keywords[["$BYTEORD"]]
  endian <- if (!is.null(byte_ord) && byte_ord == "1,2,3,4") "little" else "big"

  # Read Matrix
  seek(con, data.st)
  n.vals <- (data.en - data.st + 1) / 4

  data.mat <- matrix(
    readBin(con, "numeric", n = n.vals, size = 4, endian = endian),
    ncol = n.par, byrow = TRUE
  )

  # Column naming
  col_names <- unname(sapply(1:n.par, function(i) {
    val <- keywords[[paste0("$P", i, "N")]]
    if (is.null(val)) paste0("Channel_", i) else val
  }))

  colnames(data.mat) <- col_names

  if ( return.keywords ) {
    return(list(data = data.mat, keywords = keywords))
  } else {
    return(data.mat)
  }

}
