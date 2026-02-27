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
#' @param start.row Optional numeric specifying the row to begin reading on. Can
#' be useful for reading in just the metadata or for chunking files. Default is
#' `NULL`, which will read the whole file.
#' @param end.row Optional numeric specifying the row to end reading on. Can
#' be useful for reading in just the metadata or for chunking files. Default is
#' `NULL`, which will read the whole file.
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


readFCS <- function(
    fcs.path,
    return.keywords = FALSE,
    start.row = NULL,
    end.row = NULL
) {

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
  if (length(kv) > 0 && kv[1] == "") kv <- kv[-1]
  if (length(kv) %% 2 != 0) kv <- c(kv, "")

  # Map keys to values
  keys.names <- kv[seq(1, length(kv), by = 2)]
  keys.values <- kv[seq(2, length(kv), by = 2)]
  keywords <- stats::setNames(as.list(keys.values), keys.names)

  # Data Offsets & Metadata
  data.st <- as.numeric(keywords[["$BEGINDATA"]])
  total.events <- as.numeric(keywords[["$TOT"]])
  n.par <- as.numeric(keywords[["$PAR"]])

  # Determine Endianness
  byte.ord <- keywords[["$BYTEORD"]]
  endian <- if (!is.null(byte.ord) && byte.ord == "1,2,3,4") "little" else "big"

  # Determine Reading Range
  read.start <- if(is.null(start.row)) 1 else start.row
  read.end <- if(is.null(end.row)) total.events else end.row
  num.rows.to.read <- read.end - read.start + 1

  # Offset = Start of Data + (Skip Rows * Parameters * 4 bytes per float)
  byte.offset <- data.st + ((read.start - 1) * n.par * 4)
  seek(con, byte.offset)

  n.vals <- num.rows.to.read * n.par

  # Read Matrix
  data.mat <- matrix(
    readBin(con, "numeric", n = n.vals, size = 4, endian = endian),
    ncol = n.par, byrow = TRUE
  )

  # Column naming
  col.names <- unname(sapply(1:n.par, function(i) {
    val <- keywords[[paste0("$P", i, "N")]]
    if (is.null(val)) paste0("Channel_", i) else val
  }))

  colnames(data.mat) <- col.names

  if (return.keywords) {
    return(list(data = data.mat, keywords = keywords))
  } else {
    return(data.mat)
  }

}
