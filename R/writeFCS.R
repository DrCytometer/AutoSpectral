# writeFCS.r

#' @title Write FCS File
#'
#' @description
#' A light weight FCS file binary writer based on `flowstate`. Lower memory usage
#' and faster than `flowCore`.
#'
#' @param mat The expression data, containing both unmixed data and any
#' retained parameters such as scatter and time.
#' @param keys The keywords to be written to the file.
#' @param file.name The name of the FCS file to be written.
#' @param output.dir A character string specifying the directory to save the
#' unmixed FCS file.
#'
#' @export
#'
#' @seealso \code{\link[flowCore:write.FCS]{write.FCS}}
#'
#' @references
#' Laniewski, Nathan. \emph{flowstate}.
#' \url{https://github.com/nlaniewski/flowstate}
#'
#' Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J, Jiang M,
#' Finak G (2025). \emph{flowCore: Basic structures for flow cytometry data}.
#' \doi{10.18129/B9.bioc.flowCore}
#' \url{https://bioconductor.org/packages/flowCore}

writeFCS <- function(mat, keys, file.name, output.dir) {
  delim <- "|"

  # Ensure mandatory keys
  keys[['$TOT']] <- as.character(nrow(mat))
  keys[['$PAR']] <- as.character(ncol(mat))
  keys[['$DATATYPE']] <- 'F'
  keys[['$BYTEORD']] <- '1,2,3,4'

  # Start with 0 placeholders
  keys[['$BEGINDATA']] <- "0"
  keys[['$ENDDATA']] <- "0"

  # Create initial text segment
  text.segment <- paste0(delim, paste0(names(keys), delim, unlist(keys), delim, collapse = ""))

  TEXT.start <- 58
  TEXT.end <- nchar(text.segment, "bytes") + TEXT.start - 1
  data.stream.bytes <- nrow(mat) * ncol(mat) * 4

  # copying flowstate here since it works
  kw.len.old <- 2
  repeat {
    DATA.start <- TEXT.end + 1
    DATA.end <- DATA.start + data.stream.bytes - 1
    kw.len.new <- nchar(DATA.start) + nchar(DATA.end)
    if (kw.len.new > kw.len.old) {
      TEXT.end <- TEXT.end + kw.len.new - kw.len.old
      kw.len.old <- kw.len.new
    } else break
  }

  # Update text segment with final offsets
  text.segment <- sub("\\|\\$BEGINDATA\\|0\\|", paste0("|$BEGINDATA|", DATA.start, "|"), text.segment)
  text.segment <- sub("\\|\\$ENDDATA\\|0\\|", paste0("|$ENDDATA|", DATA.end, "|"), text.segment)

  # Prepare Header (58 bytes)
  # Large file check: if offsets > 8 chars, use 0
  h_t <- if(nchar(TEXT.start) > 8 || nchar(TEXT.end) > 8) c(0,0) else c(TEXT.start, TEXT.end)
  h_d <- if(nchar(DATA.start) > 8 || nchar(DATA.end) > 8) c(0,0) else c(DATA.start, DATA.end)

  header <- paste0(
    sprintf("%-10s", "FCS3.1"),
    sprintf("%8d%8d", h_t[1], h_t[2]),
    sprintf("%8d%8d", h_d[1], h_d[2]),
    sprintf("%8d%8d", 0, 0) # Analysis offsets
  )

  # Binary Write
  con <- file(file.path(output.dir, file.name), open = "wb")
  on.exit(close(con))

  writeChar(header, con, eos = NULL)
  writeChar(text.segment, con, eos = NULL)

  # Note: t(mat) is crucial for row-major FCS format
  writeBin(as.vector(t(mat)), con, size = 4, endian = "little")

  # FCS standard footer (8 zeros)
  writeChar("00000000", con, eos = NULL)
}
