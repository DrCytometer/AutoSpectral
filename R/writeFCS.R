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

  if (requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    return(AutoSpectralRcpp::writeFCS(
      mat        = mat,
      keys       = keys,
      file.name  = file.name,
      output.dir = output.dir
    ))
  }

  delim <- "|"

  # Ensure mandatory keys
  keys[['$TOT']] <- as.character(nrow(mat))
  keys[['$PAR']] <- as.character(ncol(mat))
  keys[['$DATATYPE']] <- 'F'
  keys[['$BYTEORD']] <- '1,2,3,4'
  keys[['$NEXTDATA']] <- '0'

  # Segment offsets: placeholders updated after layout is calculated.
  keys[['$BEGINDATA']]     <- "0"
  keys[['$ENDDATA']]       <- "0"
  keys[['$BEGINSTEXT']]    <- "0"
  keys[['$ENDSTEXT']]      <- "0"
  keys[['$BEGINANALYSIS']] <- "0"
  keys[['$ENDANALYSIS']]   <- "0"

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

  # FCS requires row-major float32: writeBin(as.vector(t(mat)), ...) would
  # allocate a full transposed double matrix + a float32 conversion buffer,
  # totalling ~20 bytes/value of transient memory on top of mat itself.
  # For large high-channel-count files (e.g. 500k events x 183 channels)
  # this can exhaust RAM.
  #
  # Instead, we transpose and write CHUNK_ROWS events at a time.
  # t(mat[rows, ]) produces a chunk of size n_col * chunk_rows (doubles),
  # and writeBin converts to float32 internally, so peak extra memory is
  # ~CHUNK_ROWS * n_col * (8 + 4) bytes = ~12 bytes/value per chunk.
  # At CHUNK_ROWS = 131072 and 183 channels that is ~288 MB.
  CHUNK_ROWS     <- 131072L
  n_row          <- nrow(mat)
  rows_remaining <- n_row
  row_offset     <- 0L

  while (rows_remaining > 0L) {
    rows_this_chunk <- min(CHUNK_ROWS, rows_remaining)
    row_idx         <- seq_len(rows_this_chunk) + row_offset
    writeBin(as.vector(t(mat[row_idx, , drop = FALSE])),
             con, size = 4, endian = "little")
    row_offset     <- row_offset + rows_this_chunk
    rows_remaining <- rows_remaining - rows_this_chunk
  }

  # FCS standard footer (8 zeros)
  writeChar("00000000", con, eos = NULL)
}
