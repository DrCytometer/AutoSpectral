# zzz.R  (AutoSpectral)
#
# Package load/unload hooks.
#
# .onLoad:  silent setup only — injects Rcpp-accelerated versions of readFCS
#           and writeFCS into this package's own namespace if AutoSpectralRcpp
#           is installed.
#
# .onAttach: emits the user-facing startup message
#
# .onUnload: restores the original pure-R implementations so the namespace is
#            clean if the package is reloaded without AutoSpectralRcpp.

# Store original pure-R implementations so .onUnload can restore them.
# Populated by .onLoad if injection occurs; NULL means no injection happened.
.original_readFCS  <- NULL
.original_writeFCS <- NULL

# Track which functions were successfully injected so .onAttach can report.
.injected_fns <- character(0)

# Helper: inject one function from AutoSpectralRcpp into this package's ns.
# Uses assignInMyNamespace(), no unlockBinding() needed.
# Returns TRUE on success, FALSE if the symbol is absent from AutoSpectralRcpp.
.inject <- function(fn_name) {
  if (!exists(fn_name,
              envir = asNamespace("AutoSpectralRcpp"),
              mode  = "function")) {
    warning(
      "AutoSpectralRcpp is installed but does not export ", fn_name, ". ",
      "Falling back to the pure-R implementation.",
      call. = FALSE
    )
    return(FALSE)
  }

  # Snapshot the original so .onUnload can restore it
  utils::assignInMyNamespace(
    paste0(".original_", fn_name),
    get(fn_name, envir = asNamespace("AutoSpectral"))
  )

  # Replace with the Rcpp version
  utils::assignInMyNamespace(
    fn_name,
    get(fn_name, envir = asNamespace("AutoSpectralRcpp"))
  )

  TRUE
}

# Helper: restore one function from its snapshot, if one exists.
.restore <- function(fn_name) {
  snapshot_name <- paste0(".original_", fn_name)
  original <- tryCatch(
    get(snapshot_name, envir = asNamespace("AutoSpectral"), inherits = FALSE),
    error = function(e) NULL
  )
  if (is.function(original)) {
    utils::assignInMyNamespace(fn_name, original)
    utils::assignInMyNamespace(snapshot_name, NULL)
  }
}

.onLoad <- function(libname, pkgname) {

  if (!requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    # AutoSpectralRcpp not installed: nothing to do.
    return(invisible(NULL))
  }

  injected <- character(0)
  if (.inject("readFCS"))  injected <- c(injected, "readFCS")
  if (.inject("writeFCS")) injected <- c(injected, "writeFCS")

  # Store for .onAttach
  utils::assignInMyNamespace(".injected_fns", injected)

  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  fns <- .injected_fns
  if (length(fns) > 0) {
    packageStartupMessage(
      "AutoSpectralRcpp detected: using Rcpp-accelerated ",
      paste(fns, collapse = " and "), "."
    )
  }
}

.onUnload <- function(libpath) {
  .restore("readFCS")
  .restore("writeFCS")
  invisible(NULL)
}
