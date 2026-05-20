# zzz.R  (AutoSpectral)
#
# Package load/unload hooks.
#
# .onLoad: if AutoSpectralRcpp is installed, replace AutoSpectral::readFCS
# with AutoSpectralRcpp::readFCS in this package's own namespace.  This is
# done once at load time so there is zero per-call overhead — callers use
# AutoSpectral::readFCS as normal and get the fast version transparently.
#
# .onUnload: reverse the injection so the namespace is clean if the package
# is reloaded without AutoSpectralRcpp (e.g. during development).

# Store original R implementations so we can restore them on unload.
# Populated by .onLoad if injection occurs.
.original_readFCS  <- NULL
.original_writeFCS <- NULL

# Helper: inject one function from AutoSpectralRcpp into the AutoSpectral ns.
# Returns TRUE if injection succeeded, FALSE otherwise.
.inject <- function(fn_name, ns) {
  if (!exists(fn_name, envir = asNamespace("AutoSpectralRcpp"), mode = "function")) {
    warning(
      "AutoSpectralRcpp is installed but does not export ", fn_name, ". ",
      "Falling back to the pure-R implementation.",
      call. = FALSE
    )
    return(FALSE)
  }
  # Snapshot original
  assign(
    paste0(".original_", fn_name),
    get(fn_name, envir = ns),
    envir = ns
  )
  unlockBinding(fn_name, ns)
  assign(fn_name, get(fn_name, envir = asNamespace("AutoSpectralRcpp")), envir = ns)
  lockBinding(fn_name, ns)
  TRUE
}

# Helper: restore one function in the AutoSpectral ns from its snapshot.
.restore <- function(fn_name, ns) {
  snapshot_name <- paste0(".original_", fn_name)
  original <- tryCatch(
    get(snapshot_name, envir = ns, inherits = FALSE),
    error = function(e) NULL
  )
  if (!is.null(original) && is.function(original)) {
    unlockBinding(fn_name, ns)
    assign(fn_name, original, envir = ns)
    lockBinding(fn_name, ns)
  }
}

.onLoad <- function(libname, pkgname) {
  
  if (!requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    # AutoSpectralRcpp not installed - use pure-R implementations, nothing to do.
    return(invisible(NULL))
  }
  
  ns        <- asNamespace("AutoSpectral")
  injected  <- character(0)
  
  if (.inject("readFCS",  ns)) injected <- c(injected, "readFCS")
  if (.inject("writeFCS", ns)) injected <- c(injected, "writeFCS")
  
  if (length(injected) > 0) {
    packageStartupMessage(
      "AutoSpectralRcpp detected: using Rcpp-accelerated ",
      paste(injected, collapse = " and "), "."
    )
  }
  
  invisible(NULL)
}

.onUnload <- function(libpath) {
  ns <- asNamespace("AutoSpectral")
  .restore("readFCS",  ns)
  .restore("writeFCS", ns)
  invisible(NULL)
}
