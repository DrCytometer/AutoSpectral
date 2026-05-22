# zzz.R  (AutoSpectral)
#
# Package load/unload hooks.
#
# .onLoad:  silent setup — injects Rcpp-accelerated readFCS and writeFCS into
#           this package's own namespace if AutoSpectralRcpp is available.
#           Works correctly under both devtools::load_all() and normal
#           install + library() workflows.
#
# .onAttach: user-facing startup message (packageStartupMessage belongs here,
#            not in .onLoad — see ?.onAttach).
#
# .onUnload: restores original pure-R implementations on detach/reload.
#
# --- Why getExportedValue() rather than exists(..., mode="function") ---
# exists() on asNamespace() is unreliable under devtools::load_all(): the
# namespace is populated but export registrations from NAMESPACE may not be
# fully reflected, causing false negatives.  getExportedValue() checks the
# exports table directly, which is the correct gate — we only want to inject
# functions that are genuinely part of AutoSpectralRcpp's public API.
#
# --- Why isNamespaceLoaded() before requireNamespace() ---
# requireNamespace() during load_all() of AutoSpectralRcpp can load the
# *installed* copy rather than the in-development one, producing a stale
# version or a namespace conflict.  Checking isNamespaceLoaded() first means
# we use whichever copy is already active — installed or load_all()'d.

# Store original pure-R implementations so .onUnload can restore them.
.original_readFCS  <- NULL
.original_writeFCS <- NULL

# Track injected functions for the .onAttach message.
.injected_fns <- character(0)

# Helper: safely retrieve an exported function from AutoSpectralRcpp.
# Returns the function if found and exported, NULL otherwise.
.get_rcpp_fn <- function(fn_name) {
  tryCatch(
    getExportedValue("AutoSpectralRcpp", fn_name),
    error = function(e) NULL
  )
}

# Helper: inject one function from AutoSpectralRcpp into this namespace.
# Returns TRUE on success.
.inject <- function(fn_name) {
  fn <- .get_rcpp_fn(fn_name)

  if (!is.function(fn)) {
    warning(
      "AutoSpectralRcpp is available but does not export ", fn_name, ". ",
      "Falling back to the pure-R implementation.",
      call. = FALSE
    )
    return(FALSE)
  }

  # Snapshot original for .onUnload restoration
  utils::assignInMyNamespace(
    paste0(".original_", fn_name),
    get(fn_name, envir = asNamespace("AutoSpectral"))
  )

  # Replace with Rcpp version
  utils::assignInMyNamespace(fn_name, fn)

  TRUE
}

# Helper: restore one function from its snapshot.
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

  # Check if AutoSpectralRcpp namespace is already loaded (e.g. via load_all())
  # before falling back to requireNamespace(), which may load the installed copy.
  rcpp_available <- isNamespaceLoaded("AutoSpectralRcpp") ||
                    requireNamespace("AutoSpectralRcpp", quietly = TRUE)

  if (!rcpp_available) return(invisible(NULL))

  injected <- character(0)
  if (.inject("readFCS"))  injected <- c(injected, "readFCS")
  if (.inject("writeFCS")) injected <- c(injected, "writeFCS")

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
