# zzz.R  (AutoSpectral)

# Mutable state lives in an environment, not as top-level namespace bindings.
.AS <- new.env(parent = emptyenv())
.AS$original_readFCS  <- NULL
.AS$original_writeFCS <- NULL
.AS$injected_fns      <- character(0)

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
  .AS[[paste0("original_", fn_name)]] <-
    get(fn_name, envir = asNamespace("AutoSpectral"), inherits = FALSE)

  # Replace with Rcpp version
  utils::assignInMyNamespace(fn_name, fn)

  TRUE
}

# Helper: restore one function from its snapshot.
.restore <- function(fn_name) {
  original <- .AS[[paste0("original_", fn_name)]]
  if (is.function(original)) {
    utils::assignInMyNamespace(fn_name, original)
    .AS[[paste0("original_", fn_name)]] <- NULL
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

  .AS$injected_fns <- injected

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
