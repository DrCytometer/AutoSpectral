# zzz.R  (AutoSpectral)

# Mutable state container — never use top-level NULL bindings mutated by
# assignInMyNamespace(); that corrupts .rdb files under R 4.6+/libdeflate.
.AS <- new.env(parent = emptyenv())

# Minimum AutoSpectralRcpp version that supports the `columns` argument to
# readFCS(). Referenced by readFCS()'s compatibility fallback and by
# .onAttach() below. Update in this one place if a future AutoSpectralRcpp
# feature requires a newer minimum.
.AS.RCPP.MIN.COLUMNS <- "1.3.0"

.onAttach <- function(libname, pkgname) {
  if (requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    packageStartupMessage(
      "AutoSpectralRcpp detected: using Rcpp-accelerated readFCS and writeFCS."
    )

    installed.version <- tryCatch(
      utils::packageVersion("AutoSpectralRcpp"),
      error = function(e) NULL
    )

    if (!is.null(installed.version) &&
        installed.version < package_version(.AS.RCPP.MIN.COLUMNS)) {
      packageStartupMessage(
        "NOTE: installed AutoSpectralRcpp (", installed.version, ") predates ",
        "AutoSpectral's `columns` argument to readFCS() (added in ",
        "AutoSpectralRcpp ", .AS.RCPP.MIN.COLUMNS, "). AutoSpectral will keep ",
        "working, but memory-optimized reads will fall back to a slower path. ",
        "Update with:\n  pak::pak(\"DrCytometer/AutoSpectralRcpp\")"
      )
    }
  }
}
