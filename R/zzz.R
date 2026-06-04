# zzz.R  (AutoSpectral)

# Mutable state container — never use top-level NULL bindings mutated by
# assignInMyNamespace(); that corrupts .rdb files under R 4.6+/libdeflate.
.AS <- new.env(parent = emptyenv())

.onAttach <- function(libname, pkgname) {
  if (requireNamespace("AutoSpectralRcpp", quietly = TRUE)) {
    packageStartupMessage(
      "AutoSpectralRcpp detected: using Rcpp-accelerated readFCS and writeFCS."
    )
  }
}
