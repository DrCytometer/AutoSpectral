# zzz.R  (AutoSpectral)

# Mutable state container — never use top-level NULL bindings mutated by
# assignInMyNamespace(); that corrupts .rdb files under R 4.6+/libdeflate.
.AS <- new.env(parent = emptyenv())

# Minimum AutoSpectralRcpp version that supports the `columns` argument to
# readFCS(). Referenced by readFCS()'s compatibility fallback and by
# .onAttach() below. Update in this one place if a future AutoSpectralRcpp
# feature requires a newer minimum.
.AS.RCPP.MIN.COLUMNS <- "1.3.1"

#' Check GitHub for a newer tagged release of a package
#'
#' Queries the GitHub tags API for `repo`, compares the newest semantic-
#' version tag against `installed.version`, and returns that version if it is
#' newer. Silently returns NULL on any failure (no network, rate limit,
#' unparseable response, non-interactive session) since this is a convenience
#' check and must never interrupt or slow down a script. Results are cached
#' on disk for `cache.days` so a session does not re-query GitHub every time
#' the package is attached.
#'
#' @param pkg Package name, used as the cache subdirectory.
#' @param repo GitHub "owner/repo" string.
#' @param installed.version A `package_version` to compare against.
#' @param cache.days Minimum days between live checks. Default 1.
#' @return A `package_version` if a newer tag exists, otherwise NULL
#'   (invisibly).
.AS.check.github.update <- function(pkg, repo, installed.version, cache.days = 1) {
  if (!interactive()) {
    return(invisible(NULL))
  }

  cache.dir <- tryCatch(
    tools::R_user_dir(pkg, which = "cache"),
    error = function(e) NULL
  )
  if (is.null(cache.dir)) {
    return(invisible(NULL))
  }

  cache.file <- file.path(cache.dir, "update_check.rds")

  if (file.exists(cache.file)) {
    cached <- tryCatch(readRDS(cache.file), error = function(e) NULL)
    if (!is.null(cached) &&
        (Sys.time() - cached$time) < as.difftime(cache.days, units = "days")) {
      if (!is.null(cached$latest) && cached$latest > installed.version) {
        return(cached$latest)
      }
      return(invisible(NULL))
    }
  }

  old.timeout <- getOption("timeout")
  on.exit(options(timeout = old.timeout), add = TRUE)
  options(timeout = 3)

  raw.json <- tryCatch({
    con <- url(paste0("https://api.github.com/repos/", repo, "/tags"), open = "rb")
    on.exit(close(con), add = TRUE)
    paste(readLines(con, warn = FALSE, n = 40), collapse = "\n")
  }, error = function(e) NULL)

  if (is.null(raw.json)) {
    return(invisible(NULL))
  }

  tag.match <- regmatches(
    raw.json,
    regexpr('"name"\\s*:\\s*"v?[0-9]+\\.[0-9]+\\.[0-9]+"', raw.json)
  )
  if (length(tag.match) == 0 || nchar(tag.match) == 0) {
    return(invisible(NULL))
  }

  latest.version <- tryCatch(
    package_version(sub('.*"v?([0-9]+\\.[0-9]+\\.[0-9]+)".*', "\\1", tag.match)),
    error = function(e) NULL
  )
  if (is.null(latest.version)) {
    return(invisible(NULL))
  }

  dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(time = Sys.time(), latest = latest.version), cache.file)

  if (latest.version > installed.version) {
    latest.version
  } else {
    invisible(NULL)
  }
}

.onAttach <- function(libname, pkgname) {
  as.installed.version <- tryCatch(
    utils::packageVersion("AutoSpectral"),
    error = function(e) NULL
  )

  if (!is.null(as.installed.version)) {
    as.latest <- .AS.check.github.update(
      pkg = "AutoSpectral",
      repo = "DrCytometer/AutoSpectral",
      installed.version = as.installed.version
    )
    if (!is.null(as.latest)) {
      packageStartupMessage(
        "A newer AutoSpectral is available on GitHub (installed: ",
        as.installed.version, ", latest: ", as.latest, "). Update with:\n",
        "  pak::pak(\"DrCytometer/AutoSpectral\")"
      )
    }
  }

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

    if (!is.null(installed.version)) {
      rcpp.latest <- .AS.check.github.update(
        pkg = "AutoSpectralRcpp",
        repo = "DrCytometer/AutoSpectralRcpp",
        installed.version = installed.version
      )
      if (!is.null(rcpp.latest)) {
        packageStartupMessage(
          "A newer AutoSpectralRcpp is available on GitHub (installed: ",
          installed.version, ", latest: ", rcpp.latest, "). Update with:\n",
          "  pak::pak(\"DrCytometer/AutoSpectralRcpp\")"
        )
      }
    }
  }
}
