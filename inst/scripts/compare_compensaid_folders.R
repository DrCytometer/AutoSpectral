# compare_compensaid_folders.R
#
# Runs CompensAID (Olsman et al., Cytometry A, 2026) on one or more
# fully-stained samples per already-unmixed folder from the unmix-comparison
# pipeline, to independently score reference-error artefacts per operator or
# unmixing method. Reuses the per-folder fluorophore -> channel maps already
# hand-verified in the setup CSVs written by setup.unmix.comparison(), so a
# fluorophore excluded there (bad file or channel match) is consistently
# excluded here too.
#
# Multiple fully-stained files per folder are treated as replicates of the
# same physical specimen(s), unmixed differently by each operator/method.
# Replicates are paired across folders by an explicit replicate identifier
# (the names of `sample.files`' per-folder entries) rather than by position,
# so the same specimen contributes a paired observation to every folder that
# unmixed it -- this is what test_compensaid_comparison.R's paired tests, and
# the error bars on plot.compensaid.comparison()'s barchart, are built on.
#
# Requires the Suggested (not Imported) packages 'CompensAID' and
# 'flowCore' -- both are only referenced via `pkg::fn()` and guarded by
# requireNamespace(), matching the AutoSpectralRcpp convention used
# elsewhere in this package. flowCore is needed solely to construct the
# flowFrame object CompensAID's public API expects; nothing here uses
# flowCore for FCS reading (readFCS() is used throughout, as usual).
#
# Depends on .resolve.df() (defined in test_unmix_comparison.R) and calls
# plot_compensaid_comparison.R's plot.compensaid.comparison() internally
# when make.plot = TRUE. .compute.compensaid.raw(), .channel.map.from.setup(),
# and .resolve.sample.files(), all defined below, are also called directly
# by plot_compensaid_dotplot_row.R's build.operator.compensaid.dotplot() to
# regenerate a single representative CompensAID run for the exemplar-pair
# figure panel, without duplicating this file's read/transform/flowFrame
# logic there. All four files (this one, plot_compensaid_comparison.R,
# test_compensaid_comparison.R, plot_compensaid_dotplot_row.R) must be
# loaded into the package namespace together.

## Internal helper. Builds a minimal flowFrame from an events x channels
## numeric matrix and sets its markernames from `channel.map` (named
## fluorophore -> channel), so CompensAID's flagged-combination output can
## be reported by fluorophore identity rather than raw detector name.
##
## @keywords internal
.build.compensaid.flowframe <- function( mat, channel.map ) {

  marker.by.channel <- stats::setNames( names( channel.map ), unname( channel.map ) )

  ff <- flowCore::flowFrame( exprs = mat )
  flowCore::markernames( ff ) <- marker.by.channel[ colnames( mat ) ]

  ff
}


## Internal helper. Resolves one folder's entry in `sample.files` to a named
## character vector (names = replicate identifiers used to pair the same
## physical specimen across folders; values = filenames within that
## folder). A bare single string is treated as one replicate. An unnamed
## vector is given default identifiers "Rep1", "Rep2", ... in order, with a
## warning, matching the "unnamed -> derive a default, warn" convention
## `setup.unmix.comparison()` already uses for `folders`. Position-based
## default names only pair correctly across folders when every folder's
## replicates are listed in the same specimen order -- name them explicitly
## whenever that isn't guaranteed.
##
## @keywords internal
.resolve.sample.files <- function( x, label ) {

  x <- stats::setNames( as.character( x ), names( x ) )

  if ( is.null( names( x ) ) || any( names( x ) == "" ) ) {
    warning(
      paste0(
        "Folder '", label, "': `sample.files` entry is unnamed (or ",
        "partially named); using \"Rep1\", \"Rep2\", ... as replicate ",
        "identifiers, in the order given. Name the entries explicitly ",
        "(e.g. c(Rep1 = \"...\", Rep2 = \"...\")) to pair replicates by ",
        "specimen identity across folders rather than by position."
      ),
      call. = FALSE
    )
    names( x ) <- paste0( "Rep", seq_along( x ) )
  }

  if ( any( duplicated( names( x ) ) ) ) {
    stop(
      paste0(
        "Folder '", label, "': `sample.files` has duplicated replicate ",
        "identifier(s): ",
        paste( unique( names( x )[ duplicated( names( x ) ) ] ), collapse = ", " )
      ),
      call. = FALSE
    )
  }

  x
}


## Internal helper. Majority-vote severity band for one (folder, pair)
## across its replicates: the most frequent `severity` value among flagged
## replicates, with ties broken toward the more severe band (`severity.labels`
## is ordered most-to-least severe, so the first tied label wins). Returns
## `NA_character_` when no replicate was flagged.
##
## @keywords internal
.majority.severity <- function( severity.vals, severity.labels ) {

  severity.vals <- severity.vals[ !is.na( severity.vals ) ]
  if ( length( severity.vals ) == 0 ) return( NA_character_ )

  counts <- table( factor( severity.vals, levels = severity.labels ) )
  best.n <- max( counts )

  severity.labels[ which( counts == best.n )[ 1 ] ]
}


## Internal helper. Reads one fully-stained sample FCS file, applies
## `transform.fun` to the channels named in `channel.map`, builds the
## flowFrame CompensAID expects (via .build.compensaid.flowframe()), and
## runs `CompensAID::CompensAID()`. Returns `list(ff = ff, res = res)`, or
## `NULL` (with a warning) on any failure -- missing file, missing channels,
## non-finite values after `transform.fun`, an error inside
## `CompensAID::CompensAID()`, or (when `timeout.seconds` is set) that call
## exceeding its time budget -- so one bad or pathological replicate doesn't
## halt or hang the rest of the comparison. Kept separate from
## .flatten.compensaid.result() so plot_compensaid_dotplot_row.R can call
## this directly to get the raw `ff`/`res` a single PlotDotSSI() panel
## needs, without re-deriving CompensAID's own SSI matrix by hand.
##
## The timeout uses base R's `setTimeLimit()` (`transient = TRUE`, reset via
## `on.exit()` either way), so it needs no extra dependency. It is
## best-effort: R only checks elapsed time at safe points between R-level
## calls, so a hang stuck entirely inside compiled C code without returning
## to R would not be interrupted. `flowDensity::deGate()`/`stats::density()`
## (CompensAID's per-channel gating step, and the most likely place for a
## pathological sample -- e.g. one with extreme spillover-corrected outliers
## from a weaker unmixing method -- to hang) are R-level with frequent
## call-outs, so this should catch the common case.
##
## @keywords internal
.compute.compensaid.raw <- function(
    fcs.path,
    channel.map,
    transform.fun,
    label,
    sample.file,
    timeout.seconds = NULL,
    ...
) {

  needed.channels <- unname( channel.map )

  mat <- tryCatch(
    readFCS( fcs.path, columns = needed.channels ),
    error = function( e ) {
      warning(
        paste0(
          "Folder '", label, "', sample '", sample.file,
          "': failed to read requested channels (", conditionMessage( e ),
          "); skipping."
        ),
        call. = FALSE
      )
      NULL
    }
  )

  if ( is.null( mat ) ) return( NULL )

  if ( !is.null( transform.fun ) ) {
    mat <- apply( mat, 2, transform.fun )
    colnames( mat ) <- needed.channels
  }

  n.nonfinite <- sum( !is.finite( mat ) )
  if ( n.nonfinite > 0 ) {
    warning(
      paste0(
        "Folder '", label, "', sample '", sample.file, "': ", n.nonfinite,
        " of ", length( mat ), " transformed value(s) are non-finite ",
        "(NA/NaN/Inf); CompensAID's per-channel density gating may hang or ",
        "behave unpredictably on this sample."
      ),
      call. = FALSE
    )
  }

  ff <- .build.compensaid.flowframe( mat, channel.map )

  res <- tryCatch(
    {
      if ( !is.null( timeout.seconds ) ) {
        setTimeLimit( elapsed = timeout.seconds, transient = TRUE )
        on.exit( setTimeLimit( elapsed = Inf, transient = TRUE ), add = TRUE )
      }
      CompensAID::CompensAID( ff = ff, ... )
    },
    error = function( e ) {
      warning(
        paste0(
          "Folder '", label, "', sample '", sample.file,
          "': CompensAID::CompensAID() failed (", conditionMessage( e ),
          "); skipping."
        ),
        call. = FALSE
      )
      NULL
    }
  )

  if ( is.null( res ) || is.null( res[[ "matrix" ]] ) ) return( NULL )

  list( ff = ff, res = res )
}


## Internal helper. Flattens a CompensAID result's SSI matrix into a
## long-format data frame with one row per off-diagonal (primary,
## secondary) channel pair, mapped back to fluorophore identity via
## `channel.map`.
##
## The (row, col) -> (secondary, primary) mapping and the assumption that
## `rownames(res$matrix)` and `colnames(res$matrix)` are identical channel
## sets both follow the CompensAID vignette's own worked example exactly
## (https://bioc.r-universe.dev/CompensAID/doc/CompensAID.Rmd); this has
## not been run against a live install, so it is worth a spot check against
## `?CompensAID::CompensAID` on the first real run.
##
## @keywords internal
.flatten.compensaid.result <- function( res, channel.map, label, sample.file ) {

  ssi.mat <- res[[ "matrix" ]]
  ch      <- rownames( ssi.mat )

  if ( is.null( ch ) ) {
    warning(
      paste0(
        "Folder '", label, "', sample '", sample.file,
        "': CompensAID output matrix has no rownames; cannot map channels ",
        "back to fluorophores. Skipping."
      ),
      call. = FALSE
    )
    return( NULL )
  }

  if ( !identical( rownames( ssi.mat ), colnames( ssi.mat ) ) ) {
    warning(
      paste0(
        "Folder '", label, "', sample '", sample.file,
        "': CompensAID output matrix rownames and colnames differ; ",
        "primary/secondary channel identity below may be unreliable. ",
        "Proceeding using rownames for both, matching the CompensAID ",
        "vignette's own convention -- verify against a live run."
      ),
      call. = FALSE
    )
  }

  n       <- length( ch )
  row.idx <- rep( seq_len( n ), times = n )
  col.idx <- rep( seq_len( n ), each  = n )
  keep    <- row.idx != col.idx

  primary.channel   <- ch[ col.idx[ keep ] ]
  secondary.channel <- ch[ row.idx[ keep ] ]
  ssi.value         <- ssi.mat[ cbind( row.idx[ keep ], col.idx[ keep ] ) ]

  ok <- is.finite( ssi.value )

  marker.by.channel <- stats::setNames( names( channel.map ), unname( channel.map ) )

  data.frame(
    primary.fluorophore   = unname( marker.by.channel[ primary.channel[ ok ] ] ),
    secondary.fluorophore = unname( marker.by.channel[ secondary.channel[ ok ] ] ),
    ssi                   = unname( ssi.value[ ok ] ),
    stringsAsFactors      = FALSE
  )
}


## Internal helper. Derives one folder's fluorophore -> channel map from its
## setup.unmix.comparison() CSV, keeping only rows CompensAID can safely use
## (excludes "Unstained" and anything flagged for exclusion or with no
## channel match). Shared by compare.compensaid.folders()'s main loop and
## plot_compensaid_dotplot_row.R's build.operator.compensaid.dotplot(), so
## the two never derive a different map for the same folder. `stop()`s
## (rather than warns) when nothing usable is found, so callers that want
## the original "warn and skip the folder" behaviour wrap this in
## `tryCatch()` themselves.
##
## @keywords internal
.channel.map.from.setup <- function( setup.csv.path, label ) {

  tb <- utils::read.csv( setup.csv.path, stringsAsFactors = FALSE, strip.white = TRUE )
  tb[ tb == "" ] <- NA

  ok.rows <- tb[
    tb$fluorophore != "Unstained" &
      tb$flag == "OK" &
      !is.na( tb$channel ) & tb$channel != "No match"
    , , drop = FALSE
  ]

  if ( nrow( ok.rows ) == 0 ) {
    stop(
      paste0( "Folder '", label, "': no usable fluorophore/channel mappings." ),
      call. = FALSE
    )
  }

  stats::setNames( ok.rows$channel, ok.rows$fluorophore )
}


## Internal helper. Thin wrapper combining .compute.compensaid.raw() and
## .flatten.compensaid.result() for one replicate FCS file, for use inside
## compare.compensaid.folders()'s main loop.
##
## @keywords internal
.run.compensaid.on.file <- function(
    fcs.path,
    channel.map,
    transform.fun,
    label,
    sample.file,
    timeout.seconds = NULL,
    ...
) {

  raw <- .compute.compensaid.raw(
    fcs.path        = fcs.path,
    channel.map     = channel.map,
    transform.fun   = transform.fun,
    label           = label,
    sample.file     = sample.file,
    timeout.seconds = timeout.seconds,
    ...
  )

  if ( is.null( raw ) ) return( NULL )

  .flatten.compensaid.result(
    res = raw$res, channel.map = channel.map, label = label, sample.file = sample.file
  )
}


#' @title Compare CompensAID-Flagged Reference Errors Across Operators or Methods
#'
#' @description
#' Runs CompensAID (Olsman et al., \emph{Cytometry A}, 2026), an independent,
#' published quality-control tool that flags marker combinations showing
#' signs of reference errors (miscalibrated single-stain controls), on one or
#' more fully-stained samples per folder already scanned by
#' `setup.unmix.comparison()`. Unlike `compare.unmix.folders()`, which scores
#' single-stained controls against the unstained reference, this function
#' scores the actual multi-color sample(s) each operator or unmixing method
#' produced, using CompensAID's own automatic (density-based) gating of the
#' positive and negative population in every channel -- no unstained control
#' is needed here.
#'
#' For every folder, the fluorophore -> channel map is taken from the
#' `flag == "OK"` rows of that folder's setup CSV (the same map
#' `compare.unmix.folders()` uses), so results are reported by fluorophore
#' identity and stay comparable across folders that use different detector
#' naming conventions. Each fully-stained sample is optionally transformed
#' (`transform.fun`, default `biexp.transform()`, matching the logicle-style
#' transform CompensAID's own pre-processing pipeline expects) before being
#' wrapped in a minimal `flowCore::flowFrame` and passed to
#' `CompensAID::CompensAID()`.
#'
#' Multiple samples per folder (`sample.files`) are treated as replicates of
#' the same physical specimen(s), paired across folders by the replicate
#' identifiers named in `sample.files` -- so the same specimen, unmixed by
#' every operator/method, contributes one paired observation per folder.
#' This pairing is what `test.compensaid.comparison()` and
#' `test.compensaid.comparison.pairs()` test on, and what the error bars on
#' `plot.compensaid.comparison()`'s barchart are computed across.
#'
#' A replicate whose transformed values contain any non-finite entries
#' (NA/NaN/Inf) triggers a warning naming the folder and sample before
#' CompensAID is run on it -- these are the replicates most likely to make
#' CompensAID's per-channel density gating misbehave (see `timeout.seconds`)
#' since they typically come from a weaker unmixing method's most extreme
#' spillover-correction outliers, exactly what CompensAID exists to flag.
#'
#' Four tables are produced, each coarser than the last:
#' \enumerate{
#'   \item `results`: one row per (folder, replicate, primary fluorophore,
#'     secondary fluorophore) -- every off-diagonal SSI value CompensAID
#'     returned, flagged when `ssi < flag.threshold` (default `-1`, matching
#'     CompensAID's own worked example) and binned into a severity band.
#'   \item `replicate.summary`: one row per (folder, replicate) -- how many
#'     marker combinations were tested and flagged in that one sample, and
#'     how many fell in each severity band. This is the table
#'     `test.compensaid.comparison()` pairs across folders.
#'   \item `pair.summary`: one row per (folder, primary, secondary) --
#'     that combination's status rolled up across the folder's replicates
#'     (majority-vote flagged/severity, plus the worst, best, and median SSI
#'     seen). This is the table `find.compensaid.example.pair()` (in
#'     `plot_compensaid_dotplot_row.R`) searches for the figure's exemplar
#'     marker combination.
#'   \item `summary`: one row per folder -- the mean and standard error (or
#'     standard deviation; see `error.type`) of `replicate.summary`'s counts
#'     across that folder's replicates. This is what
#'     `plot.compensaid.comparison()`'s grouped barchart plots.
#' }
#'
#' @param folders Named character vector of directory paths, one per
#' unmixed result set (must use the same names as `setup.unmix.comparison()`
#' and `compare.unmix.folders()`).
#' @param setup.files Named character vector of paths to the setup CSVs
#' written by `setup.unmix.comparison()`. Names must match `folders`.
#' @param sample.files Named list, names matching `folders`. Each element is
#' the fully-stained sample FCS filename(s) for that folder (paths relative
#' to that folder's directory) -- the real multi-color specimen(s) to score,
#' not the single-stain controls used by `compare.unmix.folders()`. Name the
#' elements of each per-folder vector with a shared replicate identifier
#' (e.g. `c(Rep1 = "...", Rep2 = "...")`) when more than one folder unmixed
#' the same physical specimen(s), so replicates pair correctly across
#' folders; an unnamed entry is given default identifiers "Rep1", "Rep2", ...
#' in order, with a warning, which only pairs correctly if every folder
#' lists its replicates in the same specimen order.
#' @param flag.threshold Numeric, default `-1`. A (primary, secondary) pair
#' is flagged when its SSI is strictly less than this value.
#' @param severity.breaks Numeric vector of `cut()` breaks used to bin
#' flagged pairs into severity bands, most negative first. Default `NULL`
#' derives `c(-Inf, flag.threshold - 2, flag.threshold - 1, flag.threshold)`,
#' i.e. three bands of width 1 immediately below `flag.threshold`.
#' @param severity.labels Character vector of band labels, most severe
#' first, one shorter than `severity.breaks`. Default
#' `c("Severe", "Moderate", "Mild")`.
#' @param error.type Character, one of `"sem"` (standard error of the mean,
#' the default) or `"sd"` (standard deviation). Controls what `summary`'s
#' `se.*` columns -- and hence `plot.compensaid.comparison()`'s error bars --
#' represent. A folder with only one replicate gets `NA` for every `se.*`
#' column (no error bar is drawn for it).
#' @param transform.fun A function mapping raw channel values to a
#' logicle-style display scale, applied to every channel before CompensAID
#' is run, or `NULL` to skip transformation. Default `biexp.transform()`
#' (this package's own logicle implementation) with its own defaults.
#' @param timeout.seconds Numeric or `NULL` (default). When set, each
#' `CompensAID::CompensAID()` call is given at most this many seconds
#' (`base::setTimeLimit()`); a replicate that exceeds it is warned about and
#' skipped, exactly like a read failure, rather than blocking the rest of
#' the comparison indefinitely. Best-effort: only reliably interrupts a hang
#' at an R-level call, which covers CompensAID's own density-gating step but
#' not a hang stuck entirely inside compiled C code. A sample that hangs
#' here is usually flagged by a preceding non-finite-value warning (see
#' Description) -- CompensAID's per-channel density gating is a plausible
#' place for such a sample to misbehave.
#' @param plot.dir Character. Directory for output figures. Created if
#' absent. Default `"./figure_compensaid_comparison"`.
#' @param output.csv Character. Path for the `results` table. Default
#' `"compensaid_comparison_results.csv"`.
#' @param replicate.summary.csv Character. Path for the `replicate.summary`
#' table. Default `"compensaid_comparison_replicate_summary.csv"`.
#' @param pair.summary.csv Character. Path for the `pair.summary` table.
#' Default `"compensaid_comparison_pair_summary.csv"`.
#' @param summary.csv Character. Path for the `summary` table the barchart
#' is built from. Default `"compensaid_comparison_summary.csv"`.
#' @param make.plot Logical, default `TRUE`. Calls `plot.compensaid.comparison()`
#' on the results internally.
#' @param normalize.plot Logical, default `FALSE`. Passed to
#' `plot.compensaid.comparison()`; adds a second figure of the percentage
#' (rather than raw count) of tested pairs flagged per folder, useful when
#' folders don't all have the same number of usable fluorophores.
#' @param plot.width,plot.height Numeric, defaults `7` and `5` (inches).
#' @param base.font.size Numeric, default `11`.
#' @param title.size Numeric, default `NULL`.
#' @param text.angle Numeric, default `45`. Rotation angle (degrees) of the
#' x-axis (folder) labels.
#' @param verbose Logical, default `TRUE`.
#' @param ... Additional arguments passed straight through to
#' `CompensAID::CompensAID()` (e.g. `segment.value`), for anything beyond
#' `ff` that a newer or older CompensAID release exposes.
#'
#' @return Invisibly, a named list with `results`, `replicate.summary`,
#' `pair.summary`, and `summary` (see Description), each also written to its
#' corresponding `*.csv` argument. Figures are written to `plot.dir` when
#' `make.plot = TRUE`.
#'
#' @importFrom stats aggregate setNames sd median
#'
#' @seealso [setup.unmix.comparison()], [compare.unmix.folders()],
#' [plot.compensaid.comparison()], [find.compensaid.example.pair()]
#'
#' @export

compare.compensaid.folders <- function(
    folders,
    setup.files,
    sample.files,
    flag.threshold        = -1,
    severity.breaks       = NULL,
    severity.labels       = c( "Severe", "Moderate", "Mild" ),
    error.type            = c( "sem", "sd" ),
    transform.fun         = biexp.transform(),
    timeout.seconds       = NULL,
    plot.dir              = "./figure_compensaid_comparison",
    output.csv            = "compensaid_comparison_results.csv",
    replicate.summary.csv = "compensaid_comparison_replicate_summary.csv",
    pair.summary.csv      = "compensaid_comparison_pair_summary.csv",
    summary.csv           = "compensaid_comparison_summary.csv",
    make.plot             = TRUE,
    normalize.plot        = FALSE,
    plot.width            = 7,
    plot.height           = 5,
    base.font.size        = 11,
    title.size            = NULL,
    text.angle            = 45,
    verbose               = TRUE,
    ...
) {

  error.type <- match.arg( error.type )

  error.fun <- if ( error.type == "sem" ) {
    function( x ) stats::sd( x ) / sqrt( length( x ) )
  } else {
    stats::sd
  }

  # --- optional-dependency checks ------------------------------------------

  if ( !requireNamespace( "CompensAID", quietly = TRUE ) ) {
    stop(
      paste0(
        "Package 'CompensAID' is required but not installed. Install it ",
        "with `devtools::install_github(\"Olsman/CompensAID\")`, or ",
        "`BiocManager::install(\"CompensAID\")` now that it has reached ",
        "Bioconductor release."
      ),
      call. = FALSE
    )
  }

  if ( !requireNamespace( "flowCore", quietly = TRUE ) ) {
    stop(
      paste0(
        "Package 'flowCore' is required but not installed (CompensAID's ",
        "public API expects a flowCore::flowFrame as input)."
      ),
      call. = FALSE
    )
  }

  # --- input validation ------------------------------------------------------

  if ( !setequal( names( folders ), names( setup.files ) ) ) {
    stop( "`folders` and `setup.files` must have identical names.", call. = FALSE )
  }
  if ( !setequal( names( folders ), names( sample.files ) ) ) {
    stop( "`folders` and `sample.files` must have identical names.", call. = FALSE )
  }

  missing.setup <- setup.files[ !file.exists( setup.files ) ]
  if ( length( missing.setup ) > 0 ) {
    stop(
      paste0( "Setup CSV(s) not found: ", paste( missing.setup, collapse = ", " ) ),
      call. = FALSE
    )
  }

  if ( length( severity.labels ) < 1 ) {
    stop( "`severity.labels` must have at least one entry.", call. = FALSE )
  }

  if ( is.null( severity.breaks ) ) {
    severity.breaks <- c(
      -Inf, flag.threshold - 2, flag.threshold - 1, flag.threshold
    )
  }
  if ( length( severity.breaks ) != length( severity.labels ) + 1 ) {
    stop(
      "`severity.breaks` must have exactly one more entry than `severity.labels`.",
      call. = FALSE
    )
  }

  labels <- names( folders )

  # --- process each folder ---------------------------------------------------

  all.results <- list()

  for ( label in labels ) {

    if ( verbose ) message( sprintf( "\033[34mProcessing folder: %s\033[0m", label ) )

    folder.path <- folders[[ label ]]

    channel.map <- tryCatch(
      .channel.map.from.setup( setup.files[[ label ]], label ),
      error = function( e ) {
        warning(
          paste0( "Folder '", label, "': ", conditionMessage( e ), " Skipping folder." ),
          call. = FALSE
        )
        NULL
      }
    )
    if ( is.null( channel.map ) ) next

    replicate.set <- .resolve.sample.files( sample.files[[ label ]], label )

    for ( rep.id in names( replicate.set ) ) {

      sample.file <- replicate.set[[ rep.id ]]

      if ( verbose ) message( sprintf( "  Replicate %s: %s", rep.id, sample.file ) )

      full.path <- file.path( folder.path, sample.file )

      if ( !file.exists( full.path ) ) {
        warning(
          paste0(
            "Folder '", label, "', replicate '", rep.id, "': sample file not ",
            "found: '", full.path, "'; skipping."
          ),
          call. = FALSE
        )
        next
      }

      pair.results <- .run.compensaid.on.file(
        fcs.path        = full.path,
        channel.map     = channel.map,
        transform.fun   = transform.fun,
        label           = label,
        sample.file     = sample.file,
        timeout.seconds = timeout.seconds,
        ...
      )

      if ( is.null( pair.results ) ) next

      pair.results$folder      <- label
      pair.results$replicate   <- rep.id
      pair.results$sample.file <- sample.file

      all.results[[ length( all.results ) + 1 ]] <- pair.results
    }
  }

  if ( length( all.results ) == 0 ) {
    stop( "No usable CompensAID results were produced from any folder.", call. = FALSE )
  }

  results.df <- do.call( rbind, all.results )
  rownames( results.df ) <- NULL

  results.df$flagged  <- results.df$ssi < flag.threshold
  results.df$severity <- NA_character_
  if ( any( results.df$flagged ) ) {
    results.df$severity[ results.df$flagged ] <- as.character(
      cut(
        results.df$ssi[ results.df$flagged ],
        breaks         = severity.breaks,
        labels         = severity.labels,
        include.lowest = TRUE
      )
    )
  }

  results.df <- results.df[ , c(
    "folder", "replicate", "sample.file", "primary.fluorophore",
    "secondary.fluorophore", "ssi", "flagged", "severity"
  ) ]

  utils::write.csv( results.df, output.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote results: %s\033[0m", output.csv ) )

  # --- replicate.summary: one row per (folder, replicate) -------------------

  replicate.summary.rows <- list()

  for ( label in labels ) {
    for ( rep.id in unique( results.df$replicate[ results.df$folder == label ] ) ) {

      rep.rows <- results.df[
        results.df$folder == label & results.df$replicate == rep.id,
      ]
      if ( nrow( rep.rows ) == 0 ) next

      row <- data.frame(
        folder       = label,
        replicate    = rep.id,
        sample.file  = rep.rows$sample.file[ 1 ],
        n.tested     = nrow( rep.rows ),
        n.flagged    = sum( rep.rows$flagged ),
        stringsAsFactors = FALSE
      )
      row$frac.flagged <- row$n.flagged / row$n.tested

      for ( s in severity.labels ) {
        row[[ s ]] <- sum( rep.rows$flagged & rep.rows$severity == s, na.rm = TRUE )
      }

      replicate.summary.rows[[ length( replicate.summary.rows ) + 1 ]] <- row
    }
  }

  replicate.summary.df <- do.call( rbind, replicate.summary.rows )
  rownames( replicate.summary.df ) <- NULL

  utils::write.csv( replicate.summary.df, replicate.summary.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote replicate summary: %s\033[0m", replicate.summary.csv ) )

  # --- pair.summary: one row per (folder, primary, secondary), rolled up
  # across that folder's replicates -----------------------------------------

  pair.n         <- stats::aggregate(
    ssi ~ folder + primary.fluorophore + secondary.fluorophore,
    data = results.df, FUN = length
  )
  names( pair.n )[ names( pair.n ) == "ssi" ] <- "n.replicates"

  pair.n.flagged <- stats::aggregate(
    flagged ~ folder + primary.fluorophore + secondary.fluorophore,
    data = results.df, FUN = sum
  )
  names( pair.n.flagged )[ names( pair.n.flagged ) == "flagged" ] <- "n.replicates.flagged"

  pair.median <- stats::aggregate(
    ssi ~ folder + primary.fluorophore + secondary.fluorophore,
    data = results.df, FUN = stats::median
  )
  names( pair.median )[ names( pair.median ) == "ssi" ] <- "median.ssi"

  pair.worst <- stats::aggregate(
    ssi ~ folder + primary.fluorophore + secondary.fluorophore,
    data = results.df, FUN = min
  )
  names( pair.worst )[ names( pair.worst ) == "ssi" ] <- "worst.ssi"

  pair.best <- stats::aggregate(
    ssi ~ folder + primary.fluorophore + secondary.fluorophore,
    data = results.df, FUN = max
  )
  names( pair.best )[ names( pair.best ) == "ssi" ] <- "best.ssi"

  by.cols <- c( "folder", "primary.fluorophore", "secondary.fluorophore" )

  pair.summary.df <- Reduce(
    function( a, b ) merge( a, b, by = by.cols ),
    list( pair.n, pair.n.flagged, pair.median, pair.worst, pair.best )
  )

  pair.summary.df$frac.replicates.flagged <-
    pair.summary.df$n.replicates.flagged / pair.summary.df$n.replicates
  pair.summary.df$majority.flagged <- pair.summary.df$frac.replicates.flagged > 0.5

  pair.summary.df$majority.severity <- NA_character_
  for ( i in seq_len( nrow( pair.summary.df ) ) ) {
    if ( !pair.summary.df$majority.flagged[ i ] ) next
    row.mask <- results.df$folder                 == pair.summary.df$folder[ i ] &
                results.df$primary.fluorophore     == pair.summary.df$primary.fluorophore[ i ] &
                results.df$secondary.fluorophore   == pair.summary.df$secondary.fluorophore[ i ] &
                results.df$flagged
    pair.summary.df$majority.severity[ i ] <- .majority.severity(
      results.df$severity[ row.mask ], severity.labels
    )
  }

  rownames( pair.summary.df ) <- NULL

  utils::write.csv( pair.summary.df, pair.summary.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote pair summary: %s\033[0m", pair.summary.csv ) )

  # --- summary: one row per folder, mean +/- error.type across replicates,
  # the table plot.compensaid.comparison() plots -----------------------------

  summary.rows <- list()

  for ( label in labels ) {

    lab.rows <- replicate.summary.df[ replicate.summary.df$folder == label, ]
    if ( nrow( lab.rows ) == 0 ) next

    n.rep <- nrow( lab.rows )

    row <- data.frame( folder = label, n.replicates = n.rep, stringsAsFactors = FALSE )

    row$mean.n.flagged <- mean( lab.rows$n.flagged )
    row$se.n.flagged   <- if ( n.rep > 1 ) error.fun( lab.rows$n.flagged ) else NA_real_

    row$mean.frac.flagged <- mean( lab.rows$frac.flagged )
    row$se.frac.flagged   <- if ( n.rep > 1 ) error.fun( lab.rows$frac.flagged ) else NA_real_

    for ( s in severity.labels ) {
      vals <- lab.rows[[ s ]]
      row[[ paste0( "mean.", s ) ]] <- mean( vals )
      row[[ paste0( "se.",   s ) ]] <- if ( n.rep > 1 ) error.fun( vals ) else NA_real_
    }

    summary.rows[[ length( summary.rows ) + 1 ]] <- row
  }

  summary.df <- do.call( rbind, summary.rows )
  rownames( summary.df ) <- NULL

  utils::write.csv( summary.df, summary.csv, row.names = FALSE )
  if ( verbose ) message( sprintf( "\033[32mWrote summary: %s\033[0m", summary.csv ) )

  # --- plot --------------------------------------------------------------
  # delegated to plot.compensaid.comparison() so the same plotting code can
  # be re-run standalone later, without repeating CompensAID above

  if ( make.plot ) {
    plot.compensaid.comparison(
      summary         = summary.df,
      folder.levels   = labels,
      severity.labels = severity.labels,
      error.label     = if ( error.type == "sem" ) "SEM" else "SD",
      plot.dir        = plot.dir,
      normalize       = normalize.plot,
      plot.width      = plot.width,
      plot.height     = plot.height,
      base.font.size  = base.font.size,
      title.size      = title.size,
      text.angle      = text.angle,
      verbose         = verbose
    )
  }

  return( invisible( list(
    results            = results.df,
    replicate.summary  = replicate.summary.df,
    pair.summary       = pair.summary.df,
    summary            = summary.df
  ) ) )
}
