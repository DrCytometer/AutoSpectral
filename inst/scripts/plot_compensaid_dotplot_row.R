# plot_compensaid_dotplot_row.R
#
# Finds an illustrative "exemplar" marker combination from
# compare.compensaid.folders()'s `pair.summary` table -- one flagged severe
# by a majority of non-reference operators/methods but not by the reference
# (typically AutoSpectral) -- and builds a row of CompensAID::PlotDotSSI()
# panels, one per folder, for that combination.
#
# Depends on .resolve.df() (in test_unmix_comparison.R) and
# .resolve.sample.files(), .channel.map.from.setup(), .compute.compensaid.raw()
# (all in compare_compensaid_folders.R) -- despite the filenames, these are
# production dependencies, not test-only code. All three files must be
# loaded into the package namespace together.


## Internal helper. Picks the replicate identifier from one folder's subset
## of the `results` table (already filtered to one primary/secondary pair)
## whose SSI best represents that folder's behaviour on this marker
## combination, per `replicate.select`. "median" (the default) picks the
## replicate closest to the median SSI, deliberately avoiding "worst", so a
## single unrepresentative replicate can't make a folder look better or
## worse than it typically is in a publication figure.
##
## @keywords internal
.select.compensaid.replicate <- function( pair.results, replicate.select ) {

  if ( nrow( pair.results ) == 0 ) {
    stop(
      "No `results` rows for this folder/primary/secondary combination.",
      call. = FALSE
    )
  }

  if ( replicate.select == "worst" ) {
    return( pair.results$replicate[ which.min( pair.results$ssi ) ] )
  }
  if ( replicate.select == "best" ) {
    return( pair.results$replicate[ which.max( pair.results$ssi ) ] )
  }

  med <- stats::median( pair.results$ssi )
  pair.results$replicate[ which.min( abs( pair.results$ssi - med ) ) ]
}


#' @title Find an Exemplar Marker Combination for a CompensAID Comparison Figure
#'
#' @description
#' Searches `compare.compensaid.folders()`'s `pair.summary` table for a
#' marker combination that illustrates a consistent difference between
#' `reference.folder` (typically `"AutoSpectral"`) and the other
#' operators/methods: one CompensAID flags as `severity.target` (default
#' `"Severe"`) by majority vote of replicates in more than
#' `min.other.majority.frac` of the non-reference folders, while
#' `reference.folder`'s own majority-vote severity is anything other than
#' `severity.target` (including not flagged at all).
#'
#' Eligible combinations are ranked by
#' `frac.other.majority.target - reference.frac.flagged`: the fraction of
#' non-reference folders calling it `severity.target` by majority vote,
#' minus the fraction of `reference.folder`'s own replicates that flag it
#' at all (any severity). This rewards combinations where the contrast is
#' both consistent across other folders and clean in the reference folder
#' (rarely or never flagged), rather than combinations that merely clear
#' the eligibility bar.
#'
#' @param pair.summary Either the `pair.summary` data frame returned by
#' `compare.compensaid.folders()` (or read back in from its
#' `pair.summary.csv`), or a character path to that CSV.
#' @param reference.folder Character. Name of the folder the others are
#' compared against. Default `"AutoSpectral"`.
#' @param severity.labels Character vector of severity band labels, most
#' severe first, matching `compare.compensaid.folders()`'s own
#' `severity.labels`. Default `c("Severe", "Moderate", "Mild")`.
#' @param severity.target Character, one of `severity.labels`. The band a
#' majority of non-reference folders must assign for a combination to be
#' eligible. Default `severity.labels[1]` (the most severe band).
#' @param min.other.majority.frac Numeric in `[0, 1)`, default `0.5`. A
#' combination is eligible only when strictly more than this fraction of
#' non-reference folders assign it `severity.target` by majority vote.
#' @param verbose Logical, default `TRUE`. Messages the winning combination
#' and why it was picked.
#'
#' @return Invisibly, a named list with `best` (a one-row data frame: the
#' winning `primary.fluorophore`/`secondary.fluorophore` plus its scoring
#' columns) and `candidates` (every eligible combination, ranked by score,
#' most to least illustrative).
#'
#' @importFrom stats median
#'
#' @seealso [compare.compensaid.folders()], [build.operator.compensaid.dotplot.row()]
#'
#' @export

find.compensaid.example.pair <- function(
    pair.summary,
    reference.folder        = "AutoSpectral",
    severity.labels          = c( "Severe", "Moderate", "Mild" ),
    severity.target          = severity.labels[ 1 ],
    min.other.majority.frac  = 0.5,
    verbose                  = TRUE
) {

  df <- .resolve.df( pair.summary, "pair.summary" )

  required.cols <- c(
    "folder", "primary.fluorophore", "secondary.fluorophore",
    "majority.severity", "frac.replicates.flagged"
  )
  missing.cols <- setdiff( required.cols, colnames( df ) )
  if ( length( missing.cols ) > 0 ) {
    stop(
      paste0(
        "`pair.summary` is missing required column(s): ",
        paste( missing.cols, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  if ( !severity.target %in% severity.labels ) {
    stop( "`severity.target` must be one of `severity.labels`.", call. = FALSE )
  }

  if ( !reference.folder %in% df$folder ) {
    stop(
      paste0(
        "`reference.folder` ('", reference.folder, "') not found in `pair.summary`."
      ),
      call. = FALSE
    )
  }

  other.folders <- setdiff( unique( df$folder ), reference.folder )
  if ( length( other.folders ) == 0 ) {
    stop( "No folders in `pair.summary` other than `reference.folder`.", call. = FALSE )
  }

  pair.keys <- unique( df[ , c( "primary.fluorophore", "secondary.fluorophore" ) ] )

  candidate.rows <- list()

  for ( i in seq_len( nrow( pair.keys ) ) ) {

    primary   <- pair.keys$primary.fluorophore[ i ]
    secondary <- pair.keys$secondary.fluorophore[ i ]

    row.mask <- df$primary.fluorophore == primary & df$secondary.fluorophore == secondary

    ref.row <- df[ row.mask & df$folder == reference.folder, ]
    if ( nrow( ref.row ) == 0 ) next
    ref.row <- ref.row[ 1, ]

    other.rows <- df[ row.mask & df$folder %in% other.folders, ]
    if ( nrow( other.rows ) < length( other.folders ) ) next

    frac.other.majority.target <-
      mean( other.rows$majority.severity == severity.target, na.rm = TRUE )

    ref.is.target <- !is.na( ref.row$majority.severity ) &&
      ref.row$majority.severity == severity.target

    eligible <- frac.other.majority.target > min.other.majority.frac && !ref.is.target
    if ( !eligible ) next

    score <- frac.other.majority.target - ref.row$frac.replicates.flagged

    candidate.rows[[ length( candidate.rows ) + 1 ]] <- data.frame(
      primary.fluorophore         = primary,
      secondary.fluorophore       = secondary,
      frac.other.majority.target  = frac.other.majority.target,
      reference.majority.severity = ref.row$majority.severity,
      reference.frac.flagged      = ref.row$frac.replicates.flagged,
      score                        = score,
      stringsAsFactors             = FALSE
    )
  }

  if ( length( candidate.rows ) == 0 ) {
    stop(
      paste0(
        "No marker combination is flagged '", severity.target, "' by a ",
        "majority of non-reference folders while '", reference.folder,
        "' does not also call it '", severity.target, "'. Try lowering ",
        "`min.other.majority.frac` or choosing a different `severity.target`."
      ),
      call. = FALSE
    )
  }

  candidates <- do.call( rbind, candidate.rows )
  candidates <- candidates[ order( -candidates$score ), ]
  rownames( candidates ) <- NULL

  best <- candidates[ 1, ]

  if ( verbose ) {
    ref.band <- if ( is.na( best$reference.majority.severity ) ) {
      "never flagged"
    } else {
      paste0( "worst band '", best$reference.majority.severity, "'" )
    }
    message( sprintf(
      "\033[34mExemplar pair: %s -> %s (majority '%s' in %.0f%% of non-reference folders; '%s' flags it in %.0f%% of replicates, %s)\033[0m",
      best$primary.fluorophore, best$secondary.fluorophore, severity.target,
      100 * best$frac.other.majority.target, reference.folder,
      100 * best$reference.frac.flagged, ref.band
    ) )
  }

  return( invisible( list( best = best, candidates = candidates ) ) )
}


#' @title Build One Folder's CompensAID Dot Plot for a Chosen Marker Combination
#'
#' @description
#' Re-runs CompensAID on one representative replicate of one folder and
#' returns `CompensAID::PlotDotSSI()`'s plot for the requested
#' primary/secondary marker combination -- the panel-building block behind
#' `build.operator.compensaid.dotplot.row()`. Re-running CompensAID (rather
#' than reusing `compare.compensaid.folders()`'s output) is necessary
#' because `PlotDotSSI()` needs the actual `res`/`ff` objects a run
#' produces, which the flattened `results` table does not retain.
#'
#' @param folder.path Character. Directory containing this folder's FCS
#' files, matching the corresponding entry of `compare.compensaid.folders()`'s
#' `folders` argument.
#' @param setup.csv.path Character. Path to this folder's
#' `setup.unmix.comparison()` setup CSV, matching the corresponding entry of
#' `compare.compensaid.folders()`'s `setup.files` argument.
#' @param sample.files The corresponding entry of
#' `compare.compensaid.folders()`'s `sample.files` argument for this folder
#' (a filename, or a named/unnamed vector of them).
#' @param primary.fluorophore,secondary.fluorophore Character. The marker
#' combination to plot, as named in `results`.
#' @param results Either this folder's subset of the long-format `results`
#' data frame returned by `compare.compensaid.folders()` (or read back in
#' from its `output.csv`), or a character path to that CSV; rows for other
#' folders, if present, are ignored automatically once `primary.fluorophore`/
#' `secondary.fluorophore` are matched.
#' @param transform.fun A function mapping raw channel values to a
#' logicle-style display scale, or `NULL` to skip transformation. Must
#' match what `compare.compensaid.folders()` was run with. Default
#' `biexp.transform()`.
#' @param replicate.select Character, one of `"median"` (default),
#' `"worst"`, or `"best"`. Which of this folder's replicates (by SSI on the
#' requested marker combination) to actually re-run and plot. `"median"` is
#' the safest default for a publication figure, since it can't be read as
#' cherry-picking the most dramatic replicate.
#' @param show.scores Logical, default `TRUE`. Passed to
#' `CompensAID::PlotDotSSI()`'s `showScores` argument.
#' @param subtitle Character or `NULL`. Added to the panel via
#' `ggplot2::labs(subtitle = )` when supplied -- typically the folder's
#' display label. Default `NULL`.
#' @param label Character, default `folder.path`. Used only in error/warning
#' messages, so a short folder name can be substituted for a long path.
#' @param timeout.seconds Numeric or `NULL` (default). Passed straight
#' through to `.compute.compensaid.raw()` -- see
#' `compare.compensaid.folders()`'s own `timeout.seconds` argument. Useful
#' here too, since this function re-runs CompensAID from scratch on the
#' chosen replicate rather than reusing `compare.compensaid.folders()`'s
#' output.
#' @param ... Additional arguments passed to `CompensAID::CompensAID()` via
#' `.compute.compensaid.raw()` (e.g. `segment.value`).
#'
#' @return A `ggplot` object (`CompensAID::PlotDotSSI()`'s panel for the
#' requested combination and replicate).
#'
#' @importFrom stats median
#' @importFrom ggplot2 labs
#'
#' @seealso [find.compensaid.example.pair()], [build.operator.compensaid.dotplot.row()],
#' [compare.compensaid.folders()]
#'
#' @export

build.operator.compensaid.dotplot <- function(
    folder.path,
    setup.csv.path,
    sample.files,
    primary.fluorophore,
    secondary.fluorophore,
    results,
    transform.fun     = biexp.transform(),
    replicate.select  = c( "median", "worst", "best" ),
    show.scores       = TRUE,
    subtitle          = NULL,
    label             = folder.path,
    timeout.seconds   = NULL,
    ...
) {

  replicate.select <- match.arg( replicate.select )

  results.df <- .resolve.df( results, "results" )

  pair.results <- results.df[
    results.df$primary.fluorophore   == primary.fluorophore &
      results.df$secondary.fluorophore == secondary.fluorophore,
  ]

  rep.id <- .select.compensaid.replicate( pair.results, replicate.select )

  replicate.set <- .resolve.sample.files( sample.files, label )
  sample.file   <- replicate.set[[ rep.id ]]

  full.path <- file.path( folder.path, sample.file )
  if ( !file.exists( full.path ) ) {
    stop(
      paste0(
        "'", label, "', replicate '", rep.id, "': sample file not found: '",
        full.path, "'."
      ),
      call. = FALSE
    )
  }

  channel.map <- .channel.map.from.setup( setup.csv.path, label )

  raw <- .compute.compensaid.raw(
    fcs.path        = full.path,
    channel.map     = channel.map,
    transform.fun   = transform.fun,
    label           = label,
    sample.file     = sample.file,
    timeout.seconds = timeout.seconds,
    ...
  )

  if ( is.null( raw ) ) {
    stop(
      paste0(
        "'", label, "', replicate '", rep.id, "': CompensAID could not be ",
        "re-run on this sample; see the warning above for the cause."
      ),
      call. = FALSE
    )
  }

  dot.plot <- CompensAID::PlotDotSSI(
    output     = raw$res,
    og         = raw$ff,
    primary    = primary.fluorophore,
    secondary  = secondary.fluorophore,
    showScores = show.scores
  )

  if ( !is.null( subtitle ) ) dot.plot <- dot.plot + ggplot2::labs( subtitle = subtitle )

  dot.plot
}


#' @title Build a Row of CompensAID Dot Plots, One Per Folder
#'
#' @description
#' Calls `build.operator.compensaid.dotplot()` once per folder in
#' `panel.order` for a single, shared marker combination (typically the
#' output of `find.compensaid.example.pair()`), and arranges the panels
#' side by side with `cowplot::plot_grid()`.
#'
#' @param folders Named character vector of directory paths, matching
#' `compare.compensaid.folders()`'s own `folders` argument.
#' @param setup.files Named character vector of setup CSV paths, matching
#' `compare.compensaid.folders()`'s own `setup.files` argument.
#' @param sample.files Named list of sample filename(s) per folder, matching
#' `compare.compensaid.folders()`'s own `sample.files` argument.
#' @param results The long-format `results` data frame returned by
#' `compare.compensaid.folders()` (or a path to its `output.csv`), covering
#' every folder in `panel.order`.
#' @param primary.fluorophore,secondary.fluorophore Character. The marker
#' combination to plot in every panel.
#' @param panel.order Character vector of folder names (keys of `folders`),
#' in the order the panels should appear left to right.
#' @param panel.labels Optional named list or character vector, names
#' matching `panel.order`, giving a display label per folder to use as each
#' panel's subtitle. Default `NULL` uses the `panel.order` key itself.
#' @param transform.fun,replicate.select,show.scores,... Passed through to
#' `build.operator.compensaid.dotplot()` for every panel.
#'
#' @return A `cowplot` grid object (the arranged row of panels).
#'
#' @importFrom cowplot plot_grid
#'
#' @seealso [find.compensaid.example.pair()], [build.operator.compensaid.dotplot()]
#'
#' @export

build.operator.compensaid.dotplot.row <- function(
    folders,
    setup.files,
    sample.files,
    results,
    primary.fluorophore,
    secondary.fluorophore,
    panel.order,
    panel.labels      = NULL,
    transform.fun     = biexp.transform(),
    replicate.select  = c( "median", "worst", "best" ),
    show.scores       = TRUE,
    ...
) {

  replicate.select <- match.arg( replicate.select )

  results.df <- .resolve.df( results, "results" )

  plots <- lapply( panel.order, function( key ) {

    subtitle <- key
    if ( !is.null( panel.labels ) && key %in% names( panel.labels ) ) {
      subtitle <- panel.labels[[ key ]]
    }

    build.operator.compensaid.dotplot(
      folder.path            = folders[[ key ]],
      setup.csv.path          = setup.files[[ key ]],
      sample.files            = sample.files[[ key ]],
      primary.fluorophore     = primary.fluorophore,
      secondary.fluorophore   = secondary.fluorophore,
      results                 = results.df[ results.df$folder == key, ],
      transform.fun           = transform.fun,
      replicate.select        = replicate.select,
      show.scores              = show.scores,
      subtitle                 = subtitle,
      label                    = key,
      ...
    )
  } )

  cowplot::plot_grid( plotlist = plots, nrow = 1 )
}
