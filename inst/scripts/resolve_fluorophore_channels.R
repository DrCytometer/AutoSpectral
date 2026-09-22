# resolve_fluorophore_channels.R
#
# General-purpose channel resolution for representative-plot scripts: given
# a flow cytometry dataset and one or more target fluorophore identities,
# looks up which of that dataset's own column names corresponds to each
# fluorophore, using match.fluorophores() exactly as setup.unmix.comparison()
# already does internally to build its per-folder channel map. Two datasets
# holding the same fluorophore under different column-naming conventions
# (for example, classical vs AutoSpectral-unmixed exports) can then be
# plotted on the same biplot axes without assuming their column names agree.

#' @title Punctuation/Whitespace-Insensitive Fallback Channel Match
#'
#' @description
#' Private helper called by `resolve.fluorophore.channels()` when
#' `match.fluorophores()` canonicalises a requested fluorophore successfully
#' but cannot match it to any of `data`'s own column names - typically
#' because the export spells the channel in a way `match.fluorophores()`
#' does not recognise as a synonym (for example "PE-Vio770" with no space,
#' where `fluorophore.database` holds "PE-Vio 770" and its synonyms).
#' Column names and every known spelling of each unmatched fluorophore (its
#' canonical `fluorophore` value plus all `synonym1`-`synonym5` entries) are
#' compared after stripping everything but letters and digits and a
#' trailing "-A"/"-H"/"-W" acquisition suffix, so "PE-Vio770-A",
#' "PE-Vio 770-A", and "PE Vio 770" all normalise to the same name. A
#' fluorophore is only resolved this way when exactly one column normalises
#' to a match - an ambiguous normalised match is left for
#' `resolve.fluorophore.channels()`'s own error to report, the same as an
#' outright miss.
#'
#' @param data A matrix or data frame of flow cytometry data with named
#' columns.
#' @param missing.canonical Character vector of canonical `fluorophore`
#' values that `resolve.fluorophore.channels()` could not match to a column
#' of `data` through `match.fluorophores()`.
#' @param fluorophore.database Data frame of fluorophore names/synonyms, as
#' passed to `resolve.fluorophore.channels()`.
#'
#' @return A named character vector, a subset of `missing.canonical`
#' successfully resolved this way: names are the canonical fluorophore
#' values, values are the matched column name of `data`. Entries that
#' remain unresolved (no match, or more than one) are omitted.

.resolve.fluorophore.channels.fallback <- function(
    data,
    missing.canonical,
    fluorophore.database
) {

  normalize.name <- function( x ) {
    x <- sub( "-[AHW]$", "", x )
    toupper( gsub( "[^A-Za-z0-9]", "", x ) )
  }

  normalized.columns <- normalize.name( colnames( data ) )
  syn.cols <- c( "fluorophore", paste0( "synonym", 1:5 ) )

  fallback.map <- character( 0 )

  for ( fluor in missing.canonical ) {

    db.row <- fluorophore.database[ fluorophore.database$fluorophore == fluor, ]
    if ( nrow( db.row ) == 0 ) next

    candidate.names <- unlist( db.row[ 1, intersect( syn.cols, colnames( db.row ) ) ] )
    candidate.names <- candidate.names[ !is.na( candidate.names ) & candidate.names != "" ]
    if ( length( candidate.names ) == 0 ) next

    normalized.candidates <- normalize.name( candidate.names )
    matched.idx <- which( normalized.columns %in% normalized.candidates )

    if ( length( matched.idx ) == 1 ) {
      fallback.map[ fluor ] <- colnames( data )[ matched.idx ]
    }
  }

  fallback.map
}


#' @title Resolve Fluorophore Identities to a Dataset's Own Channel Names
#'
#' @description
#' Matches `data`'s own column names to fluorophore identity with
#' `match.fluorophores()` against the full `fluorophore.database`, and
#' returns the column name matched to each requested fluorophore - the same
#' channel-map construction `setup.unmix.comparison()` uses internally for
#' each folder it scans. The match is deliberately run against the full
#' database rather than one restricted to `fluorophore.names`: tandem dyes
#' sharing a name element (for example "APC", "APC-Fire 750", and
#' "APC-Fire 810") rely on `match.fluorophores()`'s longest-match rule to
#' tell their channels apart, and that rule can only prefer the longer name
#' when the longer name is still a candidate - dropping it by restricting
#' the database first would let every tandem's channel match the shorter
#' name instead. Two datasets that name the same fluorophore's channel
#' differently (for example, classical vs AutoSpectral-unmixed exports) can
#' then be resolved independently rather than assuming a shared naming
#' convention.
#'
#' `fluorophore.names` is itself canonicalised through `match.fluorophores()`
#' before use (the same first step `setup.unmix.comparison()` takes for its
#' own `fluorophore.list`), so a synonym works exactly like the database's
#' own `fluorophore` value: `match.fluorophores()` always reports the
#' canonical `fluorophore` entry for whatever it matches, never the synonym
#' text that happened to match, so a requested name that is itself only a
#' synonym (for example "Real Blue 744", canonical "RB744") would otherwise
#' never equal what a data column resolves to.
#'
#' If `match.fluorophores()` cannot match any column of `data` to an
#' otherwise-successfully-canonicalised requested fluorophore, a
#' punctuation/whitespace-insensitive fallback comparison against
#' `fluorophore.database`'s own `fluorophore`/`synonym` columns is tried
#' before giving up, so an export that spells a tandem's channel
#' differently from every synonym on file (for example "PE-Vio770" against
#' a database entry of "PE-Vio 770") can still be resolved.
#'
#' @param data A matrix or data frame of flow cytometry data with named
#' columns.
#' @param fluorophore.names Character vector of fluorophore identities to
#' resolve - either the `fluorophore` column value itself or any recognised
#' synonym, from `fluorophore_database.csv`.
#' @param fluorophore.database Data frame of fluorophore names/synonyms
#' passed to `match.fluorophores()`. Default `NULL` loads the bundled
#' `fluorophore_database.csv`.
#'
#' @seealso [match.fluorophores()]
#'
#' @return A named character vector: names are `fluorophore.names`, exactly
#' as supplied, values are the matched column name of `data` for each.
#'
#' @export

resolve.fluorophore.channels <- function(
    data,
    fluorophore.names,
    fluorophore.database = NULL
) {

  fluorophore.names <- unique( fluorophore.names )

  if ( is.null( fluorophore.database ) ) {
    fluor.data.path <- system.file(
      "extdata", "fluorophore_database.csv", package = "AutoSpectral"
    )
    fluorophore.database <- utils::read.csv( fluor.data.path )
    fluorophore.database[ fluorophore.database == "" ] <- NA
  }

  # canonicalise the requested identities first: a requested name may itself
  # be only a synonym (e.g. "Real Blue 744" for canonical "RB744"), and
  # match.fluorophores() always reports the database's own `fluorophore`
  # value for whatever it matches, never the synonym text - so comparing
  # against the requested name verbatim would silently never match.
  canonical.requested <- suppressMessages(
    match.fluorophores( fluorophore.names, fluorophore.database )
  )
  names( canonical.requested ) <- fluorophore.names

  unmatched.requested <- fluorophore.names[ canonical.requested == "No match" ]
  if ( length( unmatched.requested ) > 0 ) {
    stop(
      paste0(
        "The following fluorophore(s) could not be matched to any entry ",
        "(or synonym) in `fluorophore.database`: ",
        paste( unmatched.requested, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  channel.matches <- suppressMessages(
    match.fluorophores( colnames( data ), fluorophore.database )
  )

  matched.idx <- which( channel.matches %in% canonical.requested )
  channel.map <- stats::setNames(
    colnames( data )[ matched.idx ], channel.matches[ matched.idx ]
  )

  dup.fluor <- unique( names( channel.map )[ duplicated( names( channel.map ) ) ] )
  if ( length( dup.fluor ) > 0 ) {
    stop(
      paste0(
        "Multiple columns matched to: ", paste( dup.fluor, collapse = ", " ),
        ". Available columns: ", paste( colnames( data ), collapse = ", " )
      ),
      call. = FALSE
    )
  }

  missing.canonical <- setdiff( canonical.requested, names( channel.map ) )

  if ( length( missing.canonical ) > 0 ) {
    fallback.map <- .resolve.fluorophore.channels.fallback(
      data = data,
      missing.canonical = missing.canonical,
      fluorophore.database = fluorophore.database
    )
    channel.map <- c( channel.map, fallback.map )
    missing.canonical <- setdiff( missing.canonical, names( fallback.map ) )
  }

  if ( length( missing.canonical ) > 0 ) {
    missing.requested <- fluorophore.names[ canonical.requested %in% missing.canonical ]
    stop(
      paste0(
        "Could not match a column to: ", paste( missing.requested, collapse = ", " ),
        ". Available columns: ", paste( colnames( data ), collapse = ", " )
      ),
      call. = FALSE
    )
  }

  resolved <- channel.map[ canonical.requested[ fluorophore.names ] ]
  names( resolved ) <- fluorophore.names

  resolved
}
