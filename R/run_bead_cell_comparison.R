# run_bead_cell_comparison.R
#
# Compares single-color control spectra measured on different particle types
# (e.g. compensation beads) against those measured on Cells, using the
# automated (non-gated) AutoSpectral extraction pipeline throughout.
#
# Expected directory layout: one folder per particle type, each containing
# single-stained FCS files named by fluorophore plus one "Unstained" file
# used as the background/AF reference for that folder.

# ---------------------------------------------------------------------------
# Control file helper
# ---------------------------------------------------------------------------

#' Finalize an auto-generated control file for bead/cell comparison
#'
#' Post-processes a control file created by \code{create.control.file()} so it
#' is usable by either the legacy or automated AutoSpectral pipeline in the
#' context of a bead-vs-cell comparison run. Sets \code{control.type} for every
#' row according to the particle type of the folder, forces the single
#' unstained/AF row to \code{control.type == "cells"} (required by both
#' pipelines' validation), and back-fills \code{universal.negative} for the
#' remaining rows with that row's filename.
#'
#' @param control.file.path Character scalar. Path to the control file CSV to
#'   finalize in place.
#' @param particle.type Character scalar. The particle type label for the
#'   folder this control file belongs to (e.g. \code{"Cells"}, \code{"Beads"}).
#' @param reference.type Character scalar. The particle type treated as the
#'   reference (cell) population; any \code{particle.type} not identical to
#'   this is treated as beads. Defaults to \code{"Cells"}.
#'
#' @details
#' Exactly one row must resolve to \code{fluorophore \%in\% c("AF", "Negative")}.
#' If zero or more than one such row is found, the function leaves
#' \code{fluorophore}/\code{universal.negative} unset for that file and emits a
#' warning instructing the user to edit the control file manually rather than
#' guessing which row is the background control.
#'
#' @return Invisibly, the finalized control file as a data frame. The same
#'   data is also written back to \code{control.file.path}.
#'
#' @keywords internal

finalize.control.file <- function( control.file.path, particle.type, reference.type = "Cells" ) {

  ctrl <- utils::read.csv(
    control.file.path, stringsAsFactors = FALSE, colClasses = "character"
  )

  # Fluorophore rows are labelled with the folder's true particle type
  # (control.type == "cells"/"beads"), which the legacy gating pipeline
  # (define.flow.control()/assign.gates()/do.gate()) needs in order to pick
  # the right gate parameters (asp$*.cells vs asp$*.beads). The single
  # unstained background row is the exception: both pipelines validate that
  # the row labelled "AF" has control.type == "cells" regardless of what
  # the particles actually are (see validate_control_file.R, af_wrong_type),
  # so that row alone is forced to "cells". This also means the automated
  # pipeline's per-file AF-PC refinement (which only runs for
  # control.type == "cells" universal negatives) is correctly skipped for
  # bead folders instead of being incorrectly applied to bead data.
  is.beads <- !identical( particle.type, reference.type )

  ctrl$control.type <- if ( is.beads ) "beads" else "cells"

  unstained.idx <- which( ctrl$fluorophore %in% c( "AF", "Negative" ) )

  if ( length( unstained.idx ) == 1 ) {
    ctrl$fluorophore[ unstained.idx ]  <- "AF"
    ctrl$control.type[ unstained.idx ] <- "cells"
    ctrl$universal.negative[ -unstained.idx ] <- ctrl$filename[ unstained.idx ]
  } else if ( length( unstained.idx ) > 1 ) {
    warning(
      "Multiple AF/Negative controls found for '", particle.type,
      "'; leaving fluorophore/universal.negative unset -- edit ",
      control.file.path, " manually.", call. = FALSE
    )
  } else {
    warning(
      "No unstained control found for '", particle.type,
      "'; leaving universal.negative unset -- edit ",
      control.file.path, " manually.", call. = FALSE
    )
  }

  utils::write.csv( ctrl, control.file.path, row.names = FALSE )
  invisible( ctrl )
}

#' Create a control file inside a particle-type directory
#'
#' Thin wrapper around \code{create.control.file()} that temporarily changes
#' the working directory to \code{p.dir} for the duration of the call, since
#' \code{create.control.file()} discovers FCS files relative to the working
#' directory. The working directory is always restored on exit, including on
#' error.
#'
#' @param p.dir Character scalar. Path to the particle-type folder containing
#'   the raw FCS files.
#' @param asp An \code{AutoSpectralParam} object.
#' @param filename Character scalar. Base filename (without extension) for the
#'   control file to be created.
#' @param legacy.pipeline Logical scalar. If \code{TRUE}, generates a control
#'   file compatible with the legacy gating pipeline (\code{fill.gate.name} and
#'   \code{legacy} both set accordingly). Defaults to \code{FALSE}.
#'
#' @return \code{NULL}, invisibly. Called for its side effect of writing the
#'   control file to \code{p.dir}.
#'
#' @noRd
.create.control.file.in.dir <- function( p.dir, asp, filename, legacy.pipeline = FALSE ) {
  old.wd <- getwd()
  on.exit( setwd( old.wd ), add = TRUE )
  setwd( p.dir )
  create.control.file(
    p.dir, asp,
    fill.gate.name = legacy.pipeline,
    filename       = filename,
    legacy         = legacy.pipeline
  )
  invisible( NULL )
}

#' Fluorophores with an insufficient spectral variant matrix for one particle type
#'
#' @param extraction.result A single element of \code{particle.results}
#'   (a list with \code{spectra} and \code{variants$variants} components,
#'   as built inside \code{run.bead.cell.comparison()}'s stage 1).
#' @param min.rows Integer scalar. Minimum number of rows a fluorophore's
#'   variant matrix must have to count as successfully extracted for this
#'   particle type.
#'
#' @return Character vector of fluorophore names (drawn from
#'   \code{rownames(extraction.result$spectra)}, excluding \code{"AF"})
#'   whose variant matrix is either absent from
#'   \code{extraction.result$variants$variants} or has fewer than
#'   \code{min.rows} rows.
#'
#' @noRd
.fluor.insufficient.variants <- function( extraction.result, min.rows ) {

  all.fluor <- rownames( extraction.result$spectra )
  all.fluor <- all.fluor[ all.fluor != "AF" ]

  variant.list <- extraction.result$variants$variants

  row.n <- vapply( all.fluor, function( fl ) {
    if ( is.null( variant.list[[ fl ]] ) ) return( 0L )
    nrow( variant.list[[ fl ]] )
  }, integer( 1 ) )

  all.fluor[ row.n < min.rows ]
}

# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

#' Compare single-color control spectra across particle types
#'
#' Runs the full AutoSpectral spectral signature extraction pipeline (legacy
#' or automated) independently on each particle-type folder in
#' \code{particle.dirs} -- typically one folder of cell-based single-stained
#' controls and one or more folders of bead-based single-stained controls --
#' then compares each non-reference particle type against the reference type
#' using cosine similarity, absolute spectral distance, and MAD-based
#' variability metrics. Produces per-particle-type diagnostic figures and a
#' consolidated results object.
#'
#' @param particle.dirs Named list of character scalars. Each name is a
#'   particle-type label (e.g. \code{"Cells"}, \code{"UltraComp Beads"}) and
#'   each value is the path to that particle type's folder of single-stained
#'   FCS files plus an unstained/AF reference file. Must contain an entry
#'   named \code{reference.type}.
#' @param asp An \code{AutoSpectralParam} object supplying cytometer,
#'   gating, and pipeline parameters. Per-folder figure and table directories
#'   are overridden internally for each particle type.
#' @param reference.type Character scalar. The name in \code{particle.dirs}
#'   to treat as the reference (typically cell-based) population that all
#'   other particle types are compared against. Defaults to \code{"Cells"}.
#' @param result.dir Character scalar. Root output directory for cached
#'   extractions, comparison figures, and the final results \code{.rds} file.
#'   Defaults to \code{"./results/bead_cell_comparison"}.
#' @param control.filename Character scalar. Base filename (without
#'   extension) used for the auto-generated control file in each particle-type
#'   folder. Defaults to \code{"fcs_control_file"}.
#' @param n.cells.brightness Integer scalar. Number of top (background-
#'   subtracted) events used per fluorophore when computing automated
#'   brightness via \code{get.brightness.automated()}. Defaults to
#'   \code{500L}.
#' @param fluor.df Optional data frame with columns \code{Fluorophore} and
#'   \code{Class}, used to annotate mismatch/variability plots by dye class.
#'   If \code{NULL} (the default), a fallback is built from the package's
#'   \code{fluorophore_database.csv}, using \code{excitation.laser} as the
#'   class.
#' @param figures Logical scalar. Whether to generate diagnostic figures
#'   during spectral extraction. Defaults to \code{TRUE}.
#' @param verbose Logical scalar. Whether to print progress messages.
#'   Defaults to \code{TRUE}.
#' @param legacy.pipeline Logical scalar. If \code{TRUE}, spectral signatures
#'   are extracted with the legacy gated pipeline
#'   (\code{define.flow.control()} / \code{clean.controls()} /
#'   \code{get.fluorophore.spectra()}); if \code{FALSE}, the automated
#'   (non-gated) pipeline (\code{get.spectra.automated()}) is used. Defaults
#'   to \code{FALSE}.
#' @param legacy.gating.system Character scalar. Gating system passed to
#'   \code{define.flow.control()} when \code{legacy.pipeline = TRUE}. Defaults
#'   to \code{"density"}.
#' @param legacy.gate Logical scalar. Whether to perform gating in
#'   \code{define.flow.control()} when \code{legacy.pipeline = TRUE}. Defaults
#'   to \code{TRUE}.
#' @param n.candidates Integer, default `1000`. Number of top-expressing
#'   candidate events selected per fluorophore before cosine-similarity
#'   filtering. Ignored in internal-negative mode, where the top 5%% of
#'   events by peak channel are used directly.
#' @param n.spectral Integer, default `50`. Number of spectral events
#'   retained after filtering for low AF cosine similarity.
#' @param k.neighbors Integer, default `2`. Number of nearest neighbours in
#'   scatter space used for per-event AF subtraction.
#' @param allow.duplicate.controls Logical, default `FALSE`. Set `TRUE` to
#' permit multiple single-stained controls for the same fluorophore
#' (diagnostic/QC use only). Each is tracked internally under a unique
#' `sample` identifier and carried through to `marker.spectra`'s rownames;
#' the true fluorophore identity is preserved as a `"fluorophore"` attribute
#' for reference-library matching and for `check.spectra.duplicates()`.
#' @param mismatch.abs.alpha Numeric scalar in `[0, 1]`. Line transparency used
#'   for the absolute-value spectral mismatch plot described below. Defaults
#'   to \code{0.5}.
#' @param cluster.peak.threshold Numeric scalar in `[0, 1]`. Passed to
#'   \code{assess.mismatch.clusters()} as \code{peak.threshold}: a
#'   fluorophore/detector pair is only used in the cluster test where that
#'   fluorophore has on-peak signal. Defaults to \code{0.1}.
#' @param cluster.min.n Integer scalar. Passed to
#'   \code{assess.mismatch.clusters()} as \code{min.n}: the minimum number of
#'   on-peak fluorophores required for a detector to receive a test
#'   statistic. Defaults to \code{3L}.
#' @param cluster.stat.threshold Numeric scalar. Passed to
#'   \code{assess.mismatch.clusters()} as \code{cluster.threshold}: the
#'   per-detector |t| cutoff used to form candidate clusters. Defaults to
#'   \code{1.96}.
#' @param cluster.n.perm Integer scalar. Passed to
#'   \code{assess.mismatch.clusters()} as \code{n.perm}: the number of
#'   sign-flip permutations used to build the null distribution of the
#'   maximum cluster mass. Defaults to \code{2000L}.
#' @param min.variant.rows Integer scalar. A fluorophore is automatically
#'   excluded from a given particle type's comparison against
#'   \code{reference.type} if its spectral variant matrix has fewer than
#'   this many rows (or is absent entirely) for \code{reference.type} or
#'   for that particle type -- evaluated independently per comparison, so
#'   a fluorophore failing this check for one bead type does not remove it
#'   from comparisons against other bead types. Defaults to \code{3L}
#'   (i.e. more than two rows required).
#' @param exclude.fluorophores Character vector, or \code{NULL} (the
#'   default). Fluorophore names to exclude from every particle type's
#'   downstream comparison/plots, in addition to the automatic exclusion
#'   from \code{min.variant.rows}. Names not present in the data are
#'   ignored.
#'
#' @details
#' Processing proceeds in four stages:
#' \enumerate{
#'   \item \strong{Per-particle-type extraction.} For each entry in
#'     \code{particle.dirs}, a control file is created if absent, validated
#'     for unresolved fluorophore matches and a resolved AF row, and then run
#'     through the selected extraction pipeline to yield spectra, automated
#'     brightness, and spectral variants. Results are cached to
#'     \code{<result.dir>/<particle.type>/<particle.type>_extraction_<pipeline>.rds}
#'     and reused on subsequent calls; switching \code{legacy.pipeline}
#'     between runs uses a distinct cache file rather than reusing a stale
#'     extraction from the other pipeline. Particle types that fail
#'     validation or processing are skipped with a warning and excluded from
#'     downstream comparison; if \code{reference.type} itself fails, the
#'     function stops.
#'   \item \strong{Fluorophore exclusion.} For each non-reference particle
#'     type, any fluorophore whose spectral variant matrix has
#'     \code{< min.variant.rows} rows (or is missing) for \strong{that}
#'     particle type or for \code{reference.type} is dropped from
#'     \strong{that particle type's} comparison only, along with any
#'     fluorophore named in \code{exclude.fluorophores}. Exclusion is
#'     evaluated independently per comparison, so a bead type that fails
#'     QC outright (e.g. from low acquisition counts) does not remove
#'     fluorophores from other, unrelated comparisons. A comparison that
#'     fails outright itself (e.g. no fluorophore survives exclusion for
#'     that pair) is skipped with a warning rather than aborting every
#'     other comparison.
#'   \item \strong{Mismatch assessment.} Each non-reference particle type's
#'     spectral variants are first renormalized from the package's native
#'     L-infinity convention (peak detector fixed to 1.0) to unit L2
#'     (Euclidean) norm via \code{l2.normalize.spectra()}, so that
#'     distance-based comparisons are not artificially biased toward zero
#'     variance at each spectrum's peak detector. The renormalized variants
#'     are then compared against the reference type's renormalized variants
#'     via \code{assess.mismatch()} (cosine similarity, unaffected by the
#'     choice of per-row normalization), \code{assess.variability()},
#'     \code{assess.variability.mad()}, and \code{bead.cell.dist()}
#'     (per-detector signed distance). \code{Extraction} in the return
#'     value retains the original L-infinity-normalized spectra/variants
#'     as produced by the extraction pipeline; only the comparison stages
#'     below use the L2-renormalized versions.
#'   \item \strong{Error-metric plots.} \code{plot.mismatch()} is called per
#'     non-reference particle type, with shared axis limits computed across
#'     all particle types for comparability.
#'   \item \strong{Spectral location of mismatch.} Per-detector MAD of the
#'     distance matrix is computed and visualized with
#'     \code{spectral.variant.plot.dens()} two ways: a signed view (raw
#'     \code{reference - other} differences, showing whether mismatch is
#'     systematically one-signed at a given detector) and a magnitude-only
#'     view (absolute differences, showing where mismatch concentrates in the
#'     spectrum irrespective of sign). Both use the same per-detector MAD of
#'     the signed distances as the reference line.
#'   \item \strong{Cluster-based mismatch testing.} For each non-reference
#'     particle type, \code{assess.mismatch.clusters()} tests whether
#'     mismatch is concentrated in a specific, contiguous region of the
#'     spectrum rather than spread evenly across detectors, using a
#'     sign-flip permutation test restricted to on-peak fluorophore/detector
#'     pairs. Detected clusters (with permutation p-values) are written to
#'     \code{<result.dir>/table_mismatch_clusters/}.
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{Extraction}{Named list (by particle type) of extraction results,
#'     each containing \code{spectra}, \code{brightness}, and \code{variants}.}
#'   \item{Comparison}{Named list (by non-reference particle type) of
#'     \code{Cosine}, \code{Variability}, \code{VariabilityMAD}, and
#'     \code{Distance} results versus \code{reference.type}.}
#'   \item{Stats}{Named list (by non-reference particle type) of
#'     \code{mismatch.plot()} outputs.}
#'   \item{MAD}{Named list (by non-reference particle type) of per-detector
#'     MAD spectra of the mismatch distance matrix.}
#'   \item{Clusters}{Named list (by non-reference particle type) of the
#'     \code{clusters} data frame from \code{assess.mismatch.clusters()}.}
#' }
#' The same list is also written to
#' \code{<result.dir>/bead_cell_comparison_results.rds}.
#'
#' @export
run.bead.cell.comparison <- function(
    particle.dirs,
    asp,
    reference.type      = "Cells",
    result.dir          = "./results/bead_cell_comparison",
    control.filename    = "fcs_control_file",
    n.cells.brightness  = 500L,
    fluor.df            = NULL,
    figures             = TRUE,
    verbose             = TRUE,
    legacy.pipeline      = FALSE,
    legacy.gating.system = "density",
    legacy.gate          = TRUE,
    n.candidates            = 1000L,
    n.spectral              = 200L,
    k.neighbors             = 2L,
    allow.duplicate.controls = TRUE,
    mismatch.abs.alpha       = 0.5,
    cluster.peak.threshold  = 0.1,
    cluster.min.n           = 3L,
    cluster.stat.threshold  = 1.96,
    cluster.n.perm          = 2000L,
    min.variant.rows        = 3L,
    exclude.fluorophores    = NULL
) {

  if ( !reference.type %in% names( particle.dirs ) )
    stop(
      "particle.dirs must include an entry named '", reference.type, "'.",
      call. = FALSE
    )

  if ( !dir.exists( result.dir ) ) dir.create( result.dir, recursive = TRUE )

  particle.types <- names( particle.dirs )
  bead.types     <- setdiff( particle.types, reference.type )

  # default dye-class annotation if none supplied --
  if ( is.null( fluor.df ) ) {
    fluor.db.path <- system.file(
      "extdata", "fluorophore_database.csv", package = "AutoSpectral"
    )
    fluor.db <- utils::read.csv( fluor.db.path, stringsAsFactors = FALSE )
    fluor.df <- data.frame(
      Fluorophore = fluor.db$fluorophore,
      Class       = fluor.db$class,
      stringsAsFactors = FALSE
    )
  }

  # -- 1. per-particle-type extraction ---------------------------------------
  particle.results <- lapply( particle.types, function( p ) {

    p.dir <- particle.dirs[[ p ]]
    p.out <- file.path( result.dir, p )
    if ( !dir.exists( p.out ) ) dir.create( p.out, recursive = TRUE )

    # cache key includes the pipeline choice -- switching legacy.pipeline
    # between runs must not silently reuse a stale extraction from the
    # other pipeline
    result.file <- file.path(
      p.out,
      paste0( p, "_extraction_", if ( legacy.pipeline ) "legacy" else "automated", ".rds" )
    )
    if ( file.exists( result.file ) ) {
      if ( verbose ) message( "Reading cached extraction for ", p )
      return( readRDS( result.file ) )
    }

    if ( verbose )
      message( "\033[36m== Processing particle type: ", p, " ==\033[0m" )

    # -- 1a. control file
    control.path <- file.path( p.dir, paste0( control.filename, ".csv" ) )

    if ( !file.exists( control.path ) ) {
      .create.control.file.in.dir( p.dir, asp, control.filename, legacy.pipeline )
      finalize.control.file( control.path, p, reference.type )
    } else if ( verbose ) {
      message( "Using existing control file for ", p )
    }

    # -- 1b. pre-flight validation (cheap check before the expensive steps)
    ctrl.check <- utils::read.csv( control.path, stringsAsFactors = FALSE )

    if ( any( grepl( "No match", ctrl.check$fluorophore ) ) ) {
      warning(
        "Skipping '", p, "': unresolved 'No match' fluorophore(s) in ",
        control.path, ". Edit the file and rerun.", call. = FALSE
      )
      return( NULL )
    }

    if ( !"AF" %in% ctrl.check$fluorophore ) {
      warning(
        "Skipping '", p, "': no AF/unstained control resolved in ",
        control.path, ". Edit the file and rerun.", call. = FALSE
      )
      return( NULL )
    }

    # -- 1c. per-folder figure/table directories (covers both pipelines)
    asp$figure.gate.dir               <- file.path( p.out, "figure_gate" )
    asp$figure.clean.control.dir      <- file.path( p.out, "figure_clean_controls" )
    asp$figure.spectral.ribbon.dir    <- file.path( p.out, "figure_spectral_ribbon" )
    asp$figure.spectra.dir            <- file.path( p.out, "figure_spectra" )
    asp$figure.similarity.heatmap.dir <- file.path( p.out, "figure_similarity_heatmap" )
    asp$figure.scatter.dir.base       <- file.path( p.out, "figure_scatter" )
    asp$table.spectra.dir             <- file.path( p.out, "table_spectra" )
    asp$variant.dir                   <- file.path( p.out, "figure_spectral_variants" )

    # -- 1d. spectral signature extraction (automated or legacy)
    if ( legacy.pipeline ) {

      legacy.ok <- tryCatch(
        {
          check.control.file(
            p.dir,
            control.path,
            asp, strict = TRUE,
            legacy = TRUE,
            allow.duplicate.controls = allow.duplicate.controls
          )
          TRUE
        },
        error = function( e ) {
          warning(
            "Skipping '", p, "': legacy control file validation failed -- ",
            conditionMessage( e ), call. = FALSE
          )
          FALSE
        }
      )

      if ( !isTRUE( legacy.ok ) ) return( NULL )

      flow.control <- define.flow.control(
        control.dir      = p.dir,
        control.def.file = control.path,
        asp              = asp,
        gate             = legacy.gate,
        gating.system    = legacy.gating.system,
        verbose          = verbose,
        allow.duplicate.controls = allow.duplicate.controls
      )

      flow.control <- clean.controls(
        flow.control = flow.control,
        asp          = asp,
        main.figures = figures,
        verbose      = verbose
      )

      spectra <- get.fluorophore.spectra(
        flow.control   = flow.control,
        asp            = asp,
        use.clean.expr = TRUE,
        title          = p,
        figures        = figures
      )

    } else {

      spectra <- get.spectra.automated(
        control.dir      = p.dir,
        control.def.file = control.path,
        asp              = asp,
        return.af        = FALSE,
        figures          = figures,
        verbose          = verbose,
        n.candidates     = n.candidates,
        n.spectral       = n.spectral,
        k.neighbors      = k.neighbors,
        allow.duplicate.controls = allow.duplicate.controls
      )
    }

    # -- 1e. brightness (top n.cells.brightness events, background-subtracted)
    # -- shared by both pipelines: it only needs the final `spectra` matrix
    # -- plus the control file/raw FCS files, not how spectra were derived
    brightness <- get.brightness.automated(
      control.dir      = p.dir,
      control.def.file = control.path,
      asp              = asp,
      spectra          = spectra,
      n.cells          = n.cells.brightness,
      verbose          = verbose
    )

    # -- 1f. spectral variants -- likewise shared by both pipelines
    # use.unmixed = FALSE: the whole point of this comparison is several
    # similar fluorophores measured across particle types, so a full-spectra
    # unmix of `spectra` against itself is unstable or unsolvable here. AF
    # extraction and variant clustering are restricted to raw detector space,
    # and the (unmixed-space) Spillover Spreading Matrix is skipped.
    variants <- get.spectral.variants(
      control.dir      = p.dir,
      control.def.file = control.path,
      asp              = asp,
      spectra          = spectra,
      figures          = figures,
      output.dir       = asp$variant.dir,
      verbose          = verbose,
      use.unmixed      = FALSE
    )

    result <- list( spectra = spectra, brightness = brightness, variants = variants )
    saveRDS( result, result.file )
    result
  } )
  names( particle.results ) <- particle.types

  ok <- !vapply( particle.results, is.null, logical( 1 ) )

  if ( !ok[ reference.type ] )
    stop(
      "Reference type '", reference.type,
      "' failed processing -- cannot continue.", call. = FALSE
    )

  particle.results <- particle.results[ ok ]
  bead.types        <- intersect( bead.types, names( particle.results ) )

  # -- 2. fluorophore exclusion (automatic + manual) -------------------------
  # Evaluated independently per comparison (reference vs one bead type), not
  # globally across every particle type in particle.dirs. A fluorophore
  # failing QC for a bead type unrelated to a given comparison must not pull
  # that fluorophore out of other, unrelated comparisons -- e.g. a third
  # bead type failing outright on low acquisition counts should not gut the
  # reference.type-vs-UltraComp comparison. Each pair still needs a
  # fluorophore's variant matrix to have >= min.variant.rows rows on BOTH
  # sides of that specific pair, since assess.variability()/
  # assess.variability.mad() slice rows 2:nrow(...) and need at least one
  # variant row beyond the reference spectrum.
  pair.exclude <- lapply( bead.types, function( b ) {
    union(
      union(
        .fluor.insufficient.variants( particle.results[[ reference.type ]], min.rows = min.variant.rows ),
        .fluor.insufficient.variants( particle.results[[ b ]],              min.rows = min.variant.rows )
      ),
      exclude.fluorophores
    )
  } )
  names( pair.exclude ) <- bead.types

  if ( verbose ) {
    for ( b in bead.types ) {
      auto.b <- setdiff( pair.exclude[[ b ]], exclude.fluorophores )
      if ( length( auto.b ) > 0 )
        message(
          "\033[33mExcluding fluorophore(s) with < ", min.variant.rows,
          " spectral variant row(s) in ", reference.type, " or ", b, ": ",
          paste( auto.b, collapse = ", " ), "\033[0m"
        )
    }
    if ( length( exclude.fluorophores ) > 0 )
      message(
        "\033[33mManually excluding fluorophore(s) from every comparison: ",
        paste( exclude.fluorophores, collapse = ", " ), "\033[0m"
      )
  }

  # -- 3. mismatch assessment vs reference -----------------------------------
  # renormalize from the package's native L-infinity convention (each
  # spectrum's peak detector fixed to exactly 1.0) to unit L2 (Euclidean)
  # norm before any distance-based comparison. L-infinity normalization
  # implicitly asserts zero measurement variance at whichever detector is
  # the peak, which pushes all of a spectrum's real cell-vs-bead disagreement
  # into the off-peak channels and biases bead.cell.dist()/the MAD traces
  # below toward reporting mismatch away from the peak even when none is
  # present. assess.mismatch()'s cosine similarity is unaffected either way.
  # cell.variants/ref.spectra are recomputed per bead type below since the
  # excluded fluorophore set (pair.exclude[[b]]) differs per comparison. A
  # comparison that fails outright (e.g. no fluorophore survives exclusion
  # for that pair) is caught and skipped so it cannot abort every other
  # comparison in this lapply.
  comparison <- lapply( bead.types, function( b ) {

    tryCatch(
      {
        if ( verbose )
          message( "\033[36m== Comparing ", b, " vs ", reference.type, " ==\033[0m" )

        fluor.exclude.b <- pair.exclude[[ b ]]

        cell.variants <- particle.results[[ reference.type ]]$variants$variants
        cell.variants <- cell.variants[ setdiff( names( cell.variants ), fluor.exclude.b ) ]
        cell.variants <- lapply( cell.variants, l2.normalize.spectra )

        # reference spectra for the on-peak mask in assess.mismatch.clusters() --
        # AF is dropped since it never appears in cell.variants/bead.variants and
        # would otherwise just trigger assess.mismatch.clusters()'s shared-rows
        # warning on every call
        ref.spectra <- particle.results[[ reference.type ]]$spectra
        ref.spectra <- ref.spectra[
          setdiff( rownames( ref.spectra ), c( "AF", fluor.exclude.b ) ), , drop = FALSE
        ]
        ref.spectra <- l2.normalize.spectra( ref.spectra )

        bead.variants <- particle.results[[ b ]]$variants$variants
        bead.variants <- bead.variants[ setdiff( names( bead.variants ), fluor.exclude.b ) ]
        bead.variants <- lapply( bead.variants, l2.normalize.spectra )
        dist.mat      <- bead.cell.dist( cell.variants, bead.variants )

        list(
          Cosine         = assess.mismatch( cell.variants, bead.variants ),
          Variability    = assess.variability( bead.variants ),
          VariabilityMAD = assess.variability.mad( bead.variants ),
          Distance       = dist.mat,
          Clusters       = assess.mismatch.clusters(
            dist.mat          = dist.mat,
            ref.spectra       = ref.spectra,
            peak.threshold    = cluster.peak.threshold,
            min.n             = cluster.min.n,
            cluster.threshold = cluster.stat.threshold,
            n.perm            = cluster.n.perm,
            bird.seed         = asp$bird.seed
          )
        )
      },
      error = function( e ) {
        warning(
          "Skipping comparison '", b, "' vs '", reference.type, "': ",
          conditionMessage( e ), call. = FALSE
        )
        NULL
      }
    )
  } )
  names( comparison ) <- bead.types

  comparison.ok <- !vapply( comparison, is.null, logical( 1 ) )

  if ( !any( comparison.ok ) )
    stop(
      "Every particle type comparison against '", reference.type,
      "' failed -- nothing to report.", call. = FALSE
    )

  comparison <- comparison[ comparison.ok ]
  bead.types <- bead.types[ comparison.ok ]

  # -- 4. error-metric plots per particle type -------------------------------
  plot.dir <- file.path( result.dir, "figures" )
  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )

  all.mismatch <- unlist( lapply(
    comparison, function( x ) rowSums( abs( x$Distance ) )
  ) )
  # the previous fixed floor of 2 was calibrated to the old L-infinity
  # normalization (peak detector == 1.0); under L2 normalization the same
  # underlying mismatch is naturally smaller in absolute terms, so the axis
  # is left fully data-driven with a small headroom margin instead
  mismatch.limits <- c( 0, max( all.mismatch, na.rm = TRUE ) * 1.05 )

  all.cosine <- unlist( lapply( comparison, function( x ) x$Cosine ) )
  sim.limits <- c( 1, min( 0.95, all.cosine, na.rm = TRUE ) )

  stats.results <- lapply( bead.types, function( b ) {
    mismatch.plot(
      cosine.data      = comparison[[ b ]]$Cosine,
      mismatch.data    = rowSums( abs( comparison[[ b ]]$Distance ) ),
      variability.data = rowSums( abs( comparison[[ b ]]$VariabilityMAD ) ),
      brightness.data  = particle.results[[ b ]]$brightness,
      fluor.df         = fluor.df,
      particle.name    = b,
      output.dir       = plot.dir,
      mismatch.limits  = mismatch.limits,
      sim.limits       = sim.limits,
      cytometer        = asp$cytometer
    )
  } )
  names( stats.results ) <- bead.types

  # -- 5. spectral location of mismatch (MAD density plots) -----------------
  mad.dir <- file.path( result.dir, "figure_mismatch_mad" )
  if ( !dir.exists( mad.dir ) ) dir.create( mad.dir, recursive = TRUE )

  abs.dir <- file.path( result.dir, "figure_mismatch_abs" )
  if ( !dir.exists( abs.dir ) ) dir.create( abs.dir, recursive = TRUE )

  mad.summary <- lapply( bead.types, function( b ) {

    dist.mat     <- comparison[[ b ]]$Distance
    mad.spectrum <- apply( dist.mat, 2, stats::mad )
    med.spectrum <- apply( abs( dist.mat ), 2, stats::median )

    spectral.variant.plot.dens(
      spectra.variants  = dist.mat,
      median.spectrum   = mad.spectrum,
      title             = paste0(
        asp$cytometer, "_", b, "_vs_", reference.type, "_mismatch_MAD"
      ),
      save              = TRUE,
      plot.dir          = mad.dir,
      variant.color     = "red",
      variant.alpha     = 0.15,
      median.line.color = "black"
    )

    # magnitude-only view: same per-fluorophore per-detector distances and
    # the same (signed) MAD reference line, but individual traces are shown
    # as |distance| so mismatch concentration reads directly off the y-axis
    # instead of cancelling visually across fluorophores of opposite sign
    spectral.variant.plot.dens(
      spectra.variants  = abs( dist.mat ),
      median.spectrum   = med.spectrum,
      title             = paste0(
        asp$cytometer, "_", b, "_vs_", reference.type, "_mismatch_abs"
      ),
      save              = TRUE,
      plot.dir          = abs.dir,
      variant.color     = "red",
      variant.alpha     = mismatch.abs.alpha,
      median.line.color = "black"
    )

    mad.spectrum
  } )
  names( mad.summary ) <- bead.types

  # -- 6. cluster-based mismatch testing (contiguous-region significance) ---
  cluster.dir <- file.path( result.dir, "table_mismatch_clusters" )
  if ( !dir.exists( cluster.dir ) ) dir.create( cluster.dir, recursive = TRUE )

  cluster.summary <- lapply( bead.types, function( b ) {

    clusters.df <- comparison[[ b ]]$Clusters$clusters

    utils::write.csv(
      clusters.df,
      file.path(
        cluster.dir,
        paste0( asp$cytometer, "_", b, "_vs_", reference.type, "_clusters.csv" )
      ),
      row.names = FALSE
    )

    if ( verbose ) {
      sig <- clusters.df[
        !is.na( clusters.df$p.value ) & clusters.df$p.value < 0.05, , drop = FALSE
      ]
      if ( nrow( sig ) > 0 ) {
        message(
          "\033[33m", b, " vs ", reference.type, ": ", nrow( sig ),
          " significant mismatch cluster(s) (p < 0.05):\033[0m"
        )
        for ( i in seq_len( nrow( sig ) ) )
          message(
            "  ", sig$start[ i ], "-", sig$end[ i ],
            " (", sig$n.detectors[ i ], " detectors, p = ",
            signif( sig$p.value[ i ], 3 ), ")"
          )
      } else {
        message(
          "\033[32m", b, " vs ", reference.type,
          ": no significant mismatch clusters.\033[0m"
        )
      }
    }

    clusters.df
  } )
  names( cluster.summary ) <- bead.types

  results <- list(
    Extraction = particle.results,
    Comparison = comparison,
    Stats      = stats.results,
    MAD        = mad.summary,
    Clusters   = cluster.summary
  )

  saveRDS( results, file.path( result.dir, "bead_cell_comparison_results.rds" ) )

  results
}
