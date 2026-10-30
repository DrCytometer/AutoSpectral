# estimate_noise_model.R

#' @title Estimate Detector Noise Model
#'
#' @description
#' Fits the per-detector noise model
#' \deqn{Var(y_d) = \sigma_{0,d}^2 + \mu_d / \kappa_d}
#' used by `unmix.gls()` and `unmix.af.gls()` from the single-stained
#' controls, and checks it on held-out events.
#'
#' Run it after `get.spectral.variants()`, which supplies the unmixed
#' positivity thresholds. The steps:
#'
#' 1. **Event selection.** Every file is read in full and its scatter core
#'    and singlets kept (`selection = "auto"`): a robust ellipse around the
#'    densest FSC-A x SSC-A mode, then an FSC-A/FSC-H singlet band. Doublets
#'    and off-core particles add variance that grows with signal, which the
#'    fit would otherwise read as shot noise.
#' 2. **Unstained cells** (the `AF` row of the control file, or
#'    `unstained.sample`) are split into training and held-out halves; their
#'    design is the leading right singular vectors of the training events.
#' 3. **Bead negatives** anchor sigma0 on beads: unstained bead samples
#'    (`Negative` rows of type `beads`), the bead controls' universal
#'    negatives, else unstained beads found inside the bead controls.
#' 4. **Single-stain controls** are split 70/30. Positives are unsaturated
#'    events above the 99.5th percentile of the control's negative in its
#'    peak detector. Spectral variants are rediscovered on the training split
#'    only (`get.fluor.variants()`), so the held-out events stay independent
#'    of the design. A control's design is its reference spectrum, its mean
#'    variant delta, its variant basis (`build.variant.basis()`) and, for
#'    cell controls, the AF rows of its negative; a bead control carries its
#'    negative's mean as an offset.
#' 5. **Joint fit.** Each group's events are projected onto its design and
#'    binned on its own abundance; the binned residual variances of every
#'    group and detector are fitted together by non-negative least squares
#'    through the mixing matrix \eqn{(I - H) \circ (I - H)}, which accounts
#'    for each detector's residual carrying the variance of the others.
#'    sigma0 is clamped above at the lowest variance a bin shows and below at
#'    `floor.frac` of its median; a detector without shot-noise information
#'    takes the median kappa. Each is flagged.
#' 6. **Calibration** on the held-out events: observed / predicted residual
#'    variance per detector and intensity tranche (1 when the model is
#'    right; the share inside 0.8-1.25 is reported, 80% or more is good),
#'    and the share of events per source whose GLS chi-square has p < 0.001.
#' 7. **Clamped data.** When no detector has 1% or more of its negative
#'    values at or below zero, the data were clamped at zero by the export
#'    and their noise is truncated there. The model is still fitted;
#'    `clamped.weights` then holds the weights recommended in its place:
#'    Poisson weights with a floor of `clamp.floor.pe` photoelectrons
#'    (`"poisson"`) or flat weights (`"flat"`), scaled to the held-out
#'    unstained cells.
#'
#' Estimator options, all off by default and intended for validation:
#' `bead.af.rows` (bead negatives and bead controls take AF rows of the bead
#' negative instead of its mean), `laser.split` (one design row per laser),
#' `per.control.quadratic` (one `c_g mu^2` per control, absorbing
#' control-specific spread; used in calibration, not for unmixing),
#' `sigma0.by.class` (one sigma0 for cells and one for beads, kappa shared)
#' and `sigma0.fallback = "laser"` (a clamped detector takes the median
#' sigma0 of its laser's unflagged detectors instead of the clamp).
#'
#' Simplifications: one unstained cell sample is used; bead negatives are
#' chosen automatically; lasers come from the cytometer database only, and
#' the laser options are skipped when any detector's laser is unknown; bead
#' controls without a universal-negative file rediscover their variants
#' against their internal negative events.
#'
#' @param control.dir Character, folder of the single-stained control FCS
#'   files.
#' @param control.def.file Character, the control definition CSV.
#' @param asp The AutoSpectral parameter list.
#' @param spectra Numeric matrix, reference spectra (fluorophores x
#'   detectors). An `AF` row is ignored.
#' @param spectral.variants List returned by `get.spectral.variants()`; its
#'   `thresholds` (unmixed positivity thresholds of the unstained cells) are
#'   used for variant discovery on the cell controls.
#' @param unstained.sample Optional path to the unstained cell FCS file, used
#'   when the control file has no `AF` row. Default `NULL`.
#' @param output.dir Folder for `noise_model.rds`. Default `asp$variant.dir`.
#' @param save Logical, save the result as `noise_model.rds`. Default `TRUE`.
#' @param selection `"auto"` (default): scatter core and singlets of every
#'   file. `"none"`: every event.
#' @param unstained.train.frac,control.train.frac Numeric, training share of
#'   the unstained cells (default `0.5`) and of each control and bead
#'   negative (default `0.7`).
#' @param max.unstained,max.test.unstained Integers, unstained events read
#'   (default `200000`) and held out at most (default `20000`).
#' @param max.control.events,max.test.control Integers, control events read
#'   (default `100000`) and held-out positives kept (default `5000`).
#' @param max.fit.events,max.cal.events Integers, events per control and bead
#'   negative in the fit (default `6000`; four times that for the unstained
#'   cells) and in the calibration (default `2000`; twice for the unstained).
#' @param n.af.rows Integer, AF design rows of a cell negative. Default `4`.
#' @param basis.var.explained,basis.max.rank Variant basis: variance
#'   explained (default `0.95`) and maximum rank (default `5`).
#' @param bead.factor Numeric, minimum ratio between the positives and an
#'   internal bead negative. Default `10`.
#' @param n.bins,min.bin.events,fit.quadratic,n.irls,n.restandardise,floor.frac,bound.sds
#'   Estimator settings: quantile bins per group (default `20`), fewest
#'   events per bin (`50`), a per-detector `c mu^2` term (`FALSE`), IRLS
#'   passes (`2`), restandardisation rounds (`2`), lower clamp on sigma0^2 as
#'   a fraction of its median (`0.01`) and sampling SDs added to its upper
#'   clamp (`1`).
#' @param bead.af.rows Integer, AF rows of a bead negative (`0` = off; `3`
#'   when on).
#' @param laser.split Logical, one design row per laser. Default `FALSE`.
#' @param laser.split.min Numeric, the part on a laser and the rest must
#'   each reach this fraction of the spectrum's peak. Default `0.01`.
#' @param per.control.quadratic,sigma0.by.class Logical options. Default
#'   `FALSE`.
#' @param sigma0.fallback `"floor"` (default) or `"laser"`.
#' @param clamp.mode `"auto"` (default) detects clamped data; `"off"` skips
#'   the check.
#' @param clamp.weights `"poisson"` (default) or `"flat"`.
#' @param clamp.floor.pe Numeric, floor of the Poisson weights in
#'   photoelectrons. Default `10`.
#' @param n.cells,som.dim,k.neighbors,sim.threshold,sim.threshold.floor,af.collinear.threshold
#'   Variant discovery settings passed to `get.fluor.variants()`. Defaults
#'   `10000`, `5`, `3`, `0.985`, `0.90`, `0.95`.
#' @param seed Integer seed of the splits and subsamples. Default
#'   `asp$bird.seed`.
#' @param parallel,threads Threading of the variant SOMs. Defaults `TRUE`,
#'   `NULL`.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A list. For `unmix.gls()` and `unmix.af.gls()`:
#' \describe{
#'   \item{`read.var`}{Named numeric vector, sigma0^2 per detector (the
#'     cells' sigma0 when fitted per class).}
#'   \item{`counts.per.unit`}{Named numeric vector, kappa per detector.}
#'   \item{`kappa.pooled`}{Median kappa.}
#' }
#' and in addition: `sigma0.sq` and `g` (= 1 / kappa), `c`, `flags` (bit
#' mask per detector: 1 sigma0 clamped high, 2 raised to the floor, 4 kappa
#' pooled, 8 sigma0 from the same laser) and `flag.text`, `classes`,
#' `sigma0.sq.by.class`, `flags.by.class` and `by.class` (a complete
#' `read.var` / `counts.per.unit` list per class, with `sigma0.by.class`),
#' `group.c`, `fit.table` (per group, bin and detector: observed and
#' predicted residual variance), `groups`, `calibration` (from held-out
#' events: `in.band`, `in.band.controls`, the tranche tables and one row per
#' source), `clamped`, `negative.share`, `flat.phi`, `clamped.weights`
#' (`kind`, `floor`, `scale` and a `noise.model` list usable in place of
#' this one), `selection` (event counts per file), `sources`, `unstained`,
#' `notes` and `settings`.
#'
#' @export

estimate.noise.model <- function(
    control.dir,
    control.def.file,
    asp,
    spectra,
    spectral.variants,
    unstained.sample      = NULL,
    output.dir            = NULL,
    save                  = TRUE,
    selection             = c( "auto", "none" ),
    unstained.train.frac  = 0.5,
    control.train.frac    = 0.7,
    max.unstained         = 200000L,
    max.test.unstained    = 20000L,
    max.control.events    = 100000L,
    max.test.control      = 5000L,
    max.fit.events        = 6000L,
    max.cal.events        = 2000L,
    n.af.rows             = 4L,
    basis.var.explained   = 0.95,
    basis.max.rank        = 5L,
    bead.factor           = 10,
    n.bins                = 20L,
    min.bin.events        = 50L,
    fit.quadratic         = FALSE,
    n.irls                = 2L,
    n.restandardise       = 2L,
    floor.frac            = 0.01,
    bound.sds             = 1,
    bead.af.rows          = 0L,
    laser.split           = FALSE,
    laser.split.min       = 0.01,
    per.control.quadratic = FALSE,
    sigma0.by.class       = FALSE,
    sigma0.fallback       = c( "floor", "laser" ),
    clamp.mode            = c( "auto", "off" ),
    clamp.weights         = c( "poisson", "flat" ),
    clamp.floor.pe        = 10,
    n.cells               = 10000L,
    som.dim               = 5L,
    k.neighbors           = 3L,
    sim.threshold         = 0.985,
    sim.threshold.floor   = 0.90,
    af.collinear.threshold = 0.95,
    seed                  = asp$bird.seed,
    parallel              = TRUE,
    threads               = NULL,
    verbose               = TRUE
) {

  selection       <- match.arg( selection )
  sigma0.fallback <- match.arg( sigma0.fallback )
  clamp.mode      <- match.arg( clamp.mode )
  clamp.weights   <- match.arg( clamp.weights )
  seed <- if ( is.null( seed ) ) 42L else as.integer( seed )

  settings <- .noise.settings(
    n.bins = n.bins, min.bin.events = min.bin.events, fit.quadratic = fit.quadratic,
    n.irls = n.irls, n.restandardise = n.restandardise, floor.frac = floor.frac,
    bound.sds = bound.sds, per.control.quadratic = per.control.quadratic,
    sigma0.by.class = sigma0.by.class, sigma0.fallback = sigma0.fallback )

  if ( is.null( spectral.variants$thresholds ) )
    stop( "`spectral.variants` must be the list returned by get.spectral.variants().",
          call. = FALSE )

  if ( is.null( output.dir ) ) output.dir <- asp$variant.dir
  if ( save && !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

  spectra <- as.matrix( spectra )
  fluorophores <- rownames( spectra )[ rownames( spectra ) != "AF" ]
  spectra   <- spectra[ fluorophores, , drop = FALSE ]
  detectors <- colnames( spectra )
  ceiling   <- asp$expr.data.max
  notes     <- character( 0 )
  note <- function( txt ) {
    notes <<- c( notes, txt )
    if ( verbose ) message( "  ", txt )
  }

  # ---------------------------------------------------------------------------
  # Control file
  # ---------------------------------------------------------------------------

  if ( !file.exists( control.def.file ) )
    stop( paste( "Unable to locate control.def.file:", control.def.file ), call. = FALSE )

  control.table <- utils::read.csv( control.def.file, stringsAsFactors = FALSE, strip.white = TRUE )
  control.table[] <- lapply( control.table, function( x ) {
    if ( is.character( x ) ) { x <- trimws( x ); x[ x == "" ] <- NA; x } else x
  } )
  control.table <- control.table[ !is.na( control.table$fluorophore ), , drop = FALSE ]

  file.of  <- stats::setNames( control.table$filename, control.table$fluorophore )
  channel  <- stats::setNames( control.table$channel, control.table$fluorophore )
  type.of  <- stats::setNames( control.table$control.type, control.table$fluorophore )
  univ.neg <- stats::setNames( control.table$universal.negative, control.table$fluorophore )
  univ.neg[ is.na( univ.neg ) ] <- "FALSE"
  is.fcs <- function( x ) !is.na( x ) && x != "FALSE" && grepl( "\\.fcs$", x, ignore.case = TRUE )

  unstained.file <- if ( "AF" %in% names( file.of ) ) file.of[[ "AF" ]] else unstained.sample
  if ( is.null( unstained.file ) )
    stop( "An unstained cell control is required: include an `AF` row in the control ",
          "file, or supply `unstained.sample`.", call. = FALSE )

  controls <- names( file.of )[ !grepl( "Negative|^AF$", names( file.of ) ) ]
  missing.fl <- setdiff( controls, fluorophores )
  if ( length( missing.fl ) )
    note( paste0( "Not in `spectra`, skipped: ", paste( missing.fl, collapse = ", " ), "." ) )
  controls <- intersect( controls, fluorophores )
  bead.negative.files <- unique( stats::na.omit( control.table$filename[
    control.table$fluorophore == "Negative" & control.table$control.type == "beads" ] ) )

  # ---------------------------------------------------------------------------
  # Loading with the automatic selection
  # ---------------------------------------------------------------------------

  scatter.channel <- read.scatter.parameter( asp )
  counts <- list()
  cache  <- new.env()

  resolve <- function( f ) {
    p <- file.path( control.dir, f )
    if ( !file.exists( p ) && file.exists( f ) ) f else p
  }

  load.events <- function( f, keep.in.cache = FALSE ) {
    if ( exists( f, envir = cache, inherits = FALSE ) ) return( get( f, envir = cache ) )
    path <- resolve( f )
    cols <- colnames( readFCS( path, start.row = 1, end.row = 1 ) )
    if ( !all( detectors %in% cols ) )
      stop( "`", f, "` does not contain every detector of `spectra`.", call. = FALSE )
    fsc.h <- if ( grepl( "-A$", scatter.channel[ 1 ] ) ) sub( "-A$", "-H", scatter.channel[ 1 ] ) else NA
    sc.cols <- intersect( scatter.channel, cols )
    dat <- readFCS( path, columns = unique( c( detectors, sc.cols, intersect( fsc.h, cols ) ) ) )
    mask <- rep( TRUE, nrow( dat ) )
    if ( selection == "auto" && all( scatter.channel[ 1:2 ] %in% cols ) ) {
      sel <- .select.noise.events(
        dat[ , scatter.channel[ 1 ] ], dat[ , scatter.channel[ 2 ] ],
        if ( !is.na( fsc.h ) && fsc.h %in% cols ) dat[ , fsc.h ] else NULL )
      mask <- sel$mask
      counts[[ f ]] <<- sel$counts
    } else {
      counts[[ f ]] <<- c( n.root = nrow( dat ), n.scatter = nrow( dat ), n.singlet = nrow( dat ) )
    }
    out <- list( fl = dat[ mask, detectors, drop = FALSE ], sc = dat[ mask, sc.cols, drop = FALSE ] )
    if ( keep.in.cache ) assign( f, out, envir = cache )
    out
  }

  if ( selection == "auto" && length( scatter.channel ) < 2L )
    note( "No FSC-A and SSC-A: every event of each file is used." )

  # ---------------------------------------------------------------------------
  # Unstained cells
  # ---------------------------------------------------------------------------

  if ( verbose ) message( "\033[34mModelling detector noise\033[0m" )
  if ( verbose ) message( "  Unstained cells: ", unstained.file )

  u <- load.events( unstained.file )
  if ( nrow( u$fl ) < 500L )
    stop( "Only ", nrow( u$fl ), " unstained events after selection.", call. = FALSE )
  keep  <- .subsample.idx( nrow( u$fl ), max.unstained, seed )
  u.fl  <- u$fl[ keep, , drop = FALSE ]
  u.sc  <- u$sc[ keep, , drop = FALSE ]
  rm( u )
  split <- .holdout.split( nrow( u.fl ), unstained.train.frac, seed, max.test.unstained )
  u.train <- u.fl[ split$train, , drop = FALSE ]
  u.test  <- u.fl[ split$test, , drop = FALSE ]
  setup.basis <- .af.basis.rows( u.train, n.af.rows, seed = seed )
  bases.af  <- list()
  bases.af[[ unstained.file ]] <- setup.basis
  negatives <- list()
  negatives[[ unstained.file ]] <- list( fl = u.train, sc = u.sc[ split$train, , drop = FALSE ] )

  # ---------------------------------------------------------------------------
  # Bead negatives
  # ---------------------------------------------------------------------------

  bead.rows <- function( neg.fl ) {
    if ( bead.af.rows <= 0L || is.null( neg.fl ) || nrow( neg.fl ) < 10L ) return( NULL )
    tr <- .holdout.split( nrow( neg.fl ), control.train.frac, seed + 7L )$train
    .af.basis.rows( neg.fl[ tr, , drop = FALSE ], bead.af.rows, seed = seed )
  }

  bead.thresholds <- function( neg.fl ) {
    neg.fl <- neg.fl[ .subsample.idx( nrow( neg.fl ), 100000L, 0L ), , drop = FALSE ]
    apply( unmix.ols( neg.fl, spectra ), 2, stats::quantile, probs = 0.995, names = FALSE )
  }

  anchors <- list()
  for ( f in bead.negative.files ) {
    if ( verbose ) message( "  Bead negative: ", f )
    b <- load.events( f, keep.in.cache = TRUE )
    if ( nrow( b$fl ) ) anchors[[ f ]] <- c( b, source = "unstained bead sample" )
  }
  bead.controls <- controls[ type.of[ controls ] == "beads" ]
  bead.own <- list()
  for ( fl.name in bead.controls ) {
    if ( !is.fcs( univ.neg[[ fl.name ]] ) ) next
    f <- univ.neg[[ fl.name ]]
    b <- load.events( f, keep.in.cache = TRUE )
    if ( !nrow( b$fl ) ) next
    bead.own[[ fl.name ]] <- c( b, source = paste( "universal negative of", fl.name ) )
    if ( is.null( anchors[[ f ]] ) ) anchors[[ f ]] <- bead.own[[ fl.name ]]
  }

  load.control <- function( n, fl.name ) {
    if ( exists( fl.name, envir = cache, inherits = FALSE ) ) return( get( fl.name, envir = cache ) )
    got  <- load.events( file.of[[ fl.name ]] )
    keep <- .subsample.idx( nrow( got$fl ), max.control.events, seed + n )
    list( fl = got$fl[ keep, , drop = FALSE ], sc = got$sc[ keep, , drop = FALSE ] )
  }

  # dark levels of every bead control, for the internal negatives
  bead.medians <- list()
  if ( any( !bead.controls %in% names( bead.own ) ) ) {
    for ( fl.name in bead.controls ) {
      n   <- match( fl.name, controls )
      got <- load.control( n, fl.name )
      assign( fl.name, got, envir = cache )
      bead.medians[[ fl.name ]] <- list( median = apply( got$fl, 2, stats::median ),
                                         fluor = match( fl.name, fluorophores ) )
    }
  }

  # ---------------------------------------------------------------------------
  # Single-stain controls
  # ---------------------------------------------------------------------------

  kept <- list()
  internal <- list()
  deferred <- character( 0 )
  sources  <- list( Unstained = unstained.file )

  process <- function( n, fl.name, ctrl, neg, background, raw.thr, unmixed.thr, b.rows, af.pcs,
                       neg.key ) {
    peak <- match( channel[[ fl.name ]], detectors )
    sp <- .holdout.split( nrow( ctrl$fl ), control.train.frac, seed + 1000L + n )
    tr.fl  <- ctrl$fl[ sp$train, , drop = FALSE ]
    pos.tr <- .positive.rows( tr.fl, peak, raw.thr, ceiling )
    pos.te <- .positive.rows( ctrl$fl[ sp$test, , drop = FALSE ], peak, raw.thr, ceiling )
    if ( length( pos.tr ) < 100L || length( pos.te ) < 50L ) {
      note( sprintf( "%s: too few positives (%d training, %d held out)", fl.name,
                     length( pos.tr ), length( pos.te ) ) )
      return( invisible( NULL ) )
    }
    test.rows <- sp$test[ pos.te ][ .subsample.idx( length( pos.te ), max.test.control, seed + n ) ]

    if ( verbose ) message( "  Variants (training split): ", fl.name )
    un <- stats::setNames( "FALSE", fl.name )
    neg.cache <- NULL
    if ( !is.null( neg ) && is.fcs( univ.neg[[ fl.name ]] ) ) {
      un[ fl.name ] <- univ.neg[[ fl.name ]]
      neg.cache <- stats::setNames( list( list( spectral = neg$fl, scatter = neg$sc ) ), univ.neg[[ fl.name ]] )
    }
    raw.thresholds <- stats::setNames( rep( Inf, length( detectors ) ), detectors )
    raw.thresholds[ peak ] <- raw.thr
    variants <- tryCatch(
      get.fluor.variants(
        fl.name,
        file.name          = file.of,
        control.dir        = control.dir,
        asp                = asp,
        spectra            = spectra,
        figures            = FALSE,
        output.dir         = output.dir,
        verbose            = FALSE,
        spectral.channel   = detectors,
        scatter.channel    = colnames( ctrl$sc ),
        universal.negative = un,
        control.type       = type.of,
        raw.thresholds     = raw.thresholds,
        unmixed.thresholds = unmixed.thr,
        flow.channel       = channel,
        af.pcs             = af.pcs,
        neg.cache          = neg.cache,
        use.unmixed        = !anyNA( unmixed.thr ),
        n.cells            = n.cells,
        som.dim            = som.dim,
        k.neighbors        = k.neighbors,
        sim.threshold      = sim.threshold,
        sim.threshold.floor    = sim.threshold.floor,
        af.collinear.threshold = af.collinear.threshold,
        parallel           = parallel,
        threads            = threads,
        control.data       = list( spectral = tr.fl, scatter = ctrl$sc[ sp$train, , drop = FALSE ] )
      ),
      error = function( e ) {
        note( paste0( fl.name, ": variant discovery failed (", conditionMessage( e ), ")" ) )
        NULL
      } )
    delta <- if ( is.null( variants ) ) NULL else
      sweep( as.matrix( variants ), 2, spectra[ fl.name, ], "-" )

    k <- length( kept )
    train.pos <- tr.fl[ pos.tr, , drop = FALSE ]
    test.ev   <- ctrl$fl[ test.rows, , drop = FALSE ]
    sources[[ fl.name ]] <<- file.of[[ fl.name ]]
    kept[[ fl.name ]] <<- list(
      fluor = match( fl.name, fluorophores ), is.cell = type.of[[ fl.name ]] == "cells",
      neg.key = neg.key, background = background, b.rows = b.rows, delta = delta,
      fit = train.pos[ .subsample.idx( nrow( train.pos ), max.fit.events, seed + k ), , drop = FALSE ],
      cal = test.ev[ .subsample.idx( nrow( test.ev ), max.cal.events, seed + k ), , drop = FALSE ] )
    invisible( NULL )
  }

  for ( n in seq_along( controls ) ) {
    fl.name <- controls[ n ]
    if ( verbose ) message( sprintf( "  Controls: %s (%d/%d)", fl.name, n, length( controls ) ) )
    peak <- match( channel[[ fl.name ]], detectors )
    if ( is.na( peak ) ) {
      note( paste0( fl.name, ": peak channel ", channel[[ fl.name ]], " is not a detector of `spectra`" ) )
      next
    }
    ctrl <- load.control( n, fl.name )
    if ( !nrow( ctrl$fl ) ) {
      note( paste0( fl.name, ": no events after selection" ) )
      next
    }

    if ( type.of[[ fl.name ]] == "cells" ) {
      neg.key <- unstained.file
      uneg <- univ.neg[[ fl.name ]]
      if ( is.fcs( uneg ) && uneg != unstained.file ) {
        if ( is.null( negatives[[ uneg ]] ) ) {
          got <- load.events( uneg )
          if ( nrow( got$fl ) ) {
            negatives[[ uneg ]] <- got
            bases.af[[ uneg ]]  <- .af.basis.rows( got$fl, n.af.rows, seed = seed )
          }
        }
        if ( !is.null( negatives[[ uneg ]] ) ) neg.key <- uneg
      }
      neg <- negatives[[ neg.key ]]
      raw.thr <- stats::quantile(
        neg$fl[ .subsample.idx( nrow( neg$fl ), 100000L, 0L ), channel[[ fl.name ]] ], 0.995, names = FALSE )
      process( n, fl.name, ctrl, if ( is.fcs( uneg ) ) neg else NULL, NULL, raw.thr,
               spectral.variants$thresholds,
               NULL, stats::setNames( list( bases.af[[ neg.key ]] ), if ( is.fcs( uneg ) ) uneg else "none" ),
               neg.key )
    } else {
      own <- bead.own[[ fl.name ]]
      if ( is.null( own ) ) {
        others <- Filter( function( m ) m$fluor != match( fl.name, fluorophores ), bead.medians )
        level  <- if ( length( others ) ) .bead.reference.levels( others, spectra )[ peak ] else NA_real_
        found  <- .internal.bead.negative( ctrl$fl, peak, level, factor = bead.factor )
        if ( !is.null( found$mask ) ) {
          own <- list( fl = ctrl$fl[ found$mask, , drop = FALSE ], sc = ctrl$sc[ found$mask, , drop = FALSE ],
                       source = paste0( "internal negative of ", fl.name, " (", found$reason, ")" ) )
          internal[[ paste( fl.name, "negative" ) ]] <- own
          sources[[ paste( fl.name, "negative" ) ]]  <- file.of[[ fl.name ]]
        } else if ( verbose ) {
          message( "  No internal bead negative in ", fl.name, " (", found$reason, ")" )
        }
      }
      if ( is.null( own ) ) {
        deferred <- c( deferred, fl.name )
        next
      }
      raw.thr <- stats::quantile( own$fl[ , peak ], 0.995, names = FALSE )
      process( n, fl.name, ctrl, own, colMeans( own$fl ), raw.thr, bead.thresholds( own$fl ),
               bead.rows( own$fl ), NULL, NULL )
    }
  }

  if ( !length( anchors ) ) anchors <- internal
  for ( fl.name in deferred ) {
    n    <- match( fl.name, controls )
    peak <- match( channel[[ fl.name ]], detectors )
    ctrl <- load.control( n, fl.name )
    if ( length( anchors ) ) {
      neg <- list( fl = do.call( rbind, lapply( anchors, `[[`, "fl" ) ),
                   sc = do.call( rbind, lapply( anchors, `[[`, "sc" ) ) )
      process( n, fl.name, ctrl, neg, colMeans( neg$fl ),
               stats::quantile( neg$fl[ , peak ], 0.995, names = FALSE ), bead.thresholds( neg$fl ),
               bead.rows( neg$fl ), NULL, NULL )
    } else {
      raw.thr <- 10^.otsu.threshold( log10( pmax( ctrl$fl[ , peak ], 1 ) ) )
      note( paste0( fl.name, ": no bead negative; positives above the Otsu split" ) )
      process( n, fl.name, ctrl, NULL, NULL, raw.thr, bead.thresholds( u.train ), NULL, NULL, NULL )
    }
  }
  rm( list = ls( envir = cache ), envir = cache )

  if ( !length( kept ) )
    stop( "No single-stain control has enough positive events.", call. = FALSE )

  # ---------------------------------------------------------------------------
  # Groups
  # ---------------------------------------------------------------------------

  delta.list <- lapply( Filter( function( k ) !is.null( k$delta ) && nrow( k$delta ) >= 3L, kept ),
                        function( k ) k$delta[ -1, , drop = FALSE ] )
  vbases <- if ( length( delta.list ) )
    build.variant.basis( delta.list, spectra = NULL, rank = basis.max.rank,
                         var.explained = basis.var.explained, min.lambda = 1e-4,
                         pooled.fallback = FALSE, verbose = FALSE ) else list()

  lasers <- .detector.lasers( detectors, asp, verbose = FALSE )
  if ( ( laser.split || sigma0.fallback == "laser" ) && is.null( lasers ) ) {
    note( "The laser of some detectors is unknown: no per-laser design rows or same-laser fallback." )
    lasers <- NULL
  }

  fit.groups <- cal.groups <- list()
  for ( fl.name in names( kept ) ) {
    k <- kept[[ fl.name ]]
    af.rows <- if ( k$is.cell ) bases.af[[ k$neg.key ]] else k$b.rows
    design  <- .control.design( spectra[ k$fluor, ], vbases[[ fl.name ]], af.rows,
                                lasers = if ( laser.split ) lasers else NULL, split.min = laser.split.min )
    offset <- if ( !is.null( k$b.rows ) ) NULL else k$background
    cls    <- if ( k$is.cell ) "cells" else "beads"
    fit.groups[[ length( fit.groups ) + 1L ]] <- list( name = fl.name, kind = "control", raw = k$fit,
      design = design, offset = offset, bin.row = 1L, cls = cls )
    cal.groups[[ length( cal.groups ) + 1L ]] <- list( name = fl.name, kind = "control", raw = k$cal,
      design = design, offset = offset, bin.row = 1L, cls = cls )
  }
  fit.groups[[ length( fit.groups ) + 1L ]] <- list( name = "Unstained", kind = "unstained",
    raw = u.train[ .subsample.idx( nrow( u.train ), 4L * max.fit.events, seed ), , drop = FALSE ],
    design = setup.basis, offset = NULL, bin.row = 1L, cls = "cells" )
  cal.groups[[ length( cal.groups ) + 1L ]] <- list( name = "Unstained", kind = "unstained",
    raw = u.test[ .subsample.idx( nrow( u.test ), 2L * max.cal.events, seed ), , drop = FALSE ],
    design = setup.basis, offset = NULL, bin.row = 1L, cls = "cells" )
  for ( a.name in names( anchors ) ) {
    bfl <- anchors[[ a.name ]]$fl
    sp  <- .holdout.split( nrow( bfl ), control.train.frac, seed + 7L )
    design <- bead.rows( bfl )
    b.tr <- bfl[ sp$train, , drop = FALSE ]
    b.te <- bfl[ sp$test, , drop = FALSE ]
    fit.groups[[ length( fit.groups ) + 1L ]] <- list( name = a.name, kind = "beads",
      raw = b.tr[ .subsample.idx( nrow( b.tr ), max.fit.events, seed ), , drop = FALSE ],
      design = design, offset = NULL, bin.row = 1L, cls = "beads" )
    cal.groups[[ length( cal.groups ) + 1L ]] <- list( name = a.name, kind = "beads",
      raw = b.te[ .subsample.idx( nrow( b.te ), max.cal.events, seed ), , drop = FALSE ],
      design = design, offset = NULL, bin.row = 1L, cls = "beads" )
    if ( is.null( sources[[ a.name ]] ) ) sources[[ a.name ]] <- a.name
  }

  # ---------------------------------------------------------------------------
  # Clamped data
  # ---------------------------------------------------------------------------

  clamp <- list( clamped = FALSE, share = NULL )
  if ( clamp.mode == "auto" ) {
    clamp <- .negatives.clamped( c( list( u.fl ), lapply( anchors, `[[`, "fl" ) ) )
    if ( clamp$clamped )
      note( paste0( "Clamped data: the negatives have no values at or below zero on any detector, so ",
                    "their noise is truncated at zero. Use `clamped.weights$noise.model` (",
                    clamp.weights, " weights) for unmixing rather than this model." ) )
  }
  u.proj   <- .noise.project( list( raw = u.test, design = setup.basis ) )
  flat.phi <- stats::median( .robust.var( u.proj$resid ) )

  # ---------------------------------------------------------------------------
  # Fit and calibration
  # ---------------------------------------------------------------------------

  if ( verbose ) message( "  Fitting the noise model" )
  model <- .fit.noise.model( fit.groups, detectors, settings, lasers )

  floor.signal <- clamp.floor.pe * stats::median( model$g )
  scale <- .noise.dispersion.scale( .poisson.noise.model( detectors, model$g, floor.signal ), u.proj )
  clamped.model <- if ( clamp.weights == "poisson" ) {
    g.s <- model$g * scale
    list( read.var = stats::setNames( g.s * floor.signal, detectors ),
          counts.per.unit = stats::setNames( 1 / g.s, detectors ), kappa.pooled = stats::median( 1 / g.s ) )
  } else {
    list( read.var = stats::setNames( rep( flat.phi, length( detectors ) ), detectors ),
          counts.per.unit = stats::setNames( rep( Inf, length( detectors ) ), detectors ), kappa.pooled = Inf )
  }

  if ( verbose ) message( "  Calibrating the noise model on held-out events" )
  calibration <- .calibrate.noise.model( model, cal.groups )

  kappa <- 1 / model$g
  by.class <- NULL
  if ( !is.null( model$sigma0.sq.by.class ) ) {
    by.class <- lapply( seq_along( model$classes ), function( k ) list(
      read.var = stats::setNames( model$sigma0.sq.by.class[ k, ], detectors ),
      counts.per.unit = kappa, kappa.pooled = stats::median( kappa ), c = model$c,
      flags = stats::setNames( model$flags.by.class[ k, ], detectors ) ) )
    names( by.class ) <- model$classes
  }

  noise.model <- c(
    list( read.var = model$sigma0.sq, counts.per.unit = kappa, kappa.pooled = stats::median( kappa ) ),
    model,
    list(
      flag.text      = stats::setNames( .noise.flag.text( model$flags ), detectors ),
      by.class       = by.class,
      calibration    = calibration,
      clamped        = clamp$clamped,
      negative.share = clamp$share,
      flat.phi       = flat.phi,
      clamped.weights = list( kind = clamp.weights, floor = floor.signal, scale = scale,
                              noise.model = clamped.model ),
      selection      = list( mode = selection, counts = do.call( rbind, counts ) ),
      sources        = unlist( sources ),
      unstained      = unstained.file,
      notes          = notes
    ) )
  noise.model$settings <- c( settings, list(
    selection = selection, bead.af.rows = bead.af.rows, laser.split = laser.split,
    laser.split.min = laser.split.min, clamp.mode = clamp.mode, clamp.weights = clamp.weights,
    clamp.floor.pe = clamp.floor.pe, seed = seed ) )

  if ( verbose ) {
    message( sprintf(
      "  Noise model: median sigma0 %.1f, median kappa %.3g, %d flagged detector(s)",
      stats::median( sqrt( model$sigma0.sq ) ), stats::median( kappa ), sum( model$flags > 0 ) ) )
    message( sprintf(
      "  Held-out calibration: %.0f%% of detector x intensity cells in 0.8-1.25 (controls %.0f%%)",
      100 * calibration$in.band, 100 * calibration$in.band.controls ) )
  }

  if ( save ) saveRDS( noise.model, file = file.path( output.dir, "noise_model.rds" ) )

  noise.model
}
