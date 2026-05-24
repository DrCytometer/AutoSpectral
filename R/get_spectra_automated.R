# get_spectra_automated.R

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

## Cosine similarity of each row of mat against a reference vector.
.cosine.sim.rows <- function( mat, ref.vec ) {
  dot.prod <- mat %*% ref.vec
  mat.norm <- sqrt( rowSums( mat^2 ) )
  ref.norm <- sqrt( sum( ref.vec^2 ) )
  as.numeric( dot.prod / ( mat.norm * ref.norm + 1e-9 ) )
}


## RLM fallback: regress each channel against peak.channel, return
## normalised [0,1] spectrum vector.
.rlm.refine.spectrum <- function(
    spectral.events,
    peak.channel,
    spectral.channels,
    rlm.iter.max
) {
  spec.vec <- stats::setNames( rep( 0, length( spectral.channels ) ), spectral.channels )
  peak.expr <- spectral.events[ , peak.channel ]

  for ( ch in spectral.channels ) {
    if ( ch == peak.channel ) {
      spec.vec[ ch ] <- 1.0
    } else {
      ch.expr <- spectral.events[ , ch ]
      coef <- fit.robust.linear.model(
        peak.expr, ch.expr, peak.channel, ch, rlm.iter.max
      )
      spec.vec[ ch ] <- coef[ 2 ]
    }
  }

  max.val <- max( spec.vec )
  if ( max.val > 0 ) spec.vec <- spec.vec / max.val
  spec.vec
}


## Load the spectral reference library for a cytometer (identified by db.col).
## Returns a matrix (fluorophores x spectral.channels) normalised to [0,1],
## or NULL when no library file exists.
.load.ref.library <- function( db.col, spectral.channels ) {
  ref.prefix <- switch(
    db.col,
    A5SE          = "Symphony",
    NorthernLights = "Aurora",
    db.col
  )

  db.extdata <- system.file( "extdata", package = "AutoSpectral" )
  ref.path   <- file.path( db.extdata,
                            paste0( ref.prefix, "_spectral_reference_library.csv" ) )

  if ( !file.exists( ref.path ) ) return( NULL )

  ref.fluors <- utils::read.csv( ref.path, row.names = 1, check.names = FALSE )
  common.cols <- intersect( spectral.channels, colnames( ref.fluors ) )
  if ( length( common.cols ) == 0 ) return( NULL )

  ref.mat <- as.matrix( ref.fluors[ , common.cols, drop = FALSE ] )
  t( apply( ref.mat, 1, function( x ) {
    mx <- max( x, na.rm = TRUE )
    if ( !is.finite( mx ) || mx <= 0 ) return( x )
    x / mx
  } ) )
}


## Map db.col to the asp$cytometer string expected by spectral.reference.plot.
.db.col.to.cytometer <- function( db.col ) {
  switch(
    db.col,
    Aurora         = "Aurora",
    NorthernLights = "Aurora",
    ID7000         = "ID7000",
    Discover       = "FACSDiscover A8",
    Opteon         = "Opteon",
    Mosaic         = "Mosaic",
    Xenith         = "Xenith",
    A5SE           = "Symphony",
    "Aurora"
  )
}


## Parse the unstained source for a single row of the control table.
## Returns list( type = "file" | "internal", file = <filename or NULL> ).
.parse.unstained.source <- function( uneg.val ) {
  if ( is.na( uneg.val ) || !nzchar( trimws( uneg.val ) ) )
    return( list( type = "internal", file = NULL ) )

  uneg.logical <- suppressWarnings( as.logical( uneg.val ) )
  if ( !is.na( uneg.logical ) && isTRUE( uneg.logical ) )
    return( list( type = "global.true", file = NULL ) )

  list( type = "file", file = trimws( uneg.val ) )
}


## Plot the cosine-similarity filter for one fluorophore using base R graphics.
.cosine.filter.single.plot <- function(
    top.mat,
    af.median.i,
    fluor.name,
    empirical.peak
) {
  cs.vals <- .cosine.sim.rows( top.mat, af.median.i )

  # non-uniform quantile bins: finer resolution at low cosine values
  bin.probs <- seq( 0.1, 1, length.out = 19 )
  bin.probs <- c( seq( 0, 0.01, length.out = 9 ), bin.probs / 10, bin.probs )
  bin.probs <- sort( unique( bin.probs ) )
  breaks    <- stats::quantile( cs.vals, probs = bin.probs )

  if ( length( unique( breaks ) ) < 3L ) {
    graphics::plot.new()
    graphics::text( 0.5, 0.5, fluor.name, adj = 0.5, cex = 0.8 )
    return( invisible( NULL ) )
  }

  groups <- as.integer( cut( cs.vals, breaks = breaks, include.lowest = TRUE ) )
  n.grp  <- max( groups, na.rm = TRUE )

  # median per bin (base R replacement for collapse::fmedian)
  res <- do.call( rbind, lapply( seq_len( n.grp ), function( g ) {
    idx <- which( groups == g )
    if ( length( idx ) == 0L ) return( NULL )
    apply( top.mat[ idx, , drop = FALSE ], 2, stats::median )
  } ) )

  if ( is.null( res ) || nrow( res ) == 0L ) {
    graphics::plot.new()
    return( invisible( NULL ) )
  }

  # edge-case: a detector appearing as peak in only one bin (artifact/unstable)
  # drop it unless it is the first bin (least AF-like) -- see BUV615::Siglec F
  if ( any( table( max.col( res ) ) == 1 ) ) {
    drop.i <- as.numeric( names( which.min( table( max.col( res ) ) ) ) )
    if ( which( max.col( res ) == drop.i ) != 1L )
      res <- res[ max.col( res ) != drop.i, , drop = FALSE ]
  }

  # peak detector from binned median matrix (as in cosine.similarity.filter)
  i.max         <- which( res == max( res, na.rm = TRUE ), arr.ind = TRUE )
  detector.peak <- colnames( res )[ i.max[ 1L, 2L ] ]

  # RdYlBu palette matching RColorBrewer::brewer.pal(11, "RdYlBu") + rev():
  # blue = least AF-like (row 1); red = most AF-like (last row)
  rdylbu <- c( "#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090",
               "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695" )
  colors <- rev( grDevices::colorRampPalette( rdylbu )( nrow( res ) ) )

  r <- range( res, na.rm = TRUE )

  for ( j in seq_len( nrow( res ) ) ) {
    x <- res[ j, ]
    if ( j == 1L ) {
      graphics::plot(
        x,
        type = "l",
        col  = colors[ j ],
        ylim = r,
        xlab = "Detectors",
        ylab = ifelse( max( x, na.rm = TRUE ) == 1,
                       "Expression [0,1]", "Expression" ),
        main = fluor.name,
        sub  = sprintf( "Detector (peak): %s", detector.peak ),
        xaxt = "n"
      )
      graphics::axis(
        1,
        at       = seq_along( x ),
        labels   = sub( "-A$", "", names( x ) ),
        las      = 2,
        cex.axis = 0.75
      )
    } else {
      graphics::lines( x, col = colors[ j ] )
    }
  }

  graphics::abline( v = which( colnames( res ) %in% detector.peak ),
                    lty = "dashed" )

  invisible( NULL )
}


## Render all per-fluorophore cosine-filter plots into a single multi-panel PDF.
.save.cosine.filter.pdf <- function( cosine.filter.data, filepath, ncol.pdf = 4L ) {
  valid <- cosine.filter.data[ !sapply( cosine.filter.data, is.null ) ]
  if ( length( valid ) == 0L ) return( invisible( NULL ) )

  nrow.pdf <- ceiling( length( valid ) / ncol.pdf )

  grDevices::pdf( filepath, width = ncol.pdf * 3, height = nrow.pdf * 2.5 )
  oldpar <- graphics::par(
    mfrow = c( nrow.pdf, ncol.pdf ),
    mar   = c( 3, 3, 2, 0.5 ),
    oma   = c( 0, 0, 0, 0 )
  )
  on.exit( { graphics::par( oldpar ); grDevices::dev.off() }, add = TRUE )

  for ( d in valid ) {
    tryCatch(
      .cosine.filter.single.plot(
        d$top.mat, d$af.median.i, d$fluor, d$empirical.peak
      ),
      error = function( e ) {
        graphics::plot.new()
        graphics::text( 0.5, 0.5, d$fluor, adj = 0.5 )
      }
    )
  }

  invisible( NULL )
}


# ---------------------------------------------------------------------------
# Exported function
# ---------------------------------------------------------------------------

#' @title Get Fluorophore Spectra - Automated Workflow
#'
#' @description
#' Single-function replacement for the three-step
#' `define.flow.control` -> `clean.controls` -> `get.fluorophore.spectra` pipeline.
#' Extracts normalised `[0, 1]` fluorophore reference spectra from a directory of
#' single-stained control FCS files without requiring scatter gating or interactive
#' input.
#'
#' The algorithm is based on the `spectracle` approach, adapted to use
#' AutoSpectral's internal FCS reader, matrix-based data layout, and the
#' control-file metadata already required by the rest of the AutoSpectral
#' pipeline:
#'
#' 1. Reads all single-stained controls with `readFCS()`.
#' 2. Removes saturating events.
#' 3. For each fluorophore, determines the AF reference from the paired
#'    universal-negative file specified in `universal.negative` column of the
#'    control table, or - if that column is empty - from the lower 25% of
#'    the control file's own events ranked by the expected peak channel
#'    (internal-negative mode).
#' 4. Performs projection-based AF orthogonalisation to identify the empirical
#'    peak detector; selects top-expressing events filtered by lowest cosine
#'    similarity to the AF vector; performs KNN scatter-matched AF subtraction.
#' 5. Applies QC: if the normalised signal in the expected peak detector is
#'    below `peak.signal.threshold` **or** the cosine similarity against the
#'    spectral reference library is below `cosine.threshold`, the spectrum is
#'    refined using the RLM approach along the expected peak channel.
#' 6. Returns the spectra matrix in AutoSpectral format (fluorophores x detectors).
#'
#' The legacy workflow (`define.flow.control` + `clean.controls` +
#' `get.fluorophore.spectra`) remains available and is unchanged.
#'
#' @param control.dir Character. Path to the directory containing the
#'   single-stained control FCS files.
#' @param control.def.file Character. Path to (or filename of) the control
#'   definition CSV. Must already exist and pass `check.control.file()`.
#' @param asp The AutoSpectral parameter list from `get.autospectral.param()`.
#' @param n.candidates Integer, default `1000`. Number of top-expressing
#'   candidate events selected per fluorophore before cosine-similarity
#'   filtering. Ignored in internal-negative mode, where the top 5%% of
#'   events by peak channel are used directly.
#' @param n.spectral Integer, default `200`. Number of spectral events
#'   retained after filtering for low AF cosine similarity.
#' @param k.neighbors Integer, default `2`. Number of nearest neighbours in
#'   scatter space used for per-event AF subtraction.
#' @param cosine.threshold Numeric, default `0.9`. Minimum cosine similarity
#'   against the spectral reference library to accept the automated spectrum;
#'   values below this threshold trigger RLM refinement.
#' @param peak.signal.threshold Numeric, default `0.5`. Minimum normalised
#'   signal in the expected peak detector to accept the automated spectrum;
#'   values below this threshold trigger RLM refinement.
#' @param top.expressing.override Named numeric vector or `NULL` (default).
#'   Override the event count for specific samples. Names should match the FCS
#'   filename (without extension) or the fluorophore name from the control file.
#' @param figures Logical, default `TRUE`. Produce spectral trace, heatmap,
#'   cosine-similarity, and reference-library QC plots.
#' @param plot.cosine.filter Logical, default `TRUE`. When `figures = TRUE`,
#'   produce a multi-panel PDF showing the per-fluorophore cosine-similarity
#'   filter traces (binned by cosine similarity to AF, coloured from least
#'   to most AF-like). Set to `FALSE` to skip.
#' @param plot.scatter.match Logical, default `TRUE`. When `figures = TRUE`,
#'   produce a multi-panel PDF showing the KNN scatter-matched AF subtraction
#'   for each fluorophore (unstained background, selected spectral events, and
#'   their matched AF events). Set to `FALSE` to skip.
#' @param verbose Logical, default `TRUE`. Print progress messages.
#' @param print.timings Logical, default `FALSE`. When `TRUE`, prints a
#'   `tictoc` timing log to the console at the end of the call, showing
#'   elapsed time for each major processing stage.
#'
#' @return A numeric matrix with fluorophores in rows and spectral detector
#'   channels in columns, values normalised to `[0, 1]` (L-infinity norm,
#'   peak = 1). Compatible with all downstream AutoSpectral functions.
#'
#' @importFrom tictoc tic toc tic.clear tic.clearlog tic.log
#'
#' @seealso [get.fluorophore.spectra()] for the legacy workflow.
#'   [get.cytometer.param()] for the cytometer auto-detection step.
#'   [spectral.reference.plot()] for the QC report produced when `figures = TRUE`.
#'
#' @export

get.spectra.automated <- function(
    control.dir,
    control.def.file,
    asp,
    n.candidates            = 1000L,
    n.spectral              = 200L,
    k.neighbors             = 2L,
    cosine.threshold        = 0.9,
    peak.signal.threshold   = 0.5,
    top.expressing.override = NULL,
    figures                 = TRUE,
    plot.cosine.filter      = TRUE,
    plot.scatter.match      = TRUE,
    verbose                 = TRUE,
    print.timings           = FALSE
) {

  tictoc::tic.clear()
  tictoc::tic.clearlog()

  # -- 0. Validate inputs
  if ( !dir.exists( control.dir ) )
    stop( "control.dir does not exist: ", control.dir, call. = FALSE )

  check.control.file( control.dir, control.def.file, asp, legacy = FALSE )

  # -- 1. Read control table
  ctrl.path <- if ( file.exists( control.def.file ) ) {
    control.def.file
  } else {
    file.path( control.dir, control.def.file )
  }

  ctrl.tbl <- utils::read.csv(
    ctrl.path, stringsAsFactors = FALSE, strip.white = TRUE
  )

  # fluorophore rows: exclude any row labelled as the universal negative itself
  # (rows where the file IS the unstained, identified by having universal.negative
  # == TRUE boolean) and exclude AF / negative rows
  uneg.col.all  <- ctrl.tbl$universal.negative
  uneg.bool.all <- suppressWarnings( as.logical( uneg.col.all ) )
  is.neg.row    <- !is.na( uneg.bool.all ) & isTRUE( uneg.bool.all )

  fluor.rows <- which(
    !is.neg.row &
    !grepl( "negative", ctrl.tbl$fluorophore, ignore.case = TRUE ) &
    !grepl( "^AF$",     ctrl.tbl$fluorophore, ignore.case = TRUE )
  )

  if ( length( fluor.rows ) == 0 )
    stop( "No fluorophore rows found in ", control.def.file, call. = FALSE )

  fluor.names    <- ctrl.tbl$fluorophore[ fluor.rows ]
  fluor.channels <- ctrl.tbl$channel[     fluor.rows ]
  fluor.files    <- ctrl.tbl$filename[    fluor.rows ]

  # parse per-fluorophore unstained source
  unstained.sources <- lapply( fluor.rows, function( i ) {
    src <- .parse.unstained.source( ctrl.tbl$universal.negative[ i ] )

    # if still global.true, resolve to the row(s) marked TRUE
    if ( src$type == "global.true" ) {
      true.rows <- which( !is.na( uneg.bool.all ) & uneg.bool.all )
      if ( length( true.rows ) > 0 ) {
        return( list( type = "file",
                      file = ctrl.tbl$filename[ true.rows[ 1L ] ] ) )
      }
      return( list( type = "internal", file = NULL ) )
    }
    src
  } )

  # collect all unique unstained filenames (non-internal)
  unique.unstained.files <- unique( vapply(
    unstained.sources,
    function( s ) if ( s$type == "file" ) s$file else NA_character_,
    character( 1L )
  ) )
  unique.unstained.files <- unique.unstained.files[ !is.na( unique.unstained.files ) ]

  n.internal <- sum( vapply( unstained.sources,
                              function( s ) s$type == "internal", logical( 1L ) ) )

  if ( length( unique.unstained.files ) == 0 && n.internal == 0 )
    stop(
      "No universal negative defined in ", control.def.file,
      ". Populate the universal.negative column or leave blank for ",
      "internal-negative mode.",
      call. = FALSE
    )

  # -- 2. Detect cytometer from first FCS file
  if ( verbose )
    message( "\033[34m-- Detecting cytometer --\033[0m" )

  first.fcs.path <- file.path( control.dir, fluor.files[ 1L ] )

  db.path <- system.file(
    "extdata", "cytometer_database.csv", package = "AutoSpectral"
  )
  if ( !nchar( db.path ) ) db.path <- NULL

  cyt.param <- tryCatch(
    get.cytometer.param( first.fcs.path, db.path = db.path, verbose = verbose ),
    error = function( e ) {
      if ( verbose )
        message(
          "\033[33m  Auto-detection failed: ", e$message,
          "\n  Falling back to asp parameters for: ", asp$cytometer, "\033[0m"
        )
      NULL
    }
  )

  if ( is.null( cyt.param ) ) {
    # derive channels from the FCS file using asp filter patterns,
    # mirroring the approach in reload.flow.control()
    all.channels <- colnames( readFCS( first.fcs.path ) )
    non.spectral.pattern <- paste0( asp$non.spectral.channel, collapse = "|" )
    spec.idx <- grep( non.spectral.pattern, all.channels, invert = TRUE )
    spectral.channels <- all.channels[ spec.idx ]
    if ( grepl( "Discover", asp$cytometer ) ) {
      spec.idx <- grep( asp$spectral.channel, spectral.channels )
      spectral.channels <- spectral.channels[ spec.idx ]
    }
    spectral.channels <- check.channels( spectral.channels, asp )
    scatter.channels  <- read.scatter.parameter( asp )
    sat.value         <- if ( !is.null( asp$saturation.ceiling ) )
                           asp$saturation.ceiling else Inf
    db.col <- if      ( grepl( "Northern", asp$cytometer, ignore.case = TRUE ) ) "NorthernLights"
              else if ( grepl( "Aurora",   asp$cytometer, ignore.case = TRUE ) ) "Aurora"
              else if ( grepl( "ID7000",   asp$cytometer, ignore.case = TRUE ) ) "ID7000"
              else if ( grepl( "Discover", asp$cytometer, ignore.case = TRUE ) ) "Discover"
              else if ( grepl( "Opteon",   asp$cytometer, ignore.case = TRUE ) ) "Opteon"
              else if ( grepl( "Mosaic",   asp$cytometer, ignore.case = TRUE ) ) "Mosaic"
              else if ( grepl( "Xenith",   asp$cytometer, ignore.case = TRUE ) ) "Xenith"
              else if ( grepl( "Symphony", asp$cytometer, ignore.case = TRUE ) ) "A5SE"
              else "Aurora"
  } else {
    spectral.channels <- cyt.param$cols.detector
    scatter.channels  <- cyt.param$scatter.param
    sat.value         <- cyt.param$sat.value
    db.col            <- cyt.param$db.col
  }

  if ( verbose )
    message( sprintf(
      "\033[32m  %d spectral channels | scatter: %s\033[0m",
      length( spectral.channels ),
      paste( scatter.channels, collapse = ", " )
    ) )

  # -- 3. Read FCS files
  if ( verbose )
    message( "\033[34m-- Reading FCS files --\033[0m" )

  tictoc::tic( "read FCS files" )

  cols.keep <- c( scatter.channels, spectral.channels )

  .read.fcs.clean <- function( path, label ) {
    mat     <- readFCS( path )
    present <- intersect( cols.keep, colnames( mat ) )
    mat     <- mat[ , present, drop = FALSE ]
    spec.present <- intersect( spectral.channels, colnames( mat ) )
    if ( length( spec.present ) > 0 ) {
      keep <- rowSums( mat[ , spec.present, drop = FALSE ] >= sat.value ) == 0
      mat  <- mat[ keep, , drop = FALSE ]
    }
    if ( verbose )
      message( sprintf( "\033[32m  %-40s  %d events\033[0m",
                         label, nrow( mat ) ) )
    mat
  }

  # pre-load unique unstained files
  unstained.cache <- list()
  for ( uf in unique.unstained.files ) {
    uf.path <- file.path( control.dir, uf )
    if ( !file.exists( uf.path ) ) {
      warning( "Unstained file not found, skipping: ", uf, call. = FALSE )
      next
    }
    unstained.cache[[ uf ]] <- .read.fcs.clean( uf.path,
                                                  paste0( "Unstained (", uf, ")" ) )
  }

  # pre-compute AF reference vectors per unique unstained file
  af.cache <- list()
  for ( uf in names( unstained.cache ) ) {
    ust.mat         <- unstained.cache[[ uf ]]
    spec.in.ust     <- intersect( spectral.channels, colnames( ust.mat ) )
    scat.in.ust     <- intersect( scatter.channels,  colnames( ust.mat ) )
    af.spec.uf      <- ust.mat[ , spec.in.ust,  drop = FALSE ]
    af.scat.uf      <- ust.mat[ , scat.in.ust,  drop = FALSE ]
    af.cache[[ uf ]] <- list(
      spectral = af.spec.uf,
      scatter  = af.scat.uf,
      mean     = colMeans( af.spec.uf ),
      median   = apply( af.spec.uf, 2, stats::median )
    )
  }

  # read fluorophore FCS files
  fluor.data <- vector( "list", length( fluor.rows ) )
  names( fluor.data ) <- fluor.names

  for ( i in seq_along( fluor.rows ) ) {
    fcs.path.i        <- file.path( control.dir, fluor.files[ i ] )
    fluor.data[[ i ]] <- .read.fcs.clean( fcs.path.i, fluor.names[ i ] )
  }

  tictoc::toc( log = TRUE, quiet = TRUE )

  # -- 4. Extract spectrum per fluorophore
  if ( verbose )
    message( "\033[34m-- Extracting spectra --\033[0m" )

  tictoc::tic( "extract spectra" )

  ref.lib <- .load.ref.library( db.col, spectral.channels )

  spectra.list         <- vector( "list", length( fluor.names ) )
  names( spectra.list ) <- fluor.names

  cosine.filter.data <- vector( "list", length( fluor.names ) )

  qc.log <- data.frame(
    Fluorophore   = fluor.names,
    EmpiricalPeak = NA_character_,
    ExpectedPeak  = fluor.channels,
    PeakSignal    = NA_real_,
    CosineSim     = NA_real_,
    Status        = "OK",
    stringsAsFactors = FALSE
  )

  for ( i in seq_along( fluor.names ) ) {
    fluor  <- fluor.names[ i ]
    expeak <- fluor.channels[ i ]
    mat.i  <- fluor.data[[ i ]]

    spec.i <- intersect( spectral.channels, colnames( mat.i ) )
    scat.i <- intersect( scatter.channels,  colnames( mat.i ) )

    if ( nrow( mat.i ) == 0 || length( spec.i ) == 0 ) {
      warning( "No usable events for '", fluor, "'.", call. = FALSE )
      spectra.list[[ i ]] <- stats::setNames(
        rep( 0, length( spectral.channels ) ), spectral.channels
      )
      qc.log$Status[ i ] <- "NO_DATA"
      next
    }

    spec.data <- mat.i[ , spec.i, drop = FALSE ]

    # -- 4a. Resolve AF reference for this fluorophore
    src.i <- unstained.sources[[ i ]]

    if ( src.i$type == "file" && src.i$file %in% names( af.cache ) ) {

      af.ref       <- af.cache[[ src.i$file ]]
      af.mean.i    <- af.ref$mean[   spec.i ]
      af.median.i  <- af.ref$median[ spec.i ]
      af.spectral.i <- af.ref$spectral[ , spec.i, drop = FALSE ]
      common.scat   <- intersect( scat.i, colnames( af.ref$scatter ) )
      af.scatter.i  <- af.ref$scatter[ , common.scat, drop = FALSE ]
      use.internal  <- FALSE

    } else {

      # -- Internal negative: lower 25% by expected peak channel
      if ( verbose )
        message( sprintf(
          "\033[33m  Internal negative for %s (no universal.negative)\033[0m",
          fluor
        ) )

      if ( expeak %in% spec.i ) {
        peak.vals <- spec.data[ , expeak ]
      } else {
        peak.vals <- rowMeans( spec.data )
      }

      n.total   <- nrow( spec.data )
      n.neg     <- max( 2L, floor( n.total * 0.25 ) )
      n.pos.int <- max( 1L, floor( n.total * 0.05 ) )

      order.peak      <- order( peak.vals )
      i.internal.neg  <- order.peak[ seq_len( n.neg ) ]
      i.internal.pos  <- utils::tail( order.peak, n.pos.int )

      neg.mat      <- spec.data[ i.internal.neg, , drop = FALSE ]
      af.mean.i    <- colMeans( neg.mat )
      af.median.i  <- apply( neg.mat, 2, stats::median )
      af.spectral.i <- neg.mat
      af.scatter.i  <- mat.i[ i.internal.neg, scat.i, drop = FALSE ]
      common.scat   <- scat.i
      use.internal  <- TRUE
    }

    v.unit.i <- af.mean.i / ( sqrt( sum( af.mean.i^2 ) ) + 1e-9 )

    # -- 4b. AF orthogonalisation -> empirical peak
    proj     <- spec.data %*% v.unit.i
    mat.orth <- spec.data - proj %*% t( v.unit.i )

    empirical.peak            <- names( which.max( colMeans( mat.orth ) ) )
    qc.log$EmpiricalPeak[ i ] <- empirical.peak

    # -- 4c. Select top candidate events
    if ( use.internal ) {

      i.top        <- i.internal.pos
      n.cand.actual <- length( i.top )

    } else {

      n.cand.i <- n.candidates
      n.spec.i <- n.spectral

      if ( !is.null( top.expressing.override ) ) {
        file.id <- tools::file_path_sans_ext( fluor.files[ i ] )
        if ( file.id %in% names( top.expressing.override ) ) {
          n.spec.i <- as.integer( top.expressing.override[[ file.id ]] )
          n.cand.i <- max( n.cand.i, n.spec.i )
        } else if ( fluor %in% names( top.expressing.override ) ) {
          n.spec.i <- as.integer( top.expressing.override[[ fluor ]] )
          n.cand.i <- max( n.cand.i, n.spec.i )
        }
      }

      n.cand.actual <- min( n.cand.i, nrow( spec.data ) )
      peak.col      <- if ( empirical.peak %in% spec.i ) empirical.peak else spec.i[ 1L ]
      i.top <- order( spec.data[ , peak.col ],
                      decreasing = TRUE )[ seq_len( n.cand.actual ) ]
    }

    # -- 4d. Cosine-similarity filter: keep least-AF-like events
    top.mat  <- spec.data[ i.top, , drop = FALSE ]
    cs.vs.af <- .cosine.sim.rows( top.mat, af.median.i )

    n.spec.actual  <- min( n.spectral, length( i.top ) )
    i.spectral     <- i.top[ order( cs.vs.af )[ seq_len( n.spec.actual ) ] ]
    spectral.events <- spec.data[ i.spectral, , drop = FALSE ]

    # -- 4e. Collect cosine filter data (rendered to PDF after loop)
    if ( figures && plot.cosine.filter ) {
      cosine.filter.data[[ i ]] <- list(
        top.mat       = top.mat,
        af.median.i   = af.median.i,
        fluor         = fluor,
        empirical.peak = empirical.peak
      )
    }

    # -- 4f. KNN scatter-matched background (AF) subtraction
    knn.idx.i <- NULL

    if ( nrow( af.scatter.i ) > 0 && length( common.scat ) >= 1L ) {

      fluor.scatter <- mat.i[ i.spectral, common.scat, drop = FALSE ]

      knn.idx.i <- FNN::knnx.index(
        data  = as.matrix( af.scatter.i[ , common.scat, drop = FALSE ] ),
        query = as.matrix( fluor.scatter ),
        k     = k.neighbors
      )

      n.ev       <- nrow( knn.idx.i )
      af.matched <- matrix( 0, n.ev, length( spec.i ) )
      colnames( af.matched ) <- spec.i

      for ( ki in seq_len( k.neighbors ) ) {
        af.matched <- af.matched +
          af.spectral.i[ knn.idx.i[ , ki ], spec.i, drop = FALSE ]
      }
      af.matched <- af.matched / k.neighbors

      spectral.sub <- spectral.events - af.matched

      # -- 4g. Scatter match plot (saved immediately per fluorophore)
      if ( figures && plot.scatter.match ) {
        matched.scatter <- af.scatter.i[
          unique( as.vector( knn.idx.i ) ), common.scat, drop = FALSE
        ]
        tryCatch(
          scatter.match.plot(
            pos.expr.data = fluor.scatter,
            neg.expr.data = matched.scatter,
            fluor.name    = fluor,
            scatter.param = common.scat,
            asp           = asp
          ),
          error = function( e )
            message( "\033[31m  scatter.match.plot error [", fluor, "]: ",
                     e$message, "\033[0m" )
        )
      }

    } else {
      spectral.sub <- spectral.events
    }

    # -- 4h. Compute spectrum: column median -> normalise [0,1]
    raw.spectrum <- apply( spectral.sub, 2, stats::median )
    max.val      <- max( raw.spectrum )

    if ( !is.finite( max.val ) || max.val <= 0 ) {
      if ( verbose )
        message( sprintf(
          "\033[31m  Non-positive maximum for %s; triggering RLM refinement\033[0m",
          fluor
        ) )
      spectrum           <- stats::setNames( rep( 0, length( spec.i ) ), spec.i )
      qc.log$Status[ i ] <- "RLM_REFINEMENT"
    } else {
      spectrum             <- raw.spectrum / max.val
      spectrum[ spectrum < 0 ] <- 0
    }

    # -- 5. QC checks
    peak.sig <- if ( expeak %in% names( spectrum ) )
      spectrum[ expeak ] else max( spectrum )

    qc.log$PeakSignal[ i ] <- round( peak.sig, 4 )

    cs.lib <- NA_real_
    if ( !is.null( ref.lib ) && fluor %in% rownames( ref.lib ) ) {
      common.ref <- intersect( spec.i, colnames( ref.lib ) )
      if ( length( common.ref ) > 0L ) {
        lib.spec <- ref.lib[ fluor, common.ref ]
        lib.spec[ is.na( lib.spec ) ] <- 0
        cs.lib <- cosine.similarity(
          rbind( spectrum[ common.ref ], lib.spec )
        )[ 1L, 2L ]
      }
    }
    qc.log$CosineSim[ i ] <- if ( !is.na( cs.lib ) ) round( cs.lib, 4 ) else NA

    needs.refinement <-
      ( is.finite( peak.sig ) && peak.sig < peak.signal.threshold ) ||
      ( !is.na( cs.lib )      && cs.lib   < cosine.threshold       ) ||
      qc.log$Status[ i ] == "RLM_REFINEMENT"

    if ( needs.refinement ) {
      if ( verbose )
        message( sprintf(
          "\033[31m  QC flag [%s]: peak=%.3f, cos=%s -> RLM refinement\033[0m",
          fluor, peak.sig,
          if ( is.na( cs.lib ) ) "N/A" else sprintf( "%.3f", cs.lib )
        ) )

      if ( expeak %in% spec.i ) {
        n.rlm    <- min( n.spectral, nrow( spec.data ) )
        i.rlm    <- order( spec.data[ , expeak ],
                           decreasing = TRUE )[ seq_len( n.rlm ) ]
        spectrum <- .rlm.refine.spectrum(
          spec.data[ i.rlm, , drop = FALSE ], expeak, spec.i, asp$rlm.iter.max
        )
        qc.log$Status[ i ] <- "RLM_REFINEMENT"
      } else {
        if ( verbose )
          message( sprintf(
            "\033[31m  Expected peak '%s' absent in file for %s; ",
            "retaining automated spectrum\033[0m",
            expeak, fluor
          ) )
        qc.log$Status[ i ] <- "QC_FAIL"
      }

    } else {
      if ( verbose )
        message( sprintf(
          "\033[32m  OK [%s]: peak=%.3f, cos=%s\033[0m",
          fluor, peak.sig,
          if ( is.na( cs.lib ) ) "N/A" else sprintf( "%.3f", cs.lib )
        ) )
    }

    # expand spectrum back to full spectral.channels vector
    full.spectrum <- stats::setNames( rep( 0, length( spectral.channels ) ),
                               spectral.channels )
    full.spectrum[ names( spectrum ) ] <- spectrum
    spectra.list[[ i ]] <- full.spectrum
  }

  tictoc::toc( log = TRUE, quiet = TRUE )

  # -- 6. Assemble matrix
  if ( verbose )
    message( "\033[34m-- Assembling spectra matrix --\033[0m" )

  marker.spectra <- do.call( rbind, spectra.list )
  rownames( marker.spectra ) <- make.unique( fluor.names )

  # -- 7. Pairwise cosine similarity warning
  sim.mat  <- cosine.similarity( marker.spectra )
  uniq.sim <- sim.mat * lower.tri( sim.mat )
  sim.idx  <- which( uniq.sim > asp$similarity.warning.n, arr.ind = TRUE )

  if ( nrow( sim.idx ) > 0 ) {
    f1   <- rownames( sim.mat )[ sim.idx[ , 1L ] ]
    f2   <- colnames( sim.mat )[ sim.idx[ , 2L ] ]
    vals <- sim.mat[ sim.idx ]
    message(
      "\033[31mHigh cosine similarity (> ", asp$similarity.warning.n,
      ") detected:\033[0m\n",
      paste( sprintf( "  %s <-> %s: %.4f", f1, f2, vals ), collapse = "\n" )
    )
  }

  # -- 8. QC summary
  if ( verbose ) {
    ok.n    <- sum( qc.log$Status == "OK" )
    rlm.n   <- sum( qc.log$Status == "RLM_REFINEMENT" )
    fail.n  <- sum( qc.log$Status == "QC_FAIL" )
    nodat.n <- sum( qc.log$Status == "NO_DATA" )
    message(
      "\n\033[34m-- QC Summary --\033[0m\n",
      sprintf(
        "  \033[32mOK: %d  \033[33mRLM refined: %d  \033[31mFailed: %d  No data: %d\033[0m",
        ok.n, rlm.n, fail.n, nodat.n
      )
    )
    print( qc.log[ , c( "Fluorophore", "EmpiricalPeak", "ExpectedPeak",
                         "PeakSignal", "CosineSim", "Status" ) ] )
  }

  # -- 9. Figures
  if ( figures ) {
    if ( verbose )
      message( "\033[34m-- Generating figures --\033[0m" )

    tictoc::tic( "generate figures" )

    if ( !dir.exists( asp$figure.spectra.dir ) )
      dir.create( asp$figure.spectra.dir, recursive = TRUE )
    if ( !dir.exists( asp$figure.similarity.heatmap.dir ) )
      dir.create( asp$figure.similarity.heatmap.dir, recursive = TRUE )
    if ( plot.scatter.match && !is.null( asp$figure.scatter.dir.base ) &&
         !dir.exists( asp$figure.scatter.dir.base ) )
      dir.create( asp$figure.scatter.dir.base, recursive = TRUE )

    tryCatch(
      expr = {
        spectral.trace(
          spectral.matrix           = marker.spectra,
          asp                       = asp,
          title                     = paste( "Automated", asp$spectra.file.name,
                                             sep = "_" ),
          plot.dir                  = asp$figure.spectra.dir,
          split.lasers              = TRUE,
          figure.spectra.line.size  = asp$figure.spectra.line.size,
          figure.spectra.point.size = asp$figure.spectra.point.size
        )

        spectral.heatmap(
          marker.spectra,
          title    = paste( "Automated", asp$spectra.file.name, sep = "_" ),
          plot.dir = asp$figure.spectra.dir
        )

        cosine.similarity.plot(
          marker.spectra,
          filename      = asp$similarity.heatmap.file.name,
          "Automated",
          output.dir    = asp$figure.similarity.heatmap.dir,
          figure.width  = asp$figure.similarity.width,
          figure.height = asp$figure.similarity.height
        )

        hotspot.matrix <- calculate.hotspot.matrix( marker.spectra )
        hotspot.max    <- max( 4, hotspot.matrix )

        create.heatmap(
          hotspot.matrix,
          title         = "Automated_Hotspot_Matrix",
          legend.label  = expression( "Hotspot Matrix"^"TM" ),
          triangular    = TRUE,
          plot.dir      = asp$figure.similarity.heatmap.dir,
          fixed.scale   = TRUE,
          scale.min     = 0,
          scale.max     = hotspot.max,
          color.palette = "inferno",
          figure.width  = asp$figure.similarity.width,
          figure.height = asp$figure.similarity.height
        )

        asp.qc           <- asp
        asp.qc$cytometer <- .db.col.to.cytometer( db.col )
        spectral.reference.plot(
          marker.spectra,
          asp.qc,
          plot.dir = asp$figure.spectra.dir,
          filename = "Automated_spectral_qc_report.pdf"
        )
      },
      error = function( e ) {
        message( "\033[31mError during figure generation: ", e$message, "\033[0m" )
      }
    )

    # cosine filter multi-panel PDF
    if ( plot.cosine.filter ) {
      tryCatch(
        .save.cosine.filter.pdf(
          cosine.filter.data,
          filepath = file.path( asp$figure.spectra.dir,
                                "Automated_cosine_filter.pdf" )
        ),
        error = function( e )
          message( "\033[31mError saving cosine filter PDF: ",
                   e$message, "\033[0m" )
      )
    }

    tictoc::toc( log = TRUE, quiet = TRUE )
  }

  # -- 10. Save CSV
  if ( !dir.exists( asp$table.spectra.dir ) )
    dir.create( asp$table.spectra.dir, recursive = TRUE )

  utils::write.csv(
    marker.spectra,
    file = file.path(
      asp$table.spectra.dir,
      paste0( "Automated_", asp$spectra.file.name, ".csv" )
    )
  )

  if ( verbose )
    message( "\033[32m-- Automated spectra extraction complete --\033[0m" )

  if ( print.timings ) {
    log <- tictoc::tic.log( format = TRUE )
    print( unlist( log ) )
  }

  return( marker.spectra )
}
