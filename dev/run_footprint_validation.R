# run_footprint_validation.R
#
# Validates two changes against the real Aurora bead/cell benchmark, found
# and tested in a Python reproduction this session (see the session report
# and corr_sig_core.py):
#
#   1. spillover.spread supplied for BOTH substrates (previously bead runs
#      had no bead_spillover_spread.csv at all; a real one now exists).
#   2. footprint.frac = 0.02, restricting the slope fit and the held-out
#      step search to each dye's own emission footprint.
#
# Calls the actual exported correct.unmixing.signatures() twice per
# direction (old defaults vs new), not a diagnostic reimplementation, so
# this is a direct test of the production function once the diff in the
# session report is applied.
#
# Run top to bottom. Requires the footprint.frac / footprint.min.channels
# diff already applied to correct_unmixing_signatures.R.

# ---------------------------------------------------------------------------
# 0. Setup - EDIT THIS SECTION (mirrors run_signature_correction_diagnostics_aurora.R)
# ---------------------------------------------------------------------------

asp.dir      <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all( asp.rcpp.dir )
devtools::load_all( asp.dir )
asp <- get.autospectral.param( cytometer = "aurora" )

cell.dir  <- "./PBMC direct stain"
bd.dir    <- "./BD beads"

cell.spectra <- read.spectra( "Cells_spectra_corrections.csv" )
bd.spectra   <- read.spectra( "BD_spectra_corrections.csv" )

concat.cells <- readFCS( "./concatenated_fcs/small_Concatenated_cells.fcs" )
concat.beads <- readFCS( "./concatenated_fcs/small_Concatenated_BD_beads.fcs" )

diag.spectra <- list( Cells = cell.spectra, Beads = bd.spectra )
diag.raw     <- list( Cells = concat.cells,  Beads = concat.beads )

diag.unstained <- list(
  Cells = readFCS( file.path( cell.dir, "A10 unstained_010_Cells.fcs" ) ),
  Beads = readFCS( file.path( bd.dir,   "A10 unstained_010_Beads_1.fcs" ) )
)

# spillover.spread per substrate. bead_spillover_spread.csv now exists
# alongside cell_spillover_spread.csv - this is the change under test in
# point 1 above, so both must be real, non-NULL matrices here.
spillover.spread <- list(
  Cells = read.csv( file.path(cell.dir, "cell_spillover_spread.csv"), row.names = 1, check.names = FALSE ),
  Beads = read.csv( file.path(bd.dir, "bead_spillover_spread.csv"), row.names = 1, check.names = FALSE )
)
spillover.spread <- lapply( spillover.spread, as.matrix )

af.name              <- "AF"
diag.threshold.probs <- 0.995
diag.seed            <- 42L

# ---------------------------------------------------------------------------
# Helpers (subset of run_signature_correction_diagnostics_aurora.R's)
# ---------------------------------------------------------------------------

.val.angle <- function( a, b ) {
  cs <- rowSums( a * b ) / ( sqrt( rowSums( a^2 ) ) * sqrt( rowSums( b^2 ) ) )
  180 / pi * acos( pmin( 1, pmax( -1, cs ) ) )
}

.val.recovered.fraction <- function( start, corrected, truth ) {
  ang.start <- .val.angle( start, truth )
  ang.end   <- .val.angle( corrected, truth )
  ( ang.start - ang.end ) / ang.start
}

b.shared.fluor <- intersect( rownames( diag.spectra$Cells ), rownames( diag.spectra$Beads ) )
b.shared.det   <- intersect( colnames( diag.spectra$Cells ), colnames( diag.spectra$Beads ) )

# ---------------------------------------------------------------------------
# Per-direction run: old defaults vs new (spillover.spread + footprint.frac)
# ---------------------------------------------------------------------------

val.results <- list()

for ( target in c( "Cells", "Beads" ) ) {

  wrong  <- if ( target == "Cells" ) "Beads" else "Cells"
  s.true <- as.matrix( diag.spectra[[ target ]][ , b.shared.det, drop = FALSE ] )
  s.wrong <- as.matrix( diag.spectra[[ wrong  ]][ , b.shared.det, drop = FALSE ] )

  scatter.det <- grep( "^(FSC|SSC)", colnames( diag.raw[[ target ]] ), value = TRUE )

  set.seed( diag.seed )
  raw.all     <- diag.raw[[ target ]][ , b.shared.det, drop = FALSE ]
  scatter.all <- diag.raw[[ target ]][ , scatter.det, drop = FALSE ]

  unst         <- diag.unstained[[ target ]][ , b.shared.det, drop = FALSE ]
  unst.um      <- unmix.ols.fast( unst, s.wrong )
  thr          <- apply( unst.um, 2, stats::quantile, probs = diag.threshold.probs )

  panel <- setdiff( rownames( s.wrong ), af.name )

  run.one <- function( spillover.spread.arg, footprint.frac.arg ) {
    set.seed( diag.seed )
    correct.unmixing.signatures(
      raw.data            = raw.all,
      spectra             = s.wrong,
      unmixed.thresholds  = thr,
      asp                 = asp,
      af.name             = af.name,
      spillover.spread    = spillover.spread.arg,
      bg.mode             = if ( target == "Cells" ) "scatter.knn" else "global.mean",
      scatter             = scatter.all,
      unstained           = unst,
      unstained.scatter   = diag.unstained[[ target ]][ , scatter.det, drop = FALSE ],
      true.spectra        = s.true,
      footprint.frac       = footprint.frac.arg,
      verbose             = FALSE
    )
  }

  fit.old <- run.one( spillover.spread.arg = NULL,                     footprint.frac.arg = 0 )
  fit.new <- run.one( spillover.spread.arg = spillover.spread[[ target ]], footprint.frac.arg = 0.02 )

  rec.old <- fit.old$recovery
  rec.new <- fit.new$recovery
  names( rec.old )[ names( rec.old ) %in% c( "recovered", "accepted" ) ] <-
    c( "recovered.old", "accepted.old" )
  names( rec.new )[ names( rec.new ) %in% c( "recovered", "accepted" ) ] <-
    c( "recovered.new", "accepted.new" )

  combo <- merge( rec.old[ , c( "fluorophore", "deg.start", "recovered.old", "accepted.old" ) ],
                   rec.new[ , c( "fluorophore", "recovered.new", "accepted.new" ) ],
                   by = "fluorophore" )
  combo <- combo[ order( combo$deg.start ), ]

  cat( sprintf( "\n===== %s (starting spectra = %s) =====\n", target, wrong ) )
  print( combo, digits = 4, row.names = FALSE )
  cat( sprintf( "\n  accepted: old %d/%d -> new %d/%d\n",
                sum( combo$accepted.old ), nrow( combo ),
                sum( combo$accepted.new ), nrow( combo ) ) )

  n.regressed <- sum( combo$recovered.new < combo$recovered.old - 0.02, na.rm = TRUE )
  cat( sprintf( "  SAFETY CHECK: %d fluorophore(s) moved meaningfully AWAY from truth ",
                n.regressed ) )
  cat( "relative to the old result (recovered.new < recovered.old - 0.02).\n" )
  if ( n.regressed > 0 )
    print( combo[ combo$recovered.new < combo$recovered.old - 0.02, ], row.names = FALSE )

  val.results[[ target ]] <- list( old = fit.old, new = fit.new, table = combo )
}

cat( "\n\nINTERPRETATION: the Python reproduction (n_per_dye=4000, idealised\n" )
cat( "IID noise) found 6/16 -> 11-12/16 accepted per direction from these two\n" )
cat( "changes, with no fluorophore moving away from truth. This script checks\n" )
cat( "whether the same pattern - more fluorophores accepted, none regressing -\n" )
cat( "holds on the real benchmark, where residual noise is structured (real\n" )
cat( "variant scatter) rather than IID, so effect sizes are not expected to\n" )
cat( "match exactly. The SAFETY CHECK line above is the one that matters most:\n" )
cat( "if it ever reports a regression, do not adopt footprint.frac = 0.02 as\n" )
cat( "the default without investigating that fluorophore specifically.\n" )
