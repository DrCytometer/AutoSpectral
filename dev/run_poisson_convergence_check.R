# run_poisson_fast_convergence_check.R
#
# Diagnoses cells that the current poisson_irls_rcpp_parallel() reverts to
# WLS, using bit-exact row equality against beta_init (no C++ rebuild
# needed). Assumes raw.data, spectra, and weights are already loaded, e.g.
# from a real FCS chunk and the panel's spectral reference library.

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all(asp.rcpp.dir)

Rcpp::compileAttributes(asp.rcpp.dir)
devtools::load_all(asp.rcpp.dir)
devtools::load_all(asp.dir)

asp <- get.autospectral.param( cytometer = "discover" )
asp$scatter.data.max.x <- 5e7
asp$scatter.data.max.y <- 5e7
a8.bd.dir <- "./Plate_BD Beads"
a8.bd.file <- "./Plate_BD Beads/fcs_control_file.csv"
a8.cell.dir <- "./Plate_Cells"
a8.cell.spectra <- read.spectra("A8 Cells_autospectral_spectra.csv")
a8.cell.variants <- readRDS(file.path(a8.cell.dir, "Spectral_variants.rds"))
a8.fc <- reload.flow.control(
  a8.bd.dir,
  a8.bd.file,
  asp
)
a8.concat.cells <- readFCS(file.path("./concatenated_fcs", "A8_concatenated_cells.fcs"))[
  , a8.fc$scatter.and.channel.spectral
]

idx <- sample(nrow(a8.concat.cells), 5000)
raw.data <- a8.concat.cells[idx, a8.fc$spectral.channel]
spectra <- a8.cell.spectra
weights <- calculate.weights(
  file.path("./concatenated_fcs", "A8_concatenated_cells.fcs"),
  a8.fc$spectral.channel
)
ols.unmix <- unmix.ols(raw.data, spectra)
wls.unmix <- unmix.wls.fast( raw.data, spectra, weights )
poisson.r.unmix <- unmix.poisson(
  raw.data = raw.data,
  spectra = spectra,
  asp = asp,
  initial.weights = weights
)
poisson.r.unmix.test <- unmix.poisson.test(
  raw.data = raw.data,
  spectra = spectra,
  asp = asp,
  initial.weights = weights
)

poisson.cpp.unmixed <- poisson_irls_rcpp_parallel(
  raw_data = raw.data,
  spectra = spectra,
  beta_init = wls.unmix,
  maxit = 100,
  tol = 1e-6,
  n_threads = 1,
  divergence_threshold = 1e4
)
colnames(poisson.cpp.unmixed) <- rownames(spectra)
poisson.cpp.unmixed.test <- poisson_irls_rcpp_parallel_test(
  raw_data = raw.data,
  spectra = spectra,
  beta_init = wls.unmix,
  maxit = 100,
  tol = 1e-6,
  n_threads = 1,
  divergence_threshold = 1e4
)
colnames(poisson.cpp.unmixed.test$beta) <- rownames(spectra)


ridge.values <- c(0, 0.01, 0.1, 0.5, 1, 2, 5)

ridge.sweep <- lapply( ridge.values, function( r ) {
  poisson_irls_rcpp_parallel_test(
    raw_data = raw.data, spectra = spectra, beta_init = wls.unmix,
    maxit = 100, tol = 1e-6, n_threads = 1,
    divergence_threshold = 1e4, noise_floor = 125, ridge = r
  )
} )

rmse.vs.glm <- sapply( ridge.sweep, function( res ) {
  sqrt( mean( ( res$beta - poisson.r.unmix.test )^2 ) )
} )
names( rmse.vs.glm ) <- ridge.values
print( rmse.vs.glm )
converged.ridge.sweep <- sapply( ridge.sweep, function( res ) {
  length(res$converged == 0)
} )
names( converged.ridge.sweep ) <- ridge.values
print( converged.ridge.sweep )

poisson.cpp.unmixed.lowfloor <- poisson_irls_rcpp_parallel_test(
  raw_data = raw.data, spectra = spectra, beta_init = wls.unmix,
  maxit = 100, tol = 1e-6, n_threads = 1,
  divergence_threshold = 1e4, noise_floor = 1e-9
)
colnames(poisson.cpp.unmixed.lowfloor$beta) <- rownames(spectra)
table( poisson.cpp.unmixed.lowfloor$converged )

create.biplot(
  poisson.r.unmix.test,
  "Alexa Fluor 647", "APC",
  asp,
  x.min = -15000,
  y.min = -15000
)
create.biplot(
  poisson.r.unmix.test,
  "BV510", "BV480",
  asp,
  x.min = -15000,
  y.min = -15000
)
result <- ridge.sweep[[6]]$beta
colnames(result) <- rownames(spectra)
create.biplot(
  poisson.cpp.unmixed.test$beta,
  "Alexa Fluor 647", "APC",
  asp,
  x.min = -15000,
  y.min = -15000
)
create.biplot(
  poisson.cpp.unmixed.lowfloor$beta,
  "BV510", "BV480",
  asp,
  x.min = -15000,
  y.min = -15000
)


# bit-exact match to WLS beta_init means that cell was reverted internally
reverted <- rowSums( abs( poisson.cpp.unmixed.test$beta - wls.unmix ) ) < 1e-12

cat( sprintf(
  "Reverted: %d / %d cells (%.2f%%)\n",
  sum( reverted ), length( reverted ), 100 * mean( reverted )
) )

# proxy for the ill-conditioning hypothesis: spread of the per-detector
# Poisson weight (1/eta) implied by the WLS starting point, per cell
eta.wls <- tcrossprod( wls.unmix, t( spectra ) )
eta.wls <- pmax( eta.wls, 1e-5 )
weight.spread <- apply( eta.wls, 1, function( e ) max( 1 / e ) / min( 1 / e ) )

diagnostic.table <- data.frame(
  cell.index     = seq_len( nrow( raw.data ) ),
  reverted       = reverted,
  weight.spread  = weight.spread,
  total.raw      = rowSums( raw.data )
)

# per-fluorophore contribution to reverted vs. converged cells
per.fluor <- data.frame(
  fluorophore = colnames( wls.unmix ),
  mean.abund.reverted  = colMeans( wls.unmix[ reverted, , drop = FALSE ] ),
  mean.abund.converged = colMeans( wls.unmix[ !reverted, , drop = FALSE ] )
)

write.csv( diagnostic.table, "poisson_fast_reverted_cells.csv", row.names = FALSE )
write.csv( per.fluor, "poisson_fast_reverted_by_fluorophore.csv", row.names = FALSE )

cat( "Weight-spread by revert status:\n" )
print( tapply( weight.spread, reverted, summary ) )

# what's failing?
table( poisson.cpp.unmixed.test$converged )
#    0    2 
# 2202 2798 
# after changes, all passing

for ( nf in c( 125, 500, 2000, 5000 ) ) {
  res <- poisson_irls_rcpp_parallel_test(
    raw_data = raw.data, spectra = spectra, beta_init = wls.unmix,
    maxit = 100, tol = 1e-6, n_threads = 1,
    divergence_threshold = 1e4, noise_floor = nf
  )
  cat( sprintf( "noise_floor=%d: %.1f%% reverted\n",
                nf, 100 * mean( rowSums( abs( res$beta - wls.unmix ) ) < 1e-12 ) ) )
}

reverted.r <- rowSums( abs( poisson.r.unmix.test - wls.unmix ) ) < 1e-12
cat( sprintf( "R glm.fit reverted: %d / %d (%.2f%%)\n",
              sum( reverted.r ), length( reverted.r ), 100 * mean( reverted.r ) ) )
