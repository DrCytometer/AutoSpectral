
asp.dir      <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"

if ( requireNamespace( "devtools", quietly = TRUE ) ) {
  try( devtools::load_all( asp.rcpp.dir ), silent = TRUE )
  devtools::load_all( asp.dir )
}
library( ggplot2 )

asp <- get.autospectral.param()

raw.lung.list <- readFCS("./FCS files/Raw stained/C1 Lung_WT_001_Samples.fcs", return.keywords = TRUE)
ols.lung.list <- readFCS("./FCS files/OLS stained/C1 Lung_WT_001_Samples.fcs", return.keywords = TRUE)
raw.lung <- raw.lung.list$data
ols.lung <- ols.lung.list$data

raw.brain.list <- readFCS("./FCS files/Raw stained/E1 Brain_WT_001_Samples.fcs", return.keywords = TRUE)
ols.brain.list <- readFCS("./FCS files/OLS stained/E1 Brain_WT_001_Samples.fcs", return.keywords = TRUE)
raw.brain <- raw.brain.list$data
ols.brain <- ols.brain.list$data

# copy $PnN/$PnR/etc for a set of column names from a source FCS keyword list
# into base.keywords at new column positions -- used below when combining
# columns from two FCS files so the appended columns keep correct metadata
# rather than falling back to generic placeholders
append.param.keywords <- function(base.keywords, source.keywords, new.cols, start.index) {
  
  field.types <- c("N", "S", "B", "E", "R", "G", "V", "DISPLAY", "TYPE")
  
  for (j in seq_along(new.cols)) {
    col.name <- new.cols[j]
    new.idx  <- start.index + j - 1L
    
    source.pN <- grep("^\\$P\\d+N$", names(source.keywords), value = TRUE)
    match.key <- source.pN[vapply(source.pN, function(k) identical(source.keywords[[k]], col.name), logical(1))]
    
    if (length(match.key) == 1) {
      old.idx <- sub("^\\$P(\\d+)N$", "\\1", match.key)
      for (ft in field.types) {
        old.key <- paste0("$P", old.idx, ft)
        if (!is.null(source.keywords[[old.key]])) {
          base.keywords[[paste0("$P", new.idx, ft)]] <- source.keywords[[old.key]]
        }
      }
    } else {
      base.keywords[[paste0("$P", new.idx, "N")]] <- col.name
      base.keywords[[paste0("$P", new.idx, "E")]] <- "0,0"
      base.keywords[[paste0("$P", new.idx, "R")]] <- as.character(asp$expr.data.max)
    }
  }
  
  base.keywords
}

# combine and write
same.cols <- intersect(colnames(raw.lung), colnames(ols.lung))
new.cols.combined <- colnames(ols.lung)[!(colnames(ols.lung) %in% same.cols)]

combined.mat <- cbind(
  raw.lung,
  ols.lung[, new.cols.combined, drop = FALSE]
)
combined.keywords <- append.param.keywords(
  raw.lung.list$keywords, ols.lung.list$keywords, new.cols.combined, ncol(raw.lung) + 1L
)
combined.keywords[["$FIL"]] <- "Combined stained lung.fcs"
writeFCS(
  mat = combined.mat, keys = combined.keywords,
  file.name = "Combined stained lung.fcs", output.dir = "./FCS files"
)

combined.brain <- cbind(
  raw.brain,
  ols.brain[, new.cols.combined, drop = FALSE]
)
combined.keywords.brain <- append.param.keywords(
  raw.brain.list$keywords, ols.brain.list$keywords, new.cols.combined, ncol(raw.brain) + 1L
)
combined.keywords.brain[["$FIL"]] <- "Combined stained brain.fcs"
writeFCS(
  mat = combined.brain, keys = combined.keywords.brain,
  file.name = "Combined stained brain.fcs", output.dir = "./FCS files"
)

# gate extensively in flowjo, export the populations of interest
# start again here
library(AutoSpectral)
library(ggplot2)
asp <- get.autospectral.param()
spectra <- read.spectra("Cells AF Removed autospectral spectra.csv")
spectral.channels <- colnames(spectra)
pbmc.af <- read.spectra("YT PBMC AF.csv")

unstained.alv.mac <- readFCS("./FCS files/Exported/alvMacUnstained.fcs")
stained.alv.mac <- readFCS("./FCS files/Exported/alvMacStained.fcs")

alv.mac.af <- colMeans(unstained.alv.mac[,spectral.channels])
alv.mac.af <- alv.mac.af/max(alv.mac.af)
alv.mac.af


unstained.microglia <- readFCS("./FCS files/Exported/MicrogliaUnstained.fcs")
stained.microglia <- readFCS("./FCS files/Exported/MicrogliaStained.fcs")
unstained.pmbc <- readFCS("./FCS files/Raw/A1 Unstained (Cells)_Controls.fcs")

microglia.af <- colMeans(unstained.microglia[,spectral.channels])
microglia.af <- microglia.af/max(microglia.af)
microglia.af
spectral.variant.plot.dens(
  pbmc.af,
  alv.mac.af,
  save = TRUE,
  title = "Alveolar Macrophage vs PBMC Autofluorescence",
  plot.height = 3
)

spectral.variant.plot.dens(
  pbmc.af,
  microglia.af,
  save = TRUE,
  title = "Microglia vs PBMC Autofluorescence",
  plot.height = 3
)


# run per cell AF on these isolated cells with full set of AF spectra to use
# only take forward cells that match to the mean same cell type AF in isolation
# this is our ground truth, pure population
microglia.ground.truth <- AutoSpectralRcpp::unmix.autospectral.rcpp(
  raw.data = unstained.microglia[,spectral.channels],
  spectra = spectra,
  af.spectra = rbind(microglia.af, pbmc.af)
)
microglia.idx <- which(microglia.ground.truth[, "AF Index"] == 1)
length(microglia.idx)/nrow(unstained.microglia) # 0.956927
alv.mac.ground.truth <- AutoSpectralRcpp::unmix.autospectral.rcpp(
  raw.data = unstained.alv.mac[,spectral.channels],
  spectra = spectra,
  af.spectra = rbind(alv.mac.af, pbmc.af)
)
alv.mac.idx <- which(alv.mac.ground.truth[, "AF Index"] == 1)
length(alv.mac.idx)/nrow(unstained.alv.mac) # 0.9177632

# import other cells to cover them up
unstained.pmbc <- readFCS("./FCS files/Raw/A1 Unstained (Cells)_Controls.fcs")

# keep the set of mac or microglia constant
# vary the number of other cells present

length(alv.mac.idx)
nrow(unstained.pmbc)


af.accuracy <- function(test.data, filler.data, percentages,
                        spectra, test.af.spectra, reference.af = NULL,
                        near.exact.threshold = 0.999, verbose = TRUE) {
  n.test <- nrow(test.data)
  if (is.null(reference.af)) reference.af <- test.af.spectra[1, ]
  
  ref.cosine <- cosine.similarity(rbind(reference.af, test.af.spectra))[1, -1]
  
  result <- lapply(percentages, function(i) {
    n.filler.use <- ceiling(n.test / (i / 100))
    if (n.filler.use > nrow(filler.data)) {
      if (verbose) {
        message(sprintf(
          "  percentage %.3g%% would need %d filler cells but only %d are available -- using all of them.",
          i, n.filler.use, nrow(filler.data)
        ))
      }
      n.filler.use <- nrow(filler.data)
    }
    data.to.unmix <- rbind(
      test.data[, spectral.channels],
      filler.data[1:n.filler.use, spectral.channels]
    )
    unmixed <- AutoSpectralRcpp::unmix.autospectral.rcpp(
      raw.data = data.to.unmix,
      spectra = spectra,
      af.spectra = test.af.spectra,
      verbose = verbose
    )
    af.idx <- unmixed[1:n.test, "AF Index"]
    matched.cosine <- ref.cosine[af.idx]
    
    exactly.correct <- mean(matched.cosine >= near.exact.threshold) * 100
    
    c(exactly.correct, mean(matched.cosine))
  })
  names(result) <- percentages
  
  return(result)
}

percentages <-c(0.1, 0.3, 1, 3, 10)
# unstained alveolar macrophages
unstained.alv.mac.res <- af.accuracy(
  test.data = unstained.alv.mac[alv.mac.idx,],
  filler.data = unstained.pmbc,
  percentages = percentages,
  spectra = spectra,
  test.af.spectra = rbind(alv.mac.af, pbmc.af)
)
unstained.alv.mac.res

# unstained microglia
unstained.microglia.res <- af.accuracy(
  test.data = unstained.microglia[microglia.idx,],
  filler.data = unstained.pmbc,
  percentages = percentages,
  spectra = spectra,
  test.af.spectra = rbind(microglia.af, pbmc.af)
)
unstained.microglia.res


# plotting results

#' Plot one or more dose-response style datasets as a side-by-side bar chart
#'
#' @param data_list A list of lists. Each inner list should look like
#'   unstained.alv.mac.res, i.e. named elements ("0.1", "0.3", "1", ...)
#'   where the first element of each entry is the y-value to plot.
#' @param set_names Character vector of legend names, one per entry in data_list.
#'   Defaults to "Set 1", "Set 2", ... if not supplied.
#' @param colors Character vector of colors, one per entry in data_list.
#'   Optional -- if NULL, ggplot's default palette is used.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param xlim Optional numeric vector of length 2, e.g. c(0.1, 3). Restricts
#'   which x-categories are shown, based on their numeric value. NULL (default)
#'   shows all categories.
#' @param ylim Optional numeric vector of length 2, e.g. c(0, 1). Sets the
#'   y-axis display range without dropping data. NULL (default) auto-scales.
#'
#' @return A ggplot object.
plot_dose_response <- function(data_list,
                               set_names = NULL,
                               colors = NULL,
                               xlab = "Concentration",
                               ylab = "Value",
                               xlim = NULL,
                               ylim = NULL) {
  
  if (is.null(set_names)) {
    set_names <- paste0("Set ", seq_along(data_list))
  }
  if (length(set_names) != length(data_list)) {
    stop("set_names must be the same length as data_list")
  }
  if (!is.null(colors) && length(colors) != length(data_list)) {
    stop("colors must be the same length as data_list")
  }
  
  # Build a tidy data frame: x = numeric version of the list names,
  # y = first element of each entry, set = which dataset it came from
  df <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    lst <- data_list[[i]]
    
    if (is.null(names(lst)) || any(names(lst) == "")) {
      stop(sprintf(
        "data_list[[%d]] (\"%s\") has no names (or has unnamed elements). ",
        i, set_names[i]
      ), "Each entry needs names like list(`0.1` = ..., `0.3` = ...). ",
      "Run names(", "data_list[[", i, "]]) to check.")
    }
    
    x_vals <- suppressWarnings(as.numeric(names(lst)))
    if (any(is.na(x_vals))) {
      stop(sprintf(
        "data_list[[%d]] (\"%s\") has names that don't convert to numbers: %s",
        i, set_names[i], paste(names(lst)[is.na(x_vals)], collapse = ", ")
      ))
    }
    
    data.frame(
      x   = x_vals,
      y   = vapply(lst, function(v) v[1], numeric(1)),
      set = set_names[i],
      stringsAsFactors = FALSE
    )
  }))
  
  df$set <- factor(df$set, levels = set_names)
  
  # Optionally restrict which x-categories are included, based on their
  # numeric value (applied before converting x to a factor).
  if (!is.null(xlim)) {
    if (length(xlim) != 2) stop("xlim must be a numeric vector of length 2, e.g. c(0.1, 3)")
    df <- df[df$x >= xlim[1] & df$x <= xlim[2], ]
    if (nrow(df) == 0) stop("No data left after applying xlim -- check the range against your x values")
  }
  
  # x is treated as a discrete axis (ordered by numeric value) so that bars
  # for each set can sit side-by-side at each x position without overlapping.
  x_levels <- sort(unique(df$x))
  df$x <- factor(df$x, levels = x_levels)
  
  p <- ggplot(df, aes(x = x, y = y, fill = set, group = set)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    labs(x = xlab, y = ylab, fill = NULL) +
    theme_classic()
  
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  if (!is.null(ylim)) {
    if (length(ylim) != 2) stop("ylim must be a numeric vector of length 2, e.g. c(0, 1)")
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  p
}

plot_dose_response(
  list(unstained.microglia.res, unstained.alv.mac.res),
  set_names = c("Microglia", "Alveolar Macrophage"),
  colors = c("black", "blue"),
  xlab = "Percentage of total cells",
  ylab = "AF assignment accuracy (%)",
  xlim = c(0.1, 10),   # optional, restricts which x-categories show
  ylim = c(0, 100)     # optional, sets the y-axis range
)


# =============================================================================
# Spike-in test 2: SOM training composition and AF library recovery
# =============================================================================
#
# Part 1 shows that per-cell AF *assignment* is deterministic and does not
# depend on how rare the target cell type is among the cells being unmixed,
# provided a matching spectrum already exists in af.spectra. This section
# asks the complementary question: if the unstained sample used to *train*
# af.spectra (via get.af.spectra()'s SOM) does not contain the target cell
# type at all, how many "correct" cells need to be spiked into the training
# sample before the SOM discovers a node that matches the true AF spectrum?
#
# Microglia (tissue macrophage) are spiked into a PBMC training pool at
# increasing fractions. At each fraction, get.af.spectra() is re-run on the
# mixed population and the best-matching resulting spectrum is compared to
# ground-truth microglia.af by cosine similarity. The same trained library is
# then used to unmix held-out microglia cells (never seen during training),
# diluted into a separate PBMC pool, to confirm SOM discovery -- not per-cell
# assignment -- is the limiting step.

training.spike.percentages <- c(0, 0.03, 0.1, 0.3, 1, 3, 10, 30)
test.reserve.n       <- min(500, floor(length(microglia.idx) * 0.3))
test.dilution.pct    <- 1
near.exact.threshold <- 0.999
spike.dir            <- "./FCS files/SpikeIn training"
if (!dir.exists(spike.dir)) dir.create(spike.dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Disjoint pools: cells used to build the spiked training file must never
# overlap with cells used afterward to test whether unmixing recovers them,
# and the same applies to the PBMC background at each stage.
# ---------------------------------------------------------------------------

set.seed(asp$bird.seed)

microglia.pure      <- unstained.microglia[microglia.idx, ]
microglia.test.rows <- sample(nrow(microglia.pure), test.reserve.n)
microglia.test.pool  <- microglia.pure[microglia.test.rows, ]
microglia.train.pool <- microglia.pure[-microglia.test.rows, ]

pmbc.split        <- sample(nrow(unstained.pmbc))
pmbc.train.pool.n <- floor(nrow(unstained.pmbc) * 0.7)
pmbc.train.pool   <- unstained.pmbc[pmbc.split[seq_len(pmbc.train.pool.n)], ]
pmbc.test.pool    <- unstained.pmbc[pmbc.split[-seq_len(pmbc.train.pool.n)], ]

n.train.total <- min(20000, nrow(pmbc.train.pool))

missing.cols <- setdiff(colnames(pmbc.train.pool), colnames(microglia.train.pool))
if (length(missing.cols) > 0) {
  stop(
    "Microglia and PBMC raw files do not share the same channel set; missing: ",
    paste(missing.cols, collapse = ", "),
    call. = FALSE
  )
}
microglia.train.pool <- microglia.train.pool[, colnames(pmbc.train.pool)]

pmbc.keywords <- as.list(readFCSheader(
  "./FCS files/Raw/A1 Unstained (Cells)_Controls.fcs"
)[[1]])

# ---------------------------------------------------------------------------
# Write a spiked training file: n.spike cells of microglia.train.pool plus
# enough PBMC filler to reach n.train.total cells. Draws with replacement
# (and warns) if a requested spike fraction exceeds the held-out microglia
# training pool.
# ---------------------------------------------------------------------------

write.spike.training.file <- function(pct, spike.pool, filler.pool, n.total,
                                      keywords, file.name, output.dir, asp) {
  
  n.spike  <- round(pct / 100 * n.total)
  n.filler <- n.total - n.spike
  
  replace.spike <- n.spike > nrow(spike.pool)
  if (replace.spike) {
    message(sprintf(
      "  spike-in %.2f%% requires %d microglia but only %d are available in the training pool -- sampling with replacement.",
      pct, n.spike, nrow(spike.pool)
    ))
  }
  set.seed(42)
  spike.rows <- if (n.spike > 0) sample(nrow(spike.pool), n.spike, replace = replace.spike) else integer(0)
  spike.mat  <- spike.pool[spike.rows, , drop = FALSE]
  
  filler.mat   <- filler.pool[seq_len(min(n.filler, nrow(filler.pool))), , drop = FALSE]
  combined.mat <- rbind(spike.mat, filler.mat)
  
  keys <- keywords
  keys[["$TOT"]] <- as.character(nrow(combined.mat))
  keys[["$FIL"]] <- file.name
  
  writeFCS(mat = combined.mat, keys = keys, file.name = file.name, output.dir = output.dir)
  
  file.path(output.dir, file.name)
}

# ---------------------------------------------------------------------------
# Main loop: train, score library recovery, test per-cell unmixing recovery
# ---------------------------------------------------------------------------

refine.values <- c(FALSE, TRUE)

run.som.recovery <- function(pct, refine.on) {
  
  message(sprintf("Training spike-in: %.2f%% (refine = %s)", pct, refine.on))
  
  file.name <- sprintf(
    "Microglia_spike_%s_refine_%s.fcs",
    gsub(".", "p", format(pct, scientific = FALSE), fixed = TRUE),
    refine.on
  )
  spike.path <- write.spike.training.file(
    pct         = pct,
    spike.pool  = microglia.train.pool,
    filler.pool = pmbc.train.pool,
    n.total     = n.train.total,
    keywords    = pmbc.keywords,
    file.name   = file.name,
    output.dir  = spike.dir,
    asp = asp
  )
  
  af.spectra.trained <- get.af.spectra(
    unstained.sample = spike.path,
    asp              = asp,
    spectra          = spectra,
    dist             = 2,
    refine           = refine.on,
    af.assign.method = "l1",
    k.neighbors = 90,
    figures          = FALSE,
    save             = FALSE,
    verbose          = FALSE,
    plot.unmixed = FALSE
  )
  print(nrow(af.spectra.trained))
  
  best.node.cosine <- max(cosine.similarity(rbind(microglia.af, af.spectra.trained))[1, -1])
  
  unmix.result <- af.accuracy(
    test.data            = microglia.test.pool,
    filler.data          = pmbc.test.pool,
    percentages          = test.dilution.pct,
    spectra              = spectra,
    test.af.spectra      = af.spectra.trained,
    reference.af         = microglia.af,
    near.exact.threshold = near.exact.threshold,
    verbose              = FALSE
  )[[1]]
  
  c(
    best.node.cosine    = best.node.cosine,
    near.exact.pct      = unmix.result[1],
    mean.matched.cosine = unmix.result[2]
  )
}

som.recovery <- lapply(refine.values, function(refine.on) {
  res <- lapply(training.spike.percentages, run.som.recovery, refine.on = refine.on)
  names(res) <- training.spike.percentages
  res
})
names(som.recovery) <- c("no.refine", "refine")

som.recovery.df <- do.call(rbind, lapply(names(som.recovery), function(tag) {
  res <- som.recovery[[tag]]
  data.frame(
    spike.percent       = training.spike.percentages,
    refine               = (tag == "refine"),
    best.node.cosine     = vapply(res, `[[`, numeric(1), "best.node.cosine"),
    near.exact.pct       = vapply(res, `[[`, numeric(1), "near.exact.pct"),
    mean.matched.cosine  = vapply(res, `[[`, numeric(1), "mean.matched.cosine")
  )
}))
write.csv(som.recovery.df, "SOM_training_spikein_recovery_l2_cos_k30.csv", row.names = FALSE)
som.recovery.df

# ---------------------------------------------------------------------------
# Plot: cosine similarity to ground-truth microglia AF vs. spike-in
# percentage during SOM training. Two series: the best matching spectrum
# anywhere in the trained library (does the SOM find the shape at all), and
# the spectrum that held-out microglia cells actually get matched to during
# per-cell unmixing (does assignment route to it).
# ---------------------------------------------------------------------------

extract.series <- function(tag, field) {
  res <- lapply(som.recovery[[tag]], `[`, field)
  names(res) <- names(som.recovery[[tag]])
  res
}

plot_dose_response(
  list(
    extract.series("no.refine", "best.node.cosine"),
    extract.series("refine",    "best.node.cosine"),
    extract.series("no.refine", "mean.matched.cosine"),
    extract.series("refine",    "mean.matched.cosine")
  ),
  set_names = c(
    "Best AF library match (no refine)",
    "Best AF library match (refine)",
    "Per-cell unmixing match (no refine)",
    "Per-cell unmixing match (refine)"
  ),
  colors = c("grey70", "grey30", "steelblue1", "steelblue4"),
  xlab = "Microglia spike-in during SOM training (%)",
  ylab = "Cosine similarity to ground-truth AF",
  ylim = c(0, 1)
)



