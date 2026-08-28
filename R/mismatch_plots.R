#' @noRd
.mismatch.plot.names <- function( x ) {
  if ( is.matrix( x ) ) rownames( x ) else names( x )
}

#' @noRd
.mismatch.plot.values <- function( x, fluors ) {
  if ( is.matrix( x ) ) as.numeric( x[ fluors, ] ) else as.numeric( x[ fluors ] )
}

#' @title Plot Spectral Mismatch, Angle, Variability, Alignment, and
#'   Brightness by Dye Class
#'
#' @description
#' Produces a standard set of diagnostic plots comparing a particle type's
#' spectral mismatch, spectral angle, variability, variability/mismatch
#' alignment, cosine similarity, and brightness against fluorophore dye
#' class, plus every pairwise correlation between those six metrics with
#' simple linear-fit statistics. Each individual plot is saved as a JPEG in
#' `output.dir`, and all plots are additionally combined into a single
#' multi-panel PDF report.
#'
#' @param cosine.data One-column numeric matrix of cosine similarity
#'   values, rownames are fluorophore names, as returned by
#'   `assess.mismatch()`.
#' @param angle.data One-column numeric matrix of spectral angle values
#'   (degrees), rownames are fluorophore names, as returned by
#'   `assess.mismatch.angle()`.
#' @param mismatch.data Named numeric vector of per-fluorophore mismatch
#'   magnitudes, e.g. `rowSums(abs(bead.cell.dist(...)))`.
#' @param variability.data Named numeric vector of per-fluorophore
#'   variability magnitudes, e.g.
#'   `rowSums(abs(assess.variability.mad(...)))`.
#' @param alignment.data One-column numeric matrix of per-fluorophore
#'   variability/mismatch alignment values, rownames are fluorophore names,
#'   as returned by `assess.variability.alignment()`.
#' @param brightness.data One-column numeric matrix of per-fluorophore
#'   brightness (MFI) values, rownames are fluorophore names, as returned
#'   by `get.brightness.automated()`.
#' @param fluor.df Data frame with at least `Fluorophore` and `Class`
#'   columns, used to annotate each fluorophore with a dye class for the
#'   violin plots.
#' @param particle.name Character. Particle type label used in plot
#'   titles and output filenames. Default `"UltraComp"`.
#' @param output.dir Character. Directory for the saved JPEGs and PDF
#'   report. Default `"./results/aurora"`.
#' @param mismatch.limits Numeric vector of length 2, y-axis limits for
#'   mismatch plots. Default `c(0, 3.5)`.
#' @param sim.limits Numeric vector of length 2, y/x-axis limits for
#'   cosine similarity plots (reversed axis). Default `c(1, 0.93)`.
#' @param angle.limits Numeric vector of length 2, y/x-axis limits for
#'   spectral angle plots. Default `c(0, 25)`.
#' @param alignment.limits Numeric vector of length 2, y/x-axis limits for
#'   variability/mismatch alignment plots. Default `c(0, 1)`.
#' @param cytometer Character. Cytometer label used in plot titles and
#'   output filenames.
#'
#' @details
#' Six metrics are compared: `Mismatch`, `Cosine`, `Angle`, `Variability`,
#' `Alignment`, and `Brightness`. For each metric, a violin/jitter plot by
#' dye class is produced. For every pair of metrics, a scatter plot with an
#' `lm()` trendline is produced, annotated with that pair's R-squared and
#' p-value; since simple linear regression's R-squared and F-test p-value
#' are symmetric in the two variables, these statistics are unaffected by
#' which metric in a pair ends up on the plot's x- versus y-axis. `Cosine`
#' is the only metric plotted on a reversed axis (cosine similarity
#' decreases with divergence, unlike the other five metrics, which all
#' increase with divergence).
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter geom_point
#' @importFrom ggplot2 geom_smooth coord_cartesian theme_classic theme
#' @importFrom ggplot2 element_text labs ggsave
#' @importFrom cowplot plot_grid ggdraw draw_label
#'
#' @return A data frame of pairwise linear-fit statistics (R-squared and
#'   p-value) for each pair of metrics, one row per comparison. Returns
#'   `invisible(NULL)` if fewer than 5 fluorophores have complete data.
#'
#' @export

mismatch.plot <- function(
    cosine.data,
    angle.data,
    mismatch.data,
    variability.data,
    alignment.data,
    brightness.data,
    fluor.df,
    particle.name = "UltraComp",
    output.dir = "./results/aurora",
    mismatch.limits = c(0, 3.5),
    sim.limits = c(1, 0.93),
    angle.limits = c(0, 25),
    alignment.limits = c(0, 1),
    cytometer
) {

  if (!dir.exists(output.dir)) dir.create(output.dir)

  metric.data <- list(
    Mismatch    = mismatch.data,
    Cosine      = cosine.data,
    Angle       = angle.data,
    Variability = variability.data,
    Alignment   = alignment.data,
    Brightness  = brightness.data
  )

  common.fluors <- Reduce(
    intersect,
    c(lapply(metric.data, .mismatch.plot.names), list(fluor.df$Fluorophore))
  )

  plot.df <- data.frame(
    Fluorophore = common.fluors,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (m in names(metric.data))
    plot.df[[m]] <- .mismatch.plot.values(metric.data[[m]], common.fluors)

  fluor.df.matched <- fluor.df[fluor.df$Fluorophore %in% common.fluors, , drop = FALSE]

  plot.df <- merge(
    plot.df, fluor.df.matched,
    by = "Fluorophore", all.x = TRUE, sort = FALSE
  )
  plot.df <- stats::na.omit(plot.df)
  plot.df$Class <- factor(plot.df$Class)

  if (nrow(plot.df) < 5) {
    warning(
      "Fewer than 5 fluorophores had complete data across all six ",
      "metrics (mismatch, cosine, angle, variability, alignment, ",
      "brightness) after merging with `fluor.df`; skipping plots for '",
      particle.name, "'.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  metric.spec <- list(
    Mismatch    = list(label = "Mismatch distance",             limits = mismatch.limits,  reverse = FALSE),
    Cosine      = list(label = "Cosine similarity",              limits = sim.limits,       reverse = TRUE),
    Angle       = list(label = "Spectral angle (deg)",           limits = angle.limits,     reverse = FALSE),
    Variability = list(label = "Variability",                    limits = NULL,             reverse = FALSE),
    Alignment   = list(label = "Variability/mismatch alignment", limits = alignment.limits, reverse = FALSE),
    Brightness  = list(label = "Fluorescence intensity",         limits = NULL,             reverse = FALSE)
  )

  metric.names <- names(metric.spec)

  get_stats <- function(fit) {
    s <- summary(fit)
    r2 <- round(s$r.squared, 3)
    f <- s$fstatistic
    p_val <- stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
    p_label <- if (p_val < 0.001) formatC(p_val, format = "e", digits = 2) else round(p_val, 4)
    paste0("R^2 = ", r2, ", p = ", p_label)
  }

  .coord <- function(x.name, y.name) {
    x.spec <- metric.spec[[x.name]]
    y.spec <- metric.spec[[y.name]]
    reverse.axes <- paste0(
      if (isTRUE(x.spec$reverse)) "x" else "",
      if (isTRUE(y.spec$reverse)) "y" else ""
    )
    ggplot2::coord_cartesian(
      xlim    = x.spec$limits,
      ylim    = y.spec$limits,
      reverse = if (nzchar(reverse.axes)) reverse.axes else "none"
    )
  }

  plot.list <- list()

  # -- violin-by-class plots, one per metric ---------------------------------
  for (m in metric.names) {
    spec <- metric.spec[[m]]

    p <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Class, y = .data[[m]], fill = Class)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.6) +
      ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.8, colour = "black") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = paste(particle.name, cytometer), x = NULL, y = spec$label) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    if (!is.null(spec$limits))
      p <- p + ggplot2::coord_cartesian(
        ylim    = spec$limits,
        reverse = if (isTRUE(spec$reverse)) "y" else "none"
      )

    plot.list[[length(plot.list) + 1]] <- p

    ggplot2::ggsave(
      file.path(
        output.dir,
        paste0("Dye class and ", tolower(m), " ", particle.name, " ", cytometer, ".jpg")
      ),
      plot = p, width = 5, height = 5
    )
  }

  # -- pairwise correlation plots, with linear-fit R^2/p-value ---------------
  metric.pairs <- utils::combn(metric.names, 2, simplify = FALSE)

  stats.comparison <- character(length(metric.pairs))
  stats.label       <- character(length(metric.pairs))

  for (i in seq_along(metric.pairs)) {
    x.name <- metric.pairs[[i]][1]
    y.name <- metric.pairs[[i]][2]

    fit   <- stats::lm(plot.df[[y.name]] ~ plot.df[[x.name]])
    label <- get_stats(fit)

    stats.comparison[i] <- paste0(x.name, "_v_", y.name)
    stats.label[i]       <- label

    p <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data[[x.name]], y = .data[[y.name]])) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
      .coord(x.name, y.name) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title    = paste(particle.name, cytometer),
        subtitle = label,
        x        = metric.spec[[x.name]]$label,
        y        = metric.spec[[y.name]]$label
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    plot.list[[length(plot.list) + 1]] <- p

    ggplot2::ggsave(
      file.path(
        output.dir,
        paste0(x.name, " and ", y.name, " ", particle.name, " ", cytometer, ".jpg")
      ),
      plot = p, width = 5, height = 5
    )
  }

  stats.df <- data.frame(
    Comparison  = stats.comparison,
    Stats_Label = stats.label,
    stringsAsFactors = FALSE
  )

  # -- combine into a single multi-panel PDF report --------------------------
  title.plot <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste("Analysis Report:", particle.name, "-", cytometer),
      fontface = "bold",
      size = 14
    )

  report.grid <- cowplot::plot_grid(plotlist = plot.list, ncol = 3)

  report.plot <- cowplot::plot_grid(
    title.plot, report.grid,
    ncol = 1, rel_heights = c(0.05, 1)
  )

  pdf.path <- file.path(output.dir, paste0("Report_", particle.name, "_", cytometer, ".pdf"))

  n.rows <- ceiling(length(plot.list) / 3)

  ggplot2::ggsave(
    pdf.path,
    plot = report.plot,
    width = 8.27,
    height = max(11.69, n.rows * 2.9 + 0.5)
  )

  return(stats.df)
}
