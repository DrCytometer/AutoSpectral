#' @title Plot Spectral Mismatch, Variability, and Brightness by Dye Class
#'
#' @description
#' Produces a standard set of diagnostic plots comparing a particle type's
#' spectral mismatch, variability, cosine similarity, and brightness
#' against fluorophore dye class, plus pairwise correlations between all
#' four metrics with simple linear-fit statistics. Each individual plot is
#' saved as a JPEG in `output.dir`, and all plots are additionally
#' combined into a single multi-panel PDF report.
#'
#' @param cosine.data One-column numeric matrix of cosine similarity
#'   values, rownames are fluorophore names, as returned by
#'   `assess.mismatch()`.
#' @param mismatch.data Named numeric vector of per-fluorophore mismatch
#'   magnitudes, e.g. `rowSums(abs(bead.cell.dist(...)))`.
#' @param variability.data Named numeric vector of per-fluorophore
#'   variability magnitudes, e.g.
#'   `rowSums(abs(assess.variability.mad(...)))`.
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
#' @param cytometer Character. Cytometer label used in plot titles and
#'   output filenames.
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
    mismatch.data,
    variability.data,
    brightness.data,
    fluor.df,
    particle.name = "UltraComp",
    output.dir = "./results/aurora",
    mismatch.limits = c(0, 3.5),
    sim.limits = c(1, 0.93),
    cytometer
) {

  if(!dir.exists(output.dir)) dir.create(output.dir)

  common.fluors <- Reduce(intersect, list(
    rownames(cosine.data),
    names(mismatch.data),
    names(variability.data),
    rownames(brightness.data),
    fluor.df$Fluorophore
  ))

  plot.df <- data.frame(
    Fluorophore = common.fluors,
    Cosine      = as.numeric(cosine.data[common.fluors,]),
    Mismatch    = as.numeric(mismatch.data[common.fluors]),
    Variability = as.numeric(variability.data[common.fluors]),
    Brightness  = as.numeric(brightness.data[common.fluors,]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  fluor.df.matched <- fluor.df[fluor.df$Fluorophore %in% common.fluors, , drop = FALSE]

  plot.df <- merge(
    plot.df, fluor.df.matched,
    by = "Fluorophore", all.x = TRUE, sort = FALSE
  )
  plot.df <- stats::na.omit(plot.df)
  plot.df$Class <- factor(plot.df$Class)

  if(nrow(plot.df) < 5){
    warning(
      "Fewer than 5 fluorophores had complete data across all four ",
      "metrics (cosine, mismatch, variability, brightness) after merging ",
      "with `fluor.df`; skipping plots for '", particle.name, "'.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  ### stats
  get_stats <- function(fit) {
    s <- summary(fit)
    r2 <- round(s$r.squared, 3)
    # Extract p-value from F-statistic
    f <- s$fstatistic
    p_val <- stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
    p_label <- if(p_val < 0.001) formatC(p_val, format = "e", digits = 2) else round(p_val, 4)
    return(paste0("R^2 = ", r2, ", p = ", p_label))
  }

  # Run models
  m_cos  <- stats::lm(Mismatch ~ Cosine, data = plot.df)
  m_var  <- stats::lm(Mismatch ~ Variability, data = plot.df)
  m_bri  <- stats::lm(Mismatch ~ Brightness, data = plot.df)
  m_c_v  <- stats::lm(Cosine ~ Variability, data = plot.df)
  m_c_b  <- stats::lm(Cosine ~ Brightness, data = plot.df)
  v_b  <- stats::lm(Variability ~ Brightness, data = plot.df)

  # Create stats dataframe for return
  stats.df <- data.frame(
    Comparison = c("Cos_v_Miss", "Var_v_Miss", "Bri_v_Miss", "Var_v_Cos", "Bri_v_Cos", "Bri_v_Var"),
    Stats_Label = c(get_stats(m_cos), get_stats(m_var), get_stats(m_bri),
                    get_stats(m_c_v), get_stats(m_c_b), get_stats(v_b))
  )

  ### plotting
  plot.list <- list()

  # plot mismatch by class
  plot.list[[1]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Class, y = Mismatch, fill = Class)) +
    ggplot2::geom_violin(
      trim = FALSE,
      alpha = 0.6
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      colour = "black"
    ) +
    ggplot2::coord_cartesian(ylim = mismatch.limits) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      x = NULL,
      y = "Mismatch distance"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Dye class and mismatch ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[1]],
    width = 5,
    height = 5
  )

  # plot variability by dye class
  plot.list[[2]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Class, y = Variability, fill = Class)) +
    ggplot2::geom_violin(
      trim = FALSE,
      alpha = 0.6
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      colour = "black"
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      x = NULL,
      y = "Variability"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Dye class and variability ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[2]],
    width = 5,
    height = 5
  )

  # plot cosine sim by dye class
  plot.list[[3]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Class, y = Cosine, fill = Class)) +
    ggplot2::geom_violin(
      trim = FALSE,
      alpha = 0.6
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      colour = "black"
    ) +
    ggplot2::coord_cartesian(ylim = sim.limits, reverse = "y") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      x = NULL,
      y = "Cosine distance"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Dye class and cosine similarity ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[3]],
    width = 5,
    height = 5
  )

  # plot brightness by dye class
  plot.list[[4]] <-  ggplot2::ggplot(plot.df, ggplot2::aes(x = Class, y = Brightness, fill = Class)) +
    ggplot2::geom_violin(
      trim = FALSE,
      alpha = 0.6
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      colour = "black"
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      x = NULL,
      y = "Fluorescence intensity"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Dye class and brightness ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[4]],
    width = 5,
    height = 5
  )

  ### plot correlations with trendline

  # plot correlation between cosine and mismatch
  plot.list[[5]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Cosine, y = Mismatch)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggplot2::coord_cartesian(xlim = sim.limits, ylim = mismatch.limits, reverse = "x") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      subtitle = stats.df$Stats_Label[1],
      x = "Cosine similarity",
      y = "Mismatch distance"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Cosine and mismatch ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[5]],
    width = 5,
    height = 5
  )

  # plot correlation between mismatch and variability
  plot.list[[6]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Variability, y = Mismatch)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggplot2::coord_cartesian(ylim = mismatch.limits) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      subtitle = stats.df$Stats_Label[2],
      x = "Variability",
      y = "Mismatch distance"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Variability and mismatch ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[6]],
    width = 5,
    height = 5
  )

  # plot correlation between mismatch and brightness
  plot.list[[7]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Brightness, y = Mismatch)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggplot2::coord_cartesian(ylim = mismatch.limits) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      subtitle = stats.df$Stats_Label[3],
      x = "Fluorescence intensity",
      y = "Mismatch distance"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Brightness and mismatch ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[7]],
    width = 5,
    height = 5
  )

  # plot correlation between cosine and variability
  plot.list[[8]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Variability, y = Cosine)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggplot2::coord_cartesian(ylim = sim.limits, reverse = "y") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      subtitle = stats.df$Stats_Label[4],
      x = "Variability",
      y = "Cosine similarity"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Variability and cosine ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[8]],
    width = 5,
    height = 5
  )

  # plot correlation between cosine and brightness
  plot.list[[9]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Brightness, y = Cosine)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggplot2::coord_cartesian(ylim = sim.limits, reverse = "y") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      subtitle = stats.df$Stats_Label[5],
      x = "Fluorescence intensity",
      y = "Cosine similarity"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Brightness and cosine ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[9]],
    width = 5,
    height = 5
  )

  # plot correlation between variability and brightness
  plot.list[[10]] <- ggplot2::ggplot(plot.df, ggplot2::aes(x = Brightness, y = Variability)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", color = "blue", se = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste(particle.name, cytometer),
      subtitle = stats.df$Stats_Label[6],
      x = "Fluorescence intensity",
      y = "Variability"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(output.dir,
              paste0("Brightness and variability ", particle.name, " ", cytometer, ".jpg")),
    plot = plot.list[[10]],
    width = 5,
    height = 5
  )

  # combine into a single multi-panel PDF report
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

  ggplot2::ggsave(
    pdf.path,
    plot = report.plot,
    width = 8.27,
    height = 11.69
  )

  return(stats.df)
}
