#' @title Handle and Plot Gating Failures
#'
#' @description
#' A helper function that dictates what occurs when a gating error is detected.
#'
#' @importFrom ggplot2 ggplot aes stat_density_2d scale_fill_viridis_c after_stat
#' @importFrom ggplot2 geom_path scale_color_manual guides guide_legend theme_bw
#' @importFrom ggplot2 coord_cartesian labs ggsave
#' @importFrom ragg agg_jpeg
#' @importFrom scattermore geom_scattermore
#'
#' @param e The error object from tryCatch
#' @param gate.id The ID of the gate that failed
#' @param files.to.gate Character vector of filenames used
#' @param scatter.coords The data matrix used for gating
#' @param samp String identifier for the sample/gate type
#' @param viability.gate Logical, is this a viability gate
#' @param control.type "beads" or "cells"
#' @param flow.scatter.and.channel.label Mapping vector for axis labels
#' @param asp The parameter list
#'
#' @return None. Stops flow.
#'
#' @export

handle.gating.error <- function(
    e,
    gate.id,
    files.to.gate,
    scatter.coords,
    samp,
    viability.gate,
    control.type,
    flow.scatter.and.channel.label,
    asp
  ) {

  # Print identifying information
  message("\n\033[31mError during gating for gate ID: ", gate.id, "\033[0m")
  message("Samples involved: ", paste(files.to.gate, collapse = ", "))

  # Identify marker names
  gate.marker <- colnames(scatter.coords)[1:2]
  x.lab <- names(which(flow.scatter.and.channel.label == gate.marker[1]))
  y.lab <- names(which(flow.scatter.and.channel.label == gate.marker[2]))

  # Downsample for plotting efficiency
  plot.coords <- scatter.coords
  if (nrow(plot.coords) > 1e5) {
    set.seed(42)
    plot.coords <- plot.coords[sample(nrow(plot.coords), 1e5), ]
  }

  # Cytometer search boundary
  cytometer.boundary.ggp <- data.frame(
    x = c(asp$scatter.data.min.x, asp$scatter.data.max.x, asp$scatter.data.max.x, asp$scatter.data.min.x, asp$scatter.data.min.x),
    y = c(asp$scatter.data.min.y, asp$scatter.data.min.y, asp$scatter.data.max.y, asp$scatter.data.max.y, asp$scatter.data.min.y)
  )

  # Select trim factors
  if (control.type == "beads") {
    tx.min <- asp$gate.data.trim.factor.x.min.beads
    tx.max <- asp$gate.data.trim.factor.x.max.beads
    ty.min <- asp$gate.data.trim.factor.y.min.beads
    ty.max <- asp$gate.data.trim.factor.y.max.beads
  } else {
    tx.min <- asp$gate.data.trim.factor.x.min.cells
    tx.max <- asp$gate.data.trim.factor.x.max.cells
    ty.min <- asp$gate.data.trim.factor.y.min.cells
    ty.max <- asp$gate.data.trim.factor.y.max.cells
  }

  # Calculate Search Region
  d.x.min <- min(plot.coords[, 1]); d.x.max <- max(plot.coords[, 1])
  d.y.min <- min(plot.coords[, 2]); d.y.max <- max(plot.coords[, 2])

  r.x.low  <- pmax(asp$scatter.data.min.x, (1 - tx.min) * d.x.min + tx.min * d.x.max)
  r.x.high <- pmin(asp$scatter.data.max.x, (1 - tx.max) * d.x.min + tx.max * d.x.max)
  r.y.low  <- pmax(asp$scatter.data.min.y, (1 - ty.min) * d.y.min + ty.min * d.y.max)
  r.y.high <- pmin(asp$scatter.data.max.y, (1 - ty.max) * d.y.min + ty.max * d.y.max)

  if (viability.gate) {
    r.x.low <- pmax(asp$scatter.data.min.x, r.x.low * asp$viability.gate.scaling)
  }

  gate.region.ggp <- data.frame(
    x = c(r.x.low, r.x.high, r.x.high, r.x.low, r.x.low),
    y = c(r.y.low, r.y.low, r.y.high, r.y.high, r.y.low)
  )

  plot.df <- as.data.frame(plot.coords)
  colnames(plot.df) <- c("x","y")

  # Generate failure plot
  diag.plot <- ggplot(plot.df, aes(x = x, y = y)) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      alpha = 0.4,
      color = "black"
    ) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      n = 60,
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(option = "mako") +
    geom_path(
      data = gate.region.ggp,
      aes(x, y, color = "Search region"),
      linewidth = asp$figure.gate.line.size,
      linetype = "dotted"
    ) +
    geom_path(
      data = cytometer.boundary.ggp,
      aes(x, y, color = "Scatter limits"),
      linewidth = asp$figure.gate.line.size
    ) +
    scale_color_manual(
      name = "Diagnostics",
      values = c("Search region" = "black", "Scatter limits" = "red")
    ) +
    guides(
      fill = "none",
      color = guide_legend(override.aes = list(fill = NA, linetype = c(1, 3), linewidth = 0.5))
    ) +
    coord_cartesian(
      xlim = range(c(plot.coords[,1], asp$scatter.data.min.x, asp$scatter.data.max.x)),
      ylim = range(c(plot.coords[,2], asp$scatter.data.min.y, asp$scatter.data.max.y))
    ) +
    theme_bw() + theme(legend.position = "right") +
    labs(title = paste("Gating Failure:", samp), x = x.lab, y = y.lab)

  fail.file <- sprintf("FAILURE_%s.jpg", samp)
  ggsave(
    file.path( asp$figure.gate.dir, fail.file ),
    plot = diag.plot,
    device = ragg::agg_jpeg,
    width = 6,
    height = 5
  )
  # Stop with original error
  stop("Gating failed. Diagnostic plot saved as: ", fail.file, "\nOriginal Error: ", e$message)
}
