# Layout-only companion to figure1_n76.R. No statistical fitting.
# Standalone panels retain their requested dimensions. This compositor is
# drawn at its actual DOCX width so labels are not shrunk after insertion.
compose_figure1_n76 <- function() {
  pt_mm <- 72.27 / 25.4
  small_theme <- theme(
    text = element_text(family = "Helvetica", size = 5.5, colour = "black"),
    plot.title = element_text(size = 7, face = "bold", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 5, margin = margin(b = 2)),
    plot.caption = element_text(size = 5, hjust = 0),
    plot.tag = element_text(size = 7, face = "bold"),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 5.5, margin = margin(1, 1, 1, 1)),
    axis.text = element_text(size = 5, colour = "black"),
    legend.title = element_text(size = 5.5), legend.text = element_text(size = 5),
    legend.key.size = unit(2.1, "mm"), legend.spacing = unit(0.4, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"),
    plot.margin = margin(2, 2, 2, 2), panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15),
    panel.border = element_rect(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.2), axis.ticks.length = unit(0.6, "mm"))
  compact <- function(p, title = NULL, tag = NULL) {
    # Clone layers: mutations must not change the standalone plots.
    p$layers <- lapply(p$layers, function(layer) ggproto(NULL, layer))
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      if (inherits(layer$geom, "GeomText") || inherits(layer$geom, "GeomTextRepel")) {
        layer$aes_params$size <- 5 / pt_mm
        if (!identical(layer$aes_params$family, "DejaVu Sans"))
          layer$aes_params$family <- "Helvetica"
        if (inherits(layer$geom, "GeomTextRepel")) {
          layer$geom_params$box.padding <- 0.15
          layer$geom_params$point.padding <- 0.1
        }
      }
      if (!is.null(layer$aes_params$linewidth))
        layer$aes_params$linewidth <- min(layer$aes_params$linewidth, 0.35)
      if (inherits(layer$geom, "GeomPoint") && !is.null(layer$aes_params$size))
        layer$aes_params$size <- min(layer$aes_params$size, 1.25)
      p$layers[[i]] <- layer
    }
    p + small_theme + labs(title = title, subtitle = NULL, caption = NULL, tag = tag)
  }

  a <- compact(pA1, "Community composition", "a") +
    labs(subtitle = sprintf("PERMANOVA R² = %.3f; p = %.3f", r2_val, p_val)) +
    scale_starshape_manual(values = COHORT_SHAPE, name = NULL,
                           labels = c("Cases", "PREDICT", "Cardiff", "Sanger")) +
    scale_size_continuous(range = c(0.9, 2.2), breaks = RICH_BRK, guide = "none") +
    scale_fill_distiller(palette = "Blues", name = "Shannon", direction = 1,
                         limits = SHANNON_LIM, breaks = c(0, 2, 4)) +
    guides(starshape = guide_legend(nrow = 1, order = 1,
                override.aes = list(size = 1.4, fill = "grey55")),
           size = "none", fill = guide_colourbar(order = 2, barwidth = unit(14, "mm"),
                                  barheight = unit(1.6, "mm"))) +
    theme(legend.box = "vertical", legend.position = "bottom")
  bm <- compact(b_main, NULL, NULL) +
    labs(x = "Shannon diversity", y = "Distance to spatial median") +
    theme(aspect.ratio = NULL, plot.margin = margin(0, 0, 0, 0))
  bt <- b_top + theme_void(base_size = 5, base_family = "Helvetica") +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  br <- b_right + theme_void(base_size = 5, base_family = "Helvetica") +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  b <- ggplotify::as.ggplot(bm |> aplot::insert_top(bt, height = 0.17) |>
                            aplot::insert_right(br, width = 0.17)) +
    theme_void(base_size = 5, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 7, face = "bold"),
          plot.subtitle = element_text(size = 5), plot.tag = element_text(size = 7, face = "bold"),
          plot.tag.position = c(0, 1),
          plot.margin = margin(2, 2, 2, 2)) +
    labs(title = "Diversity and dispersion", tag = "b",
                      subtitle = sprintf("Interaction p = %.2f", m_anova["Shannon:Group", "Pr(>F)"]))
  z <- compact(pZ, "Microbiome health", "c") +
    labs(y = "ZOE health index", subtitle = sprintf("p < 0.001; delta = %.2f", cd_new)) +
    theme(plot.margin = margin(2, 3, 2, 2))
  z$layers <- Filter(function(layer) !inherits(layer$geom, "GeomText") &&
                      !inherits(layer$geom, "GeomSegment"), z$layers)

  d <- compact(pC, "Differential abundance", "d") +
    labs(x = "MAG log fold change", y = "-log10(q)") +
    theme(legend.position = "none") + scale_size_continuous(range = c(0.4, 1.8))
  # Selected taxon labels remain in the full-size standalone volcano. At this
  # final size an unlabelled point panel preserves every MAG without collisions.
  d$layers <- Filter(function(layer) !inherits(layer$geom, "GeomTextRepel"), d$layers)
  e <- compact(pD, "DA MAGs", "e") +
    theme_void(base_size = 5, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 7, face = "bold"),
          plot.tag = element_text(size = 7, face = "bold"),
          plot.tag.position = c(0, 1),
          legend.title = element_text(size = 5.5), legend.text = element_text(size = 5),
          legend.position = "right", legend.key.size = unit(2, "mm"),
          plot.margin = margin(2, 2, 2, 2))
  f <- compact(pE, "Abundance and prevalence", "f") +
    labs(x = "MAG log fold change", y = "Prevalence difference") +
    scale_size_continuous(range = c(0.5, 2))
  g <- compact(pF, "Phylum responses", "g") +
    labs(x = "MAG log fold change") + theme(axis.text.y = element_text(size = 5))
  # Per-row full count/q strings are retained in the 4 x 6 inch standalone,
  # but cannot share a 2 mm-high row with points in the final-sized composite.
  g$layers <- Filter(function(layer) !inherits(layer$geom, "GeomText"), g$layers)

  # One strip across all 79 groups, as in the standalone panel.
  family_order <- clade_figs$family$clade
  family_strip <- function(levels, first) {
    p <- compact(pG, if (first) "Family responses" else NULL, if (first) "h" else NULL)
    p$data <- p$data[as.character(p$data$Family) %in% levels, , drop = FALSE]
    if (!first) p$layers <- Filter(function(layer) !inherits(layer$geom, "GeomText"), p$layers)
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      if (inherits(layer$geom, "GeomText") &&
          startsWith(as.character(layer$aes_params$label %||% ""), "n = "))
        layer$data$y <- max(clade_figs$data$lfc) - 2.6
      if (is.data.frame(layer$data) && "Family" %in% names(layer$data)) {
        keep <- as.character(layer$data$Family) %in% levels
        if (length(layer$aes_params$colour) == nrow(layer$data))
          layer$aes_params$colour <- layer$aes_params$colour[keep]
        layer$data <- layer$data[keep, , drop = FALSE]
      }
      p$layers[[i]] <- layer
    }
    p + scale_x_discrete(limits = levels, labels = function(x)
      as.expression(lapply(x, function(z) if (z == "Unassigned") z else bquote(italic(.(z))))),
      expand = expansion(add = 0.65)) +
      scale_size_continuous(range = c(0.35, 1.4), limits = c(0, 1), guide = "none") +
      scale_y_continuous(limits = c(min(clade_figs$data$lfc) - 0.85,
                                    max(clade_figs$data$lfc) + 1.3),
                         breaks = c(-5, 0, 5), expand = c(0, 0)) +
      labs(y = "MAG LFC", subtitle = NULL) +
      theme(legend.position = "none", axis.text.x = element_text(size = 5, angle = 90,
              hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank())
  }
  h <- family_strip(family_order, TRUE)
  stopifnot(nrow(h$data) == 584L, !anyDuplicated(h$data$taxon),
            setequal(h$data$taxon, clade_figs$data$taxon))
  row1 <- a + b + z + plot_layout(widths = c(1.28, 1.05, 0.75))
  row2 <- d + e + plot_layout(widths = c(1.65, 1))
  row3 <- f + g + plot_layout(widths = c(0.95, 1.3))
  (row1 / row2 / row3 / h) +
    plot_layout(heights = c(50, 35, 55, 62))
}
