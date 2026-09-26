# =============================================================================
# Figure 4, panels lettered in the order the Results cite them (CAMPER, then
# KEGG); the layout places d beside a to save space, so letters do not follow
# reading order.
#
#   a  CAMPER capacity by oxygen requirement and compound class (camper_n76.R)
#   b  CAMPER modules differing between groups (BH q < 0.05), fill = Cliff's delta,
#      one block per compound class
#   c  KEGG pathway ORA, all 48 pathways at BH q < 0.05, fill = signed -log10(FDR)
#   d  KEGG module ORA, fill = signed -log10(FDR)
#
# Every tile in b, c and d is the same physical square (TILE inches), fixed
# with ggh4x::force_panelsizes, so they line up when placed side by side in
# Affinity. The earlier versions let the tile stretch to the panel width, which
# drew wide slabs instead of boxes.
#
# Reads the tables written by camper_n76.R and kegg_summary_n76.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(forcats); library(ggplot2)
  library(ggh4x); library(patchwork)
})

theme_transparent <- ggplot2::theme(
  plot.background       = element_rect(fill = "transparent", color = NA),
  panel.background      = element_rect(fill = "transparent", color = NA),
  legend.background     = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA),
  legend.key            = element_rect(fill = "transparent", color = NA),
  strip.background      = element_rect(fill = "transparent", color = NA))
ggplot2::theme_set(ggplot2::theme_get() + theme_transparent)
ggsave <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

BASE <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
source(file.path(BASE, "taxon_italics_n76.R"))   # md_taxon_in(): taxon names in italics
F4   <- file.path(BASE, "result2/n76/fig4")
TILE <- 0.17   # inches, one square per row, identical across b, c and d

shorten <- function(x, n) ifelse(nchar(x) <= n, x,
  paste0(sub("\\s+\\S*$", "", substr(x, 1, n)), "…"))

# one square tile per row; optional grouping into facets drawn as left strips
square_1d <- function(df, label, value, title, subtitle, legend_name,
                      group = NULL, cap = NULL, lim = NULL, legend = TRUE, pad_to = NULL,
                      pal = c(low = "#0f4480", high = "#b35806")) {
  df <- df %>% mutate(.v = .data[[value]], .lab = .data[[label]])
  if (is.null(lim)) lim <- if (is.null(cap)) max(abs(df$.v)) else cap
  pad <- character(0)
  if (is.null(group)) {
    df <- df %>% mutate(.lab = fct_reorder(.lab, .v))
    rows <- nrow(df)
    # blank rows under the real ones, so blocks sharing a grid row are all the
    # same height and their titles line up instead of floating mid-cell
    if (!is.null(pad_to) && pad_to > rows) {
      pad <- strrep("\u200b", seq_len(pad_to - rows))
      rows <- pad_to
    }
  } else {
    df <- df %>% mutate(.g = .data[[group]]) %>%
      arrange(.g, .v) %>% mutate(.lab = factor(.lab, unique(.lab)))
    rows <- as.integer(table(factor(df$.g, unique(df$.g))))
  }
  p <- ggplot(df, aes(x = "", y = .lab, fill = .v)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_x_discrete(expand = c(0, 0)) +
    # plotmath, not ggtext: force_panelsizes below rebuilds the axis and drops
    # a markdown grob, printing the asterisks
    scale_y_discrete(expand = c(0, 0), drop = FALSE,
                     labels = function(z) expr_from_md(md_taxon_in(z)),
                     limits = c(pad, levels(droplevels(factor(df$.lab))))) +
    scale_fill_gradient2(low = pal[["low"]], mid = "grey98", high = pal[["high"]], midpoint = 0,
                         limits = c(-lim, lim), oob = scales::squish,
                         name = if (is.null(cap)) legend_name else
                           sprintf("%s\ncapped at %g", legend_name, cap)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = 6.5, colour = "black"),
          panel.grid = element_blank(), panel.spacing.y = unit(0.08, "cm"),
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6),
          legend.key.height = unit(0.45, "cm"), legend.key.width = unit(0.28, "cm"),
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 6.5, colour = "grey35"),
          plot.title.position = "plot")
  if (!legend) p <- p + theme(legend.position = "none")
  if (is.null(group)) {
    p <- p + force_panelsizes(rows = unit(rows * TILE, "in"), cols = unit(TILE, "in"))
  } else {
    p <- p + facet_grid(rows = vars(.g), scales = "free_y", space = "free_y", switch = "y") +
      theme(strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, size = 6.5, face = "bold", hjust = 1)) +
      force_panelsizes(rows = unit(rows * TILE, "in"), cols = unit(TILE, "in"))
  }
  attr(p, "n_rows") <- sum(rows)
  p
}

# c and d show the same quantity, so they share one colour scale and one legend
# (drawn on c). Capped at 5: modules reach 3.9, 34 of 48 pathways sit below 5,
# and a cap of 10 left every module tile near white.
CAP_FDR <- 5

# ---- d. KEGG modules --------------------------------------------------------
km <- read_csv(file.path(F4, "Fig4d_KEGG_module_ORA_n76.csv"), show_col_types = FALSE) %>%
  filter(q < 0.05)
pMod <- square_1d(km, "short", "signed",
                  title = "KEGG modules",
                  subtitle = sprintf("%d modules, BH q < 0.05. Colour scale as in c", nrow(km)),
                  legend_name = "Signed\n-log10(FDR)", lim = CAP_FDR, legend = FALSE)

# ---- c. KEGG pathway ORA, every pathway at BH q < 0.05 -----------------------
# The full pathway-level result, as in the published paper, rather than the
# category summary it briefly replaced (too coarse to read). Each pathway takes
# its stronger direction; the sign carries it, as on the published butterfly.
# The raw-p < 0.05 list (80 pathways) remains in Supplementary Fig. 6.
kp <- read_csv(file.path(F4, "kegg_fast/linda_results/KEGG_KO_hypergeometric_ORA_n76.csv"),
               show_col_types = FALSE) %>%
  mutate(dir = ifelse(q_case <= q_control, "Case", "Control"),
         q = pmin(q_case, q_control),
         signed = ifelse(dir == "Case", 1, -1) * -log10(q),
         lab = shorten(pathway_name, 50)) %>%
  filter(q < 0.05)
write_csv(kp, file.path(F4, "Fig4c_KEGG_pathway_ORA_q05_n76.csv"))
pPath <- square_1d(kp, "lab", "signed",
                   title = "KEGG pathway over-representation",
                   subtitle = sprintf("%d pathways, BH q < 0.05. Orange: Case, blue: Control", nrow(kp)),
                   legend_name = "Signed\n-log10(FDR)", cap = CAP_FDR)

# ---- b. CAMPER modules, one block per compound class, side by side ---------
# Each class is its own short column of square tiles; blocks sit in a 3 x 2 grid,
# the three large classes on the first row and the three small ones on the
# second so the rows balance. All blocks share one colour scale (lim_b) so a
# given shade means the same Cliff's delta everywhere; the legend is drawn once.
cm <- read_csv(file.path(F4, "camper_module_tests_n76.csv"), show_col_types = FALSE) %>%
  filter(q < 0.05) %>%
  mutate(cls = gsub("_", " ", cls), lab = shorten(module, 46))
lim_b <- max(abs(cm$delta))
cls_order <- cm %>% count(cls, sort = TRUE) %>% pull(cls)
n_in <- as.integer(table(cm$cls)[cls_order])
row_max <- c(rep(max(n_in[1:3]), 3), rep(max(n_in[4:6]), 3))
# Cliff's delta (an effect size) gets its own red-purple scale; orange-blue is
# kept for the signed -log10(FDR) of c and d, so equal colours never imply the
# same quantity. Panel a's labels report the same Cliff's delta.
PAL_DELTA <- c(low = "#542788", high = "#b2182b")
# The shared legend sits under the last (smallest) block, horizontal, so it adds
# height where there is spare room instead of width that would shift one column.
blocks <- lapply(seq_along(cls_order), function(i) {
  k <- cls_order[i]
  last <- i == length(cls_order)
  p <- square_1d(cm %>% filter(cls == k), "lab", "delta",
                 title = sprintf("%s (%d)", k, sum(cm$cls == k)), subtitle = NULL,
                 legend_name = "Cliff's delta", lim = lim_b, legend = last,
                 pad_to = if (last) NULL else row_max[i], pal = PAL_DELTA)
  if (last) p <- p + theme(legend.position = "bottom", legend.justification = "right",
                           legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.22, "cm"),
                           legend.title = element_text(size = 6.5, face = "bold", vjust = 0.8))
  p
})
blk_rows <- sapply(blocks, function(p) attr(p, "n_rows"))
row_h <- c(max(blk_rows[1:3]), max(blk_rows[4:6])) * TILE + 0.45
# A fixed-size grob is centred in its cell by default, so it floats with its
# label length and row count. pin() fixes it to a top corner instead: CAMPER
# blocks go top-right, so the two tile columns in each grid column sit on one
# vertical line; c and d go top-left, so their titles start level with a's.
pin <- function(p, side = c("left", "right")) {
  side <- match.arg(side)
  g <- ggplotGrob(p)
  g$vp <- grid::viewport(x = as.numeric(side == "right"), y = 1,
                         width = sum(g$widths), height = sum(g$heights),
                         just = c(side, "top"))
  wrap_elements(full = g)
}
pCam <- wrap_plots(lapply(blocks, pin, side = "right"), ncol = 3) +
  plot_layout(heights = row_h) +
  plot_annotation(title = "CAMPER modules differing between groups",
                  subtitle = "BH q < 0.05, grouped by compound class. Red: higher capacity in Case, purple: lower",
                  theme = theme(plot.title = element_text(face = "bold", size = 9),
                                plot.subtitle = element_text(size = 6.5, colour = "grey35"),
                                plot.margin = margin(4, 4, 4, 16)))
H_CAM <- sum(row_h) + 0.6

# ---- standalone panels ------------------------------------------------------
# canvas height = rows of tiles + room for title and margins; width fits labels
h_of <- function(p) attr(p, "n_rows") * TILE + 1.3
save_sq <- function(p, file, w) {
  ggsave(file.path(F4, file), p, width = w, height = max(h_of(p), 2.4),
         device = cairo_pdf, bg = "transparent", limitsize = FALSE)
}
save_sq(pMod, "Fig4d_KEGG_modules_n76.pdf", 3.8)
save_sq(pPath, "Fig4c_KEGG_pathways_n76.pdf", 4.0)
ggsave(file.path(F4, "Fig4b_CAMPER_modules_n76.pdf"), pCam, width = 8.3, height = H_CAM,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)
saveRDS(list(panel_b = pCam, panel_c = pPath, panel_d = pMod),
        file.path(F4, "Fig4bcd_square_panels_n76.rds"))
message(sprintf("b: %d CAMPER modules in %d blocks (%.1f in) | c %d pathways | d %d modules  (tile %.2f in)",
                nrow(cm), length(blocks), H_CAM, attr(pPath, "n_rows"), attr(pMod, "n_rows"), TILE))

# ---- Figure 4 ----------------------------------------------------------------
camper <- readRDS(file.path(F4, "Fig4ab_CAMPER_merged_panels_n76.rds"))
tag <- theme(plot.tag = element_text(face = "bold", size = 13))
# ggh4x fixes panel sizes in absolute units, which patchwork's panel alignment
# cannot handle (it fails in set_panel_dimensions). Wrapping hands those panels
# over as finished grobs, so the squares survive and nothing is realigned. The
# left margin keeps the wrapped panel's title clear of its tag.
sq <- function(p, t) pin(p + labs(tag = t) + tag + theme(plot.margin = margin(4, 4, 4, 16)))

# Panel a as saved by camper_n76.R is built for a full-width slot: boxed
# horizontal facet strips (~1.3 in of grey), a long title and 10 pt text.
# Compacted here, not in camper_n76.R, so the standalone CAMPER figure is
# unchanged. The two sets of boxplots are marked by a rotated label and a thin
# bracket line (ggh4x part-rect: right edge only) instead of a filled box.
pA <- camper$panel_a +
  facet_grid2(block ~ ., scales = "free_y", space = "free_y", switch = "y",
              labeller = labeller(block = c(`Oxygen requirement` = "Oxygen",
                                            `Compound class` = "Compound class"))) +
  labs(title = "Polyphenol-transformation capacity",
       subtitle = "Per donor: MAG abundance x module completeness, summed. Sqrt axis.") +
  theme(strip.background.y = element_part_rect(side = "r", fill = NA, colour = "grey30",
                                               linewidth = 0.5),
        strip.text.y.left = element_text(angle = 90, size = 7, face = "bold",
                                         margin = margin(0, 3, 0, 0)),
        strip.switch.pad.grid = unit(3, "pt"),
        axis.text = element_text(size = 7), axis.title.x = element_text(size = 7.5),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        plot.title.position = "plot", legend.text = element_text(size = 7),
        legend.margin = margin(0, 0, 0, 0))
pA$layers[[3]]$aes_params$size <- 2.1   # q / delta labels

# Half-page layout. Top row: a beside d (KEGG modules); b (CAMPER modules)
# underneath; the tall pathway list c as a right-hand column spanning both rows.
# Wrapped fixed-size panels cannot shrink, so every cell is sized from the real
# panel dimensions: tiles * TILE + ~1.3 in for title and margins.
W_A <- 4.3; W_B <- 3.8; W_C <- 4.0
H_TOP <- max(4.6, h_of(pMod))
H_L   <- H_TOP + H_CAM
top   <- (pA + labs(tag = "a") + tag) | sq(pMod, "d")
top   <- top + plot_layout(widths = c(W_A, W_B))
left  <- top / (wrap_elements(full = pCam) + labs(tag = "b") + tag) +
  plot_layout(heights = c(H_TOP, H_CAM))
right <- sq(pPath, "c") / plot_spacer() +
  plot_layout(heights = c(h_of(pPath), max(H_L - h_of(pPath), 0.01)))
# Panel e is the trait UpSet from fig4e_trait_upset_n76.R. The first UpSet
# counted genomes equally and so gave its largest bar to the one trait carrying
# no signal; the boxplot that replaced it fixed the weighting but dropped the
# co-occurrence structure. This one keeps the intersections and adds both
# missing dimensions, per-donor abundance share and ANCOM-BC2 direction.
# The marginal per-trait version is still built by fig4e_trait_capacity_n76.R
# and kept in this directory for the supplement.
# Given the left two-thirds of a bottom row, near its native 7.2 in width,
# rather than stretched across the full page.
# Panel f is the predicted fermentation-product formation capacity from
# fig4f_flux_capacity_n76.R, which answers the question panel e cannot: whether
# the network those genes sit in can form the compound at all on a defined diet.
# It takes the bottom-right slot that panel e leaves empty.
# Panel f is now the two-part capacity-against-use figure from
# fig4f_capacity_vs_use_n76.R, a patchwork rather than a single ggplot, so it is
# wrapped like panel e. It carries two sub-panels and needs more of the bottom
# row than the single capacity boxplot it replaces, hence the wider f slot below.
H_E <- 5.4
pE <- wrap_elements(full = readRDS(file.path(F4, "Fig4e_TraitUpSet_n76.rds")))
pF <- wrap_elements(full = readRDS(file.path(F4, "Fig4f_CapacityVsUse_n76.rds")))
H_BODY <- max(H_L, h_of(pPath))
W_E <- 6.9; W_F <- 5.2          # bottom row, f needs the extra width
FIG_W <- W_A + W_B + W_C + 0.4
FIG_H <- H_BODY + H_E + 0.7
body <- (left | right) + plot_layout(widths = c(W_A + W_B, W_C))
# f's own sub-panel titles start at its left edge, so the default tag slot sits
# on top of "(i)". Pin the tag to the far corner of the slot instead.
bottom <- ((pE + labs(tag = "e") + tag) |
           (pF + labs(tag = "f") + tag +
              theme(plot.tag.position = c(0.005, 0.995)))) +
  plot_layout(widths = c(W_E, W_F))
fig4 <- (body / bottom) + plot_layout(heights = c(H_BODY, H_E)) +
  plot_annotation(
    title = "Figure 4. Polyphenol and KEGG functional shifts, gut trait combinations, and predicted fermentation capacity against predicted use (76 donors, 584 MAGs)",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(F4, "Figure_4_n76.pdf"), fig4, width = FIG_W, height = FIG_H,
       units = "in", device = cairo_pdf, bg = "transparent", limitsize = FALSE)
message(sprintf("DONE -> Figure_4_n76.pdf (%.1f x %.1f in)", FIG_W, FIG_H))
