# =============================================================================
# Supplementary Figure 1, 76-donor catalogue
#
#   a  QC: MAG read count vs genome coverage fraction, ghost / rare / dominant
#   b  per-donor phylum composition, log-compressed, faceted by group
#
# Coverage thresholds 10% / 75%; log-compression constant k = 100.
#
# Outputs -> result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(patchwork); library(viridis)
  library(ggExtra); library(colorspace)
})

# ---- transparent canvas for Affinity ---------------------------------------
# Matches theme_transparent from Mags_annotation_Fig3.R. theme_bw/minimal/void
# are masked so every existing theme call inherits it without being edited, and
# ggsave defaults to a transparent background. theme_set covers the patchwork
# composite, whose canvas comes from the default theme rather than any panel.
theme_transparent <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA))
theme_bw      <- function(...) ggplot2::theme_bw(...)      + theme_transparent
theme_minimal <- function(...) ggplot2::theme_minimal(...) + theme_transparent
theme_void    <- function(...) ggplot2::theme_void(...)    + theme_transparent
theme_classic <- function(...) ggplot2::theme_classic(...) + theme_transparent
ggplot2::theme_set(ggplot2::theme_get() + theme_transparent)
ggsave <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "supp_fig")
source(file.path(BASE, "taxon_italics_n76.R"))   # md_taxon(): taxon names in italics
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

# Typography for a figure that has to read at A4. geom_text/annotate sizes are
# in mm: 1.8 mm is about 5 pt, which is the floor below which the panels stop
# being readable in print.
PT5 <- 1.8
small_text <- ggplot2::theme(
  axis.text    = ggplot2::element_text(size = 5.5, colour = "black"),
  axis.title   = ggplot2::element_text(size = 6.5),
  legend.text  = ggplot2::element_text(size = 5.5),
  legend.title = ggplot2::element_text(size = 6.5),
  legend.key.size = grid::unit(0.30, "cm"),
  plot.title    = ggplot2::element_text(face = "bold", size = 8),
  plot.subtitle = ggplot2::element_text(size = 6))

phylum_colors <- readRDS(file.path(BASE, "results/phylum_colors_mags.rds"))
rd <- function(f) read.delim(file.path(INDIR, f), check.names = FALSE, quote = "", comment.char = "")

# ---- a. depth vs breadth ----------------------------------------------------

cnt <- rd("MAG_read_count.tsv");      rownames(cnt) <- cnt[[1]]; cnt[[1]] <- NULL
frac <- rd("MAG_covered_fraction.tsv"); rownames(frac) <- frac[[1]]; frac[[1]] <- NULL
frac <- frac[rownames(cnt), colnames(cnt)]
mat_counts <- as.matrix(cnt); mat_frac <- as.matrix(frac)
storage.mode(mat_counts) <- "numeric"; storage.mode(mat_frac) <- "numeric"

check_df <- data.frame(Count = as.vector(mat_counts), Coverage = as.vector(mat_frac)) %>%
  filter(Count > 0)

cor_res <- suppressWarnings(cor.test(check_df$Coverage, log10(check_df$Count), method = "spearman"))
r_val <- round(cor_res$estimate, 2)
n_total <- nrow(check_df)
n_ghosts <- sum(check_df$Coverage < 0.10); pct_ghosts <- round(100 * n_ghosts / n_total, 1)
n_dominant <- sum(check_df$Coverage > 0.75); pct_dom <- round(100 * n_dominant / n_total, 1)
say(sprintf("PANEL a: %s MAG-by-donor observations with >0 reads", format(n_total, big.mark = ",")))
say(sprintf("   Spearman(coverage, log10 count) R = %.2f, p = %.3g", r_val, cor_res$p.value))
say(sprintf("   ghost hits (<10%% coverage): %s (%.1f%%) | high confidence (>75%%): %s (%.1f%%)",
            format(n_ghosts, big.mark = ","), pct_ghosts, format(n_dominant, big.mark = ","), pct_dom))

stats_label <- paste0("Global Correlation: R = ", r_val, " (p < 2.2e-16)\n",
                      "Ghost Hits (<10% cov): ", format(n_ghosts, big.mark = ","), " (", pct_ghosts, "%)\n",
                      "High Confidence (>75% cov): ", format(n_dominant, big.mark = ","), " (", pct_dom, "%)")

p_qc <- ggplot(check_df, aes(x = Coverage, y = log10(Count))) +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c(option = "magma", name = "Frequency") +
  geom_smooth(method = "gam", color = "cyan", linetype = "solid", se = FALSE, linewidth = 0.8, alpha = 0.8) +
  geom_vline(xintercept = 0.10, color = "red",  linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 0.75, color = "navy", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 0.03, y = 0.85 * max(log10(check_df$Count)),
           label = "GHOST POPULATION\n(False Positives)", color = "red", angle = 90,
           fontface = "bold", size = PT5) +
  annotate("text", x = 0.4, y = 0.4 * max(log10(check_df$Count)), label = "Diagonal Trend Curve",
           color = "#027878", fontface = "italic", size = PT5) +
  annotate("text", x = 0.9, y = 0.80 * max(log10(check_df$Count)),
           label = "DOMINANT\n(High Confidence)", color = "navy", fontface = "bold", size = PT5) +
  annotate("label", x = 0.55, y = max(log10(check_df$Count)), label = stats_label,
           hjust = 0.5, vjust = 1, fill = "white", alpha = 0.9, size = PT5, color = "black") +
  labs(title = "Quality Control: Depth vs Breadth Populations",
       subtitle = "Ghost (<10% cov) vs True Rare (Diagonal) vs Dominant (>75% cov)",
       x = "Genome Coverage Fraction (0-1)", y = "Log10 Read Count") +
  theme_minimal(base_size = 7) + small_text +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        aspect.ratio = 1)   # square panel, not a wide slab

# ---- b. per-donor phylum composition ---------------------------------------

ps <- readRDS(file.path(N76, "ps_mags_n76.rds"))
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
df_bar <- psmelt(ps_rel) %>%
  mutate(Phylum = as.character(Phylum),
         Phylum = ifelse(Phylum %in% names(phylum_colors), Phylum, "Other"))

k <- 100   # log-compression constant
df_bar_comp <- df_bar %>%
  group_by(Sample) %>%
  mutate(Abund_comp = log1p(k * Abundance),
         Abund_comp = Abund_comp / sum(Abund_comp) * 100) %>%
  ungroup()
phylum_order <- names(sort(tapply(df_bar$Abundance, df_bar$Phylum, sum), decreasing = TRUE))
df_bar_comp$Phylum <- factor(df_bar_comp$Phylum, levels = phylum_order)

# per-donor phylum total first, then average across donors (not the mean of
# per-MAG shares, which would just report a typical single MAG)
top_phyla <- df_bar %>%
  group_by(Group, Sample, Phylum) %>% summarise(tot = sum(Abundance), .groups = "drop") %>%
  group_by(Group, Phylum) %>% summarise(m = mean(tot), .groups = "drop") %>%
  group_by(Group) %>% slice_max(m, n = 3) %>%
  summarise(s = paste(sprintf("%s %.1f%%", Phylum, 100 * m), collapse = ", "), .groups = "drop")
say("PANEL b: mean relative abundance, top three phyla per group")
for (i in seq_len(nrow(top_phyla))) say(sprintf("   %-8s %s", top_phyla$Group[i], top_phyla$s[i]))

p_bar <- ggplot(df_bar_comp, aes(x = Sample, y = Abund_comp, fill = Phylum)) +
  geom_col(width = 1, position = position_stack(reverse = TRUE)) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phylum_colors, name = "Phylum (GTDB)", labels = md_taxon) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Community Composition (MAG Biomass)",
       subtitle = sprintf("Log-compressed abundance (k = %d); one bar per donor", k),
       x = NULL, y = "Compressed Relative Abundance (%)") +
  theme_minimal(base_size = 7) + small_text +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.text = ggtext::element_markdown(size = 5.5),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text = element_text(face = "bold", size = 6.5),
        panel.spacing = unit(0.2, "lines"),
        legend.position = "right")

# ---- c. MAG quality against MIMAG thresholds --------------------------------
# Catalogue-wide marginal densities, so one-MAG phyla are not drawn as
# prominently as Bacillota; phylum is encoded by the scatter-point colours.

phy_colors <- readRDS(file.path(BASE, "results/phylum_colors_mags.rds"))
qc_all <- rd("MAG_quality_taxonomy.tsv") %>%
  transmute(Genome_ID = catalog_id, Phylum = phylum, Family = family,
            Completeness = completeness, Contamination = contamination)
n_total <- nrow(qc_all)
n_high  <- sum(qc_all$Completeness >= 90 & qc_all$Contamination <= 5, na.rm = TRUE)
n_med   <- sum(qc_all$Completeness >= 50 & qc_all$Contamination <= 10, na.rm = TRUE) - n_high
say("")
say(sprintf("PANEL c: %d MAGs | high quality %d (>=90%%/<=5%%) | medium %d (>=50%%/<=10%%)",
            n_total, n_high, n_med))

stroke_colors <- colorspace::darken(phy_colors, amount = 0.3)
qc_df <- qc_all %>%
  filter(Completeness >= 50, Completeness <= 100) %>%
  mutate(Phylum = factor(ifelse(Phylum %in% names(phy_colors), Phylum, "Other")))
unmapped_qc_phyla <- setdiff(unique(as.character(qc_all$Phylum)), names(phy_colors))
unmapped_qc_phyla <- unmapped_qc_phyla[!is.na(unmapped_qc_phyla) & nzchar(unmapped_qc_phyla)]
if (length(unmapped_qc_phyla)) {
  stop("Named MAG phyla lack palette entries: ", paste(sort(unmapped_qc_phyla), collapse = ", "))
}
set.seed(42)
qc_df <- qc_df[sample(nrow(qc_df)), ]

p_qc_sc <- ggplot(qc_df, aes(Completeness, Contamination, colour = Phylum, fill = Phylum)) +
  annotate("rect", xmin = 90, xmax = 100, ymin = 0, ymax = 5, fill = "green", alpha = 0.1, color = NA) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey50") +
  geom_point(shape = 21, size = 2.5, stroke = 0.4, alpha = 0.7,
             position = position_jitter(width = 0.3, height = 0.1, seed = 42)) +
  scale_fill_manual(values = phy_colors) +
  scale_colour_manual(values = stroke_colors) +
  scale_x_continuous(breaks = seq(50, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 10, 2.5)) +
  coord_cartesian(xlim = c(50, 100), ylim = c(0, 10)) +
  labs(title = "MAG Quality Assessment (MIMAG Standards)",
       subtitle = sprintf("Total: %d | High Quality: %d | Medium Quality: %d", n_total, n_high, n_med),
       x = "Completeness (%)", y = "Contamination (%)") +
  theme_bw(base_size = 7) + small_text +
  theme(legend.position = "none", aspect.ratio = 1)
p_qc_marg <- ggMarginal(
  p_qc_sc, type = "density", groupColour = FALSE, groupFill = FALSE,
  colour = "grey35", fill = "grey70", alpha = 0.55, size = 5
)

# ---- save -------------------------------------------------------------------

ggsave(file.path(OUTDIR, "SuppFig1a_QC_Depth_vs_Breadth_n76.pdf"), p_qc, width = 7, height = 5.5,
       device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "SuppFig1b_Community_Composition_n76.pdf"), p_bar, width = 9, height = 4.5,
       device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "SuppFig1c_MAG_Quality_n76.pdf"), p_qc_marg, width = 6, height = 6,
       device = cairo_pdf, bg = "transparent")

# ---- d. ZOE reference concordance -------------------------------------------
# Drawn from the per-SGB table written by zoe_n76.R, so this panel does not
# depend on the ZOE supplementary spreadsheet being present. It gives the
# species-level concordance behind the per-donor health index of Figure 1c.
HEALTH_PAL <- c("#1B7837", "#A6DBA0", "#F7F7F7", "#C2A5CF", "#762A83")
zoe <- read.csv(file.path(N76, "ZOE_Concordance_PerSGB_n76.csv"))
zc <- suppressWarnings(cor.test(zoe$HEALTH_ranks, zoe$cor_with_case, method = "spearman"))
say(sprintf("d. ZOE concordance: Spearman rho = %.3f, p = %.3g (n = %d SGBs)",
            unname(zc$estimate), zc$p.value, nrow(zoe)))
p_zoe <- ggplot(zoe, aes(HEALTH_ranks, cor_with_case, colour = HEALTH_ranks)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.7, se = TRUE, formula = y ~ x) +
  scale_colour_gradientn(colors = HEALTH_PAL, limits = c(0, 1), name = "ZOE\nhealth-rank") +
  annotate("text", x = 0, y = max(zoe$cor_with_case), hjust = 0, vjust = 1, fontface = "bold",
           size = PT5 + 0.2, label = sprintf("Spearman rho = %.2f, p = %.2g\nn = %d SGBs",
                                       unname(zc$estimate), zc$p.value, nrow(zoe))) +
  labs(title = "Disease-associated shift aligns with the ZOE healthy-gut reference",
       subtitle = "Each point is one SGB. ZOE-favourable species are depleted in COVID; ZOE-unfavourable species are enriched.",
       x = "ZOE microbiome health-rank (0 = favourable, 1 = unfavourable)",
       y = "SGB correlation with COVID case status\n(Spearman, relative abundance)") +
  theme_bw(base_size = 7) + small_text +
  theme(legend.position = "right", aspect.ratio = 1)
ggsave(file.path(OUTDIR, "SuppFig1d_ZOE_Concordance_n76.pdf"), p_zoe, width = 7.5, height = 6.2,
       device = cairo_pdf, bg = "transparent")

tag <- theme(plot.tag = element_text(face = "bold", size = 14))
# A4 portrait: the two square panels share the top row, the donor bar chart is
# a wide strip beneath them, and the square ZOE panel closes the page. The old
# canvas was 10 x 21 in, which no journal page can carry without shrinking the
# type below legibility.
supp1 <- ((p_qc + labs(tag = "a") + tag) | (wrap_elements(full = p_qc_marg) + labs(tag = "c") + tag)) /
         (p_bar + labs(tag = "b") + tag) /
         (p_zoe + labs(tag = "d") + tag) +
  plot_layout(heights = c(1, 0.62, 1)) +
  plot_annotation(
    title = sprintf("Quality control, community composition and ZOE reference concordance (%d MAGs, %d donors)",
                    ntaxa(ps), nsamples(ps)),
    theme = theme(plot.title = element_text(face = "bold", size = 9)))
ggsave(file.path(OUTDIR, "Supplementary_Figure_1_n76.pdf"), supp1, width = 8.27, height = 11.69,
       device = cairo_pdf, bg = "transparent")

writeLines(STATS, file.path(OUTDIR, "SuppFig1_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Supplementary_Figure_1_n76.pdf"))
