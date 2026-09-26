# =============================================================================
# Figure 2, rebuilt on the 584-MAG catalogue (578 bacterial tips)
#
#   a  circular ML tree (GTDB-Tk de novo, LG + gamma, rooted on Bacteroidota)
#   b  MIMAG quality scatter with per-phylum marginal densities
#   c  phylogenetic correlogram of LFC, inset: patristic distance per GTDB rank
#   d  local Moran's I scatter (individual vs neighbour-lagged LFC)
#
# Panel code is transcribed from ~/MAG_Analysis/MAG_tree.R. Statistics come
# from figure2_n76_stats.R (fig2/figure2_stats_n76.rds); run that first.
# Tree tips equal catalog_id exactly, so the original's fuzzy name mapper is
# not needed.
#
# Outputs -> result2/n76/fig2/
# =============================================================================

suppressPackageStartupMessages({
  library(ggtree); library(treeio); library(ape); library(ggtreeExtra)
  library(ggnewscale); library(ggstar); library(phyloseq)
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(patchwork); library(ggExtra); library(colorspace)
  library(phylosignal)   # S3 plot method for the restored correlogram
  library(lme4); library(lmerTest)
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
OUTDIR <- file.path(N76, "fig2")
source(file.path(BASE, "taxon_italics_n76.R"))   # expr_from_md(): taxon names in italics
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")

STATS <- c()
say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

# Standalone panel-a label handler. The A4 composite has its own fixed label
# policy and is intentionally unaffected by this switch.
TREE_LABEL_MODE <- tolower(Sys.getenv("FIG2_TREE_LABELS", "off"))
if (!TREE_LABEL_MODE %in% c("on", "off"))
  stop("FIG2_TREE_LABELS must be 'on' or 'off'")
SHOW_STANDALONE_TREE_LABELS <- identical(TREE_LABEL_MODE, "on")
say("standalone tree labels: ", TREE_LABEL_MODE)

st <- readRDS(file.path(OUTDIR, "figure2_stats_n76.rds"))
tr <- read.tree(TREE)
phy_colors <- readRDS(file.path(BASE, "results/phylum_colors_mags.rds"))   # GTDB names
blank <- function(x) is.na(x) | x == "" | grepl("^[a-z]__$", x)

# ---- per-tip metadata -------------------------------------------------------

ps <- readRDS(file.path(N76, "ps_mags_n76.rds"))
otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
otu_rel <- sweep(otu, 2, colSums(otu), "/") * 100
cov_mat <- as.matrix(read.delim(file.path(INDIR, "MAG_covered_fraction.tsv"),
                                check.names = FALSE, row.names = 1))

meta_final <- st$meta_df %>%
  transmute(label = Genome_ID, Phylum, Family, Genus, Species, Is_Novel,
            Completeness, Contamination, lfc_GroupCase, p_GroupCase, diff_GroupCase) %>%
  mutate(
    Mean_Abundance = rowMeans(otu_rel)[label],
    Mean_Coverage  = rowMeans(cov_mat)[label] * 100,
    Phylum = ifelse(Phylum %in% names(phy_colors), Phylum, "Other"),
    diff_GroupCase = coalesce(as.logical(diff_GroupCase), FALSE),
    EnrichedStatus = factor(case_when(
      diff_GroupCase & lfc_GroupCase > 0 ~ "Enriched",
      diff_GroupCase & lfc_GroupCase < 0 ~ "Depleted",
      TRUE ~ "NS"), levels = c("Enriched", "Depleted", "NS")),
    Strong_Hit = !is.na(p_GroupCase) & p_GroupCase < 0.001 & abs(lfc_GroupCase) > 2,
    Label_Text = case_when(
      Is_Novel == "Yes" & !blank(Genus)  ~ paste0(Genus, " sp. (Novel)"),
      Is_Novel == "Yes" & !blank(Family) ~ paste0(Family, " sp. nov. (Novel)"),
      Strong_Hit & !blank(Species)       ~ Species,
      Strong_Hit & !blank(Genus)         ~ paste0(Genus, " sp."),
      TRUE ~ NA_character_))

# The full radial tree is delivered on a large standalone canvas and retains
# every strong-hit label.  At A4 size, 170 labels are not legible.  For the
# composite only, label a balanced set of five enriched and five depleted MAGs.
# A small circular separation guard prevents the strongest members of one dense
# clade from printing on top of one another. All novel MAGs remain marked by a
# gold star. This changes annotation density only; tips, rings and tests remain.
select_spaced_labels <- function(d, n_keep = 5, min_tip_gap = 10) {
  d <- d %>% arrange(p_GroupCase, desc(abs(lfc_GroupCase)))
  chosen_id <- character(); chosen_pos <- integer(); n_tip <- length(tr$tip.label)
  for (i in seq_len(nrow(d))) {
    pos <- match(d$label[i], tr$tip.label)
    gap <- if (!length(chosen_pos)) Inf else
      min(pmin(abs(pos - chosen_pos), n_tip - abs(pos - chosen_pos)))
    if (gap >= min_tip_gap) {
      chosen_id <- c(chosen_id, d$label[i]); chosen_pos <- c(chosen_pos, pos)
    }
    if (length(chosen_id) == n_keep) break
  }
  chosen_id
}
a4_top_ids <- unlist(lapply(c("Enriched", "Depleted"), function(z)
  select_spaced_labels(filter(meta_final, diff_GroupCase, EnrichedStatus == z))),
  use.names = FALSE)
meta_final <- meta_final %>%
  mutate(
    Species_Short = sub("^([[:alpha:]])[^ ]+ ", "\\1. ", Species),
    Label_Text_A4 = case_when(
      label %in% a4_top_ids & !blank(Species) ~ Species_Short,
      label %in% a4_top_ids & !blank(Genus)   ~ paste0(Genus, " sp."),
      TRUE ~ NA_character_))

stopifnot(setequal(meta_final$label, tr$tip.label))
say(sprintf("tips %d | significant %d (enriched %d, depleted %d) | strong hits %d | novel %d | labelled %d",
            nrow(meta_final), sum(meta_final$diff_GroupCase),
            sum(meta_final$EnrichedStatus == "Enriched"), sum(meta_final$EnrichedStatus == "Depleted"),
            sum(meta_final$Strong_Hit), sum(meta_final$Is_Novel == "Yes"),
            sum(!is.na(meta_final$Label_Text))))
say(sprintf("A4 tree labels %d: five spaced enriched plus five spaced depleted; all %d novel MAGs remain starred",
            sum(!is.na(meta_final$Label_Text_A4)), sum(meta_final$Is_Novel == "Yes")))

# Quality tiers over the whole 584-genome catalogue (archaea included), for the
# tree centre and panel b
qc_all <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                     quote = "", comment.char = "") %>%
  transmute(Genome_ID = catalog_id, Phylum = phylum, Species = species, Genus = genus,
            Family = family, Completeness = completeness, Contamination = contamination)
n_total <- nrow(qc_all)
n_high  <- sum(qc_all$Completeness >= 90 & qc_all$Contamination <= 5, na.rm = TRUE)
n_med   <- sum(qc_all$Completeness >= 50 & qc_all$Contamination <= 10, na.rm = TRUE) - n_high
n_named <- sum(!blank(qc_all$Species))
n_novel <- sum(blank(qc_all$Species) & (!blank(qc_all$Genus) | !blank(qc_all$Family)))
say(sprintf("catalogue: %d MAGs | high quality %d | medium quality %d | named species %d (%.0f%%) | putative novel %d",
            n_total, n_high, n_med, n_named, 100 * n_named / n_total, n_novel))

# ---- a. circular tree (MAG_tree.R lines 305-490) ----------------------------

gap_text    <- 0.25
text_space  <- 0.28   # gap from the labels to the first ring, as a fraction of
                      # tree width. 0.9 left most of the circle empty.
gap_rings   <- 0.03
ring_width  <- 0.05
base_family <- "sans"

p <- ggtree(tr, layout = "fan", open.angle = 3, size = 0.75) %<+% meta_final +
  geom_tree(aes(color = Phylum), size = 0.75) +
  scale_color_manual(values = phy_colors, name = "Phylum",
                     breaks = sort(intersect(names(phy_colors), unique(meta_final$Phylum))),
                     labels = function(x) expr_from_md(md_taxon(x)),
                     guide = guide_legend(nrow = 2, byrow = TRUE))
max_x <- max(p$data$x, na.rm = TRUE)

p_dots <- p + new_scale_color() +
  geom_tippoint(aes(color = EnrichedStatus), size = 2, na.rm = TRUE) +
  scale_color_manual(values = c("Enriched" = "#ef6548", "Depleted" = "#1c9099", "NS" = "transparent"),
                     breaks = c("Enriched", "Depleted"), name = "Differential abundance",
                     guide = guide_legend(nrow = 1))

# Radial labels without repel: each label keeps its own tip slot.  Building the
# outer layers through one helper guarantees that the standalone and A4 trees
# differ only in which labels are printed.
build_tree_plot <- function(label_col, show_labels = TRUE, ring_gap = NULL, sym = 1) {
  # rings sit just outside the tips when no labels are printed between them
  text_space <- if (!is.null(ring_gap)) ring_gap else if (show_labels) text_space else 0.03
  p_labels_local <- p_dots
  if (sym != 1) {                      # tip points are set on p_dots, before this
    for (i in seq_along(p_labels_local$layers)) {
      sz <- p_labels_local$layers[[i]]$aes_params$size
      if (!is.null(sz)) p_labels_local$layers[[i]]$aes_params$size <- sz * sym
    }
  }
  if (show_labels) {
    p_labels_local <- p_labels_local +
      geom_tiplab(
        data = function(d) dplyr::filter(d, is.na(.data[[label_col]]) | .data[[label_col]] == ""),
        aes(label = NA_character_), offset = gap_text, align = TRUE,
        linesize = 0.15, linetype = "dashed", colour = "grey85") +
      geom_tiplab(
        data = function(d) dplyr::filter(d, !is.na(.data[[label_col]]), .data[[label_col]] != ""),
        aes(label = .data[[label_col]]), offset = gap_text, align = TRUE,
        linesize = 0.2, linetype = "solid", colour = "black", family = base_family,
        fontface = "italic", size = 1.8)
  }

  p_stars_local <- p_labels_local + new_scale_fill() + new_scale("starshape") +
    geom_star(data = function(d) dplyr::filter(d, Is_Novel == "Yes"),
              aes(x = x + (max_x * 0.02), y = y, starshape = Is_Novel, fill = Is_Novel),
              size = 3 * sym, colour = "black", starstroke = 0.5 * sym) +
    scale_starshape_manual(values = c("Yes" = 1), labels = "Putative novel species", name = NULL,
                           guide = guide_legend(nrow = 1)) +
    scale_fill_manual(values = c("Yes" = "gold"), labels = "Putative novel species", name = NULL,
                      guide = guide_legend(nrow = 1))

  p_ring1_local <- p_stars_local + new_scale_fill() +
    geom_fruit(geom = geom_tile, mapping = aes(y = label, fill = Completeness),
               width = ring_width, offset = text_space, color = "white", size = 0.02) +
    scale_fill_gradient(low = "#C7E9C0", high = "#006D2C", limits = c(50, 100),
                        name = "Completeness (%)", oob = scales::squish,
                        guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                               barwidth = unit(1.4, "cm"), barheight = unit(0.12, "cm")))

  p_ring2_local <- p_ring1_local + new_scale_fill() +
    geom_fruit(geom = geom_tile, mapping = aes(y = label, fill = Contamination),
               width = ring_width, offset = gap_rings, color = "white", size = 0.02) +
    scale_fill_gradient(low = "#FFFFFF", high = "#CB181D", limits = c(0, 10),
                        name = "Contamination (%)", oob = scales::squish,
                        guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                               barwidth = unit(1.4, "cm"), barheight = unit(0.12, "cm")))

  # the abundance axis labels are unreadable at placement size, so they are
  # dropped there and kept on the large standalone tree
  abund_axis <- if (sym < 1) list(axis = "none") else
    list(axis = "x", text.size = 1.5, hjust = 0, title = "Sqrt(Abund)")
  p_ring2_local + new_scale_fill() +
    geom_fruit(geom = geom_bar,
               mapping = aes(y = label, x = sqrt(Mean_Abundance), fill = Mean_Coverage),
               stat = "identity", orientation = "y", width = 1, offset = gap_rings, pwidth = 0.2,
               axis.params = abund_axis) +
    scale_fill_gradient(low = "grey90", high = "#082d54", limits = c(0, 100),
                        name = "Genome coverage (%)", oob = scales::squish,
                        guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                               barwidth = unit(1.4, "cm"), barheight = unit(0.12, "cm"))) +
    annotate("label", x = 0, y = 0, size = 3.2 * sym, linewidth = 0, fill = "white",
             label = sprintf("%d MAGs\n%d bacterial on tree\nHQ %d | MQ %d\n%d putative novel",
                             n_total, nrow(meta_final), n_high, n_med, n_novel))
}

p_tree <- build_tree_plot("Label_Text", show_labels = SHOW_STANDALONE_TREE_LABELS)

# As published: no legends on the tree itself (they were set beside it by hand)
p_tree_clean <- p_tree + theme(legend.position = "none",
                               panel.background = element_rect(fill = "transparent", color = NA),
                               plot.background  = element_rect(fill = "transparent", color = NA))
# The published tree put 401 tips on a 400 mm canvas. Scaling the canvas with tip
# count keeps each tip's angular slot at the published size, which is what lets
# runs of adjacent labelled tips (e.g. the enriched Enterobacteriaceae) stay
# legible at label size 1.8. Scale it down when placing it, as before.
TREE_MM <- round(400 * nrow(meta_final) / 401)
say(sprintf("tree canvas %d mm (published: 400 mm for 401 tips)", TREE_MM))
ggsave(file.path(OUTDIR, "Fig2a_Tree_n76_large.pdf"), p_tree_clean,
       width = TREE_MM, height = TREE_MM, units = "mm", device = cairo_pdf, limitsize = FALSE)

# Placement size: the tree occupies about 107 mm in the Affinity composite, so
# it is drawn at that size directly. Point and star sizes then print as drawn
# instead of shrinking when the 583 mm canvas is scaled down.
# The fan fills about 66% of its canvas, so the page is sized to make the drawn
# circle 107 mm: the diameter it occupies in the Affinity composite.
TREE_CIRCLE_MM <- 107
TREE_PLACE_MM  <- round(TREE_CIRCLE_MM / 0.66)
p_tree_place <- build_tree_plot("Label_Text", show_labels = FALSE, sym = 0.62) +
  # the fan leaves a wide radial margin; removing it lets the same 107 mm carry
  # a larger tree, which is the point of drawing at placement size
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))
say(sprintf("placement tree canvas %d mm for a %d mm drawn circle; rings offset 0.03 of tree width",
            TREE_PLACE_MM, TREE_CIRCLE_MM))
ggsave(file.path(OUTDIR, "Fig2a_Tree_n76.pdf"), p_tree_place,
       width = TREE_PLACE_MM, height = TREE_PLACE_MM, units = "mm",
       device = cairo_pdf, limitsize = FALSE)

p_tree_leg <- p_tree + theme(legend.position = "right",
                             legend.text  = element_text(size = 8),
                             legend.title = element_text(size = 9, face = "bold"))

# ---- b. strain-level diversity (inStrain 1.10.0) ----------------------------
# The quality scatter that used to sit here moved to Supplementary Figure 1; it
# is a QC display and the space is better spent on a result. qc_all is still
# computed above because n_high / n_med are quoted in the text.
#
# This panel compares the SAME genome between arms, so the Enterobacteriaceae
# bloom cannot manufacture the effect. Thresholds follow inStrain's guidance:
# nucleotide diversity needs both depth and breadth before it means anything.
MIN_COV <- 5; MIN_BREADTH <- 0.5; MIN_DONORS <- 5

gi <- read_tsv(file.path(INDIR, "instrain_genome_info_76.tsv"), show_col_types = FALSE) %>%
  mutate(genome = sub("[.]fa$", "", genome))
smeta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE) %>%
  transmute(donor = sample, Group = factor(group, levels = c("Control", "Case")))
gtax <- qc_all %>% transmute(genome = Genome_ID, Phylum, Family)

ok_s <- gi %>% left_join(smeta, by = "donor") %>%
  filter(coverage >= MIN_COV, breadth_minCov >= MIN_BREADTH, !is.na(nucl_diversity))
pair_s <- ok_s %>% count(genome, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Case >= MIN_DONORS, Control >= MIN_DONORS)
cmp_s <- ok_s %>% filter(genome %in% pair_s$genome) %>% left_join(gtax, by = "genome")

m_s  <- lmer(log10(nucl_diversity + 1e-6) ~ Group + (1 | genome), data = cmp_s)
cf_s <- summary(m_s)$coefficients
cr_s <- cmp_s %>% filter(!is.na(nucl_diversity_rarefied))
m_r  <- lmer(log10(nucl_diversity_rarefied + 1e-6) ~ Group + (1 | genome), data = cr_s)
cf_r <- summary(m_r)$coefficients

pg_s <- cmp_s %>% group_by(genome, Phylum, Family) %>%
  summarise(ctrl_pi = median(nucl_diversity[Group == "Control"]),
            case_pi = median(nucl_diversity[Group == "Case"]),
            rar_ctrl = median(nucl_diversity_rarefied[Group == "Control"], na.rm = TRUE),
            rar_case = median(nucl_diversity_rarefied[Group == "Case"], na.rm = TRUE),
            .groups = "drop")
n_lower_rar <- sum(pg_s$rar_case < pg_s$rar_ctrl, na.rm = TRUE)
say(sprintf("strain panel: %d genomes with >=%d donors per arm at >=%dx and >=%.0f%% breadth",
            nrow(pg_s), MIN_DONORS, MIN_COV, 100 * MIN_BREADTH))
say(sprintf("   median nucleotide diversity: control %.5f vs case %.5f (%.1f-fold)",
            median(pg_s$ctrl_pi), median(pg_s$case_pi), median(pg_s$ctrl_pi) / median(pg_s$case_pi)))
say(sprintf("   mixed model beta = %+.3f (p = %.3g); rarefied beta = %+.3f (p = %.3g); %d of %d lower in cases (rarefied)",
            cf_s["GroupCase", 1], cf_s["GroupCase", ncol(cf_s)],
            cf_r["GroupCase", 1], cf_r["GroupCase", ncol(cf_r)], n_lower_rar, nrow(pg_s)))
say(sprintf("   coverage control %.1fx vs case %.1fx (Wilcoxon p = %.2f), so not a depth artefact",
            median(cmp_s$coverage[cmp_s$Group == "Control"]),
            median(cmp_s$coverage[cmp_s$Group == "Case"]),
            wilcox.test(coverage ~ Group, data = cmp_s, exact = FALSE)$p.value))

p_strain <- ggplot(pg_s, aes(ctrl_pi, case_pi)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linetype = "dashed") +
  geom_point(aes(fill = Phylum), shape = 21, size = 3.2, stroke = 0.4, alpha = 0.85,
             colour = "grey25") +
  scale_fill_manual(values = phy_colors, na.value = "grey70") +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x = min(pg_s$ctrl_pi), y = max(pg_s$case_pi), hjust = 0, vjust = 1,
           size = 3.4, fontface = "italic", label = "less diverse in patients") +
  labs(title = "Within-species strain diversity",
       subtitle = sprintf("%d genomes covered in both arms; %.1f-fold lower in patients (p = %.1g)",
                          nrow(pg_s), median(pg_s$ctrl_pi) / median(pg_s$case_pi),
                          cf_s["GroupCase", ncol(cf_s)]),
       x = "Nucleotide diversity, controls", y = "Nucleotide diversity, patients") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))
ggsave(file.path(OUTDIR, "Fig2b_StrainDiversity_n76.pdf"), p_strain, width = 6, height = 6,
       device = cairo_pdf, bg = "transparent")
write_csv(pg_s, file.path(OUTDIR, "Fig2b_strain_per_genome_n76.csv"))

# ---- c. correlogram + taxonomic calibration inset ---------------------------

cr <- st$corr_res
lim <- st$exact_limit
say(sprintf("correlogram: zero-crossing %.3f (lower-bound crossing, published method, %.3f); significant classes +%d / -%d of %d",
            lim, st$lower_bound_limit, sum(cr$sig == "positive"), sum(cr$sig == "negative"), nrow(cr)))

# The canonical standalone panel uses the original discrete red/blue
# significance dots.  Keep the native phylosignal rendering separately under
# the historical *_base_n76.pdf filename so the two outputs are not duplicates.
p_corr <- ggplot(cr, aes(d.mean, correlation)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(aes(y = lower, linetype = "95% confidence interval"), colour = "black", linewidth = 0.4) +
  geom_line(aes(y = upper, linetype = "95% confidence interval"), colour = "black", linewidth = 0.4) +
  geom_line(aes(linetype = "Moran's I"), colour = "black", linewidth = 1) +
  geom_point(data = dplyr::filter(cr, sig != "ns"), aes(colour = sig), size = 1.8) +
  { if (!is.na(lim)) geom_vline(xintercept = lim, colour = "red", linetype = "dashed") } +
  scale_linetype_manual(values = c("Moran's I" = "solid", "95% confidence interval" = "dashed"),
                        breaks = c("Moran's I", "95% confidence interval"), name = NULL) +
  scale_colour_manual(values = c(positive = "red", negative = "blue"),
                      labels = c(positive = "Significant positive", negative = "Significant negative"),
                      name = NULL) +
  labs(title = "Phylogenetic Correlogram of Enrichment",
       subtitle = sprintf("Global Moran's I = %.3f (p = %.4f); zero-crossing at %.2f",
                          st$global_sig$stat$I, st$global_sig$pvalue$I, lim),
       x = "Phylogenetic Distance", y = "Moran's I Autocorrelation") +
  theme_bw(base_size = 11) +
  guides(linetype = guide_legend(order = 1, keywidth = unit(0.8, "cm")),
         colour = guide_legend(order = 2, override.aes = list(size = 2.2))) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.76, 0.28),
        legend.justification.inside = c(0.5, 0.5),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_rect(fill = scales::alpha("white", 0.88), colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.key.height = unit(0.30, "cm"),
        legend.text = element_text(size = 7))

# Native S3 version retained for the explicitly named base panel and for the
# scripted A4 composite only.  It is not written over the dot-based panel.
draw_native_corr <- function(a4 = FALSE) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  if (a4) par(mar = c(3.25, 3.45, 2.15, 0.55))
  plot(st$correlogram,
       xlab = "Phylogenetic Distance", ylab = "Moran's I",
       main = "Phylogenetic Correlogram of Enrichment",
       cex.main = if (a4) 0.62 else 1.20,
       cex.lab = if (a4) 0.55 else 1.00,
       cex.axis = if (a4) 0.48 else 1.00)
  if (!is.na(lim)) abline(v = lim, col = "red", lty = 2)
  if (a4) mtext("c", side = 3, line = 1.05, adj = 0, font = 2, cex = 0.72)
  legend_labels <- if (a4)
    c("Moran's I", "95% CI", "Significant +", "Significant -") else
    c("Moran's I (Mean)", "95% Conf. Interval", "Significant Positive", "Significant Negative")
  legend(x = 1.78, y = 0.215,
         legend = legend_labels,
         col = c("black", "black", "red", "blue"),
         lty = c(1, 2, 1, 1), lwd = c(2, 1, 3, 3),
         bty = "n", cex = if (a4) 0.43 else 0.66,
         seg.len = 1.6, x.intersp = 0.65, y.intersp = 0.90)
}

p_calibration <- ggplot(st$tax_dist_df, aes(Rank, Distance)) +
  geom_boxplot(fill = "grey50", alpha = 0.7, outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
  { if (!is.na(lim)) geom_hline(yintercept = lim, colour = "red", linetype = "dashed", linewidth = 0.7) } +
  { if (!is.na(lim)) annotate("text", x = 1.6, y = lim, vjust = -0.5, size = 2.6, colour = "red",
                              fontface = "bold", label = sprintf("Signal Limit: %.2f", lim)) } +
  labs(x = NULL, y = "Patristic distance") +
  theme_bw(base_size = 7) +
  theme(plot.background = element_rect(fill = "transparent", colour = "grey50", linewidth = 0.4),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTDIR, "Fig2c_Taxonomic_Calibration_n76.pdf"), p_calibration, width = 4, height = 4)

make_corr_native_inset <- function(base, calibration = p_calibration) {
  base + inset_element(calibration, left = 0.55, bottom = 0.54, right = 0.985, top = 0.90)
}
make_corr_dot_inset <- function(base, calibration = p_calibration) {
  base + inset_element(calibration, left = 0.52, bottom = 0.53, right = 0.99, top = 0.99)
}
p_corr_native <- wrap_elements(full = ~ draw_native_corr())
p_corr_dot_inset <- make_corr_dot_inset(p_corr)
p_corr_base_inset <- make_corr_native_inset(p_corr_native)
ggsave(file.path(OUTDIR, "Fig2c_Correlogram_n76.pdf"), p_corr_dot_inset,
       width = 6.2, height = 6.2, device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "Fig2c_Correlogram_base_n76.pdf"), p_corr_base_inset,
       width = 6.2, height = 6.2, device = cairo_pdf, bg = "transparent")

# ---- d. local Moran's I (MAG_tree.R lines 1595-1690) ------------------------

cat_short <- c(
  "Phylogenetic Hotspot (Conserved Enrichment)" = "Conserved Enrichment",
  "Phylogenetic Coldspot (Conserved Depletion)" = "Conserved Depletion",
  "Evolutionary Outlier (Unique Gain)"          = "Unique Gain",
  "Evolutionary Outlier (Unique Loss)"          = "Unique Loss",
  "Non-significant"                             = "Non-significant")
phylo_colors <- c(
  "Non-significant"      = "grey90",
  "Conserved Enrichment" = "#D55E00",
  "Conserved Depletion"  = "#0072B2",
  "Unique Gain"          = "#CC79A7",
  "Unique Loss"          = "#009E73")

analysis_df <- st$analysis_df %>%
  mutate(Category = factor(cat_short[Phylo_Category], levels = names(phylo_colors)),
         Plot_Order = case_when(Category == "Non-significant" ~ 1,
                                Category %in% c("Conserved Enrichment", "Conserved Depletion") ~ 2,
                                TRUE ~ 3)) %>%
  arrange(Plot_Order)
say("local Moran categories: ", paste(names(table(analysis_df$Category)), table(analysis_df$Category),
                                      sep = " = ", collapse = "; "))

p_moran <- ggplot(analysis_df, aes(LogFC, Neighbor_Lag_LogFC, colour = Category)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, linetype = "dotdash",
              colour = "black", alpha = 0.6) +
  scale_colour_manual(values = phylo_colors, drop = TRUE, name = NULL) +
  labs(title = "Phylogenetic Clustering of Enrichment",
       subtitle = "Local Moran's I Scatterplot with Global Trend",
       x = "Log Fold Change (Individual)", y = "Log Fold Change (Phylogenetic Neighbors)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification.inside = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8))
ggsave(file.path(OUTDIR, "Fig2d_LocalMoran_n76.pdf"), p_moran, width = 6.2, height = 6.2)

p_box <- ggplot(analysis_df, aes(Category, LogFC, fill = Category)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2, seed = 42), alpha = 0.3, size = 1) +
  scale_fill_manual(values = phylo_colors) +
  labs(title = "LogFC Distribution", x = NULL, y = NULL) +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTDIR, "Fig2d_LIPA_Boxplot_n76.pdf"), p_box, width = 4, height = 5)

# ---- A4 publication composite ----------------------------------------------

# Build at final A4 size. A fixed aspect ratio lets b-d remain square; patchwork
# is allowed to leave whitespace rather than stretching those panels.
panel_a4_theme <- theme(
  aspect.ratio = 1,
  plot.title = element_text(face = "bold", size = 7.5, margin = margin(b = 1.5)),
  plot.subtitle = element_text(size = 5.5, margin = margin(b = 2)),
  axis.title = element_text(size = 6.5),
  axis.text = element_text(size = 5.5),
  legend.text = element_text(size = 5.2),
  legend.key.height = unit(0.24, "cm"),
  legend.key.width = unit(0.35, "cm"),
  plot.margin = margin(2, 2, 2, 2, unit = "pt"))
inset_a4_theme <- theme(
  axis.title = element_text(size = 5.2),
  axis.text = element_text(size = 4.7),
  plot.margin = margin(1, 1, 1, 1, unit = "pt"))
tag_a4 <- theme(plot.tag = element_text(face = "bold", size = 9),
                plot.tag.position = c(0, 1))

p_strain_a4 <- p_strain +
  labs(subtitle = sprintf("%d genomes in both arms\n%.1f-fold lower in patients (p = %.1g)",
                          nrow(pg_s), median(pg_s$ctrl_pi) / median(pg_s$case_pi),
                          cf_s["GroupCase", ncol(cf_s)]), tag = "b") +
  panel_a4_theme + tag_a4
p_calibration_a4 <- p_calibration + inset_a4_theme
p_corr_native_a4 <- wrap_elements(full = ~ draw_native_corr(a4 = TRUE))
p_corr_a4 <- make_corr_native_inset(p_corr_native_a4, p_calibration_a4)
p_moran_a4 <- p_moran +
  labs(title = "Phylogenetic clustering",
       subtitle = "Local Moran's I with global trend", tag = "d") +
  panel_a4_theme + tag_a4

# Keep the A4 tree itself unobstructed.  Its balanced 10-label annotation set
# is built separately from the full standalone tree, and its compact guide is
# extracted into the intentional whitespace above panels b-d.
p_tree_a4_data <- build_tree_plot("Label_Text_A4")
p_tree_a4 <- p_tree_a4_data +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(10, 8, 0, 8, unit = "pt"))
tree_legend_theme <- theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.box.just = "center",
        legend.text = element_text(size = 4.7),
        legend.title = element_text(size = 5.2, face = "bold"),
        legend.key.height = unit(0.18, "cm"),
        legend.key.width = unit(0.28, "cm"),
        legend.spacing.x = unit(0.06, "cm"),
        legend.spacing.y = unit(0.01, "cm"),
        legend.box.spacing = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0, 0, 0, 0))
p_tree_a4_cat_guides <- p_tree_a4_data +
  guides(fill_ggnewscale_3 = "none", fill_ggnewscale_4 = "none", fill = "none") +
  tree_legend_theme
p_tree_a4_cont_guides <- p_tree_a4_data +
  guides(colour_ggnewscale_1 = "none", colour = "none", starshape = "none",
         fill_ggnewscale_2 = "none") +
  tree_legend_theme + theme(legend.spacing.x = unit(0.35, "cm"))
tree_cat_grob <- cowplot::get_legend(p_tree_a4_cat_guides)
tree_cont_grob <- cowplot::get_legend(p_tree_a4_cont_guides)
stopifnot(!is.null(tree_cat_grob), !is.null(tree_cont_grob))
tree_legend_a4 <- (wrap_elements(full = tree_cat_grob) /
                      wrap_elements(full = tree_cont_grob)) +
  plot_layout(heights = c(1.05, 0.95))

bottom_a4 <- (p_strain_a4 | p_corr_a4 | p_moran_a4) +
  plot_layout(widths = c(1, 1, 1))
fig2_a4 <- ((wrap_elements(full = p_tree_a4) + labs(tag = "a") + tag_a4) /
              tree_legend_a4 / plot_spacer() / bottom_a4) +
  plot_layout(heights = c(2.98, 0.42, 0.03, 1))

a4_path <- file.path(OUTDIR, "Figure_2_n76_A4.pdf")
ggsave(a4_path, fig2_a4, width = 210, height = 297, units = "mm",
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)
stopifnot(file.copy(a4_path, file.path(OUTDIR, "Figure_2_n76.pdf"), overwrite = TRUE))
publication_path <- file.path(BASE, "Figures_publication", "F2 copy.pdf")
stopifnot(file.copy(a4_path, publication_path, overwrite = TRUE))

# ---- stats file -------------------------------------------------------------

# separate from Figure_2_STATS_n76.txt (written by the stats script) so re-runs never duplicate lines
writeLines(STATS, file.path(OUTDIR, "Figure_2_PLOT_NOTES_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Figure_2_n76.pdf"))
