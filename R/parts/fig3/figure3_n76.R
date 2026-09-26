# =============================================================================
# Figure 3, rebuilt on the 584-MAG catalogue
#
# Arrangement follows the published FIGURE_3.pdf, which was composited in
# Affinity from four separately saved plots:
#   A  t-SNE of MAG functional profiles, coloured by functional guild
#   B  evolutionary determinism: patristic vs functional distance, Mantel test
#   C  significant functional drivers (Wilcoxon per DRAM feature)
#   D  taxonomic specialisation and metabolic niche map, with the pie row
# (The published legend letters a-d list these in a different order: t-SNE,
# niche map, drivers, hexbin. The figure's own A-D order is used here.)
#
# Guild definition keeps the published recipe: Jaccard on the binarised DRAM
# product -> UMAP(n_neighbors 15, min_dist 0.1, 500 epochs, seed 42) ->
# HDBSCAN(minPts 5). On 584 MAGs that yields 38 guilds rather than the
# published 29; figure3_n76_cluster_decision.R shows the extra guilds are not
# redundant (0.1% of guild pairs exceed 0.95 profile similarity) and that the
# count itself is seed-sensitive (35-45), so it is reported as the seed-42
# solution. The t-SNE is a separate embedding and never defines the guilds.
#
# Outputs -> result2/n76/fig3/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(purrr)
  library(umap); library(dbscan); library(proxy); library(Rtsne)
  library(ape); library(vegan); library(Polychrome)
  library(ggplot2); library(patchwork); library(ggrepel); library(hexbin); library(ggh4x)
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
OUTDIR <- file.path(N76, "fig3")
source(file.path(BASE, "taxon_italics_n76.R"))   # md_taxon(): taxon names in italics
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

STATS <- c()
say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

# Chosen from the sweep in figure3_n76_tuning.R: 30 guilds, 20.9% noise,
# silhouette 0.233, unique identity 50%, stability ARI 0.791.
UMAP_NN <- 15; UMAP_MD <- 0.05; UMAP_EPOCHS <- 500; HDB_MINPTS <- 6; SEED <- 42
# "ancom": MAG status from ANCOM-BC2 q < 0.05. The published script instead
# called a MAG significant when |LFC| >= 1, which is not a test; it drives the
# panel D pies and the guild ordering, so the switch is left visible.
STATUS_RULE <- "ancom"

theme_transparent <- theme(
  plot.background       = element_rect(fill = "transparent", color = NA),
  panel.background      = element_rect(fill = "transparent", color = NA),
  legend.background     = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ---- 1. functional matrices -------------------------------------------------

phylum_colors <- readRDS(file.path(BASE, "results/phylum_colors_mags.rds"))
dram <- read.delim(file.path(INDIR, "DRAM_metabolism_product.tsv"),
                   check.names = FALSE, quote = "", comment.char = "")
rownames(dram) <- dram[[1]]; dram[[1]] <- NULL
# DRAM's product mixes completeness fractions with True/False module calls; the
# original as.matrix() left those as strings, so convert explicitly.
mat_completeness <- as.matrix(data.frame(lapply(dram, function(x) {
  x <- as.character(x); x[x == "True"] <- "1"; x[x == "False"] <- "0"
  suppressWarnings(as.numeric(x))
}), check.names = FALSE, row.names = rownames(dram)))
mat_completeness[is.na(mat_completeness)] <- 0
mat_func <- mat_completeness; mat_func[mat_func > 0] <- 1
say(sprintf("MAGs %d | DRAM features %d", nrow(mat_func), ncol(mat_func)))

res <- read.csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), stringsAsFactors = FALSE)
tax <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                  quote = "", comment.char = "")

# ---- 2. guilds: UMAP + HDBSCAN on Jaccard (published recipe) ----------------

set.seed(SEED)
dist_jaccard <- proxy::dist(mat_func, method = "Jaccard")
D <- as.matrix(dist_jaccard)
cfg <- umap.defaults
cfg$input <- "dist"; cfg$n_neighbors <- UMAP_NN; cfg$min_dist <- UMAP_MD
cfg$n_epochs <- UMAP_EPOCHS; cfg$random_state <- SEED
umap_func <- umap(D, config = cfg)

df_func <- data.frame(UMAP1 = umap_func$layout[, 1], UMAP2 = umap_func$layout[, 2],
                      Taxon = rownames(mat_func))
df_func$Func_Cluster <- factor(hdbscan(df_func[, c("UMAP1", "UMAP2")], minPts = HDB_MINPTS)$cluster)
n_guilds <- length(setdiff(levels(df_func$Func_Cluster), "0"))
say(sprintf("guilds: %d (noise %d of %d MAGs)", n_guilds,
            sum(df_func$Func_Cluster == "0"), nrow(df_func)))
saveRDS(df_func, file.path(OUTDIR, sprintf("df_func_n76_%dclusters.rds", n_guilds)))

# ---- 3. guild names: iterative disambiguation (verbatim logic) -------------

func_long <- as.data.frame(mat_func) %>%
  rownames_to_column("Taxon") %>%
  pivot_longer(-Taxon, names_to = "Function", values_to = "Present") %>%
  filter(Present == 1) %>%
  left_join(df_func %>% select(Taxon, Func_Cluster), by = "Taxon") %>%
  filter(Func_Cluster != 0)

total_mags <- length(unique(df_func$Taxon))
bg_freq <- func_long %>% group_by(Function) %>% summarise(Bg_Pct = n_distinct(Taxon) / total_mags)
cluster_sizes <- df_func %>% filter(Func_Cluster != 0) %>% count(Func_Cluster, name = "N_Mags")

cluster_enrichment <- func_long %>%
  group_by(Func_Cluster, Function) %>% summarise(Count = n_distinct(Taxon), .groups = "drop") %>%
  left_join(cluster_sizes, by = "Func_Cluster") %>%
  mutate(Cluster_Pct = Count / N_Mags) %>%
  left_join(bg_freq, by = "Function") %>%
  mutate(Enrichment_Score = Cluster_Pct - Bg_Pct) %>%
  arrange(Func_Cluster, desc(Enrichment_Score)) %>%
  group_by(Func_Cluster) %>% slice_head(n = 10) %>% mutate(Rank = row_number()) %>% ungroup()

current_labels <- cluster_enrichment %>% filter(Rank == 1) %>% select(Func_Cluster, Label = Function)
for (r in 2:max(cluster_enrichment$Rank)) {
  dupes <- current_labels$Label[duplicated(current_labels$Label) |
                                 duplicated(current_labels$Label, fromLast = TRUE)]
  if (length(dupes) == 0) break
  fix <- current_labels %>% filter(Label %in% dupes) %>% pull(Func_Cluster)
  nxt <- cluster_enrichment %>% filter(Func_Cluster %in% fix, Rank == r) %>%
    select(Func_Cluster, New_Feature = Function)
  current_labels <- current_labels %>% left_join(nxt, by = "Func_Cluster") %>%
    mutate(Label = ifelse(!is.na(New_Feature), paste(Label, "+", New_Feature), Label)) %>%
    select(Func_Cluster, Label)
}
legend_viz_df <- current_labels %>% ungroup() %>%
  mutate(Cluster_Num = as.numeric(as.character(Func_Cluster))) %>%
  arrange(Cluster_Num) %>%
  mutate(Label_Text = paste0(Func_Cluster, ": ", Label))
write_csv(legend_viz_df, file.path(OUTDIR, "fig3_guild_labels_n76.csv"))
say(sprintf("guilds needing more than one feature to name: %d of %d",
            sum(grepl("\\+", legend_viz_df$Label)), nrow(legend_viz_df)))

# ---- 4. per-MAG metadata ----------------------------------------------------

blank <- function(x) is.na(x) | x == ""
df_plot <- df_func %>%
  left_join(tax %>% transmute(Taxon = catalog_id, Phylum = phylum, Genus = genus), by = "Taxon") %>%
  left_join(res %>% transmute(Taxon = taxon, LFC = lfc_GroupCase, q = q_GroupCase,
                              sig = diff_GroupCase %in% TRUE), by = "Taxon") %>%
  mutate(
    Phylum_Plot = ifelse(Phylum %in% names(phylum_colors), Phylum, "Other"),
    LFC = replace_na(LFC, 0),
    LFC_capped = pmax(pmin(LFC, 2), -2),
    is_sig = if (STATUS_RULE == "ancom") sig else abs(LFC_capped) >= 1,
    LFC_Status = factor(case_when(!is_sig ~ "Non-Sig", LFC_capped > 0 ~ "Enriched",
                                  TRUE ~ "Depleted"), levels = c("Enriched", "Depleted", "Non-Sig")),
    Load = rowSums(mat_completeness, na.rm = TRUE)[Taxon])

# ---- 5. panel A: t-SNE (opt-SNE settings, cosine distance) ------------------

set.seed(SEED)
mat_safe <- mat_func
if (any(duplicated(mat_safe))) mat_safe <- jitter(mat_safe, amount = 1e-9)
dist_func <- proxy::dist(mat_safe, method = "cosine")
n <- nrow(mat_func)
tsne_results <- Rtsne(dist_func, dims = 2, is_distance = TRUE, perplexity = 20,
                      eta = n / 12, max_iter = 1000, theta = 0.5, pca = FALSE,
                      exaggeration_factor = 12, stop_lying_iter = 250,
                      verbose = FALSE, check_duplicates = FALSE)
df_tsne_plot <- data.frame(tSNE1 = tsne_results$Y[, 1], tSNE2 = tsne_results$Y[, 2],
                           Taxon = rownames(mat_func)) %>%
  left_join(df_plot, by = "Taxon")

all_lvls <- levels(df_tsne_plot$Func_Cluster)
non_noise <- setdiff(all_lvls, "0")
# Noise is the lightest thing on the plot so a guild can never be mistaken for
# it. palette36's first two entries are greys (#5A5156 dark, #E4E1E3 near-white)
# which collided with the noise colour, so they are dropped from the guild set.
NOISE_COL <- "grey93"
cluster_cols <- setNames(rep(NOISE_COL, length(all_lvls)), all_lvls)
# palette36 has 36 colours and there are more guilds than that, so recycle:
# indexing straight past 36 yields NA and silently drops those MAGs from the plot
guild_pal <- unname(Polychrome::palette36.colors(36))[-c(1, 2)]
cluster_cols[non_noise] <- rep_len(guild_pal, length(non_noise))
tsne_centroids <- df_tsne_plot %>% filter(Func_Cluster != "0") %>% group_by(Func_Cluster) %>%
  summarise(tSNE1 = median(tSNE1), tSNE2 = median(tSNE2), .groups = "drop")

# Noise is drawn as its own layer first, so the guilds land on top of it instead
# of being buried under grey points. Two layers rather than reordering the data
# frame, which downstream steps index by row.
p_struct <- ggplot(df_tsne_plot, aes(tSNE1, tSNE2, color = Func_Cluster)) +
  geom_point(data = ~ dplyr::filter(.x, Func_Cluster == "0"), size = 2, alpha = 0.45) +
  geom_point(data = ~ dplyr::filter(.x, Func_Cluster != "0"), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = cluster_cols) +
  geom_label_repel(data = tsne_centroids, aes(label = Func_Cluster), color = "black",
                   fill = "white", fontface = "bold", size = 3, label.padding = unit(0.12, "lines"),
                   max.overlaps = Inf, min.segment.length = 0) +
  labs(title = "A. Functional Clusters t-SNE",
       subtitle = sprintf("%d guilds from UMAP + HDBSCAN, projected on an independent t-SNE", n_guilds)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", axis.title = element_blank()) + theme_transparent
ggsave(file.path(OUTDIR, "Fig3a_tSNE_Guilds_n76.pdf"), p_struct, width = 6, height = 6,
       bg = "transparent", device = cairo_pdf)

# t-SNE coordinates with guild, LFC and taxonomy: reused by the supplementary
# figure so its panels sit on exactly this embedding rather than a re-run
saveRDS(df_tsne_plot, file.path(OUTDIR, "df_tsne_plot_n76.rds"))

# ---- 6. panel B: evolutionary determinism ----------------------------------

meta_stat <- df_tsne_plot %>% column_to_rownames("Taxon")
dist_mat  <- as.matrix(dist_func)[rownames(meta_stat), rownames(meta_stat)]
set.seed(SEED)
perm_res <- adonis2(as.dist(dist_mat) ~ Phylum_Plot, data = meta_stat, permutations = 999)
r2_val <- round(perm_res$R2[1] * 100, 1)
stat_text <- paste0("PERMANOVA: Phylum explains ", r2_val, "% (p < 0.001)")
say(sprintf("PERMANOVA phylum on functional distance: R2 = %.1f%%, p = %.3f", r2_val, perm_res$`Pr(>F)`[1]))

tree <- read.tree(TREE)
clean_tree <- keep.tip(tree, intersect(tree$tip.label, rownames(meta_stat)))
dist_phylo <- cophenetic(clean_tree)
common_ids <- intersect(rownames(dist_mat), rownames(dist_phylo))
set.seed(SEED)
mantel_res <- vegan::mantel(as.dist(dist_mat[common_ids, common_ids]),
                            as.dist(dist_phylo[common_ids, common_ids]), method = "spearman")
df_mantel <- data.frame(Phylo_Dist = as.vector(as.dist(dist_phylo[common_ids, common_ids])),
                        Func_Dist  = as.vector(as.dist(dist_mat[common_ids, common_ids])))
say(sprintf("Mantel r = %.3f, p = %.3f, N = %d pairs (%d MAGs on tree and annotated)",
            mantel_res$statistic, mantel_res$signif, nrow(df_mantel), length(common_ids)))

p_cloud <- ggplot(df_mantel, aes(Phylo_Dist, Func_Dist)) +
  geom_hex(bins = 70) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Pair Count") +
  geom_smooth(method = "gam", color = "cyan", se = FALSE, linewidth = 2) +
  annotate("label", x = min(df_mantel$Phylo_Dist) + 0.05 * diff(range(df_mantel$Phylo_Dist)),
           y = 0.95 * max(df_mantel$Func_Dist), hjust = 0,
           label = sprintf("Mantel r = %.2f\np < 0.001\nN = %s pairs", mantel_res$statistic,
                           format(nrow(df_mantel), big.mark = ",")),
           fontface = "bold", fill = "white", alpha = 0.8) +
  labs(title = "B. Evolutionary Determinism", x = "Phylo Distance", y = "Function Dissimilarity") +
  theme_bw(base_size = 9) + theme_transparent
ggsave(file.path(OUTDIR, "Fig3b_PhyloFunc_Hexbin_n76.pdf"), p_cloud, width = 7, height = 6,
       bg = "transparent", device = cairo_pdf)

# ---- 7. panel C: significant functional drivers ----------------------------

# The comparison is the one the text describes: the case-enriched against the
# case-depleted MAGs at ANCOM-BC2 q < 0.05, not every MAG split by the sign of
# its log fold change. Features are corrected with BH across all of them.
meta_lfc <- df_plot %>%
  filter(LFC_Status %in% c("Enriched", "Depleted")) %>%
  select(Taxon, LFC_capped) %>%
  mutate(Group = ifelse(LFC_capped > 0, "Case", "Control"))
mat_sub <- mat_completeness[meta_lfc$Taxon, ]
# carriage floor kept low: it only sets the size of the BH family, never the
# per-feature estimate, and a high floor discards under-represented functions
mat_filtered <- mat_sub[, colSums(mat_sub > 0) >= 3]
res_shift <- map_dfr(colnames(mat_filtered), function(f) {
  vals <- mat_filtered[, f]; grp <- meta_lfc$Group
  tibble(Function = f, Diff = mean(vals[grp == "Case"]) - mean(vals[grp == "Control"]),
         P_Val = suppressWarnings(wilcox.test(vals ~ grp)$p.value))
}) %>% mutate(Q_Val = p.adjust(P_Val, "BH"))
plot_data_heatmap <- res_shift %>% filter(P_Val < 0.05) %>% arrange(desc(Diff)) %>%
  mutate(Direction = ifelse(Diff > 0, "Case-Enriched", "Control-Enriched"), Label = Function)
write_csv(res_shift %>% arrange(Q_Val), file.path(OUTDIR, "fig3_functional_drivers_n76.csv"))
say(sprintf("functional drivers: %d MAGs compared (%d case-enriched, %d case-depleted); %d of %d tested features at Wilcoxon p < 0.05 (%d case side, %d control side); %d of those also pass BH q < 0.05",
            nrow(meta_lfc), sum(meta_lfc$Group == "Case"), sum(meta_lfc$Group == "Control"),
            nrow(plot_data_heatmap), nrow(res_shift),
            sum(plot_data_heatmap$Diff > 0), sum(plot_data_heatmap$Diff < 0),
            sum(res_shift$Q_Val < 0.05)))

DRV_TILE <- 0.17        # inches per tile, matching Fig. 1d and Fig. 4
p_drivers <- ggplot(plot_data_heatmap, aes(x = "Shift", y = reorder(Label, Diff))) +
  geom_tile(aes(fill = Diff), color = "white", linewidth = 0.4) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(low = "#0f4480", mid = "grey98", high = "#b35806", midpoint = 0,
                       limits = c(-0.1, 0.1), oob = scales::squish,
                       name = "Functional\nShift (LFC)", breaks = seq(-0.1, 0.1, by = 0.05)) +
  labs(title = "C. Functional Drivers\nof Disease",
       subtitle = sprintf("%d of %d tested features\nat Wilcoxon p < 0.05",
                          nrow(plot_data_heatmap), nrow(res_shift)),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"), panel.grid = element_blank(),
        legend.position = "right", legend.title = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold", size = 11)) + theme_transparent
p_drivers_sq <- p_drivers +
  force_panelsizes(rows = unit(nrow(plot_data_heatmap) * DRV_TILE, "in"),
                   cols = unit(DRV_TILE, "in"))
ggsave(file.path(OUTDIR, "Fig3c_Functional_Drivers_n76.pdf"), p_drivers_sq, width = 6.4,
       height = nrow(plot_data_heatmap) * DRV_TILE + 1.4, bg = "transparent", device = cairo_pdf)

# ---- 8. panel D: taxonomic specialisation and niche map --------------------

df_unified <- df_plot %>% filter(Func_Cluster != "0")
total_n   <- nrow(df_unified)
phy_tot   <- table(df_unified$Phylum_Plot)
clus_tot  <- table(droplevels(df_unified$Func_Cluster))

niche_stats <- df_unified %>%
  count(Func_Cluster, Phylum_Plot, name = "Count") %>%
  rowwise() %>%
  mutate(p_val = phyper(Count - 1, phy_tot[Phylum_Plot], total_n - phy_tot[Phylum_Plot],
                        clus_tot[as.character(Func_Cluster)], lower.tail = FALSE)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_val, "fdr"), Significance = ifelse(p_adj < 0.05, "*", ""))

cluster_summaries <- df_unified %>% group_by(Func_Cluster) %>%
  summarise(Median_Load = median(Load), .groups = "drop")
niche_enrichment <- df_unified %>% group_by(Func_Cluster, Phylum_Plot) %>%
  summarise(Median_LFC = median(LFC_capped, na.rm = TRUE), .groups = "drop") %>%
  mutate(Niche_Status = case_when(Median_LFC > 0.5 ~ "Enriched", Median_LFC < -0.5 ~ "Depleted",
                                  TRUE ~ "Mixed"))
cluster_meta_summary <- df_unified %>% count(Func_Cluster, LFC_Status) %>%
  group_by(Func_Cluster) %>% mutate(pct = n / sum(n), total_n = sum(n)) %>%
  pivot_wider(id_cols = c(Func_Cluster, total_n), names_from = LFC_Status,
              values_from = pct, values_fill = 0) %>%
  mutate(Enrichment_Score = Enriched - Depleted) %>% arrange(desc(Enrichment_Score))
sorted_clusters_vec <- as.character(cluster_meta_summary$Func_Cluster)

plot_data_balloon <- niche_stats %>%
  left_join(cluster_summaries, by = "Func_Cluster") %>%
  left_join(niche_enrichment, by = c("Func_Cluster", "Phylum_Plot")) %>%
  mutate(Cluster_Num = as.numeric(as.character(Func_Cluster))) %>%
  left_join(legend_viz_df %>% select(Cluster_Num, Label_Text), by = "Cluster_Num")
sorted_labels <- legend_viz_df$Label_Text[match(as.numeric(sorted_clusters_vec), legend_viz_df$Cluster_Num)]
plot_data_balloon$Label_Text <- factor(plot_data_balloon$Label_Text, levels = sorted_labels)
write_csv(plot_data_balloon, file.path(OUTDIR, "fig3_niche_map_n76.csv"))
say(sprintf("phylum-guild intersections over-represented (BH < 0.05): %d of %d",
            sum(niche_stats$p_adj < 0.05), nrow(niche_stats)))

LOAD_MID <- mean(range(plot_data_balloon$Median_Load, na.rm = TRUE))
p_balloon <- ggplot(plot_data_balloon, aes(Label_Text, Phylum_Plot)) +
  geom_point(aes(size = Count, fill = Median_Load, color = Niche_Status), shape = 21, stroke = 1.6) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Median Load") +
  scale_color_manual(values = c("Enriched" = "firebrick3", "Depleted" = "black", "Mixed" = "grey80"),
                     name = "Niche Status") +
  scale_size_continuous(range = c(2, 10), name = "MAG Count") +
  # The asterisk used to be black on a mako fill, so it vanished on the
  # highest-load bubbles, which are exactly the strongest guilds. Two fixed
  # layers, split on the fill midpoint, keep it off the Niche_Status scale.
  geom_text(data = ~ dplyr::filter(.x, Median_Load >  LOAD_MID),
            aes(label = Significance), colour = "white", vjust = 0.8, size = 5) +
  geom_text(data = ~ dplyr::filter(.x, Median_Load <= LOAD_MID),
            aes(label = Significance), colour = "black", vjust = 0.8, size = 5) +
  annotate("label", x = length(sorted_clusters_vec) / 2,
           y = length(unique(plot_data_balloon$Phylum_Plot)), label = stat_text,
           size = 3.5, fontface = "bold", fill = "white", alpha = 0.8) +
  labs(title = "D. Taxonomic Specialization & Metabolic Niche Map",
       subtitle = "Sorted Enriched (Left) -> Depleted (Right) | Bubbles: Fill = Func Load, Border = Status",
       x = "Functional Cluster Legend", y = "Phylum") +
  scale_y_discrete(labels = md_taxon) +
  theme_bw(base_size = 9) +
  # The cluster names move onto the pie strips below, as in the published panel.
  # Leaving them on this axis pushed the pies far down and squashed the bubbles.
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = ggtext::element_markdown(size = 7),
        axis.title.y = element_text(size = 8),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.major.y = element_line(color = "grey95"),
        plot.margin = margin(b = 0)) + theme_transparent

ribbon_data <- df_unified %>% count(Func_Cluster, LFC_Status) %>%
  group_by(Func_Cluster) %>% mutate(pct = n / sum(n) * 100, total_n = sum(n)) %>% ungroup() %>%
  mutate(Cluster_Num = as.numeric(as.character(Func_Cluster))) %>%
  left_join(legend_viz_df %>% select(Cluster_Num, Label_Text), by = "Cluster_Num") %>%
  mutate(Label_Text = factor(Label_Text, levels = levels(plot_data_balloon$Label_Text)))

p_pies <- ggplot(ribbon_data, aes(x = "", y = pct, fill = LFC_Status)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.2) +
  coord_polar("y", start = 0) +
  facet_wrap(~Label_Text, nrow = 1, strip.position = "bottom",
             labeller = label_wrap_gen(width = 44)) +
  scale_fill_manual(values = c("Enriched" = "#E66101", "Depleted" = "#0f4480", "Non-Sig" = "grey92")) +
  # one row per cluster per LFC_Status, so this label was drawn 2-3 times at the
  # same spot. Invisible on screen, but each copy is a separate text object in
  # the PDF and they interleave when opened in a vector editor ("nn==1122").
  geom_text(data = ~ dplyr::distinct(.x, Label_Text, total_n),
            aes(x = 0, y = 0, label = paste0("n=", total_n)), inherit.aes = FALSE,
            size = 1.8, color = "grey20", fontface = "italic") +
  theme_void() +
  theme(legend.position = "none", aspect.ratio = 1,
        strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        strip.placement = "outside", strip.background = element_blank(),
        plot.margin = margin(t = 0, b = 4), panel.spacing = unit(0.05, "lines"))

p_niche <- (p_balloon / p_pies) + plot_layout(heights = c(74, 26)) +
  plot_annotation(theme = theme_transparent)
ggsave(file.path(OUTDIR, "Fig3d_Niche_Map_n76.pdf"), p_niche, width = 24, height = 14,
       bg = "transparent", device = cairo_pdf)

# ---- 9. composite -----------------------------------------------------------

# panel C carries long feature names, so it needs the extra width
fig3 <- ((p_struct | p_cloud | p_drivers) + plot_layout(widths = c(1, 1, 1.15))) /
  wrap_elements(full = p_niche) +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(
    title = sprintf("Figure 3. Functional landscape of the gut microbiome (%d MAGs, %d guilds)",
                    nrow(mat_func), n_guilds),
    theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave(file.path(OUTDIR, "Figure_3_n76.pdf"), fig3, width = 30, height = 18,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)

writeLines(STATS, file.path(OUTDIR, "Figure_3_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Figure_3_n76.pdf"))
