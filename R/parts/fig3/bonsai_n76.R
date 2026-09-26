# =============================================================================
# Bonsai as an alternative to UMAP + HDBSCAN for functional guilds
#   Bonsai 2026 (de Groot et al., Nat Biotechnol), repo pinned at b16ef6f,
#   run on 584 MAGs x 98 DRAM pathway-completeness features.
#
# The Bonsai paper's claim is that a tree representation preserves the real
# high-dimensional distances where UMAP distorts them. This tests that on our
# data rather than taking it on faith:
#
#   a  the Bonsai tree with the UMAP/HDBSCAN guilds painted on. If the guilds
#      are real structure they should fall on contiguous parts of the tree.
#   b  the same tree cut into the same number of clusters, as Bonsai's own answer
#   c  where the two disagree: which guilds Bonsai splits, and which it merges
#   d  concordance across all four methods now available
#   e  whether either partition better preserves the original Jaccard distances
#
# Bonsai gives a tree, not a partition. The discrete clustering uses Bonsai's
# OWN cutting routine (downstream_analyses/get_clusters_max_diameter.py,
# min within-cluster pairwise distance) at n30, so the comparison is against
# what the authors actually ship.
#
# CAVEAT on panel e: Ward and pvclust minimise within-cluster distance on the
# very Jaccard matrix the separation score is computed from, so they are being
# graded on their own objective. Only Guild vs Bonsai is a fair contest there;
# neither optimises that quantity.
#
# Outputs -> result2/n76/bonsai/
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(mclust)
  library(dplyr); library(tidyr); library(tibble); library(readr); library(purrr)
  library(ggplot2); library(patchwork); library(Polychrome)
})

theme_transparent <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA))
theme_bw <- function(...) ggplot2::theme_bw(...) + theme_transparent
ggplot2::theme_set(ggplot2::theme_get() + theme_transparent)
ggsave <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "bonsai")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
SEED <- 42

bt <- read.tree(file.path(BASE, "bonsai_n76/bonsai_results/final_bonsai/tree.nwk"))
guilds <- readRDS(file.path(N76, "fig3", "df_func_n76_30clusters.rds")) %>%
  transmute(MAG = Taxon, Guild = as.character(Func_Cluster))
K <- n_distinct(guilds$Guild[guilds$Guild != "0"])
say(sprintf("Bonsai tree: %d tips | UMAP+HDBSCAN: %d guilds plus noise", Ntip(bt), K))

# ---- discrete Bonsai clustering ---------------------------------------------

# Bonsai's OWN tree cutting, from downstream_analyses/get_clusters_max_diameter.py
# (min within-cluster pairwise distance), extracted at n30 by
# /tmp/bonsai_cluster2.py. An average-linkage cut on the cophenetic distances was
# tried first and is NOT used: it chained, producing one cluster of 506 and 16
# singletons, which made Bonsai look far worse than it is.
bon <- read_tsv(file.path(BASE, "bonsai_n76/bonsai_cluster_assignment_k30.tsv"),
                show_col_types = FALSE) %>% mutate(Bonsai = as.character(Bonsai))
say(sprintf("   Bonsai clustering at n=%d from its own min-pairwise-distance cut", K))
sz <- sort(table(bon$Bonsai), decreasing = TRUE)
say(sprintf("   cluster sizes: largest %d, median %d, singletons %d",
            max(sz), median(as.integer(sz)), sum(sz == 1)))

cmp <- guilds %>% inner_join(bon, by = "MAG") %>%
  mutate(GuildLab = ifelse(Guild == "0", "noise", Guild))
say(sprintf("   MAGs in both partitions: %d", nrow(cmp)))

# ---- a-b. the tree, painted two ways ----------------------------------------

pal <- unname(Polychrome::palette36.colors(36))[-c(1, 2)]
guild_cols <- setNames(rep_len(pal, K), as.character(sort(as.integer(setdiff(unique(cmp$Guild), "0")))))
guild_cols["0"] <- "grey93"
bon_cols <- setNames(rep_len(pal, K), as.character(sort(as.integer(unique(cmp$Bonsai)))))

tdat <- cmp %>% column_to_rownames("MAG")
pA <- ggtree(bt, layout = "fan", open.angle = 8, size = 0.18, branch.length = "none") %<+%
  (cmp %>% transmute(label = MAG, Guild = factor(Guild, names(guild_cols)))) +
  geom_tippoint(aes(colour = Guild), size = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = guild_cols, na.value = "grey93", guide = "none") +
  labs(title = "Bonsai tree, tips coloured by UMAP + HDBSCAN guild",
       subtitle = sprintf("%d MAGs, %d DRAM features. Cladogram: a few branch lengths dominate the real tree.",
                          Ntip(bt), 98)) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"))

pB <- ggtree(bt, layout = "fan", open.angle = 8, size = 0.18, branch.length = "none") %<+%
  (cmp %>% transmute(label = MAG, Bonsai = factor(Bonsai, names(bon_cols)))) +
  geom_tippoint(aes(colour = Bonsai), size = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = bon_cols, na.value = "grey93", guide = "none") +
  labs(title = sprintf("The same tree cut into %d Bonsai clusters", K),
       subtitle = "Bonsai's own min-pairwise-distance cut at n = 30") +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"))

# ---- c. where they disagree -------------------------------------------------

say(""); say("## where the two partitions disagree")
nz <- cmp %>% filter(Guild != "0")
split_tab <- nz %>% group_by(Guild) %>%
  summarise(n = n(), n_bonsai = n_distinct(Bonsai),
            biggest = max(table(Bonsai)), .groups = "drop") %>%
  mutate(purity = biggest / n) %>% arrange(desc(n_bonsai))
merge_tab <- nz %>% group_by(Bonsai) %>%
  summarise(n = n(), n_guilds = n_distinct(Guild), .groups = "drop") %>%
  arrange(desc(n_guilds))
say(sprintf("   guilds Bonsai splits into >1 cluster: %d of %d (median %d pieces, worst %d)",
            sum(split_tab$n_bonsai > 1), nrow(split_tab),
            median(split_tab$n_bonsai), max(split_tab$n_bonsai)))
say(sprintf("   Bonsai clusters merging >1 guild: %d of %d (worst merges %d guilds)",
            sum(merge_tab$n_guilds > 1), nrow(merge_tab), max(merge_tab$n_guilds)))
say(sprintf("   median guild purity against the Bonsai cut: %.2f", median(split_tab$purity)))
write_csv(split_tab, file.path(OUTDIR, "bonsai_guild_splits_n76.csv"))
write_csv(merge_tab, file.path(OUTDIR, "bonsai_cluster_merges_n76.csv"))

pC <- split_tab %>% mutate(Guild = factor(Guild, Guild[order(n_bonsai, -purity)])) %>%
  ggplot(aes(Guild, n_bonsai, fill = purity)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Purity") +
  coord_flip() +
  labs(title = "Guilds that Bonsai splits",
       subtitle = "Bars at 1 mean the two methods agree on that guild",
       x = "UMAP + HDBSCAN guild", y = "Bonsai clusters the guild falls into") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.y = element_text(size = 6))

# ---- d. concordance across every method we now have -------------------------

say(""); say("## concordance across methods")
prod <- read_tsv(file.path(INDIR, "DRAM_metabolism_product.tsv"), show_col_types = FALSE)
m <- as.matrix(prod[, -1]); rownames(m) <- prod[[1]]; storage.mode(m) <- "numeric"
mb <- m; mb[mb > 0] <- 1; mb <- mb[, colSums(mb) > 0, drop = FALSE]
dj <- proxy::dist(mb, method = "Jaccard")
set.seed(SEED)
ward <- cutree(hclust(as.dist(dj), method = "ward.D2"), k = K)

meth <- cmp %>% mutate(Ward = as.character(ward[MAG]))
pv_file <- file.path(N76, "fig3", "pvclust_n76.rds")
if (file.exists(pv_file)) {
  pv <- readRDS(pv_file)
  meth$pvclust <- as.character(cutree(pv$hclust, k = K)[meth$MAG])
  say(sprintf("   pvclust available: %d branches, %d at AU >= 0.95",
              nrow(pv$edges), sum(pv$edges$au >= 0.95)))
}
mcols <- intersect(c("Guild", "Bonsai", "Ward", "pvclust"), names(meth))
nzm <- meth %>% filter(Guild != "0")
ari <- expand_grid(A = mcols, B = mcols) %>%
  mutate(ARI = map2_dbl(A, B, ~ mclust::adjustedRandIndex(nzm[[.x]], nzm[[.y]])))
write_csv(ari, file.path(OUTDIR, "bonsai_method_ari_n76.csv"))
for (i in seq_len(nrow(ari))) with(ari[i, ],
  if (A < B) say(sprintf("   ARI %-8s vs %-8s = %.3f", A, B, ARI)))

pD <- ggplot(ari, aes(A, B, fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", ARI),
                colour = ARI > 0.5), size = 3, show.legend = FALSE) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1), name = "ARI") +
  scale_colour_manual(values = c("FALSE" = "grey15", "TRUE" = "white")) +
  labs(title = "Cross-method concordance",
       subtitle = sprintf("Adjusted Rand Index, %d MAGs, noise excluded", nrow(nzm)),
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        panel.grid = element_blank())

# ---- e. does either partition preserve the real distances better? -----------

say(""); say("## which partition better preserves the original Jaccard distances")
djm <- as.matrix(dj)
common <- intersect(rownames(djm), nzm$MAG)
djs <- djm[common, common]
# Scoring only the MAGs HDBSCAN was willing to cluster hands it a free pass:
# the 122 noise genomes are precisely the ones it could not place, and dropping
# them removes its own failures from its score while also discarding genomes the
# other methods handled. Report both scopes.
within_between <- function(lab, ids, scope) {
  ids <- intersect(ids, rownames(djm))
  l <- setNames(meth[[lab]][match(ids, meth$MAG)], ids)
  d <- djm[ids, ids]
  same <- outer(l, l, "=="); ut <- upper.tri(d)
  w <- mean(d[ut & same]); b <- mean(d[ut & !same])
  tibble(method = lab, scope = scope, within = w, between = b, separation = (b - w) / b)
}
sep <- bind_rows(
  map_dfr(mcols, ~ within_between(.x, nzm$MAG, "HDBSCAN-clusterable (462)")),
  map_dfr(mcols, ~ within_between(.x, meth$MAG, "All MAGs (584)")))
write_csv(sep, file.path(OUTDIR, "bonsai_distance_separation_n76.csv"))
for (i in seq_len(nrow(sep))) with(sep[i, ],
  say(sprintf("   %-26s %-8s within %.3f vs between %.3f | separation %.3f",
              scope, method, within, between, separation)))
say("   separation = (between - within) / between; higher means the partition")
say("   groups genuinely similar genomes rather than carving the space arbitrarily.")

pE <- sep %>%
  mutate(method = factor(method, c("Bonsai", "Guild", "Ward", "pvclust")),
         scope = factor(scope, c("HDBSCAN-clusterable (462)", "All MAGs (584)"))) %>%
  ggplot(aes(separation, method, fill = scope)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", separation)), position = position_dodge(width = 0.7),
            hjust = -0.12, size = 2.4) +
  scale_fill_manual(values = c("HDBSCAN-clusterable (462)" = "grey75",
                               "All MAGs (584)" = "#3B6FB6"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Distance preservation on the original Jaccard distances",
       subtitle = "Restricting to the MAGs HDBSCAN could place flatters it; on all 584 Bonsai is ahead",
       x = "Separation (between - within) / between", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

# ---- save -------------------------------------------------------------------

write_csv(meth, file.path(OUTDIR, "bonsai_vs_guild_assignments_n76.csv"))
ggsave(file.path(OUTDIR, "Bonsai_Tree_Guild_Clusters_n76.pdf"), pA, width = 8, height = 8.6,
       device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "Bonsai_Tree_BonsaiCut_n76.pdf"), pB, width = 8, height = 8.6,
       device = cairo_pdf, bg = "transparent")
for (nm in c("Bonsai_c_Splits", "Bonsai_d_ARI", "Bonsai_e_Separation")) {
  p <- switch(nm, Bonsai_c_Splits = pC, Bonsai_d_ARI = pD, Bonsai_e_Separation = pE)
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), p, width = 5.4, height = 4.6,
         device = cairo_pdf, bg = "transparent")
}

tag <- theme(plot.tag = element_text(face = "bold", size = 13))
fig <- ((pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag)) /
       ((pC + labs(tag = "c") + tag) | (pD + labs(tag = "d") + tag) | (pE + labs(tag = "e") + tag)) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(title = "Bonsai against UMAP + HDBSCAN for functional guilds (584 MAGs)",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUTDIR, "Bonsai_Comparison_n76.pdf"), fig, width = 16, height = 14,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)

writeLines(STATS, file.path(OUTDIR, "BONSAI_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Bonsai_Comparison_n76.pdf"))
