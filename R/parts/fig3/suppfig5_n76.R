# =============================================================================
# Supplementary Figure 5, clustering-method comparison for the functional guilds
#
# Five partitions of the same 584 MAGs x 98 DRAM pathway matrix:
#   Guild    UMAP (nn 15, min_dist 0.05) + HDBSCAN (minPts 6)   -> 30 + noise
#   Bonsai   Bonsai 2026 tree, its own min-pairwise-distance cut -> 30
#   HiDeF    persistent multiscale community detection, kNN k=20 -> 30
#   Ward     hierarchical ward.D2 on Jaccard, cut at 30          -> 30
#   pvclust  the same ward.D2 with 1,000 bootstrap replicates    -> 30
#
# Two points that change the reading:
#
#   1. pvclust and Ward are the SAME partition (ARI 1.000). pvclust runs ward.D2
#      on the same distance; it contributes bootstrap support, not an
#      independent opinion. Counting them as two methods double-counts.
#
#   2. Scoring distance preservation only on the MAGs HDBSCAN chose to cluster
#      hands it a free pass, because the 122 noise genomes are exactly the ones
#      it could not place. Both scopes are reported.
#
# Outputs -> result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(mclust); library(dplyr); library(tidyr); library(tibble)
  library(readr); library(purrr); library(ggplot2); library(patchwork)
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
OUTDIR <- file.path(N76, "supp_fig")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
SEED <- 42

# ---- assemble every partition -----------------------------------------------

guilds <- readRDS(file.path(N76, "fig3", "df_func_n76_30clusters.rds")) %>%
  transmute(MAG = Taxon, Guild = as.character(Func_Cluster))
bon <- read_tsv(file.path(BASE, "bonsai_n76/bonsai_cluster_assignment_k30.tsv"),
                show_col_types = FALSE) %>% mutate(Bonsai = as.character(Bonsai))
hid <- read_tsv(file.path(BASE, "hidef_n76/hidef_cluster_assignment.tsv"),
                show_col_types = FALSE) %>% mutate(HiDeF = as.character(HiDeF))

prod <- read_tsv(file.path(INDIR, "DRAM_metabolism_product.tsv"), show_col_types = FALSE)
m <- as.matrix(prod[, -1]); rownames(m) <- prod[[1]]; storage.mode(m) <- "numeric"
mb <- m; mb[mb > 0] <- 1; mb <- mb[, colSums(mb) > 0, drop = FALSE]
dj <- proxy::dist(mb, method = "Jaccard"); djm <- as.matrix(dj)
K <- n_distinct(guilds$Guild[guilds$Guild != "0"])
set.seed(SEED)
hc <- hclust(as.dist(dj), method = "ward.D2")
ward <- cutree(hc, k = K)

meth <- guilds %>% inner_join(bon, by = "MAG") %>% inner_join(hid, by = "MAG") %>%
  mutate(Ward = as.character(ward[MAG]))
pv_file <- file.path(N76, "fig3", "pvclust_n76.rds")
pv <- if (file.exists(pv_file)) readRDS(pv_file) else NULL
if (!is.null(pv)) meth$pvclust <- as.character(cutree(pv$hclust, k = K)[meth$MAG])

mcols <- intersect(c("Guild", "Bonsai", "HiDeF", "Ward", "pvclust"), names(meth))
say(sprintf("partitions compared: %s", paste(mcols, collapse = ", ")))
for (v in mcols)
  say(sprintf("   %-8s %3d clusters%s", v, n_distinct(meth[[v]]),
              if (v == "Guild") sprintf(" (incl. %d noise MAGs)", sum(meth$Guild == "0")) else ""))
write_csv(meth, file.path(OUTDIR, "SuppFig5_all_partitions_n76.csv"))

# ---- a. cross-method concordance --------------------------------------------

nzm <- meth %>% filter(Guild != "0")
ari <- expand_grid(A = mcols, B = mcols) %>%
  mutate(ARI = map2_dbl(A, B, ~ mclust::adjustedRandIndex(nzm[[.x]], nzm[[.y]])))
write_csv(ari, file.path(OUTDIR, "SuppFig5_ari_matrix_n76.csv"))
say(""); say("## cross-method concordance (ARI, noise excluded)")
for (i in seq_len(nrow(ari))) with(ari[i, ], if (A < B) say(sprintf("   %-8s vs %-8s %.3f", A, B, ARI)))

pA <- ggplot(ari, aes(factor(A, mcols), factor(B, rev(mcols)), fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", ARI), colour = ARI > 0.55), size = 3.1, show.legend = FALSE) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1), name = "ARI") +
  scale_colour_manual(values = c("FALSE" = "grey12", "TRUE" = "white")) +
  labs(title = "A. Cross-method concordance",
       subtitle = sprintf("Adjusted Rand Index, %d MAGs, HDBSCAN noise excluded", nrow(nzm)),
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))

# ---- b. pvclust bootstrap stability -----------------------------------------

if (!is.null(pv)) {
  au <- pv$edges$au
  say(""); say(sprintf("## pvclust: %d branches, %d (%.0f%%) exceed AU 0.95 (n_boot = 1000)",
                       length(au), sum(au >= 0.95), 100 * mean(au >= 0.95)))
  pB <- ggplot(tibble(au = au), aes(au)) +
    geom_histogram(binwidth = 0.025, fill = "#3B6FB6", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.95, colour = "red", linetype = "dashed") +
    annotate("text", x = 0.95, y = Inf, label = "AU = 0.95 ", hjust = 1, vjust = 1.4,
             size = 2.6, colour = "red") +
    labs(title = "B. pvclust bootstrap stability",
         subtitle = sprintf("%d of %d branches exceed AU 0.95 (n_boot = 1000)",
                            sum(au >= 0.95), length(au)),
         x = "AU p-value (branch support under bootstrap resampling)",
         y = "Number of tree branches") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 7, colour = "grey35"))
} else {
  pB <- ggplot() + theme_void() + labs(title = "B. pvclust not available")
}

# ---- c. per-guild purity under each alternative -----------------------------

alt <- setdiff(mcols, c("Guild", "pvclust"))
purity <- map_dfr(alt, function(a) {
  nzm %>% count(Guild, .data[[a]]) %>% group_by(Guild) %>%
    summarise(purity = max(n) / sum(n), n = sum(n), .groups = "drop") %>%
    mutate(Method = a)
})
write_csv(purity, file.path(OUTDIR, "SuppFig5_guild_purity_n76.csv"))
say("")
for (a in alt) {
  p <- purity %>% filter(Method == a)
  say(sprintf("   guild purity under %-7s median %.2f | %d of %d guilds above 0.5",
              a, median(p$purity), sum(p$purity > 0.5), nrow(p)))
}
ord <- purity %>% group_by(Guild) %>% summarise(m = mean(purity), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Guild)

pC <- purity %>% mutate(Guild = factor(Guild, ord)) %>%
  ggplot(aes(Guild, purity, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Bonsai" = "#D55E00", "HiDeF" = "#5AAE61", "Ward" = "#3B6FB6"),
                    name = "Method") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "C. Per-guild purity under each alternative method",
       subtitle = "Fraction of each guild's MAGs falling into that guild's single most common alternative-method cluster",
       x = "UMAP + HDBSCAN guild", y = "Purity") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 6), legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

# ---- d. distance preservation, both scopes ----------------------------------

sep_of <- function(lab, ids, scope) {
  ids <- intersect(ids, rownames(djm))
  l <- setNames(meth[[lab]][match(ids, meth$MAG)], ids)
  d <- djm[ids, ids]; same <- outer(l, l, "=="); ut <- upper.tri(d)
  w <- mean(d[ut & same]); b <- mean(d[ut & !same])
  tibble(method = lab, scope = scope, separation = (b - w) / b)
}
sep <- bind_rows(
  map_dfr(mcols, ~ sep_of(.x, nzm$MAG,  sprintf("HDBSCAN-clusterable (%d)", nrow(nzm)))),
  map_dfr(mcols, ~ sep_of(.x, meth$MAG, sprintf("All MAGs (%d)", nrow(meth)))))
write_csv(sep, file.path(OUTDIR, "SuppFig5_distance_separation_n76.csv"))
say(""); say("## distance preservation on the original Jaccard distances")
for (i in seq_len(nrow(sep))) with(sep[i, ],
  say(sprintf("   %-26s %-8s separation %.3f", scope, method, separation)))
say("   Ward and pvclust minimise this very quantity, so they are graded on their")
say("   own objective. Guild vs Bonsai vs HiDeF is the fair contest.")

pD <- sep %>%
  mutate(method = factor(method, mcols),
         scope = factor(scope, unique(sep$scope))) %>%
  ggplot(aes(separation, method, fill = scope)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = sprintf("%.3f", separation)), position = position_dodge(width = 0.72),
            hjust = -0.12, size = 2.4) +
  scale_fill_manual(values = setNames(c("grey75", "#3B6FB6"), levels(factor(sep$scope, unique(sep$scope)))),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(title = "D. Distance preservation on the original Jaccard distances",
       subtitle = "Restricting to MAGs HDBSCAN could place flatters it; on all 584 it falls behind Bonsai and HiDeF",
       x = "Separation (between - within) / between", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6.5))

# ---- save -------------------------------------------------------------------

pl <- list(SuppFig5a_ARI = pA, SuppFig5b_pvclust_AU = pB,
           SuppFig5c_Guild_Purity = pC, SuppFig5d_Separation = pD)
for (nm in names(pl))
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), pl[[nm]],
         width = if (nm == "SuppFig5c_Guild_Purity") 9 else 5.6, height = 4.6,
         device = cairo_pdf, bg = "transparent")

tag <- theme(plot.tag = element_text(face = "bold", size = 13))
supp5 <- ((pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag)) /
         (pC + labs(tag = "c") + tag) /
         (pD + labs(tag = "d") + tag) +
  plot_layout(heights = c(1, 1, 0.9)) +
  plot_annotation(
    title = sprintf("Clustering-method comparison for the functional guilds (%d MAGs)", nrow(meth)),
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
# Scaled by 0.636 to fit A4 portrait; the aspect ratio is
# unchanged, so no panel is stretched relative to the others.
ggsave(file.path(OUTDIR, "Supplementary_Figure_5_n76.pdf"), supp5, width = 8.27, height = 8.91,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)

writeLines(STATS, file.path(OUTDIR, "SuppFig5_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Supplementary_Figure_5_n76.pdf"))
