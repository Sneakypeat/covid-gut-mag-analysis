# =============================================================================
# Figure 3 functional guilds, 584-MAG catalogue: CLUSTER TUNING
#
# How the guilds are built:
#   guild LABELS  - Jaccard distance on the binarised DRAM product matrix
#                   -> UMAP (n_neighbors 15, min_dist 0.1, n_epochs 500, seed 42)
#                   -> HDBSCAN (minPts 5)
#   plotted LAYOUT - a separate opt-SNE on cosine distance, with those labels
#                   painted on. Tuning the t-SNE therefore never changes the
#                   guilds; only the UMAP/HDBSCAN knobs below do.
#
# Run top to bottom once for the sweep, then work interactively in section 4:
# change the four TUNE_ values and re-run that section to inspect a candidate.
#
# Outputs -> result2/n76/fig3/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(umap); library(dbscan); library(proxy); library(cluster)
  library(ggplot2); library(patchwork); library(Polychrome)
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
OUTDIR <- file.path(BASE, "result2/n76/fig3")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. functional matrix (same binarisation as the original) ---------------

dram <- read.delim(file.path(INDIR, "DRAM_metabolism_product.tsv"),
                   check.names = FALSE, quote = "", comment.char = "")
rownames(dram) <- dram[[1]]; dram[[1]] <- NULL
mat_completeness <- as.matrix(data.frame(lapply(dram, function(x) {
  x <- as.character(x)
  x[x == "True"] <- "1"; x[x == "False"] <- "0"
  suppressWarnings(as.numeric(x))
}), check.names = FALSE, row.names = rownames(dram)))
mat_func <- mat_completeness
mat_func[is.na(mat_func)] <- 0
mat_func[mat_func > 0] <- 1

# Distance uses the FULL binary matrix, as the original does. Features carried
# by every MAG still count as shared presences in Jaccard, so dropping them
# would change the distances. (All-absent features are ignored by Jaccard.)
keep_feat <- colSums(mat_func) > 0 & colSums(mat_func) < nrow(mat_func)
message(sprintf("MAGs: %d | DRAM features: %d (%d vary across MAGs)",
                nrow(mat_func), ncol(mat_func), sum(keep_feat)))

set.seed(42)
dist_jaccard <- proxy::dist(mat_func, method = "Jaccard")
D <- as.matrix(dist_jaccard)

bg_freq <- colMeans(mat_func)

# ---- 2. scoring ------------------------------------------------------------
# silhouette  - separation in the ORIGINAL Jaccard space, noise excluded
# distinct    - mean over guilds of the top-3 feature enrichment
#               (share of guild MAGs carrying a feature minus catalogue share);
#               high values mean guilds are defined by features that are
#               genuinely characteristic, which is what makes them nameable
# unique_lab  - share of guilds whose single top feature is not the top
#               feature of any other guild (the original naming loop has to
#               stack features when this is low)

score_clusters <- function(cl) {
  ok <- cl != 0
  k  <- length(unique(cl[ok]))
  if (k < 2) return(tibble(n_clusters = k, noise = mean(!ok), silhouette = NA_real_,
                           distinct = NA_real_, unique_lab = NA_real_,
                           min_size = NA_integer_, median_size = NA_real_))
  sil <- mean(cluster::silhouette(cl[ok], dmatrix = D[ok, ok])[, 3])
  tops <- lapply(sort(unique(cl[ok])), function(g) {
    enr <- colMeans(mat_func[cl == g, , drop = FALSE]) - bg_freq
    sort(enr, decreasing = TRUE)[1:3]
  })
  top1 <- vapply(tops, function(x) names(x)[1], "")
  sizes <- table(cl[ok])
  tibble(n_clusters = k, noise = mean(!ok), silhouette = sil,
         distinct = mean(vapply(tops, mean, 0)),
         unique_lab = mean(!(duplicated(top1) | duplicated(top1, fromLast = TRUE))),
         min_size = as.integer(min(sizes)), median_size = as.numeric(median(sizes)))
}

embed <- function(nn, md, epochs = 500, seed = 42) {
  cfg <- umap.defaults
  cfg$input <- "dist"; cfg$n_neighbors <- nn; cfg$min_dist <- md
  cfg$n_epochs <- epochs; cfg$random_state <- seed
  umap(D, config = cfg)$layout
}

# ---- 3. nn15 md0.10 minPts 5, then sweep -------------------------------------

message("nn15 md0.10 minPts 5 ...")
lay_pub <- embed(15, 0.1)
cl_pub  <- hdbscan(lay_pub, minPts = 5)$cluster
pub     <- score_clusters(cl_pub) %>% mutate(n_neighbors = 15, min_dist = 0.1, minPts = 5)
print(pub)

GRID_NN <- c(10, 15, 20, 30, 50)
GRID_MD <- c(0.01, 0.05, 0.1)   # umap requires min_dist > 0
GRID_MP <- c(4, 5, 6, 8, 10, 15)

sweep_file <- file.path(OUTDIR, "fig3_cluster_sweep_n76.csv")
layouts <- list()
if (!file.exists(sweep_file)) {
  rows <- list()
  for (nn in GRID_NN) for (md in GRID_MD) {
    message(sprintf("  UMAP n_neighbors=%d min_dist=%.2f", nn, md))
    lay <- embed(nn, md)
    layouts[[sprintf("%d_%.2f", nn, md)]] <- lay
    for (mp in GRID_MP) {
      cl <- hdbscan(lay, minPts = mp)$cluster
      rows[[length(rows) + 1]] <- score_clusters(cl) %>%
        mutate(n_neighbors = nn, min_dist = md, minPts = mp)
    }
  }
  sweep <- bind_rows(rows)
  write_csv(sweep, sweep_file)
  saveRDS(layouts, file.path(OUTDIR, "fig3_sweep_layouts.rds"))
} else {
  sweep   <- read_csv(sweep_file, show_col_types = FALSE)
  layouts <- readRDS(file.path(OUTDIR, "fig3_sweep_layouts.rds"))
}

# Composite rank: separation, distinctness and nameability up; noise down.
ranked <- sweep %>%
  filter(!is.na(silhouette), noise <= 0.35, min_size >= 5) %>%
  mutate(score = rank(silhouette) + rank(distinct) + rank(unique_lab) + rank(-noise)) %>%
  arrange(desc(score))
write_csv(ranked, file.path(OUTDIR, "fig3_cluster_sweep_ranked_n76.csv"))
message("\nTop 12 candidates (noise <= 35%, every guild >= 5 MAGs):")
print(as.data.frame(ranked %>% select(n_neighbors, min_dist, minPts, n_clusters, noise,
                                      silhouette, distinct, unique_lab, median_size) %>%
                      head(12)), digits = 3)

# Small multiples of the top 6
pal36 <- unname(Polychrome::palette36.colors(36))
plot_candidate <- function(nn, md, mp, title = NULL) {
  lay <- layouts[[sprintf("%d_%.2f", nn, md)]]
  cl  <- hdbscan(lay, minPts = mp)$cluster
  df  <- tibble(UMAP1 = lay[, 1], UMAP2 = lay[, 2], cl = factor(cl))
  lv  <- setdiff(levels(df$cl), "0")
  cols <- setNames(c("grey80", rep(pal36, length.out = length(lv))), c("0", lv))
  cen <- df %>% filter(cl != "0") %>% group_by(cl) %>% summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
  s <- score_clusters(cl)
  ggplot(df, aes(UMAP1, UMAP2, colour = cl)) +
    geom_point(size = 0.9, alpha = 0.75) +
    geom_text(data = cen, aes(label = cl), colour = "black", size = 2.6, fontface = "bold") +
    scale_colour_manual(values = cols, guide = "none") +
    labs(title = title %||% sprintf("nn=%d  md=%.2f  minPts=%d", nn, md, mp),
         subtitle = sprintf("%d guilds | noise %.0f%% | sil %.2f | distinct %.2f | unique %.0f%%",
                            s$n_clusters, 100 * s$noise, s$silhouette, s$distinct, 100 * s$unique_lab)) +
    theme_bw(base_size = 8) + theme(axis.text = element_blank(), axis.ticks = element_blank())
}
`%||%` <- function(a, b) if (is.null(a)) b else a

top6 <- head(ranked, 6)
p_top <- wrap_plots(c(list(plot_candidate(15, 0.1, 5, "nn15 md0.10 minPts 5 on 584 MAGs")),
                      lapply(seq_len(nrow(top6)), function(i)
                        plot_candidate(top6$n_neighbors[i], top6$min_dist[i], top6$minPts[i]))),
                    ncol = 4)
ggsave(file.path(OUTDIR, "fig3_cluster_candidates_n76.pdf"), p_top,
       width = 16, height = 8.5, device = cairo_pdf, bg = "transparent")
message("\ncandidates -> ", file.path(OUTDIR, "fig3_cluster_candidates_n76.pdf"))

# ---- 4. INTERACTIVE: inspect one candidate ----------------------------------
# Edit these four and re-run from here down.
TUNE_NN     <- 15
TUNE_MD     <- 0.1
TUNE_MINPTS <- 5
TUNE_SEED   <- 42

lay_t <- if (TUNE_SEED == 42 && !is.null(layouts[[sprintf("%d_%.2f", TUNE_NN, TUNE_MD)]]))
  layouts[[sprintf("%d_%.2f", TUNE_NN, TUNE_MD)]] else embed(TUNE_NN, TUNE_MD, seed = TUNE_SEED)
cl_t <- hdbscan(lay_t, minPts = TUNE_MINPTS)$cluster
print(score_clusters(cl_t))

# What defines each guild: top 3 enriched DRAM features
guild_defs <- lapply(sort(unique(cl_t[cl_t != 0])), function(g) {
  enr <- sort(colMeans(mat_func[cl_t == g, , drop = FALSE]) - bg_freq, decreasing = TRUE)[1:3]
  tibble(guild = g, n_mags = sum(cl_t == g),
         top_features = paste(sprintf("%s (+%.0f%%)", names(enr), 100 * enr), collapse = " | "))
})
print(bind_rows(guild_defs), n = Inf, width = Inf)

# When happy, save the labels for the figure script:
# saveRDS(tibble(Taxon = rownames(mat_func), Func_Cluster = factor(cl_t),
#                UMAP1 = lay_t[, 1], UMAP2 = lay_t[, 2]),
#         file.path(OUTDIR, sprintf("df_func_n76_nn%d_md%.2f_mp%d.rds", TUNE_NN, TUNE_MD, TUNE_MINPTS)))
