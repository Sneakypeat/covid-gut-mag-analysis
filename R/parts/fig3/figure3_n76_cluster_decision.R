# =============================================================================
# Figure 3: choosing the guild resolution, with evidence
#
# Question: UMAP n_neighbors 15, min_dist 0.1 with HDBSCAN
# minPts 5 gives 38 guilds on this 584-MAG catalogue. Does 38 separate
# functions better, or does it dilute them into near-duplicate guilds?
#
# Four things decide that, none of them visual:
#   stability   mean adjusted Rand index across 6 UMAP seeds. Over-split
#               solutions are unstable: the extra splits move with the seed.
#   redundancy  share of guild PAIRS whose mean functional profiles are nearly
#               identical (cosine similarity > 0.95). This is dilution measured
#               directly: two guilds that mean the same thing.
#   identity    share of guilds whose single most enriched DRAM feature is not
#               the top feature of any other guild (nameability).
#   separation  silhouette in the original Jaccard space, noise excluded.
#
# -> result2/n76/fig3/fig3_cluster_decision_n76.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(umap); library(dbscan); library(proxy); library(cluster); library(mclust)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
OUTDIR <- file.path(BASE, "result2/n76/fig3")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ---- functional matrix (same binarisation as Mags_annotation_Fig3.R) --------

dram <- read.delim(file.path(INDIR, "DRAM_metabolism_product.tsv"),
                   check.names = FALSE, quote = "", comment.char = "")
rownames(dram) <- dram[[1]]; dram[[1]] <- NULL
mat_completeness <- as.matrix(data.frame(lapply(dram, function(x) {
  x <- as.character(x); x[x == "True"] <- "1"; x[x == "False"] <- "0"
  suppressWarnings(as.numeric(x))
}), check.names = FALSE, row.names = rownames(dram)))
mat_func <- mat_completeness
mat_func[is.na(mat_func)] <- 0
mat_func[mat_func > 0] <- 1

D <- as.matrix(proxy::dist(mat_func, method = "Jaccard"))
bg_freq <- colMeans(mat_func)
say("MAGs ", nrow(mat_func), " | DRAM features ", ncol(mat_func))

embed <- function(nn, md, seed) {
  cfg <- umap.defaults
  cfg$input <- "dist"; cfg$n_neighbors <- nn; cfg$min_dist <- md
  cfg$n_epochs <- 500; cfg$random_state <- seed
  umap(D, config = cfg)$layout
}

cos_sim <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))

evaluate <- function(cl) {
  ok <- cl != 0
  gs <- sort(unique(cl[ok]))
  k  <- length(gs)
  if (k < 3) return(NULL)
  sil <- mean(cluster::silhouette(cl[ok], dmatrix = D[ok, ok])[, 3])
  # mean completeness profile per guild -> redundancy between guilds
  cent <- t(vapply(gs, function(g) colMeans(mat_completeness[cl == g, , drop = FALSE], na.rm = TRUE),
                   numeric(ncol(mat_completeness))))
  cent[is.na(cent)] <- 0
  pairs <- combn(seq_len(k), 2)
  sims <- apply(pairs, 2, function(ij) cos_sim(cent[ij[1], ], cent[ij[2], ]))
  top1 <- vapply(gs, function(g) {
    enr <- colMeans(mat_func[cl == g, , drop = FALSE]) - bg_freq
    names(sort(enr, decreasing = TRUE))[1]
  }, "")
  tibble(n_guilds = k, noise = mean(!ok), silhouette = sil,
         redundant_pairs = mean(sims > 0.95),
         max_pair_similarity = max(sims),
         unique_identity = mean(!(duplicated(top1) | duplicated(top1, fromLast = TRUE))),
         min_size = min(table(cl[ok])), median_size = median(table(cl[ok])))
}

# ---- candidates ------------------------------------------------------------

cands <- tribble(
  ~label,                      ~nn, ~md,  ~mp,
  "nn15 md0.10 minPts 5",       15, 0.10,  5,
  "nn15 md0.10 minPts 6",       15, 0.10,  6,
  "nn15 md0.10 minPts 8",       15, 0.10,  8,
  "nn15 md0.05 minPts 6",       15, 0.05,  6,
  "nn20 md0.05 minPts 6",       20, 0.05,  6,
  "nn10 md0.05 minPts 6",       10, 0.05,  6,
  "nn20 md0.05 minPts 10",      20, 0.05, 10)

SEEDS <- c(42, 7, 123, 2024, 5, 99)

rows <- list()
for (i in seq_len(nrow(cands))) {
  cd <- cands[i, ]
  say("evaluating: ", cd$label)
  labs <- lapply(SEEDS, function(s) hdbscan(embed(cd$nn, cd$md, s), minPts = cd$mp)$cluster)
  # stability: ARI between every pair of seeds, over MAGs clustered in both
  aris <- combn(length(SEEDS), 2, function(ij) {
    a <- labs[[ij[1]]]; b <- labs[[ij[2]]]
    keep <- a != 0 & b != 0
    if (sum(keep) < 20) return(NA_real_)
    mclust::adjustedRandIndex(a[keep], b[keep])
  })
  ks <- vapply(labs, function(x) length(unique(x[x != 0])), 1L)
  ev <- evaluate(labs[[1]])                       # seed 42, the reported solution
  rows[[i]] <- bind_cols(cd, ev,
                         tibble(stability_ARI = mean(aris, na.rm = TRUE),
                                guilds_across_seeds = paste(ks, collapse = "/")))
}

res <- bind_rows(rows) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  arrange(desc(stability_ARI))
write_csv(res, file.path(OUTDIR, "fig3_cluster_decision_n76.csv"))
say("")
print(as.data.frame(res %>% select(label, n_guilds, guilds_across_seeds, stability_ARI,
                                   redundant_pairs, max_pair_similarity, unique_identity,
                                   silhouette, noise, min_size)), row.names = FALSE)
