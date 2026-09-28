# =============================================================================
# Figure 2 statistics, 584-MAG catalogue (578 bacterial tips on the bac120 tree)
#
# Global signal, correlogram, taxonomic distance calibration, picante
# clustering, Blomberg's K, PGLS lambda, phylogenetic logistic regression,
# node-level Fisher enrichment and local Moran's I. 9,999 permutations;
# correlogram 1,000 bootstraps over 80 distance classes.
#
# LFC input: the ANCOM-BC2 estimate for EVERY MAG, NA -> 0.
#
# This is slow (the Lambda permutations refit by ML each time), so it writes
# everything the figure needs to fig2/figure2_stats_n76.rds and figure2_n76.R
# only plots.
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(phylobase); library(phylosignal); library(adephylo)
  library(picante); library(phytools); library(caper); library(phylolm); library(phangorn)
  library(dplyr); library(tidyr); library(tibble); library(readr)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "fig2")
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
say <- function(...) cat(sprintf("[%s | +%.1f min] ", format(Sys.time(), "%H:%M:%S"),
                                 as.numeric(difftime(Sys.time(), t0, units = "mins"))), ..., "\n", sep = "")
REPS <- as.integer(Sys.getenv("FIG2_REPS", "9999"))
say("permutations: ", REPS)

# ---- 1. tree + per-MAG metadata ---------------------------------------------

tree <- read.tree(TREE)
tree$node.label <- NULL                     # GTDB-Tk decorations break phylobase

res <- read.csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), stringsAsFactors = FALSE)
tax <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
blank <- function(x) is.na(x) | x == "" | grepl("^[a-z]__$", x)

meta_df <- tax %>%
  transmute(Genome_ID = catalog_id, Phylum = phylum, Class = class, Order = order,
            Family = family, Genus = genus, Species = species,
            Completeness = completeness, Contamination = contamination) %>%
  left_join(res %>% select(Genome_ID = taxon, lfc_GroupCase, p_GroupCase, q_GroupCase, diff_GroupCase),
            by = "Genome_ID") %>%
  mutate(
    Tier = case_when(
      !is.na(p_GroupCase) & p_GroupCase < 0.001 & abs(lfc_GroupCase) > 2 ~ "Strong",
      !is.na(p_GroupCase) & p_GroupCase < 0.05  & abs(lfc_GroupCase) > 1 ~ "Moderate",
      TRUE ~ "NS"),
    ANCOM_Status = case_when(
      Tier != "NS" & lfc_GroupCase > 0 ~ "Enriched in Case",
      Tier != "NS" & lfc_GroupCase < 0 ~ "Enriched in Control",
      TRUE ~ "NS"),
    Is_Novel = ifelse(blank(Species) & (!blank(Genus) | !blank(Family)), "Yes", "No"),
    Log_Fold_Change = ifelse(is.na(lfc_GroupCase), 0, lfc_GroupCase),
    enriched_case   = as.integer(ANCOM_Status == "Enriched in Case"),
    Is_Novel_bin    = as.integer(Is_Novel == "Yes")
  )

missing_tips <- setdiff(tree$tip.label, meta_df$Genome_ID)
if (length(missing_tips)) stop("FATAL: tree tips without metadata: ", paste(head(missing_tips), collapse = ", "))
meta_df <- meta_df %>% filter(Genome_ID %in% tree$tip.label)
rownames(meta_df) <- meta_df$Genome_ID
meta_df$Completeness[is.na(meta_df$Completeness)]   <- mean(meta_df$Completeness, na.rm = TRUE)
meta_df$Contamination[is.na(meta_df$Contamination)] <- mean(meta_df$Contamination, na.rm = TRUE)

say(sprintf("tips: %d | metadata rows: %d | LFC NA set to 0: %d | novel: %d | enriched_case: %d",
            Ntip(tree), nrow(meta_df), sum(is.na(meta_df$lfc_GroupCase)),
            sum(meta_df$Is_Novel == "Yes"), sum(meta_df$enriched_case)))

clean_df <- as.matrix(meta_df[tree$tip.label, "Log_Fold_Change", drop = FALSE])
p4d <- phylo4d(tree, clean_df)

# ---- 2. global signal + correlogram -----------------------------------------

say("phyloSignal (Cmean, I, Lambda) ...")
set.seed(42)
global_sig <- phyloSignal(p4d, methods = c("Cmean", "I", "Lambda"), rep = REPS)
say(sprintf("  Cmean = %.4f (p = %.4g) | Moran's I = %.4f (p = %.4g) | Lambda = %.4f (p = %.4g)",
            global_sig$stat$Cmean, global_sig$pvalue$Cmean, global_sig$stat$I,
            global_sig$pvalue$I, global_sig$stat$Lambda, global_sig$pvalue$Lambda))
write.csv(global_sig, file.path(OUTDIR, "Table_S3_Global_Phylo_Stats_n76.csv"))

say("phyloCorrelogram (80 classes, 1000 bootstraps) ...")
set.seed(42)
correlogram <- phyloCorrelogram(p4d, trait = "Log_Fold_Change", ci.bs = 1000, n.points = 80)

# phyloCorrelogram returns an UNNAMED matrix. From the phylosignal 1.3.1 source:
#   res[,1] distance, res[,2:3] bootstrap CI (lower, upper), res[,4] observed I.
# MAG_tree.R renamed column 2 as "correlation", so its "Signal Limit" was the
# zero-crossing of the LOWER bound, not of the correlogram itself. Both are
# computed here; the mean-curve crossing is the one used in the figure.
corr_res <- data.frame(d.mean = correlogram$res[, 1], lower = correlogram$res[, 2],
                       upper = correlogram$res[, 3], correlation = correlogram$res[, 4])
corr_res$sig <- with(corr_res, ifelse(lower > 0 & upper > 0, "positive",
                                ifelse(lower < 0 & upper < 0, "negative", "ns")))

zero_cross <- function(x, y) {
  i <- suppressWarnings(min(which(y < 0)))
  if (is.infinite(i) || i == 1) return(NA_real_)
  x[i - 1] - y[i - 1] * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
}
exact_limit       <- zero_cross(corr_res$d.mean, corr_res$correlation)
lower_bound_limit <- zero_cross(corr_res$d.mean, corr_res$lower)
say(sprintf("  zero-crossing of mean Moran's I = %.3f | of lower CI bound = %.3f",
            exact_limit, lower_bound_limit))
say(sprintf("  distance classes significantly positive: %d, negative: %d (of %d)",
            sum(corr_res$sig == "positive"), sum(corr_res$sig == "negative"), nrow(corr_res)))

# ---- 3. taxonomic distance calibration --------------------------------------

full_dist_mat <- cophenetic(tree)
ranks_to_test <- c("Genus", "Family", "Order", "Class", "Phylum")
tax_dist_df <- bind_rows(lapply(ranks_to_test, function(r) {
  groups <- split(rownames(meta_df), meta_df[[r]])
  groups <- groups[names(groups) != "" & !is.na(names(groups))]
  d <- unlist(lapply(groups, function(m) {
    m <- intersect(m, tree$tip.label)
    if (length(m) < 2) return(NULL)
    s <- full_dist_mat[m, m]; s[lower.tri(s)]
  }))
  if (is.null(d)) NULL else data.frame(Rank = r, Distance = d)
}))
tax_dist_df$Rank <- factor(tax_dist_df$Rank, levels = ranks_to_test)
rank_medians <- tax_dist_df %>% group_by(Rank) %>% summarise(median = median(Distance), .groups = "drop")
say("  within-rank median patristic distance: ",
    paste(sprintf("%s %.3f", rank_medians$Rank, rank_medians$median), collapse = " | "))

# ---- 4. picante clustering of case-enriched MAGs ----------------------------

say("picante ses.mpd / ses.mntd ...")
valid_enriched <- intersect(meta_df$Genome_ID[meta_df$enriched_case == 1], tree$tip.label)
comm_mat <- matrix(0, nrow = 2, ncol = Ntip(tree), dimnames = list(c("Target", "Dummy"), tree$tip.label))
comm_mat[, valid_enriched] <- 1
set.seed(42)
mpd_res  <- ses.mpd(comm_mat,  full_dist_mat, null.model = "taxa.labels", runs = REPS)
mntd_res <- ses.mntd(comm_mat, full_dist_mat, null.model = "taxa.labels", runs = REPS)
picante_res <- bind_rows(
  mpd_res[1, ]  %>% transmute(Metric = "MPD (deep clustering)",  ntaxa, obs = mpd.obs,  rand.mean = mpd.rand.mean,  z = mpd.obs.z,  p = mpd.obs.p),
  mntd_res[1, ] %>% transmute(Metric = "MNTD (tip clustering)", ntaxa, obs = mntd.obs, rand.mean = mntd.rand.mean, z = mntd.obs.z, p = mntd.obs.p))
write.csv(picante_res, file.path(OUTDIR, "Table_S4_Picante_Clustering_n76.csv"), row.names = FALSE)
print(picante_res)

# ---- 5. Blomberg's K, PGLS lambda, phylogenetic logistic regression ---------

say("Blomberg's K ...")
trait_vec <- setNames(meta_df$Log_Fold_Change, meta_df$Genome_ID)
set.seed(42)
k_res <- phylosig(tree, trait_vec, method = "K", test = TRUE, nsim = REPS)
say(sprintf("  K = %.4f, p = %.4g", k_res$K, k_res$P))

say("PGLS lambda (ML) ...")
comp_dat <- comparative.data(phy = tree, data = meta_df, names.col = "Genome_ID",
                             vcv = TRUE, na.omit = FALSE)
lambda_fit <- pgls(Log_Fold_Change ~ 1, data = comp_dat, lambda = "ML")
lam    <- lambda_fit$param["lambda"]
lam_ci <- tryCatch(lambda_fit$param.CI$lambda$ci.val, error = function(e) c(NA, NA))
say(sprintf("  Pagel's lambda = %.3f, 95%% CI [%.3f, %.3f]", lam, lam_ci[1], lam_ci[2]))

say("phyloglm ...")
fit_phyloglm <- phyloglm(enriched_case ~ Is_Novel_bin + Completeness + Contamination,
                         phy = tree, data = meta_df, method = "logistic_MPLE")
glm_coef <- summary(fit_phyloglm)$coefficients
print(glm_coef)

capture.output(
  list(Blombergs_K = k_res, PGLS_Lambda = summary(lambda_fit), PhyloGLM = summary(fit_phyloglm)),
  file = file.path(OUTDIR, "Comparative_Models_Summary_n76.txt"))

# ---- 6. node-level clade enrichment -----------------------------------------

say("node Fisher tests ...")
enriched_vec <- setNames(meta_df$enriched_case, meta_df$Genome_ID)
ntip <- Ntip(tree); x_bg <- sum(enriched_vec == 1); n_bg <- length(enriched_vec)
node_results <- bind_rows(lapply((ntip + 1):(ntip + tree$Nnode), function(nd) {
  mags <- tree$tip.label[Descendants(tree, nd, type = "tips")[[1]]]
  x <- sum(enriched_vec[mags] == 1); n <- length(mags)
  tbl <- matrix(c(x, n - x, x_bg - x, (n_bg - n) - (x_bg - x)), nrow = 2, byrow = TRUE)
  data.frame(node = nd, n_tips = n, n_enriched = x,
             p_value = fisher.test(tbl, alternative = "greater")$p.value)
})) %>% mutate(p_adj = p.adjust(p_value, "BH")) %>% arrange(p_adj)
write.csv(node_results, file.path(OUTDIR, "Table_S5_Node_Enrichment_n76.csv"), row.names = FALSE)
say("  nodes with BH q < 0.05: ", sum(node_results$p_adj < 0.05))

# ---- 7. local Moran's I -----------------------------------------------------

say("lipaMoran ...")
set.seed(42)
lipa_results <- lipaMoran(p4d, rep = REPS, as.p4d = FALSE)

inv_dist_mat <- 1 / full_dist_mat
diag(inv_dist_mat) <- 0
W <- inv_dist_mat / rowSums(inv_dist_mat)
common_tips <- rownames(W)
lag_vector <- W %*% clean_df[common_tips, 1]

analysis_df <- data.frame(
  Genome_ID = common_tips,
  LogFC = clean_df[common_tips, 1],
  Neighbor_Lag_LogFC = as.numeric(lag_vector),
  LIPA_Index = lipa_results$lipa[common_tips, 1],
  P_Value = lipa_results$p.value[common_tips, 1]
) %>%
  mutate(Phylo_Category = case_when(
    P_Value > 0.05 ~ "Non-significant",
    LogFC > 0 & Neighbor_Lag_LogFC > 0 ~ "Phylogenetic Hotspot (Conserved Enrichment)",
    LogFC < 0 & Neighbor_Lag_LogFC < 0 ~ "Phylogenetic Coldspot (Conserved Depletion)",
    LogFC > 0 & Neighbor_Lag_LogFC < 0 ~ "Evolutionary Outlier (Unique Gain)",
    LogFC < 0 & Neighbor_Lag_LogFC > 0 ~ "Evolutionary Outlier (Unique Loss)",
    TRUE ~ "Non-significant"))
print(table(analysis_df$Phylo_Category))

write.csv(meta_df %>% left_join(analysis_df %>% select(Genome_ID, Neighbor_Lag_LogFC, LIPA_Index,
                                                       LIPA_P_Value = P_Value, Phylo_Category),
                                by = "Genome_ID"),
          file.path(OUTDIR, "Table_S2_Phylogenetic_Categorization_n76.csv"), row.names = FALSE)

# ---- 8. hand-off ------------------------------------------------------------

saveRDS(list(meta_df = meta_df, global_sig = global_sig, correlogram = correlogram,
             corr_res = corr_res, exact_limit = exact_limit, lower_bound_limit = lower_bound_limit,
             tax_dist_df = tax_dist_df,
             picante_res = picante_res, k_res = k_res, lambda = lam, lambda_ci = lam_ci,
             glm_coef = glm_coef, node_results = node_results, analysis_df = analysis_df,
             reps = REPS),
        file.path(OUTDIR, "figure2_stats_n76.rds"))

writeLines(c(
  sprintf("tips: %d (bacterial MAGs on the bac120 tree); permutations: %d", Ntip(tree), REPS),
  sprintf("Cmean = %.4f, p = %.4g", global_sig$stat$Cmean, global_sig$pvalue$Cmean),
  sprintf("Global Moran's I = %.4f, p = %.4g", global_sig$stat$I, global_sig$pvalue$I),
  sprintf("Lambda (phyloSignal) = %.4f, p = %.4g", global_sig$stat$Lambda, global_sig$pvalue$Lambda),
  sprintf("Pagel's lambda (PGLS ML) = %.3f, 95%% CI [%.3f, %.3f]", lam, lam_ci[1], lam_ci[2]),
  sprintf("Blomberg's K = %.4f, p = %.4g", k_res$K, k_res$P),
  sprintf("Correlogram zero-crossing distance (mean Moran's I) = %.3f", exact_limit),
  sprintf("  (zero-crossing of the lower CI bound, as MAG_tree.R computed it = %.3f)", lower_bound_limit),
  sprintf("Distance classes significantly positive: %d, negative: %d, of %d",
          sum(corr_res$sig == "positive"), sum(corr_res$sig == "negative"), nrow(corr_res)),
  sprintf("Within-rank median patristic distance: %s",
          paste(sprintf("%s %.3f", rank_medians$Rank, rank_medians$median), collapse = "; ")),
  sprintf("phyloglm: %s", paste(sprintf("%s beta=%.4f p=%.3g", rownames(glm_coef), glm_coef[, 1], glm_coef[, ncol(glm_coef)]), collapse = "; ")),
  sprintf("Picante: %s", paste(sprintf("%s z=%.2f p=%.4g", picante_res$Metric, picante_res$z, picante_res$p), collapse = "; ")),
  sprintf("Nodes enriched for case MAGs (BH q < 0.05): %d of %d", sum(node_results$p_adj < 0.05), nrow(node_results)),
  sprintf("LIPA categories: %s", paste(names(table(analysis_df$Phylo_Category)), table(analysis_df$Phylo_Category), sep = " = ", collapse = "; "))
), file.path(OUTDIR, "Figure_2_STATS_n76.txt"))
say("done")
