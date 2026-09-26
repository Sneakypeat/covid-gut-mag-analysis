# =============================================================================
# Rarefaction seed stability: how much of the ANCOM-BC2 call set is one random
# draw? Companion to depth_rarefaction_n76.R (which uses seed 42).
# Outputs -> result2/n76/supp_fig/depth/rarefaction_seed_stability_n76.csv
# =============================================================================
suppressPackageStartupMessages({library(phyloseq); library(vegan); library(ANCOMBC); library(dplyr); library(readr)})
N76 <- file.path(Sys.getenv("COVID_MAG_BASE", unset = getwd()), "result2/n76")
ps  <- readRDS(file.path(N76, "ps_mags_n76.rds"))
otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
prim <- read_csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), show_col_types = FALSE) %>%
  transmute(taxon, lfc_prim = lfc_GroupCase, sig_prim = diff_GroupCase %in% TRUE)
scale_to <- min(colSums(otu))
out <- list()
for (sd in c(1, 7, 99, 2024)) {
  set.seed(sd)
  rar <- t(rrarefy(t(otu), sample = scale_to))
  f <- ancombc2(data = phyloseq(otu_table(rar, taxa_are_rows = TRUE), sample_data(ps)),
                fix_formula = "Group", p_adj_method = "BH", pseudo_sens = TRUE, prv_cut = 0.01,
                lib_cut = 1000, group = "Group", struc_zero = FALSE, neg_lb = TRUE, alpha = 0.05, global = FALSE)
  r <- f$res %>% transmute(taxon, lfc = lfc_GroupCase, sig = diff_GroupCase %in% TRUE)
  j <- inner_join(prim, r, by = "taxon")
  ok <- complete.cases(j$lfc_prim, j$lfc)
  out[[as.character(sd)]] <- tibble(seed = sd, n_sig = sum(j$sig),
    up = sum(j$sig & j$lfc > 0), down = sum(j$sig & j$lfc < 0),
    kept = sum(j$sig_prim & j$sig), dir_held = sum(j$sig_prim & sign(j$lfc_prim) == sign(j$lfc)),
    rho = cor(j$lfc_prim[ok], j$lfc[ok], method = "spearman"))
  message(sd, " done")
}
res <- bind_rows(out)
print(as.data.frame(res))
write_csv(res, file.path(N76, "supp_fig/depth/rarefaction_seed_stability_n76.csv"))
