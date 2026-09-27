# =============================================================================
# Sequencing depth: audit and test of consequence (76 donors, 584 MAGs)
#
# The four cohorts were sequenced to different depths and the arms are not
# matched, so this asks whether depth, rather than infection, could produce the
# community and MAG-level results. Everything is recomputed after rarefying
# every donor to the shallowest library:
#   1. depth per cohort and per arm
#   2. depth as a covariate on Shannon, richness and distance to the group spatial median
#   3. the same endpoints on the rarefied table
#   4. ANCOM-BC2 re-run on the rarefied table, against the primary result
#
# Outputs -> result2/n76/supp_fig/depth/
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(ANCOMBC)
  library(dplyr); library(readr); library(tibble)
})
set.seed(42)

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "supp_fig/depth")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

STATS <- character()
say <- function(s) { STATS <<- c(STATS, s); message(s) }   # callers format their own text

ps  <- readRDS(file.path(N76, "ps_mags_n76.rds"))
aud <- read_csv(file.path(N76, "supp_fig/SuppFig4_sequencing_audit_n76.csv"), show_col_types = FALSE) %>%
  group_by(donor, Group, cohort) %>%
  summarise(depth_M = sum(depth_M), mapped = sum(mapped_reads), .groups = "drop") %>%
  mutate(cohort = gsub("\n", " ", cohort))

otu  <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
sd_df <- data.frame(sample_data(ps))
sd_df$donor <- NULL                       # the sample table already carries a donor column
meta <- sd_df %>% rownames_to_column("donor") %>%
  left_join(aud %>% select(donor, depth_M, cohort), by = "donor")
stopifnot(!any(is.na(meta$depth_M)), identical(meta$donor, colnames(otu)))
meta$ld <- log10(meta$depth_M)

# ---- 1. what the depths are -------------------------------------------------
say("## 1. sequencing depth by cohort (analysis-ready read pairs, millions)")
by_coh <- aud %>% group_by(cohort) %>%
  summarise(n = n(), median = median(depth_M), min = min(depth_M), max = max(depth_M), .groups = "drop")
for (i in seq_len(nrow(by_coh))) with(by_coh[i, ],
  say(sprintf("   %-24s n = %2d | median %6.1f M | range %5.1f to %6.1f M", cohort, n, median, min, max)))
say(sprintf("   spread across donors: %.1f-fold | Kruskal-Wallis across cohorts p = %.3g",
            max(aud$depth_M) / min(aud$depth_M), kruskal.test(depth_M ~ cohort, data = aud)$p.value))
say(sprintf("   case %.1f M against control %.1f M, Wilcoxon p = %.3g (the arms are not depth-matched)",
            median(aud$depth_M[aud$Group == "Case"]), median(aud$depth_M[aud$Group == "Control"]),
            wilcox.test(depth_M ~ Group, data = aud)$p.value))
write_csv(aud, file.path(OUTDIR, "depth_per_donor_n76.csv"))

# ---- 2. depth as a covariate ------------------------------------------------
meta$shannon  <- diversity(t(otu), "shannon")[meta$donor]
meta$richness <- colSums(otu > 0)[meta$donor]
clr <- function(x) { x <- x + 0.5; t(apply(t(x), 1, function(v) log(v) - mean(log(v)))) }
d_raw <- vegdist(clr(otu), "euclidean")
meta$dist_centroid <- as.numeric(betadisper(d_raw, meta$Group)$distances)

say("")
say("## 2. depth as a covariate on the community endpoints")
for (v in c("shannon", "richness", "dist_centroid")) {
  f0 <- lm(reformulate("Group", v), data = meta)
  f1 <- lm(reformulate(c("Group", "ld"), v), data = meta)
  s1 <- summary(f1)$coefficients
  say(sprintf("   %-13s group alone %+8.3f (p = %.2g) | with log10 depth %+8.3f (p = %.2g) | depth term p = %.2g",
              v, coef(f0)[2], summary(f0)$coefficients[2, 4], coef(f1)[2], s1[2, 4], s1[3, 4]))
}
for (g in c("Case", "Control")) {
  sub <- meta[meta$Group == g, ]
  ct <- suppressWarnings(cor.test(sub$shannon, sub$ld, method = "spearman"))
  say(sprintf("   within %-7s depth vs Shannon: rho = %+.2f (p = %.2g)", g, ct$estimate, ct$p.value))
}

# ---- 3. rarefy to the shallowest library ------------------------------------
scale_to <- min(colSums(otu))
rar <- t(rrarefy(t(otu), sample = scale_to))
sh_r <- diversity(t(rar), "shannon"); rich_r <- colSums(rar > 0)
g <- meta$Group
say("")
say(sprintf("## 3. every donor rarefied to the shallowest library (%s mapped reads, seed 42)",
            format(scale_to, big.mark = ",")))
say(sprintf("   Shannon  raw case %.2f vs control %.2f, p = %.2g | rarefied %.2f vs %.2f, p = %.2g",
            median(meta$shannon[g == "Case"]), median(meta$shannon[g == "Control"]),
            wilcox.test(meta$shannon ~ g)$p.value,
            median(sh_r[g == "Case"]), median(sh_r[g == "Control"]), wilcox.test(sh_r ~ g)$p.value))
say(sprintf("   richness raw case %.0f vs control %.0f, p = %.2g | rarefied %.0f vs %.0f, p = %.2g",
            median(meta$richness[g == "Case"]), median(meta$richness[g == "Control"]),
            wilcox.test(meta$richness ~ g)$p.value,
            median(rich_r[g == "Case"]), median(rich_r[g == "Control"]), wilcox.test(rich_r ~ g)$p.value))
say(sprintf("   Shannon per donor, raw against rarefied: Spearman rho = %.3f", cor(meta$shannon, sh_r, method = "spearman")))
d_rar <- vegdist(clr(rar), "euclidean")
a1 <- adonis2(d_raw ~ g, permutations = 999); a2 <- adonis2(d_rar ~ g, permutations = 999)
say(sprintf("   Aitchison PERMANOVA raw R2 = %.3f (p = %.3f) | rarefied R2 = %.3f (p = %.3f)",
            a1$R2[1], a1$`Pr(>F)`[1], a2$R2[1], a2$`Pr(>F)`[1]))

# ---- 4. differential abundance on the rarefied table ------------------------
# same settings as the primary run in rebuild76_ancombc2.R
ps_rar <- phyloseq(otu_table(rar, taxa_are_rows = TRUE), sample_data(ps))
say("")
say("## 4. ANCOM-BC2 re-run on the rarefied table (prv_cut 0.01, lib_cut 1000, BH)")
fit <- ancombc2(data = ps_rar, fix_formula = "Group", p_adj_method = "BH", pseudo_sens = TRUE,
                prv_cut = 0.01, lib_cut = 1000, group = "Group", struc_zero = FALSE,
                neg_lb = TRUE, alpha = 0.05, global = FALSE)
rar_res <- fit$res %>% transmute(taxon, lfc_rar = lfc_GroupCase, q_rar = q_GroupCase,
                                 sig_rar = diff_GroupCase %in% TRUE)
prim <- read_csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), show_col_types = FALSE) %>%
  transmute(taxon, lfc_prim = lfc_GroupCase, q_prim = q_GroupCase, sig_prim = diff_GroupCase %in% TRUE)
cmp <- inner_join(prim, rar_res, by = "taxon")
write_csv(cmp, file.path(OUTDIR, "ancombc2_raw_vs_rarefied_n76.csv"))

both <- cmp %>% filter(sig_prim)
say(sprintf("   significant: primary %d | rarefied %d | overlap %d of %d primary calls (%.0f%%)",
            sum(cmp$sig_prim), sum(cmp$sig_rar), sum(both$sig_rar), nrow(both),
            100 * mean(both$sig_rar)))
say(sprintf("   direction held for %d of %d primary calls (%.1f%%)",
            sum(sign(both$lfc_prim) == sign(both$lfc_rar)), nrow(both),
            100 * mean(sign(both$lfc_prim) == sign(both$lfc_rar))))
ok <- complete.cases(cmp$lfc_prim, cmp$lfc_rar)   # a MAG can fall under prv_cut in one run
say(sprintf("   log fold change agreement over the %d MAGs tested in both: Spearman rho = %.3f, Pearson r = %.3f",
            sum(ok), cor(cmp$lfc_prim[ok], cmp$lfc_rar[ok], method = "spearman"),
            cor(cmp$lfc_prim[ok], cmp$lfc_rar[ok])))
say(sprintf("   enriched/depleted split: primary %d/%d | rarefied %d/%d",
            sum(cmp$sig_prim & cmp$lfc_prim > 0), sum(cmp$sig_prim & cmp$lfc_prim < 0),
            sum(cmp$sig_rar  & cmp$lfc_rar  > 0), sum(cmp$sig_rar  & cmp$lfc_rar  < 0)))

writeLines(STATS, file.path(OUTDIR, "DEPTH_RAREFACTION_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "DEPTH_RAREFACTION_STATS_n76.txt"))
