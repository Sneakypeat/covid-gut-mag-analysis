# =============================================================================
# Is the patient-side KEGG signal separable from the Enterobacteriaceae bloom?
#
# The same question that was already asked of CAMPER (where all 58 module
# effects vanished once Enterobacteriaceae abundance entered the model) and of
# the mobilome (where the difference vanished once those genomes were removed).
# Here it is asked of the gene catalogue, with the published LinDA settings and
# one extra term:
#
#   published   ~ Group
#   adjusted    ~ Group + Enterobacteriaceae_log10_z
#
# then the same hypergeometric ORA is re-run on the adjusted KO sets, so the 48
# pathways of Figure 4c can be compared directly.
#
# This is a conditional association, not a mediation analysis. Enterobacteriaceae
# abundance is itself disease-associated, so conditioning on it changes the
# estimand: it asks which KO shifts are visible between patients and controls of
# comparable Enterobacteriaceae load, not what the bloom "causes". The bloom and
# the case label are strongly collinear, and the variance inflation is reported
# for that reason.
#
# Outputs -> result2/n76/fig4/kegg_fast/entero_adjust/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(LinDA)
})

BASE <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
N76  <- file.path(BASE, "result2/n76")
KF   <- file.path(N76, "fig4/kegg_fast")
OUT  <- file.path(KF, "entero_adjust")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

mat <- readRDS(file.path(KF, "feature_matrices/KO_abundance_fractional_RPKM_n76.rds"))
ent <- read_tsv(file.path(N76, "fig4/followup/enterobacteriaceae_metaphlan_n76.tsv"),
                show_col_types = FALSE)
meta <- ent %>%
  transmute(donor, Group = factor(Group, c("Control", "Case")), source,
            entero_z = Enterobacteriaceae_log10_z,
            entero_pct = Enterobacteriaceae_percent) %>%
  as.data.frame()
rownames(meta) <- meta$donor
stopifnot(setequal(colnames(mat), meta$donor))
meta <- meta[colnames(mat), ]

# how far apart can the two terms be pulled at all
r <- cor(as.numeric(meta$Group == "Case"), meta$entero_z)
say(sprintf("collinearity: r(Group, Enterobacteriaceae) = %.3f, VIF = %.2f", r, 1 / (1 - r^2)))
say(sprintf("Enterobacteriaceae median %.3f%% in cases against %.4f%% in controls",
            median(meta$entero_pct[meta$Group == "Case"]),
            median(meta$entero_pct[meta$Group == "Control"])))

run_linda <- function(formula) {
  fit <- linda(otu.tab = mat, meta = meta, formula = formula, type = "count",
               adaptive = TRUE, imputation = FALSE, pseudo.cnt = 0.5, corr.cut = 0.1,
               p.adj.method = "BH", alpha = 0.05, prev.cut = 0.10, lib.cut = 1,
               winsor.quan = 0.97, n.cores = 4)
  nm <- grep("Group", names(fit$output), value = TRUE)[1]
  fit$output[[nm]] %>% tibble::rownames_to_column("KO") %>%
    transmute(KO, lfc = log2FoldChange, p = pvalue, q = padj,
              Status = case_when(padj < 0.05 & log2FoldChange > 1 ~ "Higher in Case",
                                 padj < 0.05 & log2FoldChange < -1 ~ "Higher in Control",
                                 TRUE ~ "Not significant"))
}

say(""); say("## LinDA")
pub <- run_linda("~ Group")
adj <- run_linda("~ Group + entero_z")
say(sprintf("   published ~Group          : %d tested | %d higher in Case | %d higher in Control",
            nrow(pub), sum(pub$Status == "Higher in Case"), sum(pub$Status == "Higher in Control")))
say(sprintf("   adjusted  ~Group+entero   : %d tested | %d higher in Case | %d higher in Control",
            nrow(adj), sum(adj$Status == "Higher in Case"), sum(adj$Status == "Higher in Control")))

both <- pub %>% select(KO, lfc_pub = lfc, q_pub = q, st_pub = Status) %>%
  inner_join(adj %>% select(KO, lfc_adj = lfc, q_adj = q, st_adj = Status), by = "KO")
write_csv(both, file.path(OUT, "KEGG_KO_entero_adjusted_n76.csv"))
for (d in c("Higher in Case", "Higher in Control")) {
  sub <- both %>% filter(st_pub == d)
  say(sprintf("   %-18s kept %d of %d (%.0f%%) | same sign %d | median |log2FC| %.2f -> %.2f",
              d, sum(sub$st_adj == d), nrow(sub), 100 * mean(sub$st_adj == d),
              sum(sign(sub$lfc_adj) == sign(sub$lfc_pub)),
              median(abs(sub$lfc_pub)), median(abs(sub$lfc_adj))))
}
say(sprintf("   log2FC correlation between models: Spearman rho = %.3f",
            cor(both$lfc_pub, both$lfc_adj, method = "spearman")))

# ---- the same ORA, on each KO set --------------------------------------------
ko_pw <- read_tsv(file.path(N76, "gene_catalogue/KEGG_KO_pathway_map_20260917.tsv"),
                  show_col_types = FALSE) %>%
  transmute(KO = KO_ID, pid = sub("^ko", "", pathway_id), pathway_name)

ora <- function(res) {
  tested <- res$KO
  up <- res$KO[res$Status == "Higher in Case"]; dn <- res$KO[res$Status == "Higher in Control"]
  univ <- intersect(tested, ko_pw$KO)
  up <- intersect(up, univ); dn <- intersect(dn, univ); N <- length(univ)
  ko_pw %>% filter(KO %in% univ) %>% group_by(pid, pathway_name) %>%
    summarise(M = n_distinct(KO), k_case = sum(unique(KO) %in% up),
              k_control = sum(unique(KO) %in% dn), .groups = "drop") %>%
    mutate(p_case    = phyper(k_case - 1,    M, N - M, length(up), lower.tail = FALSE),
           p_control = phyper(k_control - 1, M, N - M, length(dn), lower.tail = FALSE),
           q_case = p.adjust(p_case, "BH"), q_control = p.adjust(p_control, "BH"),
           dir = ifelse(q_case <= q_control, "Case", "Control"),
           q = pmin(q_case, q_control))
}
o_pub <- ora(pub); o_adj <- ora(adj)
cmp <- o_pub %>% filter(q < 0.05) %>%
  select(pid, pathway_name, dir_pub = dir, q_pub = q) %>%
  left_join(o_adj %>% select(pid, dir_adj = dir, q_adj = q), by = "pid") %>%
  mutate(kept = q_adj < 0.05 & dir_adj == dir_pub) %>%
  arrange(desc(kept), q_pub)
write_csv(cmp, file.path(OUT, "KEGG_pathway_ORA_entero_adjusted_n76.csv"))
write_csv(o_adj, file.path(OUT, "KEGG_pathway_ORA_adjusted_full_n76.csv"))

say(""); say("## pathway ORA (BH q < 0.05)")
say(sprintf("   published model: %d pathways (%d Case, %d Control)", nrow(cmp),
            sum(cmp$dir_pub == "Case"), sum(cmp$dir_pub == "Control")))
say(sprintf("   adjusted  model: %d pathways total; of the published %d, %d survive with the same direction",
            sum(o_adj$q < 0.05), nrow(cmp), sum(cmp$kept)))
for (d in c("Case", "Control")) {
  sub <- cmp %>% filter(dir_pub == d)
  say(sprintf("   %-8s side: %d of %d survive", d, sum(sub$kept), nrow(sub)))
}
say(""); say("   published pathways that do NOT survive adjustment:")
for (i in which(!cmp$kept)) say(sprintf("      %-52s q %.2g -> %.2g (%s)", cmp$pathway_name[i],
                                        cmp$q_pub[i], cmp$q_adj[i], cmp$dir_adj[i]))
say(""); say("   top surviving pathways by adjusted q:")
top <- cmp %>% filter(kept) %>% arrange(q_adj) %>% head(12)
for (i in seq_len(nrow(top))) say(sprintf("      %-52s %s  q %.2g -> %.2g", top$pathway_name[i],
                                          top$dir_pub[i], top$q_pub[i], top$q_adj[i]))

writeLines(STATS, file.path(OUT, "KEGG_ENTERO_ADJUST_STATS_n76.txt"))
message("DONE -> ", OUT)
