# =============================================================================
# Audit item 3: the phylum and family BLUPs are random-effect deviations fitted
# on the 290 significant MAGs, so they describe how a clade's MAGs deviate from
# the average significant MAG, not whether the clade itself changed in the gut.
#
# This runs the matching absolute test: MAG counts summed to phylum and to
# family per donor, then ANCOM-BC2 on all 76 donors with no pre-selection.
#
# Outputs -> result2/n76/taxon_level/
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(ANCOMBC); library(dplyr); library(tidyr); library(readr)
})
set.seed(42)

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "taxon_level")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

STATS <- character()
say <- function(s) { STATS <<- c(STATS, s); message(s) }

ps  <- readRDS(file.path(N76, "ps_mags_n76.rds"))
otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)

tax <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                  quote = "", comment.char = "") %>%
  transmute(taxon = catalog_id,
            Phylum = sub("^p__", "", sapply(strsplit(gtdbtk_classification, ";"), `[`, 2)),
            Family = sub("^f__", "", sapply(strsplit(gtdbtk_classification, ";"), `[`, 5))) %>%
  filter(taxon %in% rownames(otu))
tax$Phylum[!nzchar(tax$Phylum)] <- "Unassigned"
tax$Family[!nzchar(tax$Family)] <- "Unassigned"

run_level <- function(rank) {
  key <- setNames(tax[[rank]], tax$taxon)[rownames(otu)]
  agg <- rowsum(otu, group = key)                       # counts summed within clade
  agg <- agg[rownames(agg) != "Unassigned", , drop = FALSE]
  ps_l <- phyloseq(otu_table(agg, taxa_are_rows = TRUE), sample_data(ps))
  fit <- ancombc2(data = ps_l, fix_formula = "Group", p_adj_method = "BH", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, group = "Group", struc_zero = FALSE,
                  neg_lb = TRUE, alpha = 0.05, global = FALSE)
  out <- fit$res %>%
    transmute(clade = taxon, lfc = lfc_GroupCase, se = se_GroupCase,
              q = q_GroupCase, sig = diff_GroupCase %in% TRUE,
              ci_lo = lfc_GroupCase - 1.96 * se_GroupCase,
              ci_hi = lfc_GroupCase + 1.96 * se_GroupCase) %>%
    arrange(lfc)
  n_mags <- tibble(clade = names(table(key)), n_mags = as.integer(table(key)))
  out <- left_join(out, n_mags, by = "clade")
  write_csv(out, file.path(OUTDIR, sprintf("ancombc2_%s_level_n76.csv", tolower(rank))))
  say("")
  say(sprintf("## %s level, all 76 donors, no MAG pre-selection (prevalence >= 10%%)", rank))
  say(sprintf("   %d %s tested | %d differ at BH q < 0.05 (%d higher in cases, %d lower)",
              nrow(out), tolower(rank), sum(out$sig), sum(out$sig & out$lfc > 0), sum(out$sig & out$lfc < 0)))
  for (i in which(out$sig)) with(out[i, ],
    say(sprintf("   %-24s LFC %+6.2f (95%% CI %+6.2f to %+6.2f), q = %.3g, %d MAGs",
                clade, lfc, ci_lo, ci_hi, q, n_mags)))
  out
}

phy <- run_level("Phylum")
fam <- run_level("Family")

# do the donor-level results agree in direction with the BLUPs reported in Fig. 1g,h?
blup_p <- read_tsv(file.path(N76, "phylum_glmm_blups_n76.tsv"), show_col_types = FALSE)
if ("Phylum" %in% names(blup_p)) {
  j <- inner_join(phy, blup_p, by = c("clade" = "Phylum"))
  if (nrow(j) > 2) {
    say("")
    say(sprintf("phylum BLUP against donor-level LFC: Spearman rho = %.2f over %d phyla; direction agrees for %d",
                cor(j$lfc, j$blup, method = "spearman"), nrow(j), sum(sign(j$lfc) == sign(j$blup))))
  }
}

writeLines(STATS, file.path(OUTDIR, "TAXON_LEVEL_DA_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "TAXON_LEVEL_DA_STATS_n76.txt"))
