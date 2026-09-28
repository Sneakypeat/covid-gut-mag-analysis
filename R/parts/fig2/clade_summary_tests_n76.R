# =============================================================================
# Which clade-level summary of a MAG-level result is
# defensible? Three summaries of the same 584 genomes, reported side by side.
#
#   B. over-representation: does the clade hold more depleted (or enriched)
#      MAGs than the catalogue rate, Fisher exact, BH-corrected
#   C. clade as one unit: ANCOM-BC2 on summed counts (taxon_level_da_n76.R)
#
# The directional two-sided Fisher panels in clade_fisher_figure_n76.R are the
# version used in Figure 1g,h. B summarises the species-level result; C asks a
# different question and is bounded by the compositional constraint.
#
# Outputs -> result2/n76/taxon_level/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "taxon_level")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

STATS <- character()
say <- function(s) { STATS <<- c(STATS, s); message(s) }

anc <- read_csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), show_col_types = FALSE) %>%
  transmute(taxon, lfc = lfc_GroupCase, se = se_GroupCase, sig = diff_GroupCase %in% TRUE)
tax <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                  quote = "", comment.char = "") %>%
  transmute(taxon = catalog_id,
            Phylum = sub("^p__", "", sapply(strsplit(gtdbtk_classification, ";"), `[`, 2)),
            Family = sub("^f__", "", sapply(strsplit(gtdbtk_classification, ";"), `[`, 5)))
d <- inner_join(anc, tax, by = "taxon") %>% filter(!is.na(lfc), !is.na(se), se > 0)
say(sprintf("MAGs with an effect size and a standard error: %d", nrow(d)))

# ---- B. over-representation among the differentially abundant MAGs ----------
ora <- function(rank, direction) {
  hits <- d %>% mutate(hit = sig & if (direction == "depleted") lfc < 0 else lfc > 0)
  tot_hit <- sum(hits$hit); tot <- nrow(hits)
  hits %>% group_by(clade = .data[[rank]]) %>%
    summarise(n_mags = n(), k = sum(hit), .groups = "drop") %>%
    filter(n_mags >= 3) %>%
    rowwise() %>%
    mutate(p = fisher.test(matrix(c(k, n_mags - k, tot_hit - k, tot - n_mags - tot_hit + k), 2),
                           alternative = "greater")$p.value) %>%
    ungroup() %>% mutate(direction = direction, rate = k / n_mags, q = p.adjust(p, "BH")) %>%
    arrange(q)
}
for (rk in c("Phylum", "Family")) {
  both <- bind_rows(ora(rk, "depleted"), ora(rk, "enriched"))
  write_csv(both, file.path(OUTDIR, sprintf("clade_overrepresentation_%s_n76.csv", tolower(rk))))
  say("")
  say(sprintf("## B. %s: clades over-represented among the differentially abundant MAGs (Fisher, BH)", rk))
  sigs <- both %>% filter(q < 0.05)
  say(sprintf("   %d clade-direction combinations at q < 0.05", nrow(sigs)))
  for (i in seq_len(min(nrow(sigs), 14))) with(sigs[i, ],
    say(sprintf("   %-24s %-9s %2d of %2d MAGs (%.0f%%), q = %.3g", clade, direction, k, n_mags, 100 * rate, q)))
}

writeLines(STATS, file.path(OUTDIR, "CLADE_SUMMARY_TESTS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "CLADE_SUMMARY_TESTS_n76.txt"))
