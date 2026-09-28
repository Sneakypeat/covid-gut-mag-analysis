# =============================================================================
# ZOE gut-microbiome health rank, recomputed on the 76-donor cohort
#
# Two analyses:
#
#   A  per-SGB concordance  -> correlation figure, SUPPLEMENTARY
#       For each SGB detected here, its ZOE health-rank against its Spearman
#       correlation with case status. Tests whether the species ZOE calls
#       unfavourable are the ones that rise in COVID.
#
#   B  per-donor health index -> box plot, MAIN Figure 1 panel c
#       health_index = sum(rel_ab * (1 - 2*HEALTH_ranks)) / sum(rel_ab),
#       bounded [-1, 1], higher = healthier.
#
# Every one of the 76 donors is profiled at donor level, so there is no
# run-level averaging step.
#
# Outputs -> result2/n76/
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(readr); library(ggplot2)
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

source("R/config.R", local = TRUE)
ZOE_XLSX <- Sys.getenv("COVID_MAG_ZOE_XLSX", unset = file.path(BASE, "references", "41586_2025_9854_MOESM3_ESM.xlsx"))
dir.create(file.path(OUTDIR, "fig1"), showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

HEALTH_PAL <- c("#1B7837", "#A6DBA0", "#F7F7F7", "#C2A5CF", "#762A83")
cols <- c("Control" = "#00BFC4", "Case" = "#F8766D")

# ---- inputs -----------------------------------------------------------------

rankings <- read_excel(ZOE_XLSX, sheet = "S5", skip = 3) %>%
  select(SGB, HEALTH_ranks) %>%
  filter(!is.na(SGB), !is.na(HEALTH_ranks)) %>%
  mutate(SGB = as.character(SGB), HEALTH_ranks = as.numeric(HEALTH_ranks))
say("ZOE ranked SGBs: ", nrow(rankings))

abundance <- read_tsv(file.path(INDIR, "merged_abundance_SGB_76.tsv"),
                      comment = "#", show_col_types = FALSE) %>% rename(Taxon = 1)
sgb_abundance <- abundance %>%
  filter(str_detect(Taxon, "t__SGB")) %>%
  mutate(SGB = str_extract(Taxon, "SGB[0-9]+(?:_group)?")) %>%
  select(-Taxon)
say("SGB rows in cohort profile: ", nrow(sgb_abundance))

meta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), stringsAsFactors = FALSE)
meta$cohort <- ifelse(meta$group == "Case", "COVID", "Control")

# merge_metaphlan_tables suffixes every column with ".profile"
norm_id <- function(x) x %>% str_remove("\\.profile(\\.txt)?$") %>% trimws()
abund_cols <- setdiff(colnames(sgb_abundance), "SGB")
col_map <- tibble(orig_col = abund_cols, join_id = norm_id(abund_cols))
keep <- col_map %>% inner_join(meta %>% select(sample, cohort), by = c("join_id" = "sample"))
keep$group_bin <- ifelse(keep$cohort == "COVID", 1, 0)
say(sprintf("donors matched: %d (%d COVID, %d Control)",
            nrow(keep), sum(keep$group_bin), sum(1 - keep$group_bin)))
if (nrow(keep) != 76) stop("FATAL: expected 76 donors after matching, got ", nrow(keep))

long <- sgb_abundance %>%
  select(SGB, all_of(keep$orig_col)) %>%
  pivot_longer(-SGB, names_to = "orig_col", values_to = "rel_abundance") %>%
  left_join(keep, by = "orig_col") %>%
  inner_join(rankings, by = "SGB")
say("SGBs carrying a ZOE rank and detected here: ", length(unique(long$SGB)))

# ---- A. per-SGB concordance  (SUPPLEMENTARY correlation figure) ------------

per_sgb_cor <- long %>%
  group_by(SGB) %>%
  filter(sum(rel_abundance > 0) >= 5) %>%
  summarise(n_samples = n(), prevalence = mean(rel_abundance > 0),
            cor_with_case = suppressWarnings(cor(rel_abundance, group_bin, method = "spearman")),
            HEALTH_ranks = first(HEALTH_ranks), .groups = "drop") %>%
  filter(!is.na(cor_with_case))

conc <- suppressWarnings(cor.test(per_sgb_cor$HEALTH_ranks, per_sgb_cor$cor_with_case,
                                  method = "spearman"))
rho <- unname(conc$estimate); pval <- conc$p.value
say(sprintf("A per-SGB concordance: Spearman rho = %.3f, p = %.3g (n = %d SGBs)",
            rho, pval, nrow(per_sgb_cor)))
write_csv(per_sgb_cor, file.path(OUTDIR, "ZOE_Concordance_PerSGB_n76.csv"))

pCorr <- ggplot(per_sgb_cor, aes(x = HEALTH_ranks, y = cor_with_case, color = HEALTH_ranks)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.7, se = TRUE) +
  scale_color_gradientn(colors = HEALTH_PAL, limits = c(0, 1), name = "ZOE\nhealth-rank") +
  annotate("text", x = 0, y = max(per_sgb_cor$cor_with_case), hjust = 0, vjust = 1,
           label = sprintf("Spearman rho = %.2f, p = %.2g\nn = %d SGBs", rho, pval,
                           nrow(per_sgb_cor)),
           fontface = "bold", size = 3.8) +
  labs(title = "Disease-associated shift aligns with the ZOE healthy-gut reference",
       subtitle = "Each point is one SGB. ZOE-favourable species are depleted in COVID; ZOE-unfavourable species are enriched.",
       x = "ZOE microbiome health-rank (0 = favourable, 1 = unfavourable)",
       y = "SGB correlation with COVID case status\n(Spearman, relative abundance)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", plot.subtitle = element_text(size = 8.5), aspect.ratio = 1)

ggsave(file.path(OUTDIR, "fig1", "SuppFig_ZOE_Concordance_n76.pdf"), pCorr,
       width = 7.5, height = 6.2, device = cairo_pdf, bg = "transparent")

# ---- B. per-donor health index  (MAIN figure box plot) --------------------

per_sample <- long %>%
  group_by(orig_col, cohort) %>%
  summarise(health_index     = sum(rel_abundance * (1 - 2 * HEALTH_ranks)) / sum(rel_abundance),
            mb_health_score  = sum(rel_abundance * (1 - HEALTH_ranks)),
            ranked_abundance = sum(rel_abundance), .groups = "drop") %>%
  mutate(donor = norm_id(orig_col))

cliffs_delta <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))
ci_ctrl  <- per_sample$health_index[per_sample$cohort == "Control"]
ci_covid <- per_sample$health_index[per_sample$cohort == "COVID"]
wt <- wilcox.test(ci_covid, ci_ctrl, exact = FALSE)
cd <- cliffs_delta(ci_covid, ci_ctrl)
say(sprintf("B health index (n=%d ctrl vs %d case): Wilcoxon W = %.0f, p = %.3g, Cliff's delta = %.2f",
            length(ci_ctrl), length(ci_covid), wt$statistic, wt$p.value, cd))
say(sprintf("   medians: Control %.3f (range %.3f to %.3f) ; COVID %.3f (range %.3f to %.3f)",
            median(ci_ctrl), min(ci_ctrl), max(ci_ctrl),
            median(ci_covid), min(ci_covid), max(ci_covid)))

# Written in the shape figure1_n76.R expects for panel c
write_csv(per_sample %>% select(donor, cohort, health_index, mb_health_score, ranked_abundance),
          file.path(OUTDIR, "ZOE_HealthIndex_PerDonor_n76.csv"))

writeLines(c(
  sprintf("ZOE per-SGB concordance: Spearman rho = %.3f, p = %.3g, n = %d SGBs", rho, pval, nrow(per_sgb_cor)),
  sprintf("ZOE health index: Wilcoxon W = %.0f, p = %.3g, Cliff's delta = %.2f", wt$statistic, wt$p.value, cd),
  sprintf("  Control n = %d, median %.3f [%.3f, %.3f]", length(ci_ctrl), median(ci_ctrl), min(ci_ctrl), max(ci_ctrl)),
  sprintf("  COVID   n = %d, median %.3f [%.3f, %.3f]", length(ci_covid), median(ci_covid), min(ci_covid), max(ci_covid))
), file.path(OUTDIR, "fig1", "ZOE_STATS_n76.txt"))

say("done -> ", OUTDIR)
