# =============================================================================
# Strain-diversity model. Each donor contributes up to 20 genomes, so the
# residuals are nested in donor as well as in genome.
#
# Fit with crossed random intercepts for genome and donor, plus the
# coverage-adjusted and depth-rarefied sensitivities. Effect, 95% CI and p are
# reported for each.
#
# Outputs -> result2/n76/fig2/strain/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(lme4); library(lmerTest)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "fig2/strain")

MIN_COV <- 5; MIN_BREADTH <- 0.5; MIN_DONORS <- 5     # as in figure2_n76.R

STATS <- character()
say <- function(s) { STATS <<- c(STATS, s); message(s) }

gi <- read_tsv(file.path(INDIR, "instrain_genome_info_76.tsv"), show_col_types = FALSE) %>%
  mutate(genome = sub("[.]fa$", "", genome))
smeta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE) %>%
  transmute(donor = sample, Group = factor(group, levels = c("Control", "Case")))

ok <- gi %>% left_join(smeta, by = "donor") %>%
  filter(coverage >= MIN_COV, breadth_minCov >= MIN_BREADTH, !is.na(nucl_diversity))
keep <- ok %>% count(genome, Group) %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  filter(Case >= MIN_DONORS, Control >= MIN_DONORS)
cmp <- ok %>% filter(genome %in% keep$genome)

say(sprintf("observations %d | genomes %d | donors %d (case %d, control %d)",
            nrow(cmp), n_distinct(cmp$genome), n_distinct(cmp$donor),
            n_distinct(cmp$donor[cmp$Group == "Case"]), n_distinct(cmp$donor[cmp$Group == "Control"])))
per_donor <- cmp %>% count(donor)
say(sprintf("genomes measured per donor: median %d, range %d to %d",
            median(per_donor$n), min(per_donor$n), max(per_donor$n)))

fit_one <- function(form, data, label) {
  m  <- lmer(form, data = data)
  cf <- summary(m)$coefficients["GroupCase", ]
  ci <- confint(m, parm = "GroupCase", method = "Wald")
  say(sprintf("%-46s beta %+.3f (95%% CI %+.3f to %+.3f), p = %.3g, %.1f-fold lower in cases",
              label, cf[["Estimate"]], ci[1], ci[2], cf[["Pr(>|t|)"]], 10^(-cf[["Estimate"]])))
  tibble(model = label, beta = cf[["Estimate"]], ci_lo = ci[1], ci_hi = ci[2],
         se = cf[["Std. Error"]], df = cf[["df"]], p = cf[["Pr(>|t|)"]],
         fold = 10^(-cf[["Estimate"]]), n_obs = nobs(m))
}

say("")
say("## nucleotide diversity, log10 scale")
res <- bind_rows(
  fit_one(log10(nucl_diversity + 1e-6) ~ Group + (1 | genome), cmp,
          "genome only"),
  fit_one(log10(nucl_diversity + 1e-6) ~ Group + (1 | genome) + (1 | donor), cmp,
          "genome + donor"),
  fit_one(log10(nucl_diversity + 1e-6) ~ Group + log10(coverage) + (1 | genome) + (1 | donor), cmp,
          "genome + donor, coverage-adjusted"),
  fit_one(log10(nucl_diversity_rarefied + 1e-6) ~ Group + (1 | genome) + (1 | donor),
          cmp %>% filter(!is.na(nucl_diversity_rarefied)), "genome + donor, depth-rarefied")
)

# how much of the variance sits at the donor level
m_full <- lmer(log10(nucl_diversity + 1e-6) ~ Group + (1 | genome) + (1 | donor), data = cmp)
vc <- as.data.frame(VarCorr(m_full))
say("")
say(sprintf("variance components: genome %.4f | donor %.4f | residual %.4f (donor holds %.0f%% of the random variance)",
            vc$vcov[vc$grp == "genome"], vc$vcov[vc$grp == "donor"], vc$vcov[vc$grp == "Residual"],
            100 * vc$vcov[vc$grp == "donor"] / sum(vc$vcov)))

# a donor-level summary carries one value per donor, so no pseudo-replication at all
pd <- cmp %>% group_by(donor, Group) %>% summarise(pi = median(nucl_diversity), .groups = "drop")
w  <- wilcox.test(pi ~ Group, data = pd, exact = FALSE)
say(sprintf("donor-level median diversity: control %.2e vs case %.2e, Wilcoxon p = %.3g (n = %d donors)",
            median(pd$pi[pd$Group == "Control"]), median(pd$pi[pd$Group == "Case"]), w$p.value, nrow(pd)))

pg <- cmp %>% group_by(genome) %>%
  summarise(ctrl = median(nucl_diversity[Group == "Control"]),
            case = median(nucl_diversity[Group == "Case"]), .groups = "drop")
say(sprintf("per-genome medians: %d of %d genomes lower in cases (sign test p = %.3g)",
            sum(pg$case < pg$ctrl), nrow(pg),
            binom.test(sum(pg$case < pg$ctrl), nrow(pg))$p.value))

write_csv(res, file.path(OUTDIR, "strain_model_refit_n76.csv"))
writeLines(STATS, file.path(OUTDIR, "STRAIN_MODEL_REFIT_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "STRAIN_MODEL_REFIT_n76.txt"))
