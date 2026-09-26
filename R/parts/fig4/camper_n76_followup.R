# =============================================================================
# CAMPER follow-up audit: cohort robustness, quality sensitivity, duplicated
# module profiles, and MAG/taxon drivers for the 76-donor COVID comparison.
#
# This script does not replace camper_n76.R. It interrogates the claims made by
# that analysis and writes audit tables to result2/n76/fig4/followup/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
CAMPER <- file.path(N76, "fig4")
OUT    <- file.path(CAMPER, "followup")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rd <- function(path) read_tsv(path, show_col_types = FALSE, progress = FALSE)
cliffs <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))
safe_wilcox <- function(a, b) {
  if (!length(a) || !length(b) || length(unique(c(a, b))) < 2L) return(1)
  wilcox.test(a, b, exact = FALSE)$p.value
}

cap_primary <- rd(file.path(CAMPER, "camper_donor_module_capacity_n76.tsv"))
comp <- rd(file.path(CAMPER, "camper_module_completeness_n76.tsv"))
primary <- read_csv(file.path(CAMPER, "camper_module_tests_n76.csv"),
                    show_col_types = FALSE, progress = FALSE) %>%
  rename(primary_delta = delta, primary_p = p, primary_q = q)
tested_modules <- primary$module

meta <- rd(file.path(INDIR, "sample_metadata_76.tsv")) %>%
  transmute(donor = sample, Group = group, source)
quality <- rd(file.path(INDIR, "MAG_quality_taxonomy.tsv"))

ab <- rd(file.path(INDIR, "MAG_relative_abundance_percent.tsv"))
names(ab)[1] <- "MAG"
ab_long <- ab %>%
  pivot_longer(-MAG, names_to = "donor", values_to = "raw_percent")

da <- read_csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"),
               show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    MAG = taxon,
    da_direction = case_when(
      diff_GroupCase %in% TRUE & lfc_GroupCase > 0 ~ "Enriched",
      diff_GroupCase %in% TRUE & lfc_GroupCase < 0 ~ "Depleted",
      TRUE ~ "NS"
    ),
    da_lfc = lfc_GroupCase,
    da_q = q_GroupCase
  )

# ---- 1. Repeat each module comparison against each control source -----------

control_sources <- sort(unique(meta$source[meta$Group == "Control"]))
source_tests <- bind_rows(lapply(control_sources, function(control_source) {
  cap_primary %>%
    filter(module %in% tested_modules,
           Group == "Case" | source == control_source) %>%
    group_by(module, oxygen, cls, n_steps) %>%
    summarise(
      control_source = control_source,
      n_case = sum(Group == "Case"),
      n_control = sum(Group == "Control"),
      med_case = median(capacity[Group == "Case"]),
      med_control = median(capacity[Group == "Control"]),
      delta = cliffs(capacity[Group == "Case"], capacity[Group == "Control"]),
      p = safe_wilcox(capacity[Group == "Case"], capacity[Group == "Control"]),
      .groups = "drop"
    ) %>%
    mutate(q_within_source = p.adjust(p, "BH"))
})) %>%
  left_join(primary %>% select(module, primary_delta, primary_q), by = "module") %>%
  mutate(same_direction_as_primary =
           sign(delta) == sign(primary_delta) | (delta == 0 & primary_delta == 0))
write_csv(source_tests, file.path(OUT, "module_control_source_sensitivity_n76.csv"))

source_summary <- source_tests %>%
  group_by(module, primary_delta, primary_q) %>%
  summarise(
    direction_concordant_all_sources = all(same_direction_as_primary),
    significant_all_sources = all(q_within_source < 0.05),
    min_abs_source_delta = min(abs(delta)),
    max_abs_source_delta = max(abs(delta)),
    .groups = "drop"
  )

# ---- 2. High-quality-MAG and total-community-abundance sensitivities --------

module_meta <- primary %>% select(module, oxygen, cls, n_steps)

capacity_tests <- function(weights, analysis) {
  caps <- weights %>%
    inner_join(comp %>% select(MAG, module, completeness),
               by = "MAG", relationship = "many-to-many") %>%
    filter(module %in% tested_modules) %>%
    group_by(donor, module) %>%
    summarise(capacity = sum(weight * completeness), .groups = "drop") %>%
    complete(donor = meta$donor, module = tested_modules,
             fill = list(capacity = 0)) %>%
    left_join(meta, by = "donor")

  caps %>%
    group_by(module) %>%
    summarise(
      med_control = median(capacity[Group == "Control"]),
      med_case = median(capacity[Group == "Case"]),
      delta = cliffs(capacity[Group == "Case"], capacity[Group == "Control"]),
      p = safe_wilcox(capacity[Group == "Case"], capacity[Group == "Control"]),
      .groups = "drop"
    ) %>%
    mutate(q = p.adjust(p, "BH"), analysis = analysis)
}

weights_primary <- ab_long %>%
  group_by(donor) %>% mutate(weight = raw_percent / sum(raw_percent)) %>% ungroup()

hq_ids <- quality %>%
  filter(completeness >= 90, contamination <= 5) %>%
  pull(catalog_id)
weights_hq <- ab_long %>%
  filter(MAG %in% hq_ids) %>%
  group_by(donor) %>% mutate(weight = raw_percent / sum(raw_percent)) %>% ungroup()

# MAG_relative_abundance_percent.tsv is relative to all reads. Dividing by 100
# retains the unmapped fraction instead of renormalising within the MAG set.
weights_raw <- ab_long %>% mutate(weight = raw_percent / 100)

sens_long <- bind_rows(
  capacity_tests(weights_primary, "primary_MAG_resolved_renormalized"),
  capacity_tests(weights_hq, "HQ_MAGs_renormalized"),
  capacity_tests(weights_raw, "raw_total_community_fraction")
) %>% left_join(module_meta, by = "module")
write_csv(sens_long, file.path(OUT, "module_abundance_quality_sensitivity_long_n76.csv"))

sens_wide <- sens_long %>%
  select(module, analysis, delta, q) %>%
  pivot_wider(names_from = analysis, values_from = c(delta, q)) %>%
  left_join(primary %>% select(module, oxygen, cls, n_steps, primary_delta, primary_q),
            by = "module") %>%
  left_join(source_summary, by = c("module", "primary_delta", "primary_q"))
write_csv(sens_wide, file.path(OUT, "module_robustness_summary_n76.csv"))

# ---- 3. Exact duplicate module profiles -------------------------------------

# Exact equality is assessed at the MAG-completeness level. Such modules cannot
# be treated as independent observations in this catalogue even if CAMPER gives
# them different substrate labels.
profile_matrix <- comp %>%
  filter(module %in% tested_modules) %>%
  select(MAG, module, completeness) %>%
  pivot_wider(names_from = MAG, values_from = completeness, values_fill = 0)

profile_key <- apply(as.matrix(profile_matrix[-1]), 1, function(x)
  paste(format(x, digits = 17, scientific = TRUE, trim = TRUE), collapse = "|"))
profile_groups <- split(profile_matrix$module, profile_key)
profile_groups <- profile_groups[lengths(profile_groups) > 1]
duplicate_profiles <- imap_dfr(profile_groups, function(modules, key) {
  tibble(profile_group = sprintf("duplicate_%02d", match(key, names(profile_groups))),
         n_modules = length(modules), module = sort(modules))
}) %>%
  left_join(primary %>% select(module, oxygen, cls, primary_delta, primary_q), by = "module")
write_csv(duplicate_profiles, file.path(OUT, "exact_duplicate_module_profiles_n76.csv"))

# ---- 4. MAG and taxon drivers of significant module differences -------------

sig_modules <- primary %>% filter(primary_q < 0.05) %>% pull(module)
drivers <- weights_primary %>%
  inner_join(comp %>% filter(module %in% sig_modules, completeness > 0) %>%
               select(MAG, module, module_completeness = completeness),
             by = "MAG", relationship = "many-to-many") %>%
  left_join(meta, by = "donor") %>%
  mutate(contribution = weight * module_completeness) %>%
  group_by(module, MAG, module_completeness) %>%
  summarise(
    mean_case_contribution = mean(contribution[Group == "Case"]),
    mean_control_contribution = mean(contribution[Group == "Control"]),
    contribution_difference = mean_case_contribution - mean_control_contribution,
    .groups = "drop"
  ) %>%
  group_by(module) %>%
  mutate(
    module_net_difference = sum(contribution_difference),
    percent_of_net_difference = if_else(
      module_net_difference == 0, NA_real_,
      100 * contribution_difference / module_net_difference
    ),
    absolute_difference_rank = min_rank(desc(abs(contribution_difference)))
  ) %>%
  ungroup() %>%
  left_join(quality %>%
              select(catalog_id, catalogue_source_group = group, source_binner,
                     genome_completeness = completeness,
                     genome_contamination = contamination,
                     quality_tier, phylum, class, order, family, genus, species),
            by = c("MAG" = "catalog_id")) %>%
  left_join(da, by = "MAG") %>%
  left_join(primary %>% select(module, oxygen, cls, primary_delta, primary_q), by = "module") %>%
  arrange(module, absolute_difference_rank)
write_csv(drivers, file.path(OUT, "significant_module_MAG_drivers_n76.csv"))

genus_drivers <- drivers %>%
  mutate(genus = coalesce(na_if(genus, ""), "Unclassified")) %>%
  group_by(module, genus) %>%
  summarise(contribution_difference = sum(contribution_difference),
            .groups = "drop") %>%
  group_by(module) %>%
  mutate(module_net_difference = sum(contribution_difference),
         percent_of_net_difference = if_else(
           module_net_difference == 0, NA_real_,
           100 * contribution_difference / module_net_difference),
         absolute_difference_rank = min_rank(desc(abs(contribution_difference)))) %>%
  ungroup() %>%
  left_join(primary %>% select(module, oxygen, cls, primary_delta, primary_q), by = "module") %>%
  arrange(module, absolute_difference_rank)
write_csv(genus_drivers, file.path(OUT, "significant_module_genus_drivers_n76.csv"))

# ---- 5. Compact, machine-verifiable audit summary ----------------------------

n_sig <- sum(primary$primary_q < 0.05)
n_source_concordant <- source_summary %>%
  filter(primary_q < 0.05, direction_concordant_all_sources) %>% nrow()
n_source_sig <- source_summary %>%
  filter(primary_q < 0.05, significant_all_sources) %>% nrow()
n_hq_sig_same <- sens_wide %>%
  filter(primary_q < 0.05,
         q_HQ_MAGs_renormalized < 0.05,
         sign(delta_HQ_MAGs_renormalized) == sign(primary_delta)) %>% nrow()
n_raw_sig_same <- sens_wide %>%
  filter(primary_q < 0.05,
         q_raw_total_community_fraction < 0.05,
         sign(delta_raw_total_community_fraction) == sign(primary_delta)) %>% nrow()

summary_lines <- c(
  sprintf("tested_modules\t%d", nrow(primary)),
  sprintf("primary_q_lt_0.05\t%d", n_sig),
  sprintf("primary_significant_same_direction_all_3_control_sources\t%d", n_source_concordant),
  sprintf("primary_significant_q_lt_0.05_in_all_3_source_specific_tests\t%d", n_source_sig),
  sprintf("primary_significant_retained_HQ_MAG_sensitivity\t%d", n_hq_sig_same),
  sprintf("primary_significant_retained_raw_abundance_sensitivity\t%d", n_raw_sig_same),
  sprintf("HQ_MAGs\t%d", length(hq_ids)),
  sprintf("exact_duplicate_profile_groups\t%d", length(profile_groups)),
  sprintf("modules_in_exact_duplicate_profile_groups\t%d", nrow(duplicate_profiles))
)
writeLines(summary_lines, file.path(OUT, "CAMPER_FOLLOWUP_STATS_n76.txt"))
message(paste(summary_lines, collapse = "\n"))
message("DONE -> ", OUT)
