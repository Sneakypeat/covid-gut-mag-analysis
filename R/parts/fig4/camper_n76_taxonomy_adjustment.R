# CAMPER follow-up requested for the 76-donor COVID analysis:
#   1) control-source robustness and taxonomic-contribution panels;
#   2) case-status association conditional on Enterobacteriaceae abundance.
#
# The adjustment is a conditional association analysis, not a causal mediation
# model. Enterobacteriaceae abundance may lie on the disease-associated path.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(sandwich)
  library(lmtest)
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

BASE <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
IN <- file.path(BASE, "result2/inputfiles/new")
CAMPER <- file.path(BASE, "result2/n76/fig4")
FOLLOW <- file.path(CAMPER, "followup")
FIG4 <- file.path(BASE, "result2/n76/fig4")
dir.create(FOLLOW, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG4, recursive = TRUE, showWarnings = FALSE)

cols <- c("Control" = "#00A6A6", "Case" = "#E85D5D")
theme_pub <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7, face = "bold"),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6, colour = "grey30"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )

primary <- read_csv(file.path(CAMPER, "camper_module_tests_n76.csv"),
                    show_col_types = FALSE)
source_tests <- read_csv(file.path(FOLLOW, "module_control_source_sensitivity_n76.csv"),
                         show_col_types = FALSE)
drivers <- read_csv(file.path(FOLLOW, "significant_module_genus_drivers_n76.csv"),
                    show_col_types = FALSE)
duplicates <- read_csv(file.path(FOLLOW, "exact_duplicate_module_profiles_n76.csv"),
                       show_col_types = FALSE)

# ---- Panel 1a: every primary-significant module against each control source --

source_n <- source_tests %>%
  distinct(control_source, n_control) %>%
  mutate(source_label = sprintf("%s (n = %d)", control_source, n_control))

source_plot <- source_tests %>%
  filter(primary_q < 0.05) %>%
  left_join(source_n, by = c("control_source", "n_control")) %>%
  mutate(
    source_label = factor(source_label, levels = source_n$source_label),
    source_FDR = if_else(q_within_source < 0.05, "q < 0.05", "q >= 0.05")
  )

p_source <- ggplot(source_plot, aes(primary_delta, delta, colour = source_FDR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_point(size = 1.35, alpha = 0.82) +
  facet_wrap(~source_label, nrow = 1) +
  scale_colour_manual(values = c("q < 0.05" = "#7B3294", "q >= 0.05" = "#BDBDBD")) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    title = "a  Control-cohort robustness",
    subtitle = "Each point is one primary-significant CAMPER module; dashed line is equal effect",
    x = "Primary Cliff's delta (all controls)",
    y = "Source-specific Cliff's delta",
    colour = NULL
  ) + theme_pub + theme(legend.position = "bottom")

# ---- Panel 1b: taxonomic drivers without duplicate-module overweighting -----

module_profile <- primary %>%
  filter(q < 0.05) %>%
  select(module) %>%
  left_join(duplicates %>% select(module, profile_group), by = "module") %>%
  mutate(profile_id = coalesce(profile_group, paste0("unique::", module)))

representatives <- module_profile %>%
  group_by(profile_id) %>%
  arrange(module, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()
n_profiles <- nrow(representatives)

taxon_summary <- drivers %>%
  semi_join(representatives, by = "module") %>%
  group_by(module) %>%
  mutate(signed_share = contribution_difference / sum(abs(contribution_difference))) %>%
  ungroup() %>%
  group_by(genus) %>%
  summarise(
    positive_mean_share = sum(pmax(signed_share, 0), na.rm = TRUE) / n_profiles,
    negative_mean_share = sum(pmin(signed_share, 0), na.rm = TRUE) / n_profiles,
    mean_absolute_share = sum(abs(signed_share), na.rm = TRUE) / n_profiles,
    profiles_present = n_distinct(module),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_absolute_share))
write_csv(taxon_summary,
          file.path(FOLLOW, "camper_taxonomic_contribution_unique_profiles_n76.csv"))

top_genera <- taxon_summary %>% slice_head(n = 12) %>% pull(genus)
taxon_plot <- taxon_summary %>%
  filter(genus %in% top_genera) %>%
  mutate(genus = fct_reorder(genus, mean_absolute_share)) %>%
  pivot_longer(c(positive_mean_share, negative_mean_share),
               names_to = "direction", values_to = "mean_share") %>%
  mutate(direction = recode(direction,
                            positive_mean_share = "Higher in Case",
                            negative_mean_share = "Lower in Case"))

p_taxon <- ggplot(taxon_plot, aes(mean_share, genus, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("Higher in Case" = cols[["Case"]],
                               "Lower in Case" = cols[["Control"]])) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "b  Taxa contributing to module shifts",
    subtitle = sprintf("Mean signed share across %d distinct module profiles; duplicate profiles counted once", n_profiles),
    x = "Mean contribution to absolute module difference",
    y = NULL,
    fill = NULL
  ) + theme_pub + theme(legend.position = "bottom")

robustness_figure <- p_source / p_taxon +
  plot_layout(heights = c(1.0, 1.2))
ggsave(file.path(FIG4, "Fig4ab_CAMPER_control_source_taxon_panel_n76.pdf"),
       robustness_figure, width = 7.2, height = 7.0, units = "in", device = cairo_pdf)
ggsave(file.path(FIG4, "Fig4ab_CAMPER_control_source_taxon_panel_n76.png"),
       robustness_figure, width = 7.2, height = 7.0, units = "in", dpi = 400)
saveRDS(
  list(panel_a = p_source, panel_b = p_taxon, combined_ab = robustness_figure),
  file.path(FIG4, "Fig4ab_CAMPER_plot_objects_n76.rds"),
  compress = "xz"
)

# ---- Test 2: conditional case effect after Enterobacteriaceae adjustment ----

meta <- read_tsv(file.path(IN, "sample_metadata_76.tsv"), show_col_types = FALSE) %>%
  transmute(donor = sample,
            Group = factor(group, levels = c("Control", "Case")), source)
capacity <- read_tsv(file.path(CAMPER, "camper_donor_module_capacity_n76.tsv"),
                     show_col_types = FALSE)

mp <- read_tsv(file.path(IN, "merged_abundance_SGB_76.tsv"),
               comment = "#", show_col_types = FALSE)
family_row <- mp %>%
  filter(grepl("\\|f__Enterobacteriaceae$", clade_name))
stopifnot(nrow(family_row) == 1L)

entero <- family_row %>%
  select(-clade_name) %>%
  pivot_longer(everything(), names_to = "profile", values_to = "Enterobacteriaceae_percent") %>%
  mutate(donor = sub("\\.profile$", "", profile)) %>%
  select(donor, Enterobacteriaceae_percent)
stopifnot(setequal(meta$donor, entero$donor), nrow(entero) == 76L)

min_positive <- min(entero$Enterobacteriaceae_percent[entero$Enterobacteriaceae_percent > 0])
pseudo <- min_positive / 2
entero <- entero %>%
  mutate(
    Enterobacteriaceae_log10 = log10(Enterobacteriaceae_percent + pseudo),
    Enterobacteriaceae_log10_z = as.numeric(scale(Enterobacteriaceae_log10))
  ) %>%
  left_join(meta, by = "donor")
write_tsv(entero, file.path(FOLLOW, "enterobacteriaceae_metaphlan_n76.tsv"))

robust_term <- function(fit, term) {
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  stopifnot(term %in% rownames(ct))
  tibble(
    beta = unname(ct[term, 1]),
    SE_HC3 = unname(ct[term, 2]),
    t_HC3 = unname(ct[term, 3]),
    p_HC3 = unname(ct[term, 4]),
    conf_low = beta - qt(0.975, df.residual(fit)) * SE_HC3,
    conf_high = beta + qt(0.975, df.residual(fit)) * SE_HC3
  )
}

model_data <- capacity %>%
  semi_join(primary %>% select(module), by = "module") %>%
  select(donor, module, oxygen, cls, n_steps, capacity) %>%
  inner_join(entero %>% select(donor, Group, source,
                               Enterobacteriaceae_percent,
                               Enterobacteriaceae_log10_z), by = "donor") %>%
  mutate(capacity_asin_sqrt = asin(sqrt(pmin(pmax(capacity, 0), 1))))

adjusted <- model_data %>%
  group_by(module, oxygen, cls, n_steps) %>%
  group_modify(~{
    unadj_fit <- lm(capacity_asin_sqrt ~ Group, data = .x)
    full_fit <- lm(capacity_asin_sqrt ~ Group + Enterobacteriaceae_log10_z, data = .x)
    tax_fit <- lm(capacity_asin_sqrt ~ Enterobacteriaceae_log10_z, data = .x)
    unadj <- robust_term(unadj_fit, "GroupCase")
    adj <- robust_term(full_fit, "GroupCase")
    tax <- robust_term(full_fit, "Enterobacteriaceae_log10_z")
    delta_r2 <- summary(full_fit)$r.squared - summary(tax_fit)$r.squared
    tibble(
      unadjusted_beta = unadj$beta,
      unadjusted_SE_HC3 = unadj$SE_HC3,
      unadjusted_p_HC3 = unadj$p_HC3,
      adjusted_beta = adj$beta,
      adjusted_SE_HC3 = adj$SE_HC3,
      adjusted_t_HC3 = adj$t_HC3,
      adjusted_p_HC3 = adj$p_HC3,
      adjusted_conf_low = adj$conf_low,
      adjusted_conf_high = adj$conf_high,
      enterobacteriaceae_beta = tax$beta,
      enterobacteriaceae_p_HC3 = tax$p_HC3,
      partial_R2_Group = delta_r2 / (1 - summary(tax_fit)$r.squared),
      full_model_R2 = summary(full_fit)$r.squared
    )
  }) %>%
  ungroup() %>%
  mutate(
    unadjusted_q_HC3 = p.adjust(unadjusted_p_HC3, "BH"),
    adjusted_q_HC3 = p.adjust(adjusted_p_HC3, "BH"),
    enterobacteriaceae_q_HC3 = p.adjust(enterobacteriaceae_p_HC3, "BH"),
    retained_after_adjustment = adjusted_q_HC3 < 0.05 & sign(adjusted_beta) == sign(unadjusted_beta)
  ) %>%
  arrange(adjusted_q_HC3)
write_csv(adjusted,
          file.path(FOLLOW, "camper_enterobacteriaceae_adjusted_models_n76.csv"))

# With two predictors, the predictor VIF is 1/(1-r^2). This quantifies whether
# the taxonomic bloom and case label leave enough overlap for conditional tests.
group_binary <- as.numeric(entero$Group == "Case")
r_group_tax <- cor(group_binary, entero$Enterobacteriaceae_log10_z)
vif_group <- 1 / (1 - r_group_tax^2)

family_w <- wilcox.test(Enterobacteriaceae_percent ~ Group, data = entero, exact = FALSE)
cliffs_delta <- with(entero,
  mean(outer(Enterobacteriaceae_percent[Group == "Case"],
             Enterobacteriaceae_percent[Group == "Control"],
             function(x, y) sign(x - y))))

plot_adjusted <- adjusted %>%
  mutate(
    status = case_when(
      adjusted_q_HC3 < 0.05 & adjusted_beta > 0 ~ "Case effect retained (+)",
      adjusted_q_HC3 < 0.05 & adjusted_beta < 0 ~ "Case effect retained (-)",
      TRUE ~ "Not retained"
    ),
    label = if_else(adjusted_q_HC3 < 0.05, module, NA_character_)
  )

p_adjust <- ggplot(plot_adjusted,
                   aes(adjusted_beta, -log10(pmax(adjusted_q_HC3, .Machine$double.xmin)),
                       colour = status)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3,
             colour = "grey45") +
  geom_point(size = 1.7, alpha = 0.85) +
  scale_colour_manual(values = c(
    "Case effect retained (+)" = cols[["Case"]],
    "Case effect retained (-)" = cols[["Control"]],
    "Not retained" = "#BDBDBD"
  )) +
  labs(
    title = "CAMPER association conditional on Enterobacteriaceae abundance",
    subtitle = "HC3-robust linear model on asin-sqrt capacity; BH across all tested modules",
    x = "Adjusted case coefficient",
    y = expression(-log[10](italic(q))),
    colour = NULL
  ) + theme_pub + theme(legend.position = "bottom")
# The conditional Enterobacteriaceae model is retained as an audit table only.
# It is not a Figure 4 panel because conditioning on a likely disease-linked
# mediator changes the estimand and does not provide a valid functional
# denominator comparison.

summary_lines <- c(
  sprintf("controls_by_source\t%s", paste(sprintf("%s:%d", source_n$control_source,
                                                   source_n$n_control), collapse = ";")),
  sprintf("primary_significant_modules\t%d", sum(primary$q < 0.05)),
  sprintf("primary_significant_same_direction_all_control_sources\t%d",
          source_tests %>% filter(primary_q < 0.05) %>% group_by(module) %>%
            summarise(ok = all(same_direction_as_primary), .groups = "drop") %>%
            summarise(n = sum(ok)) %>% pull(n)),
  sprintf("primary_significant_q_lt_0.05_all_control_sources\t%d",
          source_tests %>% filter(primary_q < 0.05) %>% group_by(module) %>%
            summarise(ok = all(q_within_source < 0.05), .groups = "drop") %>%
            summarise(n = sum(ok)) %>% pull(n)),
  sprintf("distinct_significant_module_profiles\t%d", n_profiles),
  sprintf("enterobacteriaceae_half_minimum_pseudocount_percent\t%.10g", pseudo),
  sprintf("enterobacteriaceae_control_median_percent\t%.6g",
          median(entero$Enterobacteriaceae_percent[entero$Group == "Control"])),
  sprintf("enterobacteriaceae_case_median_percent\t%.6g",
          median(entero$Enterobacteriaceae_percent[entero$Group == "Case"])),
  sprintf("enterobacteriaceae_case_vs_control_wilcoxon_p\t%.10g", family_w$p.value),
  sprintf("enterobacteriaceae_case_vs_control_cliffs_delta\t%.6f", cliffs_delta),
  sprintf("group_enterobacteriaceae_point_biserial_r\t%.6f", r_group_tax),
  sprintf("group_enterobacteriaceae_VIF\t%.6f", vif_group),
  sprintf("modules_tested_adjusted\t%d", nrow(adjusted)),
  sprintf("case_effect_retained_adjusted_q_lt_0.05\t%d",
          sum(adjusted$retained_after_adjustment)),
  sprintf("case_effect_not_retained_adjusted_q_ge_0.05\t%d",
          sum(!adjusted$retained_after_adjustment))
)
writeLines(summary_lines,
           file.path(FOLLOW, "CAMPER_TAXONOMY_ADJUSTMENT_STATS_n76.txt"))
message(paste(summary_lines, collapse = "\n"))
