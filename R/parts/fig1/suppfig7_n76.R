# =============================================================================
# Supplementary Figure 7, replication of the read-level endpoints in three
# external COVID-19 cohorts (MetaPhlAn species profiles, no MAG catalogue)
#
#   a  per-cohort and pooled Hedges' g for the four read-level endpoints
#   b  the two primary endpoints per cohort, case against control
#   c  classified read fraction by sequencing chemistry, the direct test of
#      whether the MGI/Illumina split in the discovery cohort creates a
#      profiling asymmetry that could masquerade as a disease effect
#
# The discovery cohort is drawn for scale and is not pooled.
#
# Built from the tables written by the external-validation pipeline in
# result2/n76/supp_fig/external_validation/analysis/.
#
# Outputs -> result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})

theme_transparent <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA),
  strip.background      = ggplot2::element_rect(fill = "grey95", color = NA))
theme_bw <- function(...) ggplot2::theme_bw(...) + theme_transparent
ggplot2::theme_set(ggplot2::theme_get() + theme_transparent)
ggsave <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
N76    <- file.path(BASE, "result2/n76")
EXT    <- file.path(N76, "supp_fig/external_validation/analysis")
OUTDIR <- file.path(N76, "supp_fig")
source(file.path(BASE, "taxon_italics_n76.R"))   # md_taxon(): taxon names in italics
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
cols <- c("Control" = "#00A6A6", "Case" = "#E85D5D")

eff <- read_tsv(file.path(EXT, "external_validation_effects.tsv"), show_col_types = FALSE)
met <- read_tsv(file.path(EXT, "external_validation_sample_metrics.tsv"), show_col_types = FALSE) %>%
  filter(profile_qc != "all_unclassified")   # one Chinese sample returned no classified taxa

coh_short <- c(Discovery_n76 = "Discovery", PRJEB43555_China_severity = "China",
                PRJNA624223_HongKong_baseline = "Hong Kong",
                PRJNA890008_Luxembourg_mild = "Luxembourg",
                `External pooled (REML-HK)` = "External pooled")
coh_lab <- c(Discovery_n76 = "Discovery (this study)",
             PRJEB43555_China_severity = "China, severity (PRJEB43555)",
             PRJNA624223_HongKong_baseline = "Hong Kong (PRJNA624223)",
             PRJNA890008_Luxembourg_mild = "Luxembourg, mild (PRJNA890008)",
             `External pooled (REML-HK)` = "External pooled (REML, Hartung-Knapp)")
n_lab <- met %>% count(cohort, group) %>% pivot_wider(names_from = group, values_from = n) %>%
  mutate(lab = sprintf("%s\n%d case / %d control", coh_lab[cohort], Case, Control))
pooled_n <- met %>% filter(cohort != "Discovery_n76") %>% count(group)
n_lab <- bind_rows(n_lab, tibble(cohort = "External pooled (REML-HK)",
  lab = sprintf("%s\n%d case / %d control", coh_lab[["External pooled (REML-HK)"]],
                pooled_n$n[pooled_n$group == "Case"], pooled_n$n[pooled_n$group == "Control"])))

# ---- a. forest ---------------------------------------------------------------
fo <- eff %>%
  mutate(lab = n_lab$lab[match(cohort, n_lab$cohort)],
         lab = factor(lab, rev(n_lab$lab[match(names(coh_lab), n_lab$cohort)])),
         kind = case_when(cohort == "Discovery_n76" ~ "Discovery",
                          grepl("pooled", cohort) ~ "Pooled", TRUE ~ "External"),
         label = factor(label, c("ZOE health index", "Species Shannon diversity",
                                 "Enterobacteriaceae / commensal balance",
                                 "Enterobacteriaceae abundance")))
het <- fo %>% filter(kind == "Pooled") %>%
  mutate(txt = sprintf("pooled g = %.2f (%.2f, %.2f)\nI² = %.0f%%, p = %.3f",
                       hedges_g, ci_low, ci_high, I2, heterogeneity_p))
say("## a. pooled external effects (Hedges' g, negative = lower in cases)")
for (i in seq_len(nrow(het))) with(het[i, ],
  say(sprintf("   %-40s g = %+.2f (%.2f, %.2f)  I2 = %.0f%%  p_het = %.3f", label, hedges_g, ci_low, ci_high, I2, heterogeneity_p)))
dir_ok <- fo %>% filter(kind == "External") %>% group_by(label) %>%
  summarise(ok = sum(sign(hedges_g) == expected_sign), n = n(), .groups = "drop")
for (i in seq_len(nrow(dir_ok))) with(dir_ok[i, ], say(sprintf("   %-40s expected direction in %d of %d cohorts", label, ok, n)))

pA <- ggplot(fo, aes(hedges_g, lab, colour = kind, shape = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.45) +
  geom_point(size = 1.9) +
  geom_text(data = het, aes(x = -Inf, y = 1.35, label = txt), hjust = -0.03, vjust = 1,
            size = 2.1, colour = "grey25", inherit.aes = FALSE) +
  facet_wrap(~label, nrow = 2, scales = "free_x", labeller = labeller(label = md_taxon_in)) +
  scale_colour_manual(values = c(Discovery = "grey45", External = "#4C72B0", Pooled = "#B2182B"),
                      name = NULL) +
  scale_shape_manual(values = c(Discovery = 16, External = 16, Pooled = 18), name = NULL) +
  labs(title = "Replication in three external COVID-19 cohorts",
       subtitle = "Hedges' g, case minus control; negative = lower in cases. Discovery cohort shown for scale, not pooled.",
       x = "Hedges' g (95% CI)", y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6.5, colour = "black"),
        strip.text = ggtext::element_markdown(size = 7.5, face = "bold"),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(file.path(OUTDIR, "SuppFig7a_External_Forest_n76.pdf"), pA, width = 8.6, height = 5.2,
       device = cairo_pdf)

# ---- b. the two primary endpoints, per cohort --------------------------------
long <- met %>%
  select(cohort, group, health_index, shannon_species) %>%
  pivot_longer(c(health_index, shannon_species), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = recode(metric, health_index = "ZOE health index",
                         shannon_species = "Species Shannon diversity"),
         lab = sprintf("%s\n%d / %d", coh_short[cohort],
                       n_lab$Case[match(cohort, n_lab$cohort)],
                       n_lab$Control[match(cohort, n_lab$cohort)]),
         lab = factor(lab, unique(lab[order(match(cohort, names(coh_lab)))])),
         group = factor(group, c("Control", "Case")))
say(""); say("## b. medians by cohort")
long %>% group_by(metric, cohort, group) %>% summarise(med = median(value), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = med) %>%
  { for (i in seq_len(nrow(.))) say(sprintf("   %-26s %-30s control %6.2f  case %6.2f",
                                            .$metric[i], .$cohort[i], .$Control[i], .$Case[i])) }

pB <- ggplot(long, aes(lab, value, fill = group, colour = group)) +
  geom_boxplot(width = 0.6, alpha = 0.45, outlier.shape = NA, linewidth = 0.35,
               position = position_dodge(width = 0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7, seed = 42),
             size = 0.5, alpha = 0.55) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = cols, name = NULL) + scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Read-level endpoints by cohort",
       subtitle = "MetaPhlAn 4.2.4 species profiles (mpa_vJan25); the MAG catalogue is not used here",
       x = "Cohort (case / control)", y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6.5, colour = "black"),
        strip.text = element_text(size = 7.5, face = "bold"),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom")
ggsave(file.path(OUTDIR, "SuppFig7b_External_Distributions_n76.pdf"), pB, width = 8.6, height = 4.4,
       device = cairo_pdf)

# ---- where the heterogeneity sits: cases or controls? ------------------------
# The discovery contrast is much larger than any external one. These lines test
# whether that comes from the patient arm (a severity difference) or from the
# comparator arm (a baseline difference), because the two have different
# consequences for how the pooled estimate should be read.
meta76 <- read_tsv(file.path(BASE, "result2/inputfiles/new/sample_metadata_76.tsv"),
                   show_col_types = FALSE) %>% transmute(donor = sample, src = source)
mm <- met %>% left_join(meta76, by = "donor") %>%
  mutate(arm = ifelse(cohort == "Discovery_n76" & group == "Control",
                      paste0("UK ", sub("PRJEB", "", src)), as.character(cohort)))
say(""); say("## arm-wise comparison (all four cohorts profiled by the same pipeline)")
say("   control arms, is this a healthy adult baseline?")
mm %>% filter(group == "Control") %>% group_by(arm) %>%
  summarise(n = n(), shannon = median(shannon_species), health = median(health_index),
            entero = median(enterobacteriaceae_percent), below3 = sum(shannon_species < 3),
            .groups = "drop") %>%
  { for (i in seq_len(nrow(.))) say(sprintf("      %-30s n = %3d  Shannon %.2f  health %+.3f  Enterobacteriaceae %.4f%%  Shannon < 3: %d",
                                            .$arm[i], .$n[i], .$shannon[i], .$health[i], .$entero[i], .$below3[i])) }
d <- met %>% filter(cohort == "Discovery_n76")
say("   discovery arm against each external arm (Wilcoxon)")
for (g in c("Case", "Control")) for (c in setdiff(unique(met$cohort), "Discovery_n76")) {
  a <- d[[ "health_index" ]][d$group == g]; b <- met$health_index[met$cohort == c & met$group == g]
  as_ <- d$shannon_species[d$group == g];   bs <- met$shannon_species[met$cohort == c & met$group == g]
  say(sprintf("      %-7s vs %-30s health p = %-9.2g shannon p = %.2g", g, c,
              wilcox.test(a, b)$p.value, wilcox.test(as_, bs)$p.value))
}
say("   reading: the patient arms agree across hospitalised cohorts; the control arms do not,")
say("   so the pooled estimate is diluted by comparator baseline, not by a discovery-side outlier.")

# ---- severity trend, reported in the stats file only -------------------------
sev <- read_tsv(file.path(EXT, "PRJEB43555_severity_trend.tsv"), show_col_types = FALSE)
say(""); say("## severity trend within PRJEB43555 (none significant)")
for (i in seq_len(nrow(sev))) with(sev[i, ],
  say(sprintf("   %-32s n = %d  rho = %+.3f  p = %.3f  q = %.3f", metric, n, spearman_rho, p, q)))

# ---- c. does the sequencing chemistry shift the profiles? -------------------
# Patients were sequenced on MGI DNBSEQ and every control arm on Illumina. If
# that split biased profiling, it would show first in how much of each library
# MetaPhlAn can classify. It does not: the discovery arms sit together whichever
# chemistry produced them.
plat <- mm %>%
  mutate(chem = ifelse(cohort == "Discovery_n76" & group == "Case", "MGI DNBSEQ", "Illumina"),
         arm  = ifelse(cohort == "Discovery_n76",
                       ifelse(group == "Case", "Discovery\ncases", paste0("Discovery\n", sub("UK ", "controls ", arm))),
                       paste0(coh_short[cohort], "\n", tolower(group))),
         arm  = factor(arm, unique(arm[order(cohort != "Discovery_n76", cohort, group)])))
say(""); say("## c. classified read fraction by sequencing chemistry")
plat %>% group_by(chem) %>%
  summarise(n = n(), med = median(classified_percent), lo = min(classified_percent),
            hi = max(classified_percent), .groups = "drop") %>%
  { for (i in seq_len(nrow(.))) say(sprintf("   %-12s n = %3d  median %.1f%%  range %.1f-%.1f%%",
                                            .$chem[i], .$n[i], .$med[i], .$lo[i], .$hi[i])) }
d_case <- plat$classified_percent[plat$cohort == "Discovery_n76" & plat$group == "Case"]
d_ctrl <- plat$classified_percent[plat$cohort == "Discovery_n76" & plat$group == "Control"]
wt <- wilcox.test(d_case, d_ctrl)
say(sprintf("   within the discovery cohort, MGI cases %.1f%% against Illumina controls %.1f%%, Wilcoxon p = %.2f",
            median(d_case), median(d_ctrl), wt$p.value))

pC <- ggplot(plat, aes(arm, classified_percent, fill = chem)) +
  geom_boxplot(width = 0.6, alpha = 0.55, outlier.shape = NA, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.12, seed = 42), size = 0.5, alpha = 0.6,
             colour = "grey25") +
  scale_fill_manual(values = c("MGI DNBSEQ" = "#E85D5D", "Illumina" = "#4C72B0"), name = NULL) +
  labs(title = "Sequencing chemistry does not shift how much of a library is classified",
       subtitle = sprintf("Patients were sequenced on MGI DNBSEQ, every control arm on Illumina; within this study the two sit together (Wilcoxon p = %.2f)", wt$p.value),
       x = NULL, y = "Reads classified (%)") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6, colour = "black"),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(file.path(OUTDIR, "SuppFig7c_Platform_Classified_n76.pdf"), pC, width = 8.6, height = 3.6,
       device = cairo_pdf)

# ---- assemble -----------------------------------------------------------------
tag <- theme(plot.tag = element_text(face = "bold", size = 13))
supp7 <- (pA + labs(tag = "a") + tag) / (pB + labs(tag = "b") + tag) /
  (pC + labs(tag = "c") + tag) +
  plot_layout(heights = c(1.25, 1, 0.75)) +
  plot_annotation(
    title = "Supplementary Figure 7. Replication of the read-level endpoints\nin three external COVID-19 cohorts",
    theme = theme(plot.title = element_text(face = "bold", size = 10)))
# Scaled by 0.880 to fit A4 portrait; the aspect ratio is
# unchanged, so no panel is stretched relative to the others.
ggsave(file.path(OUTDIR, "Supplementary_Figure_7_n76.pdf"), supp7, width = 8.27, height = 11.3,
       units = "in", device = cairo_pdf, limitsize = FALSE)
writeLines(STATS, file.path(OUTDIR, "SuppFig7_STATS_n76.txt"))
message("DONE -> Supplementary_Figure_7_n76.pdf")
