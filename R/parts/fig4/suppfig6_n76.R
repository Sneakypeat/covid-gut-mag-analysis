# =============================================================================
# Supplementary Figure 6, CAMPER module robustness, taxa and carriage, plus the
# KO-level shifts behind the KEGG over-representation
#
#   a  every Figure 4b module (BH q < 0.05) re-tested against each control study
#   b  genera contributing to the module shifts (distinct module profiles only)
#   c  modules carried more often by enriched than depleted MAGs (Fisher, q < 0.05)
#   d  KO-level shifts, LinDA volcano (the input to the Figure 4c-d ORA)
#
# a-c back the CAMPER statements that lost their panel when Figure 4 was
# reduced to CAMPER a-b plus KEGG c-d. d moved here from the old Supplementary
# Figure 6, whose other panel (the raw-p < 0.05 pathway butterfly) was dropped
# as redundant with Figure 4c. Panels a and b are the plot objects
# written by camper_n76_taxonomy_adjustment.R; c is drawn here from
# camper_carriage_by_direction_n76.csv (camper_n76.R). The earlier carriage
# bubble showed only the enriched-minus-depleted skew; c shows both carriage
# fractions, which are the numbers the text quotes.
#
# The Enterobacteriaceae-conditional model from the same script stays an audit
# table, not a panel (conditioning on a likely mediator changes the estimand).
#
# Outputs -> result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(forcats)
  library(ggplot2); library(patchwork)
})

theme_transparent <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA))
theme_bw <- function(...) ggplot2::theme_bw(...) + theme_transparent
ggplot2::theme_set(ggplot2::theme_get() + theme_transparent)
ggsave <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
F4     <- file.path(N76, "fig4")
FOLLOW <- file.path(F4, "followup")
OUTDIR <- file.path(N76, "supp_fig")
source(file.path(BASE, "taxon_italics_n76.R"))   # md_taxon(): taxon names in italics
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
CARRY_THRESHOLD <- 0.5   # as in camper_n76.R
dir_cols <- c("Enriched" = "#E76F51", "Depleted" = "#457B9D", "NS" = "grey65")   # Figure 1e
tag <- theme(plot.tag = element_text(face = "bold", size = 13))

primary <- read_csv(file.path(F4, "camper_module_tests_n76.csv"), show_col_types = FALSE)
obj     <- readRDS(file.path(F4, "Fig4ab_CAMPER_plot_objects_n76.rds"))
tax     <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv")) %>%
  transmute(MAG = catalog_id, family, genus)

# ---- a. control-cohort robustness -------------------------------------------
src <- read_csv(file.path(FOLLOW, "module_control_source_sensitivity_n76.csv"), show_col_types = FALSE)
# the follow-up tables predate the last camper_n76.R run; stop if they drifted
chk <- src %>% distinct(module, primary_delta, primary_q) %>%
  inner_join(primary %>% select(module, delta, q), by = "module")
stopifnot(nrow(chk) == nrow(primary), max(abs(chk$primary_delta - chk$delta)) < 1e-12)
sig_src <- src %>% filter(primary_q < 0.05)
by_mod  <- sig_src %>% group_by(module) %>%
  summarise(same = all(same_direction_as_primary), all_q = all(q_within_source < 0.05), .groups = "drop")
say("## a. control-cohort robustness")
say(sprintf("   %d Figure 4b modules; same direction against all three control studies: %d; q < 0.05 in all three: %d",
            nrow(by_mod), sum(by_mod$same), sum(by_mod$all_q)))
sig_src %>% group_by(control_source, n_control) %>% summarise(n = sum(q_within_source < 0.05), .groups = "drop") %>%
  { for (i in seq_len(nrow(.))) say(sprintf("   %-10s (n = %2d controls) significant: %d", .$control_source[i], .$n_control[i], .$n[i])) }
# camper_n76_checks.R (CAMPER_CHECKS_n76.txt, the numbers in RESULTS §7) runs BH
# over the 55 significant modules only; camper_n76_followup.R, which colours
# this panel, runs it over all tested modules. Same p values, different family.
bh55 <- sig_src %>% group_by(control_source) %>% mutate(q55 = p.adjust(p, "BH")) %>% ungroup()
say(sprintf("   BH family: all %d tested modules -> %d in all three; the %d significant only -> %d in all three (per source: %s)",
            n_distinct(src$module), sum(by_mod$all_q), nrow(by_mod),
            bh55 %>% group_by(module) %>% summarise(ok = all(q55 < 0.05)) %>% pull(ok) %>% sum(),
            bh55 %>% group_by(control_source) %>% summarise(n = sum(q55 < 0.05)) %>%
              { paste(sprintf("%s %d", .$control_source, .$n), collapse = ", ") }))

pA <- obj$panel_a +
  labs(title = "Control-cohort robustness",
       subtitle = sprintf("Each point is one of the %d Figure 4b modules; dashed line is equal effect", nrow(by_mod)))
ggsave(file.path(OUTDIR, "SuppFig6a_Control_Source_n76.pdf"), pA, width = 7.5, height = 3.4,
       device = cairo_pdf)

# ---- b. taxa contributing to the shifts --------------------------------------
# Genus shares from the driver partition; Enterobacteriaceae genera are marked
# from the GTDB family of the catalogue MAGs, which ties the panel to the
# family-level partition quoted in the text.
ent <- read_csv(file.path(F4, "camper_entero_contribution_n76.csv"), show_col_types = FALSE) %>%
  inner_join(primary %>% filter(q < 0.05, delta > 0) %>% select(module), by = "module")
say(""); say("## b. Enterobacteriaceae share of the net case increase (camper_entero_contribution_n76.csv)")
say(sprintf("   higher-in-case modules: %d | median share %.1f%% | share > 50%% in %d",
            nrow(ent), 100 * median(ent$share), sum(ent$share > 0.5)))
ent_genera <- tax %>% filter(family == "Enterobacteriaceae") %>% distinct(genus) %>% pull(genus)
pB <- obj$panel_b +
  # genus in italics; the family marker stays a literal asterisk, escaped so the
  # markdown renderer does not read it as emphasis
  scale_y_discrete(labels = function(g) ifelse(g %in% ent_genera,
                                               paste0("*", g, "* \\*"), paste0("*", g, "*"))) +
  labs(title = "Taxa contributing to module shifts",
       subtitle = paste0(obj$panel_b$labels$subtitle, "<br>\\* *Enterobacteriaceae* (GTDB)")) +
  theme(axis.text.y = ggtext::element_markdown(),
        plot.subtitle = ggtext::element_markdown())
say(sprintf("   Enterobacteriaceae genera among the 12 plotted: %s",
            paste(intersect(levels(obj$panel_b$data$genus), ent_genera), collapse = ", ")))
ggsave(file.path(OUTDIR, "SuppFig6b_Taxon_Drivers_n76.pdf"), pB, width = 4.6, height = 4.2,
       device = cairo_pdf)

# ---- c. carriage by MAG direction -------------------------------------------
cw <- read_csv(file.path(F4, "camper_carriage_by_direction_n76.csv"), show_col_types = FALSE)
dir_df <- read.csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv")) %>%
  transmute(MAG = taxon,
            direction = case_when(diff_GroupCase %in% TRUE & lfc_GroupCase > 0 ~ "Enriched",
                                  diff_GroupCase %in% TRUE & lfc_GroupCase < 0 ~ "Depleted",
                                  TRUE ~ "NS"))
n_dir <- table(dir_df$direction)
sig_c <- cw %>% filter(q < 0.05)
carriers <- read_tsv(file.path(F4, "camper_module_completeness_n76.tsv"), show_col_types = FALSE) %>%
  filter(module %in% sig_c$module, completeness >= CARRY_THRESHOLD) %>%
  distinct(MAG) %>% inner_join(dir_df, by = "MAG") %>% left_join(tax, by = "MAG")
enr_ent <- dir_df %>% left_join(tax, by = "MAG") %>%
  filter(direction == "Enriched", family == "Enterobacteriaceae")
say(""); say("## c. carriage (completeness >= 0.5), Fisher enriched vs depleted MAGs")
say(sprintf("   MAGs: %d enriched, %d depleted, %d NS | modules tested %d, at q < 0.05: %d (all more often in enriched: %s)",
            n_dir[["Enriched"]], n_dir[["Depleted"]], n_dir[["NS"]], sum(!is.na(cw$q)), nrow(sig_c),
            all(sig_c$Enriched > sig_c$Depleted)))
say(sprintf("   enriched MAGs carrying >= 1 of these modules: %d of %d; enriched Enterobacteriaceae among them: %d of %d",
            sum(carriers$direction == "Enriched"), n_dir[["Enriched"]],
            sum(carriers$direction == "Enriched" & carriers$family %in% "Enterobacteriaceae"), nrow(enr_ent)))
carriers %>% filter(direction == "Enriched") %>% count(family, sort = TRUE) %>%
  { say(sprintf("   carrier families (enriched): %s", paste(sprintf("%s %d", .$family, .$n), collapse = "; "))) }
for (i in seq_len(nrow(sig_c))) with(sig_c[i, ],
  say(sprintf("   %-62s enriched %5.1f%%  depleted %4.1f%%  q = %.2g", module, 100 * Enriched, 100 * Depleted, q)))
write_csv(sig_c, file.path(OUTDIR, "SuppFig6c_carriage_q05_n76.csv"))

cdat <- sig_c %>%
  mutate(module = fct_reorder(module, Enriched - Depleted),
         lab = sprintf("%d/%d vs %d/%d, q = %s",
                       round(Enriched * n_dir[["Enriched"]]), n_dir[["Enriched"]],
                       round(Depleted * n_dir[["Depleted"]]), n_dir[["Depleted"]],
                       formatC(q, format = "e", digits = 1)))
clong <- cdat %>% select(module, Enriched, Depleted, NS) %>%
  pivot_longer(-module, names_to = "direction", values_to = "frac") %>%
  mutate(direction = factor(direction, c("NS", "Depleted", "Enriched")))
pC <- ggplot(cdat, aes(y = module)) +
  geom_segment(aes(x = 100 * Depleted, xend = 100 * Enriched, yend = module),
               colour = "grey70", linewidth = 0.6) +
  geom_point(data = clong %>% arrange(direction),
             aes(x = 100 * frac, colour = direction, size = direction)) +
  geom_text(aes(x = Inf, label = lab), hjust = 1.02, size = 2, colour = "grey25") +
  scale_colour_manual(values = dir_cols, name = NULL,
                      breaks = c("Enriched", "Depleted", "NS"),
                      labels = sprintf("%s (n = %d)", c("Enriched", "Depleted", "NS"),
                                       as.integer(n_dir[c("Enriched", "Depleted", "NS")]))) +
  scale_size_manual(values = c(Enriched = 2.2, Depleted = 2.2, NS = 1.3), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.62))) +
  labs(title = "Modules carried more often by enriched MAGs",
       subtitle = sprintf("%d of %d modules, Fisher's exact test (enriched vs depleted), BH q < 0.05; carried = completeness >= %.1f",
                          nrow(sig_c), sum(!is.na(cw$q)), CARRY_THRESHOLD),
       x = "MAGs carrying the module (%)", y = NULL) +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black", size = 6),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6, colour = "grey30"),
        legend.position = "bottom", legend.text = element_text(size = 6))
ggsave(file.path(OUTDIR, "SuppFig6c_Carriage_n76.pdf"), pC, width = 6.2, height = 4.2,
       device = cairo_pdf)

# ---- d. KO-level shifts (LinDA volcano) --------------------------------------
# Plot object from the HPC KEGG pipeline (gene_catalogue_linda_ora_n76.R); only
# the title is changed (its "c" prefix belonged to the old Figure 4 layout).
volc <- readRDS(file.path(F4, "kegg_fast/fig4/Fig4cd_KEGG_plot_objects_n76.rds"))
pD <- volc$panel_c + labs(title = "KEGG orthologue shifts") + theme_transparent
say(""); say(sprintf("## d. volcano: %s", pD$labels$subtitle))
ggsave(file.path(OUTDIR, "SuppFig6d_KO_Volcano_n76.pdf"), pD, width = 5, height = 4.2,
       device = cairo_pdf)

# ---- assemble -----------------------------------------------------------------
top    <- (pA + labs(tag = "a") + tag) | (pD + labs(tag = "d") + tag)
bottom <- (pB + labs(tag = "b") + tag) | (pC + labs(tag = "c") + tag)
supp6 <- (top + plot_layout(widths = c(1.7, 1))) /
  (bottom + plot_layout(widths = c(1, 1.45))) +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(
    title = "CAMPER module robustness, contributing taxa and carriage; KEGG orthologue shifts (n = 76)",
    theme = theme(plot.title = element_text(face = "bold", size = 12)))
# Scaled by 0.899 to fit A4 landscape; the aspect ratio is
# unchanged, so no panel is stretched relative to the others.
ggsave(file.path(OUTDIR, "Supplementary_Figure_6_n76.pdf"), supp6, width = 11.69, height = 7.91,
       units = "in", device = cairo_pdf, limitsize = FALSE)
writeLines(STATS, file.path(OUTDIR, "SuppFig6_STATS_n76.txt"))
message("DONE -> Supplementary_Figure_6_n76.pdf")
