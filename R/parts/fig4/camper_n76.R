# =============================================================================
# CAMPER polyphenol-transformation capacity, 76-donor cohort (new figure)
#
# CAMPER 1.0.13 module counts per MAG -> module completeness per MAG ->
# abundance-weighted capacity per donor -> Case vs Control.
#
# Completeness respects CAMPER's step structure. The subheader vocabulary is
# inconsistent (Step/step, Subunit/subunit, subunit-only rows, "route N step M",
# "step 1 v2", "step 1a", "enzyme (+)/(-)"), so each row is parsed into:
#   route   - "route N" marks an alternative full pathway; otherwise one route
#   step    - "step N"; rows with no step number belong to step 1
#   variant - "vN" or "enzyme (+/-)": interchangeable enzymes for one step
#   part    - "subunit X" or the letter in "step 1a": all parts are required
# A step is satisfied when some variant has every part. Route completeness is
# satisfied steps / steps in that route; module completeness is the best route.
# Counting genes instead would let one subunit of a multi-subunit enzyme score
# the step, and would let alternative routes double-count it.
#
# Donor capacity is threshold-free: sum over MAGs of relative abundance x
# completeness. The >= 0.5 carriage threshold is used ONLY for the panel that
# compares enriched with depleted MAGs, where a binary call is needed.
#
# Outputs -> result2/n76/fig4/  (CAMPER panels of the new Figure 4)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(forcats)
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

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "fig4")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

STATS <- c()
say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
rd  <- function(f) read.delim(file.path(INDIR, f), check.names = FALSE, quote = "", comment.char = "")

cols    <- c("Control" = "#00BFC4", "Case" = "#F8766D")
OXY_COL <- c("anoxic" = "#2C5F8A", "oxic" = "#E08A1E", "both" = "#8C8C8C")
CARRY_THRESHOLD <- 0.5
cliffs <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))

# ---- 1. module completeness per MAG -----------------------------------------

cm <- rd("CAMPER_module_counts.tsv")
def_cols <- c("gene_id", "gene_description", "module", "header", "subheader",
              "specific_reaction", "oxygen", "EC", "Notes")
mag_cols <- setdiff(names(cm), def_cols)
say("CAMPER definition rows: ", nrow(cm), " ; MAGs: ", length(mag_cols))

defs <- cm %>%
  select(all_of(def_cols)) %>%
  mutate(row_id  = row_number(),
         sh      = str_squish(tolower(subheader)),
         route   = coalesce(str_match(sh, "route\\s*(\\d+)")[, 2], "all"),
         step    = coalesce(as.integer(str_match(sh, "step\\s*(\\d+)")[, 2]), 1L),
         variant = coalesce(str_match(sh, "\\bv(\\d+)\\b")[, 2],
                            str_match(sh, "enzyme\\s*\\(([+-])\\)")[, 2], "1"),
         part    = coalesce(str_match(sh, "subunit\\s*([a-z0-9]+)")[, 2],
                            str_match(sh, "step\\s*\\d+([a-z])\\b")[, 2], "-"),
         cls     = str_split_fixed(header, ";", 4)[, 2])
stopifnot(!anyNA(defs$step))
say("routed modules: ", n_distinct(defs$module[defs$route != "all"]))

hits <- cm %>%
  select(all_of(mag_cols)) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = "MAG", values_to = "n") %>%
  mutate(hit = as.numeric(n) > 0) %>%
  select(-n)

step_tbl <- hits %>%
  left_join(defs %>% select(row_id, module, route, step, variant, part), by = "row_id") %>%
  group_by(MAG, module, route, step, variant, part) %>% summarise(ok = any(hit), .groups = "drop") %>%
  group_by(MAG, module, route, step, variant)       %>% summarise(ok = all(ok), .groups = "drop") %>%
  group_by(MAG, module, route, step)                %>% summarise(ok = any(ok), .groups = "drop")

completeness <- step_tbl %>%
  group_by(MAG, module, route) %>%
  summarise(n_steps = n(), steps_ok = sum(ok), .groups = "drop") %>%
  mutate(route_c = steps_ok / n_steps) %>%
  group_by(MAG, module) %>%
  slice_max(route_c, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(MAG, module, best_route = route, n_steps, steps_ok, completeness = route_c)

mod_meta <- defs %>%
  group_by(module) %>%
  summarise(oxygen = first(oxygen), cls = first(cls),
            n_steps = n_distinct(paste(route, step)), .groups = "drop")
write_tsv(completeness, file.path(OUTDIR, "camper_module_completeness_n76.tsv"))
say("modules: ", nrow(mod_meta), " ; MAG x module pairs with completeness > 0: ",
    sum(completeness$completeness > 0))

# ---- 2. abundance-weighted capacity per donor ------------------------------

ab <- rd("MAG_relative_abundance_percent.tsv")
names(ab)[1] <- "MAG"
# Renormalised within each donor to the MAG-resolved community. Catalogue
# mapping rate ranges 27-97% across case donors, so raw percentages would let
# mapping rate masquerade as capacity. Capacity is therefore the
# abundance-weighted completeness of the community the catalogue resolves.
ab_long <- ab %>% pivot_longer(-MAG, names_to = "donor", values_to = "relab") %>%
  mutate(relab = as.numeric(relab)) %>%
  group_by(donor) %>% mutate(relab = relab / sum(relab)) %>% ungroup()

meta <- rd("sample_metadata_76.tsv") %>%
  transmute(donor = sample, Group = factor(group, levels = c("Control", "Case")), source)
stopifnot(nrow(meta) == 76)

cap_module <- ab_long %>%
  inner_join(completeness %>% select(MAG, module, completeness), by = "MAG",
             relationship = "many-to-many") %>%
  group_by(donor, module) %>%
  summarise(capacity = sum(relab * completeness), .groups = "drop") %>%
  complete(donor, module, fill = list(capacity = 0)) %>%
  left_join(meta, by = "donor") %>%
  left_join(mod_meta, by = "module")
write_tsv(cap_module, file.path(OUTDIR, "camper_donor_module_capacity_n76.tsv"))

# Per-MAG mean completeness within an oxygen class, then abundance-weighted:
# "what fraction of the community's polyphenol routes of this kind are intact".
mag_by <- function(key) {
  completeness %>%
    left_join(mod_meta %>% select(module, key = all_of(key)), by = "module") %>%
    group_by(MAG, key) %>% summarise(mean_c = mean(completeness), .groups = "drop") %>%
    complete(MAG = mag_cols, key, fill = list(mean_c = 0))
}
agg_capacity <- function(key) {
  ab_long %>%
    inner_join(mag_by(key), by = "MAG", relationship = "many-to-many") %>%
    group_by(donor, key) %>% summarise(capacity = sum(relab * mean_c), .groups = "drop") %>%
    left_join(meta, by = "donor")
}
cap_oxy <- agg_capacity("oxygen")
cap_cls <- agg_capacity("cls")
# Drop classes no MAG in this catalogue carries (stilbenes): an all-zero facet
# has no test and would print NaN.
empty_cls <- cap_cls %>% group_by(key) %>% summarise(tot = sum(capacity)) %>% filter(tot == 0) %>% pull(key)
if (length(empty_cls)) say("compound classes with zero capacity in every donor, dropped: ", paste(empty_cls, collapse = ", "))
cap_cls <- cap_cls %>% filter(!key %in% empty_cls)

test_by <- function(df) {
  df %>% group_by(key) %>%
    summarise(med_ctrl = median(capacity[Group == "Control"]),
              med_case = median(capacity[Group == "Case"]),
              W = unname(wilcox.test(capacity[Group == "Case"], capacity[Group == "Control"], exact = FALSE)$statistic),
              p = wilcox.test(capacity[Group == "Case"], capacity[Group == "Control"], exact = FALSE)$p.value,
              delta = cliffs(capacity[Group == "Case"], capacity[Group == "Control"]),
              .groups = "drop") %>%
    mutate(q = p.adjust(p, "BH"))
}
t_oxy <- test_by(cap_oxy); t_cls <- test_by(cap_cls)
write_csv(t_oxy, file.path(OUTDIR, "camper_oxygen_tests_n76.csv"))
write_csv(t_cls, file.path(OUTDIR, "camper_class_tests_n76.csv"))
say("")
say("-- abundance-weighted capacity by OXYGEN requirement --")
for (i in seq_len(nrow(t_oxy))) with(t_oxy[i, ], say(sprintf(
  "   %-7s Control %.4f | Case %.4f | W=%.0f p=%.3g q=%.3g delta=%+.2f", key, med_ctrl, med_case, W, p, q, delta)))
say("-- by COMPOUND CLASS --")
for (i in seq_len(nrow(t_cls))) with(t_cls[i, ], say(sprintf(
  "   %-20s Control %.4f | Case %.4f | W=%.0f p=%.3g q=%.3g delta=%+.2f", key, med_ctrl, med_case, W, p, q, delta)))

# ---- 3. per-module tests ----------------------------------------------------

t_mod <- cap_module %>%
  group_by(module, oxygen, cls, n_steps) %>%
  filter(sum(capacity > 0) >= 10) %>%
  summarise(med_ctrl = median(capacity[Group == "Control"]),
            med_case = median(capacity[Group == "Case"]),
            p = wilcox.test(capacity[Group == "Case"], capacity[Group == "Control"], exact = FALSE)$p.value,
            delta = cliffs(capacity[Group == "Case"], capacity[Group == "Control"]),
            .groups = "drop") %>%
  mutate(q = p.adjust(p, "BH")) %>%
  arrange(delta)
write_csv(t_mod, file.path(OUTDIR, "camper_module_tests_n76.csv"))
say("")
say(sprintf("modules tested (non-zero capacity in >= 10 donors): %d ; q < 0.05: %d (%d lower in Case, %d higher)",
            nrow(t_mod), sum(t_mod$q < 0.05),
            sum(t_mod$q < 0.05 & t_mod$delta < 0), sum(t_mod$q < 0.05 & t_mod$delta > 0)))
sig_by_oxy <- t_mod %>% filter(q < 0.05) %>% count(oxygen, dir = ifelse(delta < 0, "lower_in_case", "higher_in_case"))
for (i in seq_len(nrow(sig_by_oxy))) with(sig_by_oxy[i, ], say(sprintf("   %-7s %-15s %d", oxygen, dir, n)))

# ---- 4. carriage by MAG direction ------------------------------------------

res <- read.csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), stringsAsFactors = FALSE)
dir_df <- res %>% transmute(MAG = taxon,
                            direction = case_when(diff_GroupCase %in% TRUE & lfc_GroupCase > 0 ~ "Enriched",
                                                  diff_GroupCase %in% TRUE & lfc_GroupCase < 0 ~ "Depleted",
                                                  TRUE ~ "NS"))
carry <- completeness %>%
  complete(MAG = mag_cols, module, fill = list(completeness = 0)) %>%
  mutate(carries = completeness >= CARRY_THRESHOLD) %>%
  inner_join(dir_df, by = "MAG")

carry_tab <- carry %>%
  group_by(module, direction) %>% summarise(frac = mean(carries), k = sum(carries), n = n(), .groups = "drop")
fisher_mod <- carry %>% filter(direction != "NS") %>%
  group_by(module) %>%
  summarise(p = tryCatch(fisher.test(table(factor(carries, c(FALSE, TRUE)),
                                           factor(direction, c("Depleted", "Enriched"))))$p.value,
                         error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(q = p.adjust(p, "BH"))
carry_n <- carry_tab %>% group_by(module) %>% summarise(n_carrying = sum(k), .groups = "drop")
carry_wide <- carry_tab %>%
  select(module, direction, frac) %>%
  pivot_wider(names_from = direction, values_from = frac) %>%
  left_join(carry_n, by = "module") %>%
  left_join(fisher_mod, by = "module") %>%
  left_join(mod_meta, by = "module")
write_csv(carry_wide, file.path(OUTDIR, "camper_carriage_by_direction_n76.csv"))
say("")
say(sprintf("carriage (completeness >= %.1f): modules differing Enriched vs Depleted MAGs at q < 0.05: %d",
            CARRY_THRESHOLD, sum(carry_wide$q < 0.05, na.rm = TRUE)))

# ---- 5. figure --------------------------------------------------------------

theme_c <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey40", size = 8), strip.background = element_rect(fill = "grey95"),
        axis.text = element_text(color = "black"))

lab_q <- function(q) ifelse(q < 0.001, "q < 0.001", sprintf("q = %.3f", q))

# a. capacity by oxygen requirement and compound class, in one panel
# Previously nine separate boxplot facets (three oxygen, six class). Merged to a
# single horizontal panel: the categories share an axis, so they can be read
# against each other instead of nine independent y-scales. sqrt x because
# aromatic hydrocarbon capacity is two orders below non-specific.
cap_all <- bind_rows(
  cap_oxy %>% mutate(block = "Oxygen requirement"),
  cap_cls %>% mutate(block = "Compound class")) %>%
  mutate(key = gsub("_", " ", key),
         block = factor(block, c("Oxygen requirement", "Compound class")))
lab_all <- bind_rows(
  t_oxy %>% mutate(block = "Oxygen requirement"),
  t_cls %>% mutate(block = "Compound class")) %>%
  mutate(key = gsub("_", " ", key),
         block = factor(block, c("Oxygen requirement", "Compound class")),
         lab = sprintf("%s, %+.2f", lab_q(q), delta))
key_order <- cap_all %>% group_by(key) %>% summarise(m = median(capacity), .groups = "drop") %>%
  arrange(m) %>% pull(key)
cap_all$key <- factor(cap_all$key, key_order)
lab_all$key <- factor(lab_all$key, key_order)

pA <- ggplot(cap_all, aes(capacity, key, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.62, alpha = 0.45, outlier.shape = NA, linewidth = 0.4,
               position = position_dodge(width = 0.72)) +
  geom_point(position = position_jitterdodge(jitter.height = 0.14, dodge.width = 0.72, seed = 42),
             size = 0.75, alpha = 0.65) +
  geom_text(data = lab_all, aes(x = Inf, y = key, label = lab), inherit.aes = FALSE,
            hjust = 1.03, size = 2.5, colour = "grey25") +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = cols, name = NULL) + scale_colour_manual(values = cols, name = NULL) +
  scale_x_sqrt(expand = expansion(mult = c(0.02, 0.22))) +
  labs(title = "Polyphenol-transformation capacity by oxygen requirement and compound class",
       subtitle = "Per donor: MAG-resolved relative abundance x mean module completeness, summed. Square-root axis.",
       x = "Abundance-weighted completeness", y = NULL) +
  theme_c + theme(legend.position = "bottom", strip.placement = "outside",
                  strip.text.y.left = element_text(angle = 0, size = 8, face = "bold"),
                  panel.grid.major.y = element_blank())

# b. modules: effect size, carriage breadth and carriage skew in one panel
# Previously a lollipop (c) and a carriage heatmap (d) sharing the same module
# axis. Merged: position is the effect size, bubble area is how many MAGs carry
# the module, fill is the enriched-minus-depleted carriage difference.
fc <- t_mod %>% filter(q < 0.05) %>%
  left_join(carry_wide %>% select(module, n_carrying, Enriched, Depleted), by = "module") %>%
  mutate(carry_skew = Enriched - Depleted,
         module = fct_reorder(module, delta), cls = gsub("_", " ", cls))
say(sprintf("merged module panel: %d modules, carriage %d-%d MAGs, skew %+.2f to %+.2f",
            nrow(fc), min(fc$n_carrying, na.rm = TRUE), max(fc$n_carrying, na.rm = TRUE),
            min(fc$carry_skew, na.rm = TRUE), max(fc$carry_skew, na.rm = TRUE)))

pB <- ggplot(fc, aes(delta, module)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = delta, yend = module), linewidth = 0.4, colour = "grey65") +
  geom_point(aes(size = n_carrying, fill = carry_skew), shape = 21, colour = "grey20", stroke = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey96", high = "#B2182B", midpoint = 0,
                       name = "Carriage skew\n(enriched - depleted)", labels = scales::percent) +
  scale_size_continuous(range = c(1.5, 7), name = "MAGs carrying") +
  facet_grid(cls ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Modules differing between groups (BH q < 0.05)",
       subtitle = "Position: Cliff's delta, negative = lower capacity in Case. Size: MAGs carrying. Fill: carriage skew across MAG response.",
       x = "Cliff's delta", y = NULL) +
  theme_c + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 7),
                  axis.text.y = element_text(size = 7), legend.position = "right")

tag <- theme(plot.tag = element_text(face = "bold", size = 13))
# free() stops panel c's outside facet strips from widening the left margin of a and b
figC <- (free(pA + labs(tag = "a") + tag) / (pB + labs(tag = "b") + tag)) +
  plot_layout(heights = c(1, 1.9)) +
  plot_annotation(title = "CAMPER polyphenol-transformation capacity (76 donors, 584 MAGs)",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUTDIR, "Figure_CAMPER_n76.pdf"), figC, width = 11, height = 13,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)
ggsave(file.path(OUTDIR, "FigCAMPER_a_n76.pdf"), pA, width = 9, height = 5.2,
       device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "FigCAMPER_b_n76.pdf"), pB, width = 9, height = 9,
       device = cairo_pdf, bg = "transparent")
# handed to figure4_n76.R as Figure 4 panels a and b
saveRDS(list(panel_a = pA, panel_b = pB),
        file.path(OUTDIR, "Fig4ab_CAMPER_merged_panels_n76.rds"))

writeLines(STATS, file.path(OUTDIR, "CAMPER_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Figure_CAMPER_n76.pdf"))
