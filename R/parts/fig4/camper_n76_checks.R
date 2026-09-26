# =============================================================================
# Verification of the CAMPER claims before they go into the Results section:
#   1. direction and significance against each of the three control cohorts
#   2. how much of each module's case-control difference Enterobacteriaceae carry
#   3. redundancy: modules whose per-MAG carriage profiles are exactly identical
#   4. sensitivity: unrenormalised abundance, and high-quality MAGs only
# -> result2/n76/fig4/CAMPER_CHECKS_n76.txt
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
OUTDIR <- file.path(BASE, "result2/n76/fig4")
OUT <- c(); say <- function(...) { l <- paste0(...); cat(l, "\n"); OUT <<- c(OUT, l) }
rd <- function(f) read.delim(file.path(INDIR, f), check.names = FALSE, quote = "", comment.char = "")
cliffs <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))

comp <- read_tsv(file.path(OUTDIR, "camper_module_completeness_n76.tsv"), show_col_types = FALSE)
tests <- read_csv(file.path(OUTDIR, "camper_module_tests_n76.csv"), show_col_types = FALSE)
sig_mods <- tests %>% filter(q < 0.05)
say(sprintf("significant modules carried forward: %d (%d higher, %d lower in Case)",
            nrow(sig_mods), sum(sig_mods$delta > 0), sum(sig_mods$delta < 0)))

tax <- rd("MAG_quality_taxonomy.tsv")
meta <- rd("sample_metadata_76.tsv") %>%
  transmute(donor = sample, Group = factor(group, levels = c("Control", "Case")), source)
ab <- rd("MAG_relative_abundance_percent.tsv"); names(ab)[1] <- "MAG"
ab_raw <- ab %>% pivot_longer(-MAG, names_to = "donor", values_to = "pct") %>%
  mutate(pct = as.numeric(pct) / 100)
ab_norm <- ab_raw %>% group_by(donor) %>% mutate(relab = pct / sum(pct)) %>% ungroup()

capacity <- function(abund, keep_mags = NULL) {
  cc <- comp %>% select(MAG, module, completeness)
  if (!is.null(keep_mags)) cc <- cc %>% filter(MAG %in% keep_mags)
  abund %>% inner_join(cc, by = "MAG", relationship = "many-to-many") %>%
    group_by(donor, module) %>% summarise(capacity = sum(relab * completeness), .groups = "drop") %>%
    complete(donor, module, fill = list(capacity = 0)) %>% left_join(meta, by = "donor")
}
test_modules <- function(cap, mods, ctrl_rows = NULL) {
  cap %>% filter(module %in% mods) %>%
    { if (is.null(ctrl_rows)) . else filter(., Group == "Case" | source %in% ctrl_rows) } %>%
    group_by(module) %>%
    summarise(delta = cliffs(capacity[Group == "Case"], capacity[Group == "Control"]),
              p = tryCatch(wilcox.test(capacity[Group == "Case"], capacity[Group == "Control"],
                                       exact = FALSE)$p.value, error = function(e) NA_real_),
              .groups = "drop") %>%
    mutate(q = p.adjust(p, "BH"))
}

cap_main <- capacity(ab_norm)
mods <- sig_mods$module
all_mods <- tests$module   # BH family = every module tested, as in camper_n76_followup.R
base_dir <- setNames(sign(sig_mods$delta), sig_mods$module)

# ---- 1. per control cohort --------------------------------------------------
say("")
say("## 1. Against each control cohort separately")
per <- list()
for (s in c("PRJEB7331", "PRJEB7949", "PRJEB39223")) {
  # corrected across all tested modules, then read off the significant ones, so
  # this count matches module_control_source_sensitivity_n76.csv exactly
  t <- test_modules(cap_main, all_mods, ctrl_rows = s) %>% filter(module %in% mods)
  same <- sum(sign(t$delta) == base_dir[t$module], na.rm = TRUE)
  sig  <- sum(t$q < 0.05 & sign(t$delta) == base_dir[t$module], na.rm = TRUE)
  n_ctrl <- sum(meta$source == s)
  say(sprintf("   %-11s (n = %2d controls): same direction %d/%d | still q < 0.05 %d",
              s, n_ctrl, same, nrow(t), sig))
  per[[s]] <- t %>% transmute(module, dir_ok = sign(delta) == base_dir[module], sig_ok = q < 0.05 & dir_ok)
}
all3_dir <- Reduce(`+`, lapply(per, function(x) as.integer(x$dir_ok[match(mods, x$module)])))
all3_sig <- Reduce(`+`, lapply(per, function(x) as.integer(x$sig_ok[match(mods, x$module)])))
say(sprintf("   consistent direction in all three cohorts: %d of %d", sum(all3_dir == 3, na.rm = TRUE), length(mods)))
say(sprintf("   individually significant in all three:      %d of %d", sum(all3_sig == 3, na.rm = TRUE), length(mods)))

# ---- 2. Enterobacteriaceae share of the difference --------------------------
say("")
say("## 2. Share of each module's case-control difference carried by Enterobacteriaceae")
fam <- tax %>% transmute(MAG = catalog_id, family, phylum)
contrib <- ab_norm %>%
  inner_join(comp %>% select(MAG, module, completeness), by = "MAG", relationship = "many-to-many") %>%
  filter(module %in% mods) %>%
  left_join(fam, by = "MAG") %>% left_join(meta, by = "donor") %>%
  mutate(entero = family == "Enterobacteriaceae", contrib = relab * completeness) %>%
  group_by(module, entero, Group) %>% summarise(mean_c = mean(contrib) * 76 / n_distinct(meta$donor), .groups = "drop") %>%
  group_by(module, entero) %>%
  summarise(diff = mean_c[Group == "Case"] - mean_c[Group == "Control"], .groups = "drop") %>%
  pivot_wider(names_from = entero, values_from = diff, names_prefix = "e_", values_fill = 0) %>%
  mutate(total = e_TRUE + e_FALSE, share = ifelse(total > 0, e_TRUE / total, NA_real_))
up <- contrib %>% filter(module %in% sig_mods$module[sig_mods$delta > 0])
say(sprintf("   modules higher in Case: %d | median Enterobacteriaceae share of the net increase: %.1f%%",
            nrow(up), 100 * median(up$share, na.rm = TRUE)))
say(sprintf("   share > 50%% in %d of %d", sum(up$share > 0.5, na.rm = TRUE), nrow(up)))
write_csv(contrib, file.path(OUTDIR, "camper_entero_contribution_n76.csv"))

# ---- 3. redundancy ----------------------------------------------------------
say("")
say("## 3. Redundancy among the significant modules")
wide <- comp %>% filter(module %in% mods) %>%
  select(MAG, module, completeness) %>%
  pivot_wider(names_from = module, values_from = completeness, values_fill = 0)
prof <- as.matrix(wide[, -1]); colnames(prof) <- names(wide)[-1]
keys <- apply(prof, 2, function(x) paste(round(x, 6), collapse = "|"))
grp <- split(names(keys), keys)
dupe_groups <- grp[lengths(grp) > 1]
say(sprintf("   %d significant modules collapse to %d distinct per-MAG profiles",
            length(keys), length(unique(keys))))
say(sprintf("   %d modules sit in %d groups of exactly identical profiles",
            sum(lengths(dupe_groups)), length(dupe_groups)))
for (g in dupe_groups) {
  n_mags <- sum(prof[, g[1]] > 0)
  say(sprintf("      %d modules carried by %d MAGs: %s", length(g), n_mags,
              paste(substr(g, 1, 46), collapse = "; ")))
}

# ---- 4. sensitivity ---------------------------------------------------------
say("")
say("## 4. Sensitivity")
t_raw <- test_modules(capacity(ab_raw %>% rename(relab = pct)), mods)
say(sprintf("   unrenormalised (total-community) abundance: %d of %d keep q < 0.05 and direction",
            sum(t_raw$q < 0.05 & sign(t_raw$delta) == base_dir[t_raw$module], na.rm = TRUE), length(mods)))
hq <- tax$catalog_id[tax$completeness >= 90 & tax$contamination <= 5]
t_hq <- test_modules(capacity(ab_norm, keep_mags = hq), mods)
say(sprintf("   high-quality MAGs only (n = %d genomes): %d of %d keep q < 0.05 and direction",
            length(hq), sum(t_hq$q < 0.05 & sign(t_hq$delta) == base_dir[t_hq$module], na.rm = TRUE), length(mods)))
lost_hq <- t_hq %>% filter(!(q < 0.05 & sign(delta) == base_dir[module]))
if (nrow(lost_hq)) {
  say("   modules losing significance under the high-quality restriction:")
  for (i in seq_len(nrow(lost_hq)))
    say(sprintf("      %-52s q = %.3f", substr(lost_hq$module[i], 1, 52), lost_hq$q[i]))
}
write_csv(t_raw, file.path(OUTDIR, "camper_sensitivity_unnormalised_n76.csv"))
write_csv(t_hq,  file.path(OUTDIR, "camper_sensitivity_highquality_n76.csv"))

writeLines(OUT, file.path(OUTDIR, "CAMPER_CHECKS_n76.txt"))
