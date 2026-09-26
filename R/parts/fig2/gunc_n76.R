# =============================================================================
# GUNC chimerism screen of the 584-MAG catalogue
#   GUNC 1.1.1 against proGenomes 2.1 (the reference its CSS thresholds were
#   calibrated on; gunc_db_v111 also holds progenomes_3 and gtdb_214)
#
# Two claims in the manuscript turn on this.
#
#   Changed claim 6. The 21 MAGs without a GTDB species are the worst genomes in
#   the catalogue. A bin merged from two unrelated organisms has no GTDB species
#   for the same reason it has no coherent taxonomy: it is not one genome.
#   CheckM2 cannot see this, because single-copy markers from a chimera still
#   look complete and uncontaminated.
#
#   Changed claim 7. Contamination predicted case enrichment and nothing
#   explained it. If chimeric bins are what CheckM2 contamination was partly
#   tracking, removing them should settle whether the association is biology or
#   assembly artefact.
#
# Outputs -> result2/n76/supp_fig/gunc/
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(phylolm)
  library(dplyr); library(tidyr); library(tibble); library(readr); library(purrr)
  library(ggplot2); library(patchwork)
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
OUTDIR <- file.path(N76, "supp_fig/gunc")
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
rd <- function(f) read.delim(file.path(INDIR, f), check.names = FALSE, quote = "", comment.char = "")

g    <- read_tsv(file.path(INDIR, "gunc_maxCSS_584.tsv"), show_col_types = FALSE)
tax  <- rd("MAG_quality_taxonomy.tsv")
st   <- readRDS(file.path(N76, "fig2", "figure2_stats_n76.rds"))
tree <- read.tree(TREE); tree$node.label <- NULL

# GUNC names each genome by splitting the FILENAME at the first literal ".fa".
# That truncates any catalog_id containing ".fa_sub_" while leaving ".fna_"
# untouched, because ".fna" does not contain the substring ".fa". Joining on the
# raw id silently matches only 206 of 584, so rebuild the key and assert.
gunc_key <- function(x) vapply(strsplit(paste0(x, ".fa"), ".fa", fixed = TRUE), `[`, "", 1)
tax$gkey <- gunc_key(tax$catalog_id)
stopifnot(length(unique(tax$gkey)) == nrow(tax),
          sum(tax$gkey %in% g$genome) == nrow(g))

d <- tax %>%
  transmute(catalog_id, gkey, arm = group, family, phylum, species,
            completeness, contamination) %>%
  inner_join(g, by = c("gkey" = "genome")) %>%
  mutate(novel = is.na(species) | trimws(species) == "",
         hq    = completeness > 90 & contamination < 5)
stopifnot(nrow(d) == 584)

say("GUNC 1.1.1, proGenomes 2.1, 584 MAGs")
say(sprintf("   fail GUNC: %d of 584 (%.1f%%); clean set %d",
            sum(!d$pass.GUNC), 100 * mean(!d$pass.GUNC), sum(d$pass.GUNC)))
say(sprintf("   median clade-separation score %.3f, median contamination portion %.3f",
            median(d$clade_separation_score), median(d$contamination_portion)))
arm <- d %>% group_by(arm) %>% summarise(n = n(), fail = sum(!pass.GUNC),
                                         pct = 100 * mean(!pass.GUNC), .groups = "drop")
for (i in seq_len(nrow(arm))) with(arm[i, ],
  say(sprintf("   %-8s-derived: %d fail of %d (%.1f%%)", arm, fail, n, pct)))
say("   an even split across arms means this is a property of the binning, not of a cohort")

# ---- changed claim 6 --------------------------------------------------------

say(""); say("## novel MAGs and chimerism")
nv <- d %>% filter(novel)
tt <- d %>% group_by(novel) %>% summarise(n = n(), fail = sum(!pass.GUNC),
                                          pct = 100 * mean(!pass.GUNC),
                                          css = median(clade_separation_score), .groups = "drop")
for (i in seq_len(nrow(tt))) with(tt[i, ],
  say(sprintf("   %-14s %3d MAGs | fail %3d (%.1f%%) | median CSS %.3f",
              ifelse(novel, "no GTDB sp.", "known sp."), n, fail, pct, css)))
say(sprintf("   Fisher p = %.3g", fisher.test(table(d$novel, d$pass.GUNC))$p.value))
say(sprintf("   of the %d novel MAGs: %d pass GUNC, %d are >90%% complete and <5%% contaminated, %d are both",
            nrow(nv), sum(nv$pass.GUNC), sum(nv$hq), sum(nv$pass.GUNC & nv$hq)))
keep <- nv %>% filter(pass.GUNC, hq)
if (nrow(keep)) {
  say("   the defensible novel-taxa set:")
  for (i in seq_len(nrow(keep))) with(keep[i, ],
    say(sprintf("      %-22s %-16s completeness %.1f%%, CheckM2 contamination %.2f%%, CSS %.3f",
                family, phylum, completeness, contamination, clade_separation_score)))
}
write_csv(d, file.path(OUTDIR, "gunc_584_annotated_n76.csv"))

# ---- changed claim 7 --------------------------------------------------------

say(""); say("## contamination, novelty and chimerism")
for (v in c("clade_separation_score", "contamination_portion", "n_effective_surplus_clades")) {
  r <- suppressWarnings(cor.test(d$contamination, d[[v]], method = "spearman"))
  say(sprintf("   Spearman(CheckM2 contamination, %-27s) rho = %+.3f, p = %.3g",
              v, r$estimate, r$p.value))
}

md <- st$meta_df %>%
  left_join(d %>% select(Genome_ID = catalog_id, CSS = clade_separation_score,
                         cont_portion = contamination_portion, passG = pass.GUNC),
            by = "Genome_ID") %>%
  filter(!is.na(CSS))
say(sprintf("   tree-placed MAGs with a GUNC result: %d (%d GUNC-clean)",
            nrow(md), sum(md$passG)))

fit_one <- function(dat, label, formula) {
  dat <- as.data.frame(dat); rownames(dat) <- dat$Genome_ID
  tr <- keep.tip(tree, intersect(tree$tip.label, dat$Genome_ID))
  dd <- dat[tr$tip.label, , drop = FALSE]
  f <- tryCatch(phyloglm(formula, phy = tr, data = dd, method = "logistic_MPLE"),
                error = function(e) { message("  phyloglm failed: ", conditionMessage(e)); NULL })
  if (is.null(f)) return(NULL)
  cf <- summary(f)$coefficients
  pc <- if ("p.value" %in% colnames(cf)) "p.value" else colnames(cf)[ncol(cf)]
  tibble(model = label, term = rownames(cf), beta = cf[, "Estimate"],
         se = cf[, "StdErr"], p = cf[, pc], n = nrow(dd))
}
models <- bind_rows(
  fit_one(md, "All MAGs", enriched_case ~ Is_Novel_bin + Completeness + Contamination),
  fit_one(md, "Adjusted for CSS", enriched_case ~ Is_Novel_bin + Completeness + Contamination + CSS),
  fit_one(md %>% filter(passG), "GUNC-clean only",
          enriched_case ~ Is_Novel_bin + Completeness + Contamination))
write_csv(models, file.path(OUTDIR, "gunc_phyloglm_models_n76.csv"))
say("")
for (i in seq_len(nrow(models))) with(models[i, ],
  if (term != "(Intercept)")
    say(sprintf("   %-18s %-14s beta %+.4f  p = %.3g  (n = %d)", model, term, beta, p, n)))
gb <- function(m, t) models %>% filter(model == m, term == t) %>% pull(beta)
gp <- function(m, t) models %>% filter(model == m, term == t) %>% pull(p)
say("")
say(sprintf("   contamination: %+.4f (p = %.3g) -> %+.4f (p = %.3g) once chimeric genomes are dropped",
            gb("All MAGs", "Contamination"), gp("All MAGs", "Contamination"),
            gb("GUNC-clean only", "Contamination"), gp("GUNC-clean only", "Contamination")))
say(sprintf("   novelty:       %+.4f (p = %.3g) -> %+.4f (p = %.3g)",
            gb("All MAGs", "Is_Novel_bin"), gp("All MAGs", "Is_Novel_bin"),
            gb("GUNC-clean only", "Is_Novel_bin"), gp("GUNC-clean only", "Is_Novel_bin")))
say("   note: adjusting FOR the clade-separation score does not do this. Chimerism is")
say("   not a linear nuisance term; those genomes are wrong, so exclusion beats adjustment.")

# ---- figure -----------------------------------------------------------------

pA <- d %>% count(arm, pass.GUNC) %>%
  ggplot(aes(arm, n, fill = pass.GUNC)) +
  geom_col(position = "fill", width = 0.6) +
  scale_fill_manual(values = c("FALSE" = "#D55E00", "TRUE" = "grey75"),
                    labels = c("Fails GUNC", "Passes"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Chimerism by discovery arm",
       subtitle = sprintf("%d of 584 genomes fail (%.1f%%)", sum(!d$pass.GUNC), 100 * mean(!d$pass.GUNC)),
       x = "MAG discovered in", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "GUNC_a_arm_n76.pdf"), pA, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pB <- d %>% mutate(grp = ifelse(novel, "No GTDB species", "Known species")) %>%
  ggplot(aes(grp, clade_separation_score, fill = grp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.55, linewidth = 0.35, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.15, seed = 42), size = 0.6, alpha = 0.4) +
  scale_fill_manual(values = c("Known species" = "grey75", "No GTDB species" = "#D55E00")) +
  labs(title = "Novelty is largely chimerism",
       subtitle = sprintf("%d of %d novel MAGs fail GUNC (%.0f%%), against %.0f%% of known species",
                          sum(nv$pass.GUNC == FALSE), nrow(nv),
                          100 * mean(!nv$pass.GUNC),
                          100 * mean(!d$pass.GUNC[!d$novel])),
       x = NULL, y = "Clade separation score") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "GUNC_b_novel_n76.pdf"), pB, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pC <- ggplot(d, aes(clade_separation_score, contamination)) +
  geom_point(aes(colour = pass.GUNC), size = 1.2, alpha = 0.7) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = c("FALSE" = "#D55E00", "TRUE" = "grey70"),
                      labels = c("Fails GUNC", "Passes"), name = NULL) +
  labs(title = "CheckM2 contamination tracks chimerism",
       subtitle = sprintf("Spearman rho = %+.3f",
                          cor(d$contamination, d$clade_separation_score, method = "spearman")),
       x = "Clade separation score", y = "CheckM2 contamination (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "GUNC_c_contamination_n76.pdf"), pC, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pD <- models %>% filter(term %in% c("Contamination", "Is_Novel_bin")) %>%
  mutate(term = recode(term, Is_Novel_bin = "Novel species", Contamination = "Contamination"),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         model = factor(model, levels = unique(models$model))) %>%
  ggplot(aes(beta, model, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.35) +
  facet_wrap(~ term, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = c("grey25", "#3B6FB6", "#5AAE61"), name = NULL) +
  scale_y_discrete(labels = NULL) +
  labs(title = "Both awkward findings are chimera artefacts",
       subtitle = "Phylogenetic logistic regression for case enrichment",
       x = "Coefficient (log odds)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(nrow = 3))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "GUNC_d_phyloglm_n76.pdf"), pD, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pl <- list(GUNC_a_arm = pA, GUNC_b_novel = pB, GUNC_c_contamination = pC, GUNC_d_phyloglm = pD)
for (nm in names(pl))
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), pl[[nm]], width = 5.2, height = 4.6,
         device = cairo_pdf, bg = "transparent")

tag <- theme(plot.tag = element_text(face = "bold", size = 13))
fig <- ((pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag)) /
       ((pC + labs(tag = "c") + tag) | (pD + labs(tag = "d") + tag)) +
  plot_annotation(title = "GUNC chimerism screen of the 584-MAG catalogue",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUTDIR, "GUNC_Figure_n76.pdf"), fig, width = 11, height = 9,
       device = cairo_pdf, bg = "transparent")

writeLines(STATS, file.path(OUTDIR, "GUNC_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "GUNC_Figure_n76.pdf"))
