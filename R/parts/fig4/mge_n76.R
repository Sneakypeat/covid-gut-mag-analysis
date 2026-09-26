# =============================================================================
# Mobilome of the 76-donor MAG catalogue (geNomad 1.12.0)
#
# Asks three things, in order of how much they would change the manuscript:
#
#   Q1  Does mobile-element content explain the unexplained contamination
#       result of changed-claim 7? CheckM2 scores duplicated and foreign
#       sequence as contamination, and a mis-binned plasmid or free phage
#       contig is exactly that. If MGE load absorbs the contamination
#       coefficient, the "assembly artefact" reading becomes a biological one.
#
#   Q2  Do case-enriched genomes carry a heavier mobilome than depleted ones,
#       at matched completeness and phylogeny?
#
#   Q3  At community level, weighting each MAG by its abundance, is the COVID
#       gut mobilome-richer than the control gut? This is the only one of the
#       three that is a finding rather than a control.
#
# Three MGE classes are kept apart throughout: provirus (phage integrated into
# a host contig), free_virus (a whole contig called viral but binned into the
# MAG) and plasmid. Only the first is prophage burden in the usual sense; the
# other two are also the direct candidates for inflating contamination.
#
# Outputs -> result2/n76/mge/
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
OUTDIR <- file.path(N76, "mge")
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
rd <- function(f) read.delim(file.path(INDIR, f), check.names = FALSE, quote = "", comment.char = "")
cols <- c("Control" = "#00BFC4", "Case" = "#F8766D")

mge  <- read_tsv(file.path(INDIR, "mge_per_mag_584.tsv"), show_col_types = FALSE)
tax  <- rd("MAG_quality_taxonomy.tsv")
meta <- rd("sample_metadata_76.tsv") %>%
  transmute(donor = sample, Group = factor(group, levels = c("Control", "Case")), source)
st   <- readRDS(file.path(N76, "fig2", "figure2_stats_n76.rds"))
tree <- read.tree(TREE); tree$node.label <- NULL

ab <- rd("MAG_relative_abundance_percent.tsv"); names(ab)[1] <- "MAG"
ab_long <- ab %>% pivot_longer(-MAG, names_to = "donor", values_to = "pct") %>%
  mutate(pct = as.numeric(pct) / 100) %>%
  group_by(donor) %>% mutate(relab = pct / sum(pct)) %>% ungroup()

d <- tax %>%
  transmute(catalog_id, source_arm = group, family, phylum, genus,
            completeness, contamination) %>%
  inner_join(mge, by = "catalog_id") %>%
  mutate(genome_mb        = genome_bp / 1e6,
         provirus_per_mb  = n_provirus / genome_mb,
         mge_per_mb       = n_mge / genome_mb,
         mge_pct          = 100 * mge_frac,
         entero           = family == "Enterobacteriaceae")
stopifnot(nrow(d) == 584)

say(sprintf("584 MAGs, %.2f Gbp, geNomad 1.12.0", sum(d$genome_bp) / 1e9))
say("")
say("## mobilome of the catalogue")
for (k in c("provirus", "free_virus", "plasmid")) {
  n  <- d[[paste0("n_", k)]]; bp <- d[[paste0(k, "_bp")]]
  say(sprintf("   %-11s %5d regions in %3d of 584 MAGs (%.1f%%); median in carriers %d; %.1f Mbp total",
              k, sum(n), sum(n > 0), 100 * mean(n > 0), median(n[n > 0]), sum(bp) / 1e6))
}
say(sprintf("   any MGE: %d of 584 MAGs (%.1f%%); median mobilome %.2f%% of genome (IQR %.2f-%.2f)",
            sum(d$n_mge > 0), 100 * mean(d$n_mge > 0), 100 * median(d$mge_frac),
            100 * quantile(d$mge_frac, .25), 100 * quantile(d$mge_frac, .75)))

# =============================================================================
# Q1. does the mobilome explain CheckM2 contamination?
# =============================================================================

say(""); say("## Q1. mobilome vs CheckM2 contamination")
ct <- function(x, y, lab) {
  r <- suppressWarnings(cor.test(x, y, method = "spearman"))
  say(sprintf("   %-34s Spearman rho = %+.3f, p = %.3g", lab, r$estimate, r$p.value))
  r
}
ct(d$mge_frac,        d$contamination, "mobilome fraction")
ct(d$plasmid_frac,    d$contamination, "plasmid fraction")
ct(d$free_virus_frac, d$contamination, "free-virus fraction")
ct(d$provirus_frac,   d$contamination, "provirus fraction")

m_cont <- lm(contamination ~ completeness + plasmid_frac + free_virus_frac + provirus_frac, data = d)
say(sprintf("   lm(contamination ~ completeness + 3 MGE classes): adjusted R2 = %.3f",
            summary(m_cont)$adj.r.squared))
co <- summary(m_cont)$coefficients
for (i in 2:nrow(co))
  say(sprintf("      %-16s beta %+.4g  p = %.3g", rownames(co)[i], co[i, 1], co[i, 4]))

# =============================================================================
# Q2. mobilome and case enrichment, phylogeny-aware
# =============================================================================

say(""); say("## Q2. mobilome and case enrichment")
md <- st$meta_df %>%
  left_join(d %>% select(catalog_id, n_provirus, n_free_virus, n_plasmid,
                         provirus_frac, plasmid_frac, free_virus_frac,
                         mge_frac, mge_pct, mge_per_mb, provirus_per_mb, entero),
            by = c("Genome_ID" = "catalog_id"))
stopifnot(sum(is.na(md$mge_frac)) == 0)

for (v in c("mge_pct", "provirus_per_mb")) {
  w <- wilcox.test(md[[v]] ~ md$enriched_case, exact = FALSE)
  say(sprintf("   %-16s case-enriched median %.3f vs rest %.3f; Wilcoxon p = %.3g", v,
              median(md[[v]][md$enriched_case == 1]), median(md[[v]][md$enriched_case == 0]), w$p.value))
}
say(sprintf("   Enterobacteriaceae mobilome %.2f%% of genome vs %.2f%% elsewhere; Wilcoxon p = %.3g",
            100 * median(md$mge_frac[md$entero]), 100 * median(md$mge_frac[!md$entero]),
            wilcox.test(mge_frac ~ entero, data = md, exact = FALSE)$p.value))

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
  fit_one(md, "Published-style (no mobilome)",
          enriched_case ~ Is_Novel_bin + Completeness + Contamination),
  fit_one(md, "Plus mobilome",
          enriched_case ~ Is_Novel_bin + Completeness + Contamination + mge_pct),
  fit_one(md, "Mobilome, no contamination",
          enriched_case ~ Is_Novel_bin + Completeness + mge_pct))
write_csv(models, file.path(OUTDIR, "mge_phyloglm_models_n76.csv"))
say("")
for (i in seq_len(nrow(models))) with(models[i, ],
  if (term != "(Intercept)")
    say(sprintf("   %-30s %-14s beta %+.4f  p = %.3g", model, term, beta, p)))
b1 <- models %>% filter(model == "Published-style (no mobilome)", term == "Contamination") %>% pull(beta)
b2 <- models %>% filter(model == "Plus mobilome", term == "Contamination") %>% pull(beta)
p2 <- models %>% filter(model == "Plus mobilome", term == "Contamination") %>% pull(p)
say(sprintf("   contamination coefficient %+.4f -> %+.4f on adding the mobilome (%.0f%% change, p = %.3g)",
            b1, b2, 100 * (b2 - b1) / abs(b1), p2))

# =============================================================================
# Q3. community-level mobilome load per donor
# =============================================================================

say(""); say("## Q3. abundance-weighted community mobilome")
donor <- ab_long %>%
  inner_join(d %>% select(MAG = catalog_id, mge_frac, provirus_frac, plasmid_frac,
                          free_virus_frac, provirus_per_mb), by = "MAG") %>%
  group_by(donor) %>%
  summarise(w_mge      = sum(relab * mge_frac),
            w_provirus = sum(relab * provirus_frac),
            w_plasmid  = sum(relab * plasmid_frac),
            w_free     = sum(relab * free_virus_frac),
            w_prov_mb  = sum(relab * provirus_per_mb), .groups = "drop") %>%
  left_join(meta, by = "donor")

comm <- map_dfr(c("w_mge", "w_provirus", "w_plasmid", "w_free", "w_prov_mb"), function(v) {
  a <- donor[[v]][donor$Group == "Case"]; b <- donor[[v]][donor$Group == "Control"]
  tibble(metric = v, case = median(a), control = median(b),
         cliffs = mean(outer(a, b, function(x, y) sign(x - y))),
         p = wilcox.test(a, b, exact = FALSE)$p.value)
}) %>% mutate(q = p.adjust(p, "BH"))
write_csv(comm, file.path(OUTDIR, "mge_community_tests_n76.csv"))
for (i in seq_len(nrow(comm))) with(comm[i, ],
  say(sprintf("   %-11s case %.4f vs control %.4f | Cliff's delta %+.3f | p = %.3g, q = %.3g",
              metric, case, control, cliffs, p, q)))

# THE decisive control. Enterobacteriaceae are mobilome-rich by nature and they
# are what blooms in the cases, so an abundance-weighted community metric will
# rise mechanically whether or not the rest of the community changed at all.
# Recompute with them, and with all of Pseudomonadota, removed.
say("")
say("   confound test: recompute after dropping the blooming clade")
weighted <- function(keep, lab) {
  sub <- ab_long %>% filter(MAG %in% keep) %>%
    group_by(donor) %>% mutate(relab = pct / sum(pct)) %>% ungroup() %>%
    inner_join(d %>% select(MAG = catalog_id, mge_frac, plasmid_frac, provirus_frac,
                            n_amr, genome_mb), by = "MAG") %>%
    group_by(donor) %>%
    summarise(w_mge = sum(relab * mge_frac), w_plasmid = sum(relab * plasmid_frac),
              w_provirus = sum(relab * provirus_frac),
              w_amr = sum(relab * n_amr / genome_mb), .groups = "drop") %>%
    left_join(meta, by = "donor")
  map_dfr(c("w_mge", "w_plasmid", "w_provirus", "w_amr"), function(v) {
    a <- sub[[v]][sub$Group == "Case"]; b <- sub[[v]][sub$Group == "Control"]
    tibble(set = lab, metric = v, case = median(a), control = median(b),
           cliffs = mean(outer(a, b, function(x, y) sign(x - y))),
           p = wilcox.test(a, b, exact = FALSE)$p.value)
  })
}
confound <- bind_rows(
  weighted(d$catalog_id, "All 584 MAGs"),
  weighted(d$catalog_id[d$family != "Enterobacteriaceae"], "Enterobacteriaceae removed"),
  weighted(d$catalog_id[d$phylum != "Pseudomonadota"], "Pseudomonadota removed"))
write_csv(confound, file.path(OUTDIR, "mge_confound_test_n76.csv"))
for (i in seq_len(nrow(confound))) with(confound[i, ],
  say(sprintf("   %-27s %-11s case %.4f vs control %.4f | Cliff's delta %+.3f | p = %.3g",
              set, metric, case, control, cliffs, p)))

# each control cohort separately, because the controls are three studies
say("")
for (s in c("PRJEB7331", "PRJEB7949", "PRJEB39223")) {
  a <- donor$w_mge[donor$Group == "Case"]; b <- donor$w_mge[donor$source == s]
  say(sprintf("   w_mge, cases vs %-11s (n = %2d): Cliff's delta %+.3f, p = %.3g", s, length(b),
              mean(outer(a, b, function(x, y) sign(x - y))),
              wilcox.test(a, b, exact = FALSE)$p.value))
}

# =============================================================================
# figure
# =============================================================================

pA <- d %>% select(catalog_id, provirus_frac, plasmid_frac, free_virus_frac) %>%
  pivot_longer(-catalog_id, names_to = "class", values_to = "frac") %>%
  mutate(class = recode(class, provirus_frac = "Provirus", plasmid_frac = "Plasmid",
                        free_virus_frac = "Free virus")) %>%
  ggplot(aes(100 * frac, class, fill = class)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.8, linewidth = 0.35, show.legend = FALSE) +
  scale_x_sqrt() +
  scale_fill_manual(values = c("Provirus" = "#8073AC", "Plasmid" = "#5AAE61", "Free virus" = "#D55E00")) +
  labs(title = "Mobilome content of the 584 MAGs",
       subtitle = "Square-root axis; each point is one genome",
       x = "Percent of genome", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "MGE_a_classes_n76.pdf"), pA, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pB <- ggplot(d, aes(100 * mge_frac, contamination)) +
  geom_point(aes(colour = entero), size = 1.3, alpha = 0.7) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = c("FALSE" = "grey65", "TRUE" = "#D55E00"),
                      labels = c("Other families", "Enterobacteriaceae"), name = NULL) +
  labs(title = "Mobilome against CheckM2 contamination",
       subtitle = "Does mis-binned mobile DNA account for the contamination signal?",
       x = "Mobilome (% of genome)", y = "CheckM2 contamination (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"), legend.position = "bottom")

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "MGE_b_contamination_n76.pdf"), pB, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pC <- md %>% mutate(status = ifelse(enriched_case == 1, "Case-enriched",
                             ifelse(diff_GroupCase %in% TRUE, "Depleted", "Not significant"))) %>%
  ggplot(aes(status, mge_pct, fill = status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.35, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.15, seed = 42), size = 0.6, alpha = 0.45) +
  scale_y_sqrt() +
  scale_fill_manual(values = c("Case-enriched" = "#F8766D", "Depleted" = "#00BFC4",
                               "Not significant" = "grey80")) +
  labs(title = "Mobilome by differential-abundance status",
       x = NULL, y = "Mobilome (% of genome, sqrt axis)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 7.5))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "MGE_c_status_n76.pdf"), pC, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pD <- models %>% filter(term != "(Intercept)") %>%
  mutate(term = recode(term, Is_Novel_bin = "Novel species", mge_pct = "Mobilome (%)"),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         model = factor(model, levels = unique(models$model))) %>%
  ggplot(aes(beta, model, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.35) +
  facet_wrap(~ term, ncol = 1, scales = "free") +
  scale_colour_manual(values = c("grey25", "#3B6FB6", "#5AAE61"), name = NULL) +
  scale_y_discrete(labels = NULL) +
  labs(title = "Case enrichment with and without the mobilome",
       subtitle = "Phylogenetic logistic regression",
       x = "Coefficient (log odds)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(nrow = 3))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "MGE_d_phyloglm_n76.pdf"), pD, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pE <- confound %>%
  filter(metric %in% c("w_mge", "w_plasmid")) %>%
  mutate(metric = recode(metric, w_mge = "Whole mobilome", w_plasmid = "Plasmid only"),
         set = factor(set, levels = c("All 584 MAGs", "Enterobacteriaceae removed",
                                      "Pseudomonadota removed"))) %>%
  ggplot(aes(set, cliffs, fill = metric)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = ifelse(p < 0.05, sprintf("p=%.2g", p), "ns")),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("Whole mobilome" = "#3B6FB6", "Plasmid only" = "#5AAE61"), name = NULL) +
  scale_y_continuous(limits = c(-0.1, 1), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "The community signal is the bloom, restated",
       subtitle = "Abundance-weighted case-control effect collapses once Enterobacteriaceae are dropped",
       x = NULL, y = "Cliff's delta, case vs control") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 7), legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "MGE_e_community_n76.pdf"), pE, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pF <- d %>% group_by(phylum) %>%
  filter(n() >= 10) %>%
  summarise(med = 100 * median(mge_frac), n = n(), .groups = "drop") %>%
  ggplot(aes(reorder(phylum, med), med)) +
  geom_col(fill = "#7BA7C7", width = 0.7) +
  geom_text(aes(label = paste0("n=", n)), hjust = -0.15, size = 2.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Median mobilome by phylum", subtitle = "Phyla with at least 10 MAGs",
       x = NULL, y = "Mobilome (% of genome)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "MGE_f_phylum_n76.pdf"), pF, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pl <- list(MGE_a_classes = pA, MGE_b_contamination = pB, MGE_c_status = pC,
           MGE_d_phyloglm = pD, MGE_e_community = pE, MGE_f_phylum = pF)
for (nm in names(pl))
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), pl[[nm]], width = 5.2, height = 4.6,
         device = cairo_pdf, bg = "transparent")

tag <- theme(plot.tag = element_text(face = "bold", size = 13))
fig <- ((pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag) | (pC + labs(tag = "c") + tag)) /
       ((pD + labs(tag = "d") + tag) | (pE + labs(tag = "e") + tag) | (pF + labs(tag = "f") + tag)) +
  plot_annotation(title = "Mobilome of the COVID gut MAG catalogue (geNomad) - negative result",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUTDIR, "MGE_Figure_n76.pdf"), fig, width = 16, height = 10,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)

write_tsv(d, file.path(OUTDIR, "mge_per_mag_annotated_n76.tsv"))
write_tsv(donor, file.path(OUTDIR, "mge_per_donor_n76.tsv"))
writeLines(STATS, file.path(OUTDIR, "MGE_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "MGE_Figure_n76.pdf"))
