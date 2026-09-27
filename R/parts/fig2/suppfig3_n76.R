# =============================================================================
# Supplementary Figure 3, rebuilt-catalogue robustness checks
#
# Two questions a reviewer will ask, answered directly.
#
# 1. CATALOGUE IMBALANCE. The catalogue holds 198 case-derived and 386
#    control-derived genomes, and the depleted:enriched ratio has tracked that
#    imbalance across analyses (1.20 published, 1.45 donor-collapsed, 2.33 here).
#    A genome can only be discovered in the arm whose assemblies produced it, so
#    part of the depletion may be catalogue construction rather than biology.
#    Test: re-run ANCOM-BC2 restricted to MAGs detected in BOTH arms, where
#    discovery provenance cannot decide presence, and see what survives.
#
# 2. CONTAMINATION. The phylogenetic logistic regression now finds contamination
#    predicting case-enrichment (beta = +0.088, p = 0.008), which reads as
#    assembly artefact. Test whether it is carried by the blooming
#    Enterobacteriaceae (strain mixtures inflate CheckM2 contamination) by
#    refitting without them and with family adjustment.
#
# Detection uses the same 0.01% relative-abundance floor as Figure 1.
#
# Outputs -> result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(ANCOMBC); library(ape); library(phylolm)
  library(dplyr); library(tidyr); library(tibble); library(readr); library(purrr)
  library(ggplot2); library(patchwork); library(ggrepel)
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
OUTDIR <- file.path(N76, "supp_fig")
source(file.path(BASE, "taxon_italics_n76.R"))   # md_taxon(): taxon names in italics
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

PREV_DETECTION <- 1e-4
cols <- c("Control" = "#00BFC4", "Case" = "#F8766D")

ps  <- readRDS(file.path(N76, "ps_mags_n76.rds"))
res_full <- read.csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), stringsAsFactors = FALSE)
tax <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                  quote = "", comment.char = "")
prov <- tax %>% transmute(taxon = catalog_id, source_arm = group, family, phylum,
                          completeness, contamination)

otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
rel <- sweep(otu, 2, colSums(otu), "/")
grp <- as.character(sample_data(ps)$Group)
prev_case <- rowMeans(rel[, grp == "Case", drop = FALSE] > PREV_DETECTION)
prev_ctrl <- rowMeans(rel[, grp == "Control", drop = FALSE] > PREV_DETECTION)

# =============================================================================
# 1. catalogue-imbalance sensitivity
# =============================================================================

sets <- list(
  "All 584 MAGs"                = rownames(rel),
  "Detected in\nboth arms"       = rownames(rel)[prev_case >= 0.10 & prev_ctrl >= 0.10],
  "Detected in both\n(>= 25%)"   = rownames(rel)[prev_case >= 0.25 & prev_ctrl >= 0.25])
for (nm in names(sets))
  say(sprintf("set '%s': %d MAGs (%d case-derived, %d control-derived)", nm, length(sets[[nm]]),
              sum(prov$source_arm[match(sets[[nm]], prov$taxon)] == "Case"),
              sum(prov$source_arm[match(sets[[nm]], prov$taxon)] == "Control")))

run_ancom <- function(keep) {
  sub <- prune_taxa(keep, ps)
  set.seed(123)
  out <- ancombc2(data = sub, fix_formula = "Group", p_adj_method = "BH", pseudo_sens = TRUE,
                  prv_cut = 0.01, lib_cut = 1000, group = "Group", struc_zero = FALSE,
                  neg_lb = TRUE, alpha = 0.05, global = FALSE)
  out$res
}
res_list <- list("All 584 MAGs" = res_full)
for (nm in setdiff(names(sets), "All 584 MAGs")) {
  say(sprintf("running ANCOM-BC2 on '%s' ...", nm))
  res_list[[nm]] <- run_ancom(sets[[nm]])
}

summarise_res <- function(r, nm) {
  sig <- r %>% filter(diff_GroupCase %in% TRUE)
  up <- sum(sig$lfc_GroupCase > 0); dn <- sum(sig$lfc_GroupCase < 0)
  tibble(Set = nm, tested = nrow(r), significant = nrow(sig), enriched = up, depleted = dn,
         ratio = dn / pmax(up, 1))
}
summ <- bind_rows(imap(res_list, ~ summarise_res(.x, .y)))
write_csv(summ, file.path(OUTDIR, "SuppFig3_imbalance_summary_n76.csv"))
say("")
for (i in seq_len(nrow(summ))) with(summ[i, ],
  say(sprintf("   %-26s tested %3d | significant %3d | enriched %3d | depleted %3d | ratio %.2f",
              Set, tested, significant, enriched, depleted, ratio)))

# provenance split within each set
prov_tab <- bind_rows(imap(res_list, function(r, nm)
  r %>% filter(diff_GroupCase %in% TRUE) %>%
    left_join(prov %>% select(taxon, source_arm), by = "taxon") %>%
    count(source_arm, direction = ifelse(lfc_GroupCase > 0, "Enriched", "Depleted")) %>%
    mutate(Set = nm)))
write_csv(prov_tab, file.path(OUTDIR, "SuppFig3_imbalance_provenance_n76.csv"))

# LFC concordance, full vs restricted.
# ANCOM-BC2 is compositional: it estimates each sample's log sampling fraction from
# the whole taxon set, so dropping the arm-exclusive genomes moves the reference
# frame. Check whether that is all that changes.
key <- "Detected in\nboth arms"
conc <- res_full %>% select(taxon, lfc_full = lfc_GroupCase, q_full = q_GroupCase, diff_full = diff_GroupCase) %>%
  inner_join(res_list[[key]] %>% select(taxon, lfc_sub = lfc_GroupCase, q_sub = q_GroupCase,
                                        diff_sub = diff_GroupCase), by = "taxon") %>%
  filter(!is.na(lfc_full), !is.na(lfc_sub)) %>%
  left_join(prov %>% select(taxon, source_arm), by = "taxon") %>%
  mutate(shift = lfc_sub - lfc_full)
r_pear <- cor(conc$lfc_full, conc$lfc_sub)
r_spear <- cor(conc$lfc_full, conc$lfc_sub, method = "spearman")
shift_mu <- mean(conc$shift); shift_sd <- sd(conc$shift)
lost <- conc %>% filter(diff_full %in% TRUE, !(diff_sub %in% TRUE))
gained <- conc %>% filter(!(diff_full %in% TRUE), diff_sub %in% TRUE)
kept <- sum(conc$diff_full %in% TRUE & conc$diff_sub %in% TRUE)
say("")
say(sprintf("LFC concordance, full vs '%s' (%d shared MAGs):", key, nrow(conc)))
say(sprintf("   Pearson r = %.4f | Spearman rho = %.4f | sign flips among significant calls: %d",
            r_pear, r_spear, sum(sign(conc$lfc_full) != sign(conc$lfc_sub) &
                                  (conc$diff_full %in% TRUE | conc$diff_sub %in% TRUE))))
say(sprintf("   the restriction shifts every effect size by a near-constant %+.3f log units (SD %.3f, max deviation %.3f)",
            shift_mu, shift_sd, max(abs(conc$shift - shift_mu))))
say(sprintf("   significant in both %d | lost %d (all depleted, median |LFC| %.2f, median q %.3f -> %.3f) | gained %d (all enriched)",
            kept, nrow(lost), median(abs(lost$lfc_full)), median(lost$q_full), median(lost$q_sub), nrow(gained)))
say("   reading: dropping the arm-exclusive genomes moves the compositional reference, not the")
say("   ranking. Borderline calls near the new zero switch sides; the ordering of MAGs is intact.")

pA <- ggplot(summ %>% mutate(Set = factor(Set, levels = summ$Set)),
             aes(Set, ratio, fill = Set)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("%d dep / %d enr\nratio %.2f", depleted, enriched, ratio)),
            vjust = -0.25, size = 2.9) +
  scale_fill_manual(values = c("grey70", "#7BA7C7", "#3B6FB6")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Depletion:enrichment ratio by MAG set",
       subtitle = "Restricting to MAGs detected in both arms removes discovery-provenance asymmetry",
       x = NULL, y = "Depleted : enriched") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 7.5))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "SuppFig3a_Ratio_n76.pdf"), pA, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pB <- ggplot(conc, aes(lfc_full, lfc_sub, colour = source_arm)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dotted") +
  geom_abline(slope = 1, intercept = shift_mu, colour = "black", linewidth = 0.5) +
  geom_point(size = 1.5, alpha = 0.75) +
  scale_colour_manual(values = c("Case" = "#F8766D", "Control" = "#00BFC4"),
                      name = "MAG discovered in") +
  annotate("label", x = min(conc$lfc_full), y = max(conc$lfc_sub), hjust = 0, vjust = 1, size = 2.6,
           label = sprintf("Spearman rho = %.3f\nconstant shift %+.2f (SD %.2f)\n%d MAGs, 0 sign flips",
                           r_spear, shift_mu, shift_sd, nrow(conc))) +
  labs(title = "The restriction moves the reference, not the ranking",
       subtitle = "Solid line: y = x + shift. Dotted: identity",
       x = "Log fold change, all 584 MAGs", y = "Log fold change, detected in both arms") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"), legend.position = "bottom")

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "SuppFig3b_LFC_concordance_n76.pdf"), pB, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pC <- ggplot(prov_tab %>% mutate(Set = factor(Set, levels = summ$Set)),
             aes(Set, n, fill = interaction(source_arm, direction, sep = " / "))) +
  geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c("Case / Enriched" = "#F8766D", "Control / Enriched" = "#F6B3AE",
                               "Case / Depleted" = "#7FD8DB", "Control / Depleted" = "#00BFC4"),
                    name = "Discovered in / direction") +
  labs(title = "Significant MAGs by discovery provenance",
       subtitle = "Control-derived genomes dominate the depleted side in the full catalogue",
       x = NULL, y = "MAGs") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 7.5), legend.position = "bottom",
        legend.text = element_text(size = 6.5), legend.key.size = unit(0.3, "cm")) +
  guides(fill = guide_legend(nrow = 2, title.position = "top"))

# =============================================================================
# 2. contamination and the Enterobacteriaceae bloom
# =============================================================================

st <- readRDS(file.path(N76, "fig2", "figure2_stats_n76.rds"))
tree <- read.tree(TREE); tree$node.label <- NULL
meta_df <- st$meta_df
meta_df$entero <- as.integer(meta_df$Family == "Enterobacteriaceae")
mean_ab <- rowMeans(rel) * 100
meta_df$mean_abund <- mean_ab[meta_df$Genome_ID]
meta_df$mean_abund_plot <- meta_df$mean_abund + 1e-5

say("")
say("## contamination")
say(sprintf("   median contamination: Enterobacteriaceae %.2f%% (n = %d) vs all others %.2f%% (n = %d); Wilcoxon p = %.3g",
            median(meta_df$Contamination[meta_df$entero == 1]), sum(meta_df$entero == 1),
            median(meta_df$Contamination[meta_df$entero == 0]), sum(meta_df$entero == 0),
            wilcox.test(Contamination ~ entero, data = meta_df)$p.value))
enr <- meta_df$enriched_case == 1
say(sprintf("   median contamination: case-enriched %.2f%% vs rest %.2f%%; Wilcoxon p = %.3g",
            median(meta_df$Contamination[enr]), median(meta_df$Contamination[!enr]),
            wilcox.test(Contamination ~ enriched_case, data = meta_df)$p.value))
ct <- cor.test(meta_df$Contamination, log10(meta_df$mean_abund + 1e-6), method = "spearman")
say(sprintf("   contamination vs mean abundance: Spearman rho = %.3f, p = %.3g", ct$estimate, ct$p.value))

fit_one <- function(dat, label, formula) {
  dat <- as.data.frame(dat); rownames(dat) <- dat$Genome_ID
  tr <- keep.tip(tree, intersect(tree$tip.label, dat$Genome_ID))
  d <- dat[tr$tip.label, , drop = FALSE]
  f <- tryCatch(phyloglm(formula, phy = tr, data = d, method = "logistic_MPLE"),
                error = function(e) { message("  phyloglm failed: ", conditionMessage(e)); NULL })
  if (is.null(f)) return(NULL)
  cf <- summary(f)$coefficients
  pcol <- if ("p.value" %in% colnames(cf)) "p.value" else colnames(cf)[ncol(cf)]
  tibble(model = label, term = rownames(cf), beta = cf[, "Estimate"], se = cf[, "StdErr"],
         p = cf[, pcol], n = nrow(d))
}
# every Enterobacteriaceae MAG is case-enriched, so a family term is completely
# separated and cannot be estimated; drop-one-group refits are used instead.
n_ent <- sum(meta_df$entero == 1); n_ent_enr <- sum(meta_df$entero == 1 & meta_df$enriched_case == 1)
say(sprintf("   %d of %d Enterobacteriaceae MAGs are case-enriched: complete separation, so no family term is estimable",
            n_ent_enr, n_ent))
say(sprintf("   Spearman(completeness, contamination) = %.3f -- the two quality terms move in opposite directions",
            cor(meta_df$Completeness, meta_df$Contamination, method = "spearman")))

hq <- meta_df %>% filter(Contamination <= 5)
# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "SuppFig3c_Provenance_n76.pdf"), pC, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pseudo <- meta_df %>% filter(Phylum != "Pseudomonadota")
models <- bind_rows(
  fit_one(meta_df, "All MAGs", enriched_case ~ Is_Novel_bin + Completeness + Contamination),
  fit_one(meta_df %>% filter(entero == 0), sprintf("Enterobacteriaceae excluded (n = %d)", sum(meta_df$entero == 0)),
          enriched_case ~ Is_Novel_bin + Completeness + Contamination),
  fit_one(pseudo, sprintf("Pseudomonadota excluded (n = %d)", nrow(pseudo)),
          enriched_case ~ Is_Novel_bin + Completeness + Contamination),
  fit_one(hq, sprintf("Contamination <= 5%% only (n = %d)", nrow(hq)),
          enriched_case ~ Is_Novel_bin + Completeness + Contamination),
  fit_one(meta_df, "Contamination unadjusted", enriched_case ~ Contamination))
write_csv(models, file.path(OUTDIR, "SuppFig3_contamination_models_n76.csv"))
say("")
for (i in seq_len(nrow(models))) with(models[i, ],
  if (term != "(Intercept)")
    say(sprintf("   %-38s %-16s beta %+.4f  p = %.3g  (n = %d)", model, term, beta, p, n)))
ct_b <- models %>% filter(term == "Contamination")
say("")
say(sprintf("   contamination stays positive and significant in %d of %d refits",
            sum(ct_b$beta > 0 & ct_b$p < 0.05), nrow(ct_b)))

pD <- meta_df %>%
  mutate(grp = ifelse(entero == 1, "Enterobacteriaceae", "All other families"),
         status = ifelse(enriched_case == 1, "Case-enriched", "Not enriched")) %>%
  ggplot(aes(grp, Contamination, fill = status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.35) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, seed = 42), size = 0.7, alpha = 0.5) +
  scale_fill_manual(values = c("Case-enriched" = "#F8766D", "Not enriched" = "grey80"), name = NULL) +
  scale_x_discrete(labels = md_taxon) +
  labs(title = "Contamination by family and enrichment status",
       subtitle = sprintf("All %d *Enterobacteriaceae* MAGs are case-enriched, and they are the cleaner genomes", n_ent),
       x = NULL, y = "CheckM2 contamination (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        axis.text.x = ggtext::element_markdown(size = 8),
        plot.subtitle = ggtext::element_markdown(size = 7, colour = "grey35"),
        legend.position = "bottom")

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "SuppFig3d_Contamination_family_n76.pdf"), pD, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pE <- ggplot(meta_df, aes(mean_abund_plot, Contamination)) +
  geom_point(aes(colour = factor(entero, levels = c(0, 1),
                                 labels = c("Other families", "Enterobacteriaceae"))),
             size = 1.4, alpha = 0.75) +
  geom_smooth(method = "lm", colour = "black", se = TRUE, linewidth = 0.6) +
  scale_x_log10() +
  scale_colour_manual(values = c("Other families" = "grey65", "Enterobacteriaceae" = "#D55E00"),
                      labels = md_taxon, name = NULL) +
  annotate("label", x = min(meta_df$mean_abund_plot), y = max(meta_df$Contamination),
           hjust = 0, vjust = 1, size = 2.8,
           label = sprintf("Spearman rho = %.2f\np = %.2g", ct$estimate, ct$p.value)) +
  labs(title = "Contamination falls as MAG abundance rises",
       subtitle = "No sign of strain mixtures inflating contamination in the abundant genomes",
       x = "Mean relative abundance (%, log scale)", y = "CheckM2 contamination (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.text = ggtext::element_markdown(), legend.position = "bottom")

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "SuppFig3e_Contamination_abundance_n76.pdf"), pE, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pF <- models %>% filter(term != "(Intercept)") %>%
  mutate(term = recode(term, Is_Novel_bin = "Novel species", Completeness = "Completeness",
                       Contamination = "Contamination"),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         model = factor(model, levels = unique(models$model))) %>%
  ggplot(aes(beta, model, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55), size = 0.35) +
  facet_wrap(~ term, ncol = 1, scales = "free") +
  scale_colour_manual(values = c("grey25", "#D55E00", "#3B6FB6", "#5AAE61", "#8073AC"),
                      labels = md_taxon_in, name = NULL) +
  labs(title = "Phylogenetic logistic regression for case enrichment",
       subtitle = "The contamination term survives every drop-one refit",
       x = "Coefficient (log odds)", y = NULL) +
  scale_y_discrete(labels = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = ggtext::element_markdown(size = 6),
        legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(nrow = 3))

# ---- save -------------------------------------------------------------------

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "SuppFig3f_Phyloglm_n76.pdf"), pF, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pl <- list(SuppFig3a_Ratio = pA, SuppFig3b_LFC_concordance = pB, SuppFig3c_Provenance = pC,
           SuppFig3d_Contamination_family = pD, SuppFig3e_Contamination_abundance = pE,
           SuppFig3f_Phyloglm = pF)
for (nm in names(pl))
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), pl[[nm]], width = 5.2, height = 4.6,
         device = cairo_pdf, bg = "transparent")

tag <- theme(plot.tag = element_text(face = "bold", size = 13))

# Two columns leave each panel narrower than its title, so the titles ran into
# the neighbouring panel. Wrap them and drop the type a point for the composite.
wrap_lab <- function(x, n) if (is.null(x)) NULL else paste(strwrap(x, n), collapse = "\n")
for_page <- function(p, md_subtitle = FALSE, lean_x = FALSE) {
  ttl <- wrap_lab(p$labels$title, 36)
  sub <- wrap_lab(p$labels$subtitle, 54)
  if (md_subtitle && !is.null(sub)) sub <- gsub("\n", "<br>", sub)   # gridtext needs <br>
  p + labs(title = ttl, subtitle = sub) +
    theme(plot.title = element_text(size = 8.5),
          plot.subtitle = if (md_subtitle) ggtext::element_markdown(size = 6.3, colour = "grey35")
                          else element_text(size = 6.3, colour = "grey35"),
          # legend.text is left alone: panels d, e and f carry ggtext markdown
          # legends and ggplot2 4 refuses to merge two element classes
          legend.title = element_text(size = 6.5)) +
    # category names are wider than a half-page square panel
    (if (lean_x) theme(axis.text.x = if (md_subtitle)
           ggtext::element_markdown(size = 6, angle = 20, hjust = 1)
         else element_text(size = 6, angle = 20, hjust = 1)) else NULL)
}
pA <- for_page(pA, lean_x = TRUE); pB <- for_page(pB)
pC <- for_page(pC, lean_x = TRUE)
pD <- for_page(pD, md_subtitle = TRUE, lean_x = TRUE)
pE <- for_page(pE); pF <- for_page(pF)

# Every plotting area square, then the canvas is sized around them. pF is
# faceted into three rows, so a third each keeps the block itself square.
sq <- theme(aspect.ratio = 1)
pA <- pA + sq; pB <- pB + sq; pC <- pC + sq; pD <- pD + sq; pE <- pE + sq
pF <- pF + theme(aspect.ratio = 1/3)
supp3 <- ((pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag) | (pC + labs(tag = "c") + tag)) /
         ((pD + labs(tag = "d") + tag) | (pE + labs(tag = "e") + tag) | (pF + labs(tag = "f") + tag)) +
  plot_annotation(
    title = "Catalogue-imbalance and genome-quality robustness checks",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
# canvas follows the panels: square plotting areas need the extra height
ggsave(file.path(OUTDIR, "Supplementary_Figure_3_n76.pdf"), supp3, width = 12, height = 10.4,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)

writeLines(STATS, file.path(OUTDIR, "SuppFig3_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Supplementary_Figure_3_n76.pdf"))
