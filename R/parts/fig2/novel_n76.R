# =============================================================================
# Novel taxa in the 76-donor catalogue
#
# Evidence for the candidate novel-species paragraph:
#
#   1. how many MAGs are novel, and at what rank
#   2. whether they are good enough genomes to name (MIMAG)
#   3. where they sit taxonomically and which arm discovered them
#   4. whether they are differentially abundant, and in which direction
#   5. whether they are phylogenetically clustered
#   6. whether they are rare or genuinely prevalent
#
# Outputs -> result2/n76/novel/
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(phyloseq)
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
OUTDIR <- file.path(N76, "novel")
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }
rd <- function(f) read.delim(file.path(INDIR, f), check.names = FALSE, quote = "", comment.char = "")

PREV_DETECTION <- 1e-4   # same floor as Figure 1

tax  <- rd("MAG_quality_taxonomy.tsv")
res  <- read.csv(file.path(N76, "ANCOMBC2_MAG_results_n76.csv"), stringsAsFactors = FALSE)
ps   <- readRDS(file.path(N76, "ps_mags_n76.rds"))
tree <- read.tree(TREE); tree$node.label <- NULL

blank <- function(x) is.na(x) | trimws(x) == ""
d <- tax %>%
  transmute(catalog_id, source_arm = group, phylum, class, order, family, genus, species,
            completeness, contamination,
            n_contigs = checkm2_Total_Contigs, n50 = checkm2_Contig_N50,
            genome_size = checkm2_Genome_Size, gc = checkm2_GC_Content,
            coding_density = checkm2_Coding_Density) %>%
  mutate(novel_species = blank(species),
         novel_genus   = blank(genus),
         novel_family  = blank(family),
         rank_novel = case_when(novel_family ~ "New family",
                                novel_genus  ~ "New genus",
                                novel_species ~ "New species",
                                TRUE ~ "Known species"),
         # MIMAG high quality needs >90% complete, <5% contaminated. rRNA and
         # tRNA counts are not in this table, so this is the genome-quality part
         # of the standard only, and is labelled as such throughout.
         hq = completeness > 90 & contamination < 5)

otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
rel <- sweep(otu, 2, colSums(otu), "/")
grp <- as.character(sample_data(ps)$Group)
prev <- tibble(catalog_id = rownames(rel),
               prev_case = rowMeans(rel[, grp == "Case",    drop = FALSE] > PREV_DETECTION),
               prev_ctrl = rowMeans(rel[, grp == "Control", drop = FALSE] > PREV_DETECTION),
               mean_abund = 100 * rowMeans(rel))

d <- d %>%
  left_join(prev, by = "catalog_id") %>%
  left_join(res %>% select(catalog_id = taxon, lfc = lfc_GroupCase,
                           q = q_GroupCase, diff = diff_GroupCase), by = "catalog_id") %>%
  mutate(status = case_when(diff %in% TRUE & lfc > 0 ~ "Enriched in case",
                            diff %in% TRUE & lfc < 0 ~ "Depleted in case",
                            TRUE ~ "Not significant"))

# ---- 1-3. how much novelty, of what quality, from where ---------------------

say(sprintf("catalogue: %d MAGs", nrow(d)))
say("")
say("## 1. novelty by rank")
for (r in c("New family", "New genus", "New species", "Known species"))
  say(sprintf("   %-14s %3d MAGs", r, sum(d$rank_novel == r)))
nov <- d %>% filter(novel_species)
say(sprintf("   any novelty at species level or above: %d of %d (%.1f%%)",
            nrow(nov), nrow(d), 100 * nrow(nov) / nrow(d)))

say("")
say("## 2. are they good enough to name")
say(sprintf("   novel MAGs meeting >90%% complete and <5%% contaminated: %d of %d",
            sum(nov$hq), nrow(nov)))
say(sprintf("   median completeness %.1f%% (range %.1f-%.1f), median contamination %.2f%%",
            median(nov$completeness), min(nov$completeness), max(nov$completeness),
            median(nov$contamination)))
# the CheckM2 assembly columns exist for the 386 control-derived genomes; the
# 198 case-derived MAGs have completeness and contamination but no contig
# statistics
n_stat <- sum(!is.na(nov$n50))
say(sprintf("   contig statistics available for %d of %d novel MAGs (the carried-over case genomes have none):",
            n_stat, nrow(nov)))
say(sprintf("      median contigs %.0f, median N50 %s bp",
            median(nov$n_contigs, na.rm = TRUE),
            format(median(nov$n50, na.rm = TRUE), big.mark = ",")))
say("   note: rRNA and tRNA counts are not in this table, so full MIMAG high-quality")
say("   status cannot be asserted from these numbers alone.")

say("")
say("## 3. placement and discovery arm")
for (i in seq_len(nrow(count(nov, phylum, sort = TRUE))))
  NULL
ph <- nov %>% count(phylum, sort = TRUE)
say(sprintf("   phyla: %s", paste(sprintf("%s (%d)", ph$phylum, ph$n), collapse = ", ")))
fa <- nov %>% count(family, sort = TRUE) %>% filter(n > 1)
if (nrow(fa))
  say(sprintf("   families with more than one novel MAG: %s",
              paste(sprintf("%s (%d)", fa$family, fa$n), collapse = ", ")))
arm <- nov %>% count(source_arm)
say(sprintf("   discovered in: %s", paste(sprintf("%s %d", arm$source_arm, arm$n), collapse = ", ")))
tst <- fisher.test(matrix(c(sum(nov$source_arm == "Case"), sum(nov$source_arm == "Control"),
                            sum(d$source_arm == "Case") - sum(nov$source_arm == "Case"),
                            sum(d$source_arm == "Control") - sum(nov$source_arm == "Control")),
                          nrow = 2))
say(sprintf("   novelty rate: case-derived %.1f%% vs control-derived %.1f%%; Fisher p = %.3g",
            100 * sum(nov$source_arm == "Case") / sum(d$source_arm == "Case"),
            100 * sum(nov$source_arm == "Control") / sum(d$source_arm == "Control"), tst$p.value))

# ---- 4. differential abundance ----------------------------------------------

say("")
say("## 4. differential abundance of the novel MAGs")
st <- nov %>% count(status)
for (i in seq_len(nrow(st))) say(sprintf("   %-18s %d", st$status[i], st$n[i]))
sig <- nov %>% filter(diff %in% TRUE) %>% arrange(lfc)
if (nrow(sig)) {
  say("   significant novel MAGs:")
  for (i in seq_len(nrow(sig))) with(sig[i, ],
    say(sprintf("      %-28s %-22s LFC %+.2f  q = %.3g",
                substr(paste(genus, family), 1, 28), substr(catalog_id, 1, 22), lfc, q)))
}
# is novelty itself associated with direction, against the rest of the catalogue?
ct <- table(d$novel_species, d$status)
say(sprintf("   novel vs known, enrichment direction: Fisher p = %.3g",
            fisher.test(ct[, c("Enriched in case", "Depleted in case")])$p.value))

# ---- 5. are the novel MAGs phylogenetically clustered? ----------------------

say("")
say("## 5. phylogenetic clustering of novel MAGs")
tips <- intersect(tree$tip.label, nov$catalog_id)
say(sprintf("   %d of %d novel MAGs are on the tree", length(tips), nrow(nov)))
cd <- cophenetic(tree)

# picante's ses.mpd expects a multi-row community table; with a single set it is
# clearer to draw the null directly. 999 random tip sets of the same size.
mpd_of  <- function(tp) mean(cd[tp, tp][lower.tri(diag(length(tp)))])
mntd_of <- function(tp) { m <- cd[tp, tp]; diag(m) <- Inf; mean(apply(m, 1, min)) }
set.seed(123)
null <- replicate(999, {
  s_ <- sample(tree$tip.label, length(tips))
  c(mpd = mpd_of(s_), mntd = mntd_of(s_))
})
zp <- function(obs, nullv, lab) {
  z <- (obs - mean(nullv)) / sd(nullv)
  p <- (sum(nullv <= obs) + 1) / (length(nullv) + 1)
  say(sprintf("   %-9s observed %.3f vs null %.3f; z = %+.3f, p = %.3f  (z < 0 means clustered)",
              lab, obs, mean(nullv), z, p))
}
zp(mpd_of(tips),  null["mpd", ],  "mean pairwise")
zp(mntd_of(tips), null["mntd", ], "nearest-taxon")

say(sprintf("   %d novel MAGs are Lachnospiraceae",
            sum(nov$family == "Lachnospiraceae", na.rm = TRUE)))

# ---- 6. prevalence ----------------------------------------------------------

say("")
say("## 6. prevalence")
say(sprintf("   novel MAGs: median prevalence %.1f%% of controls, %.1f%% of cases",
            100 * median(nov$prev_ctrl), 100 * median(nov$prev_case)))
say(sprintf("   known MAGs: median prevalence %.1f%% of controls, %.1f%% of cases",
            100 * median(d$prev_ctrl[!d$novel_species]), 100 * median(d$prev_case[!d$novel_species])))
say(sprintf("   novel MAGs present in at least a quarter of donors on either side: %d of %d",
            sum(nov$prev_ctrl >= 0.25 | nov$prev_case >= 0.25), nrow(nov)))
say(sprintf("   median abundance: novel %.4f%% vs known %.4f%%; Wilcoxon p = %.3g",
            median(nov$mean_abund), median(d$mean_abund[!d$novel_species]),
            wilcox.test(mean_abund ~ novel_species, data = d, exact = FALSE)$p.value))

# A control-derived genome is guaranteed present in at least one control, so
# comparing novel to known across the whole catalogue confounds novelty with
# discovery arm. Repeat the comparison inside the control-derived genomes only.
say("")
say("   fair comparison, control-derived genomes only (removes the discovery guarantee):")
cd_only <- d %>% filter(source_arm == "Control")
for (v in c("prev_ctrl", "prev_case")) {
  a <- cd_only[[v]][cd_only$novel_species]; b <- cd_only[[v]][!cd_only$novel_species]
  say(sprintf("      %-10s novel %.1f%% vs known %.1f%%; Wilcoxon p = %.3g", v,
              100 * median(a), 100 * median(b), wilcox.test(a, b, exact = FALSE)$p.value))
}
drop_nov <- median(cd_only$prev_ctrl[cd_only$novel_species] - cd_only$prev_case[cd_only$novel_species])
drop_kno <- median(cd_only$prev_ctrl[!cd_only$novel_species] - cd_only$prev_case[!cd_only$novel_species])
say(sprintf("      prevalence drop control->case: novel %.1f points vs known %.1f points; Wilcoxon p = %.3g",
            100 * drop_nov, 100 * drop_kno,
            wilcox.test((cd_only$prev_ctrl - cd_only$prev_case) ~ cd_only$novel_species,
                        exact = FALSE)$p.value))

write_csv(nov %>% arrange(phylum, family, catalog_id), file.path(OUTDIR, "novel_MAGs_n76.csv"))
write_csv(d, file.path(OUTDIR, "all_MAGs_novelty_n76.csv"))

# ---- figure -----------------------------------------------------------------

sc <- c("Enriched in case" = "#F8766D", "Depleted in case" = "#00BFC4",
        "Not significant" = "grey75")

pA <- d %>% count(rank_novel) %>%
  mutate(rank_novel = factor(rank_novel,
           levels = c("Known species", "New species", "New genus", "New family"))) %>%
  ggplot(aes(rank_novel, n, fill = rank_novel)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("grey75", "#7BA7C7", "#3B6FB6", "#20419A")) +
  scale_y_sqrt(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Novelty in the 584-MAG catalogue",
       subtitle = "GTDB-Tk assignment; square-root axis", x = NULL, y = "MAGs") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 7.5))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "Novel_a_rank_n76.pdf"), pA, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pB <- ggplot(nov, aes(completeness, contamination, colour = status)) +
  geom_rect(aes(xmin = 90, xmax = 100, ymin = -Inf, ymax = 5),
            fill = "grey92", colour = NA, inherit.aes = FALSE) +
  geom_point(size = 2, alpha = 0.85) +
  geom_text_repel(aes(label = ifelse(diff %in% TRUE, genus, "")), size = 2.3,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = sc, name = NULL) +
  labs(title = "Genome quality of the novel MAGs",
       subtitle = "Shaded: >90% complete and <5% contaminated",
       x = "Completeness (%)", y = "Contamination (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "Novel_b_quality_n76.pdf"), pB, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pC <- nov %>% count(phylum, status) %>%
  ggplot(aes(reorder(phylum, n, sum), n, fill = status)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = sc, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Novel MAGs by phylum and direction", x = NULL, y = "MAGs") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "Novel_c_phylum_n76.pdf"), pC, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pD <- ggplot(d, aes(100 * prev_ctrl, 100 * prev_case)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
  geom_point(data = filter(d, !novel_species), colour = "grey85", size = 1) +
  geom_point(data = nov, aes(colour = status), size = 2) +
  geom_text_repel(data = filter(nov, prev_ctrl >= 0.25 | prev_case >= 0.25),
                  aes(label = genus), size = 2.2, max.overlaps = 15) +
  scale_colour_manual(values = sc, name = NULL) +
  labs(title = "Prevalence of the novel MAGs",
       subtitle = sprintf("Detection floor %.2f%% relative abundance; grey = known species",
                          100 * PREV_DETECTION),
       x = "Prevalence in controls (%)", y = "Prevalence in cases (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

# inline save, so each panel can be inspected where it is built
ggsave(file.path(OUTDIR, "Novel_d_prevalence_n76.pdf"), pD, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")

pl <- list(Novel_a_rank = pA, Novel_b_quality = pB, Novel_c_phylum = pC, Novel_d_prevalence = pD)
for (nm in names(pl))
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), pl[[nm]], width = 5.2, height = 4.6,
         device = cairo_pdf, bg = "transparent")

tag <- theme(plot.tag = element_text(face = "bold", size = 13))
fig <- ((pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag)) /
       ((pC + labs(tag = "c") + tag) | (pD + labs(tag = "d") + tag)) +
  plot_annotation(title = sprintf("Novel taxa in the 76-donor catalogue (%d MAGs without a GTDB species)",
                                  nrow(nov)),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUTDIR, "Novel_Figure_n76.pdf"), fig, width = 11, height = 9,
       device = cairo_pdf, bg = "transparent")

writeLines(STATS, file.path(OUTDIR, "NOVEL_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Novel_Figure_n76.pdf"))
