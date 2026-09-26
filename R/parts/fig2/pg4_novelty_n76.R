# =============================================================================
# Protein-level novelty against proGenomes4
#   DIAMOND blastp, 1,365,648 prodigal-gv proteins from the 584-MAG catalogue
#   vs 139,409,217 proGenomes4 representative proteins
#   --sensitive, e-value 1e-5, --max-target-seqs 1 (existence of a homolog only)
#
# The point of running the whole catalogue rather than the 21 no-GTDB-species
# genomes is the baseline: "X% of proteins have no homolog" means nothing
# without knowing what a known species scores.
#
# The reading changed once GUNC came back. 13 of the 21 novel MAGs are chimeric,
# and a chimera could carry orphan proteins for a boring reason, so the only
# interpretable comparison is within the GUNC-clean set.
#
# RESULT: null, in both directions. Inside the GUNC-clean set, genomes without a
# GTDB species score 4.16% orphan against 3.63% for known species (p = 0.67).
# Chimerism does not raise it either (3.22% vs 3.63%, p = 0.58; rho = +0.05
# against the clade-separation score). The median hit identity is 98.8%, so the
# catalogue is already well represented in proGenomes4 at near-identity. These
# genomes are unplaced in GTDB, not novel in gene content, and no novel-biology
# claim can rest on them.
#
# Outputs -> result2/n76/supp_fig/gunc/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
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
OUTDIR <- file.path(BASE, "result2/n76/supp_fig/gunc")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

pg <- read_tsv(file.path(INDIR, "pg4_protein_novelty_per_mag.tsv"), show_col_types = FALSE)
gu <- read_csv(file.path(OUTDIR, "gunc_584_annotated_n76.csv"), show_col_types = FALSE)

d <- gu %>%
  select(catalog_id, arm, family, phylum, species, completeness, contamination,
         novel, pass.GUNC, clade_separation_score) %>%
  inner_join(pg, by = "catalog_id") %>%
  mutate(orphan_pct = 100 * orphan_fraction,
         grp = case_when(!novel ~ "Known species",
                         novel & pass.GUNC ~ "No GTDB species, GUNC-clean",
                         TRUE ~ "No GTDB species, chimeric"))
stopifnot(nrow(d) == 584)

say("proGenomes4 protein search: 1,365,648 proteins, 1,290,546 with a homolog")
say(sprintf("   catalogue orphan fraction: median %.2f%%, IQR %.2f-%.2f%%, max %.1f%%",
            100 * median(d$orphan_fraction),
            100 * quantile(d$orphan_fraction, .25), 100 * quantile(d$orphan_fraction, .75),
            100 * max(d$orphan_fraction)))
say(sprintf("   median identity of the hits that do land: %.1f%%",
            median(d$median_pident_of_hits, na.rm = TRUE)))

say(""); say("## orphan fraction by novelty and chimerism")
tab <- d %>% group_by(grp) %>%
  summarise(n = n(), med = 100 * median(orphan_fraction),
            iqr_lo = 100 * quantile(orphan_fraction, .25),
            iqr_hi = 100 * quantile(orphan_fraction, .75), .groups = "drop") %>%
  arrange(desc(med))
for (i in seq_len(nrow(tab))) with(tab[i, ],
  say(sprintf("   %-30s n = %3d | median %.2f%% (IQR %.2f-%.2f)", grp, n, med, iqr_lo, iqr_hi)))

# the decisive comparison: inside the GUNC-clean set only
clean <- d %>% filter(pass.GUNC)
w <- wilcox.test(orphan_fraction ~ novel, data = clean, exact = FALSE)
say("")
say(sprintf("   GUNC-clean genomes only (n = %d): novel %.2f%% vs known %.2f%%; Wilcoxon p = %.3g",
            nrow(clean),
            100 * median(clean$orphan_fraction[clean$novel]),
            100 * median(clean$orphan_fraction[!clean$novel]), w$p.value))
say("   this is the only comparison where a raised orphan fraction can mean")
say("   uncharacterised biology rather than a bin-merging artefact.")

# and confirm chimerism itself inflates it, which is why the split is needed
wc <- wilcox.test(orphan_fraction ~ pass.GUNC, data = d, exact = FALSE)
say(sprintf("   chimeric vs clean across the whole catalogue: %.2f%% vs %.2f%%; Wilcoxon p = %.3g",
            100 * median(d$orphan_fraction[!d$pass.GUNC]),
            100 * median(d$orphan_fraction[d$pass.GUNC]), wc$p.value))
r <- suppressWarnings(cor.test(d$orphan_fraction, d$clade_separation_score, method = "spearman"))
say(sprintf("   orphan fraction vs clade-separation score: Spearman rho = %+.3f, p = %.3g",
            r$estimate, r$p.value))

# the five genomes that survive both filters
say(""); say("## the five defensible candidate taxa")
five <- d %>% filter(novel, pass.GUNC, completeness > 90, contamination < 5) %>%
  arrange(desc(orphan_fraction))
for (i in seq_len(nrow(five))) with(five[i, ],
  say(sprintf("   %-22s %-16s %5d proteins | orphan %.2f%% | completeness %.1f%%",
              family, phylum, n_proteins, orphan_pct, completeness)))
if (nrow(five))
  say(sprintf("   their median orphan fraction %.2f%% against %.2f%% for known species",
              100 * median(five$orphan_fraction), 100 * median(d$orphan_fraction[!d$novel])))

write_csv(d, file.path(OUTDIR, "pg4_novelty_annotated_n76.csv"))

pA <- ggplot(d, aes(reorder(grp, orphan_pct, median), orphan_pct, fill = grp)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85, linewidth = 0.35, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c("Known species" = "grey75",
                               "No GTDB species, GUNC-clean" = "#3B6FB6",
                               "No GTDB species, chimeric" = "#D55E00")) +
  labs(title = "Orphan protein fraction separates nothing",
       subtitle = "No GTDB species vs known, GUNC-clean only: p = 0.67. Chimeric vs clean: p = 0.58",
       x = NULL, y = "Orphan proteins (% of genome's proteins)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.y = element_text(size = 7))

pB <- ggplot(d, aes(clade_separation_score, orphan_pct)) +
  geom_point(aes(colour = grp), size = 1.3, alpha = 0.75) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = c("Known species" = "grey70",
                                 "No GTDB species, GUNC-clean" = "#3B6FB6",
                                 "No GTDB species, chimeric" = "#D55E00"), name = NULL) +
  labs(title = "Chimerism does not raise the orphan fraction either",
       subtitle = sprintf("Spearman rho = %+.3f, p = %.2f -- flat", r$estimate, r$p.value),
       x = "GUNC clade separation score", y = "Orphan proteins (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(OUTDIR, "PG4_a_orphan_by_group_n76.pdf"), pA, width = 5.6, height = 4.2,
       device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "PG4_b_orphan_vs_css_n76.pdf"), pB, width = 5.2, height = 4.6,
       device = cairo_pdf, bg = "transparent")
tag <- theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(OUTDIR, "PG4_Figure_n76.pdf"),
       (pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag),
       width = 11, height = 4.6, device = cairo_pdf, bg = "transparent")

writeLines(STATS, file.path(OUTDIR, "PG4_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "PG4_Figure_n76.pdf"))
