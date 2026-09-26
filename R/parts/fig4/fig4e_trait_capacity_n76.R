#!/usr/bin/env Rscript
# Figure 4e: abundance-weighted gut functional trait capacity per donor.
#
# The UpSet this replaces counted genomes equally, so a trait carried by 527 of
# 584 MAGs dominated the display while carrying no signal, and the butyrate
# terminal pair, carried by 17 genomes but by 2.6% of a healthy community and
# 0.1% of a patient one, was invisible. Presence of the machinery in a genome
# says nothing about whether that genome is abundant in a donor.
#
# Each donor's value is the share of its detected community carried by genomes
# positive for the trait. The denominator is restricted to genomes whose trait
# status is known, which matters only for the model-derived acid requirement:
# it is defined for the 294 modelled MAGs and unknown for the rest.
#
# This is encoded capacity weighted by abundance, not expression or flux.
#
# Outputs -> result2/n76/fig4/Fig4e_TraitCapacity_n76.{pdf,png,rds}
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "fig4")
TRAITS <- file.path(N76, "metabolic_support/gut_traits_v1")
STATS  <- character(); say <- function(s) { message(s); STATS <<- c(STATS, s) }

LAB <- c(butyrate_terminal_pair  = "Butyrate kinase-route\nterminal pair (ptb/buk)",
         acetate_terminal_pair   = "Acetate terminal\npair (pta/ackA)",
         mucin_gh_repertoire     = "Putative mucin-GH\nmarker combination",
         plant_backbone_cazyme   = "Plant-associated\nCAZyme-family hits",
         strict_acid_uptake      = "Model-derived acid\nrequirement/uptake")
GCOL <- c(Control = "#00BFC4", Case = "#F8766D")

wide <- read_tsv(file.path(TRAITS, "genome_traits_wide.tsv"), show_col_types = FALSE)
ab   <- read.delim(file.path(INDIR, "MAG_relative_abundance_percent.tsv"),
                   check.names = FALSE, row.names = 1)
meta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE)
stopifnot(setequal(wide$catalog_id, rownames(ab)), nrow(meta) == 76)
ab <- as.matrix(ab[wide$catalog_id, meta$sample])
grp <- factor(meta$group, levels = c("Control", "Case"))

cliff <- function(a, b) { m <- outer(a, b, "-"); (sum(m > 0) - sum(m < 0)) / length(m) }

per_donor <- lapply(names(LAB), function(tr) {
  state <- wide[[tr]]
  known <- state != "unknown"
  pos   <- known & state == "positive"
  share <- 100 * colSums(ab[pos, , drop = FALSE]) / colSums(ab[known, , drop = FALSE])
  tibble(trait = tr, sample = meta$sample, group = grp, share = share,
         n_genomes = sum(pos), n_known = sum(known))
}) %>% bind_rows()

tests <- per_donor %>% group_by(trait, n_genomes, n_known) %>%
  summarise(case_median = median(share[group == "Case"]),
            control_median = median(share[group == "Control"]),
            delta = cliff(share[group == "Case"], share[group == "Control"]),
            p = wilcox.test(share ~ group, exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(q = p.adjust(p, "BH")) %>% arrange(delta)

say("Abundance-weighted trait capacity, share of the detected community per donor")
for (i in seq_len(nrow(tests))) with(tests[i, ],
  say(sprintf("  %-24s %3d/%3d genomes | case %6.2f%% vs control %6.2f%% | delta %+.2f | q = %.3g",
              trait, n_genomes, n_known, case_median, control_median, delta, q)))

ord <- tests$trait
per_donor$trait <- factor(per_donor$trait, ord)
tests$trait <- factor(tests$trait, ord)
# Precomputed, not called inside aes(): the saved plot object is reloaded by
# fig4_square_panels_n76.R, where a helper from this script would not exist.
tests <- tests %>% mutate(
  label = sprintf("%s, δ = %+.2f",
                  ifelse(q < 1e-4, "q < 1e-4", sprintf("q = %.3g", q)), delta))

p4e <- ggplot(per_donor, aes(share, trait, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = .62, linewidth = .22,
               position = position_dodge(width = .72), alpha = .85) +
  geom_point(aes(colour = group), position = position_jitterdodge(
    jitter.width = .16, dodge.width = .72, seed = 42), size = .55, alpha = .6, show.legend = FALSE) +
  geom_text(data = tests, aes(x = 104, y = trait, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.1, colour = "grey25") +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_colour_manual(values = GCOL) +
  scale_y_discrete(labels = LAB) +
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(.01, 0))) +
  coord_cartesian(clip = "off") +
  labs(title = "Abundance-weighted gut functional trait capacity",
       subtitle = "Per donor: share of the detected community carried by trait-positive genomes",
       x = "Share of detected community (%)", y = NULL,
       caption = paste("Encoded capacity weighted by abundance, not expression or flux.",
                       "Denominator excludes genomes of unknown trait status.")) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(lineheight = .95),
        axis.title.x = element_text(size = 7.5),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        plot.caption = element_text(size = 5.5, colour = "grey35", hjust = 0),
        plot.title.position = "plot",
        legend.position = c(.88, .12), legend.background = element_blank(),
        legend.text = element_text(size = 7), legend.key.size = grid::unit(8, "pt"),
        plot.margin = margin(4, 40, 4, 4))

ggsave(file.path(OUTDIR, "Fig4e_TraitCapacity_n76.pdf"), p4e, width = 6.6, height = 3.4,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(OUTDIR, "Fig4e_TraitCapacity_n76.png"), p4e, width = 6.6, height = 3.4,
       units = "in", dpi = 600, bg = "white", type = "cairo-png")
saveRDS(p4e, file.path(OUTDIR, "Fig4e_TraitCapacity_n76.rds"))
write_csv(tests, file.path(OUTDIR, "Fig4e_trait_capacity_tests_n76.csv"))
write_csv(per_donor, file.path(OUTDIR, "Fig4e_trait_capacity_per_donor_n76.csv"))
writeLines(STATS, file.path(OUTDIR, "FIG4E_TRAIT_CAPACITY_n76.txt"))
message("DONE -> Fig4e_TraitCapacity_n76.pdf")
