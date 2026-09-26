# =============================================================================
# Standalone prevalence legend: one gradient strip per phylum, light (low
# prevalence) to saturated (high), matching the shading used on the MAG points
# in Figure 1. Same rule as shade_by_prev() in figure1_n76.R:
#   colour = white * (1 - p) + phylum colour * p,  p = scaled prevalence^gamma
#
# Outputs -> result2/n76/fig1/Fig1_PrevalenceLegend_n76.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(phyloseq)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "fig1")
source(file.path(BASE, "taxon_italics_n76.R"))   # expr_from_md(): taxon names in italics

PREV_GAMMA <- 0.5     # as in figure1_n76.R
STEPS      <- 120     # gradient resolution

pal <- read_tsv(file.path(N76, "phylum_palette_n76.tsv"), show_col_types = FALSE)
cols <- setNames(pal[[2]], pal[[1]])

# phyla present in the catalogue, ordered by how many MAGs they carry
tax <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                  quote = "", comment.char = "") %>%
  transmute(taxon = catalog_id,
            Phylum = sub("^p__", "", sapply(strsplit(gtdbtk_classification, ";"), `[`, 2))) %>%
  filter(nzchar(Phylum))

ps  <- readRDS(file.path(N76, "ps_mags_n76.rds"))
otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
rel <- sweep(otu, 2, colSums(otu), "/")
prev <- rowMeans(rel > 1e-4)                     # the 0.01% detection floor used in Fig. 1f

# same order as the phylum panel (Fig. 1g), taken from the file that panel
# writes, so the two stay in sync if its ordering ever changes
fisher_order <- read_csv(file.path(N76, "taxon_level/directional_fisher_phylum_n76.csv"),
                         show_col_types = FALSE)$clade
phy <- tax %>% filter(taxon %in% names(prev)) %>%
  mutate(prev = prev[taxon]) %>%
  group_by(Phylum) %>%
  summarise(n_mags = n(), prev_max = max(prev), prev_med = median(prev), .groups = "drop") %>%
  filter(Phylum %in% names(cols), Phylum %in% fisher_order) %>%
  arrange(match(Phylum, fisher_order))               # file order is bottom-to-top, as in Fig. 1g

shade <- function(base_hex, p) {
  p <- pmin(pmax(p, 0), 1); m <- col2rgb(base_hex) / 255
  grDevices::rgb((1 - p) + p * m[1, ], (1 - p) + p * m[2, ], (1 - p) + p * m[3, ])
}

grid <- phy %>%
  mutate(Phylum = factor(Phylum, phy$Phylum)) %>%
  tidyr::crossing(step = seq_len(STEPS)) %>%
  mutate(p = (step / STEPS)^PREV_GAMMA,
         fill = shade(cols[as.character(Phylum)], p))

BAR_H <- 0.62          # strip height in row units; the rest is the gap between strips

p_leg <- ggplot(grid, aes(step, Phylum, fill = fill)) +
  geom_tile(height = BAR_H, width = 1) + scale_fill_identity() +
  # one black frame per strip, drawn over the gradient
  geom_tile(data = phy %>% mutate(Phylum = factor(Phylum, phy$Phylum)),
            aes(x = (STEPS + 1) / 2, y = Phylum), inherit.aes = FALSE,
            width = STEPS, height = BAR_H, fill = NA, colour = "black", linewidth = 0.35) +
  scale_x_continuous(expand = expansion(mult = 0.012),
                     breaks = c(1, STEPS), labels = c("low", "high")) +
  scale_y_discrete(expand = expansion(add = 0.45),
                   labels = function(x) expr_from_md(md_taxon(x))) +
  labs(title = "Prevalence", x = NULL, y = NULL,
       caption = paste("Shading of the MAG points in Fig. 1: white to phylum colour",
                       "with donor prevalence (0.01% detection floor, sqrt scale).",
                       sep = "\n")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.text.x = element_text(size = 6.5, colour = "grey25"),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold", size = 9),
        plot.caption = element_text(size = 5.5, colour = "grey45", hjust = 0),
        plot.margin = margin(4, 10, 4, 4, "pt"),
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

ggsave(file.path(OUTDIR, "Fig1_PrevalenceLegend_n76.pdf"), p_leg,
       width = 3.6, height = 0.20 * nrow(phy) + 1.0, device = cairo_pdf, bg = "transparent")
write_csv(phy, file.path(OUTDIR, "Fig1_PrevalenceLegend_phyla_n76.csv"))
message(sprintf("DONE -> %s (%d phyla)", file.path(OUTDIR, "Fig1_PrevalenceLegend_n76.pdf"), nrow(phy)))
