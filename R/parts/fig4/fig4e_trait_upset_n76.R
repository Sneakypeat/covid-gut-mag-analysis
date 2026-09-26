#!/usr/bin/env Rscript
# Figure 4e: gut functional trait UpSet, carrying ANCOM-BC2 status and abundance.
#
# The first UpSet counted genomes equally, so plant CAZyme hits (527 of 584 MAGs)
# dominated the display while carrying no case-control signal, and the butyrate
# terminal pair (17 genomes) was invisible despite a 21-fold abundance collapse.
# The boxplot that replaced it fixed the weighting but lost the co-occurrence
# structure: which traits travel together in the same genome.
#
# This keeps the UpSet and adds the two missing dimensions to it.
#   top    per-donor share of the community carried by each trait combination,
#          cases against controls, which is the abundance weighting
#   middle genome counts split by ANCOM-BC2 direction, which is the differential
#          -abundance status of the genomes making up each combination
#   bottom the intersection matrix, with an acid-requirement strip
#
# Intersections are defined over the four traits called across all 584 genomes.
# Model-derived acid requirement is reported separately from the four genomic
# traits, as the positive fraction among modelled members of each column.
# Unknown model-derived calls are excluded from that strip's denominator.
#
# This is encoded capacity weighted by abundance. It is not expression, and it is
# not flux: see FLUX_FEASIBILITY_n76.md for why the draft reconstructions cannot
# supply a flux answer.
#
# Outputs -> result2/n76/fig4/Fig4e_TraitUpSet_n76.{pdf,png,rds}
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- Sys.getenv("FIG4E_OUTDIR", file.path(N76, "fig4"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
TRAITS <- file.path(N76, "metabolic_support/gut_traits_v2")
STATS  <- character(); say <- function(s) { message(s); STATS <<- c(STATS, s) }
PT <- 2.845276  # points per mm, for size = pt/PT

SETS <- c(plant_backbone_cazyme  = "Plant-associated CAZyme hits",
          acetate_terminal_pair  = "Acetate pair (pta/ackA)",
          mucin_gh_repertoire    = "Putative mucin-GH combination",
          butyrate_terminal_pair = "Butyrate pair (ptb/buk)")
GCOL <- c(Control = "#00BFC4", Case = "#F8766D")
ACOL <- c(Depleted = "#0F7C86", `Not significant` = "#C9C9C9", Enriched = "#C2453B")

wide <- read_tsv(file.path(TRAITS, "genome_traits_wide.tsv"), show_col_types = FALSE)
ab   <- read.delim(file.path(INDIR, "MAG_relative_abundance_percent.tsv"),
                   check.names = FALSE, row.names = 1)
meta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE)
stopifnot(setequal(wide$catalog_id, rownames(ab)), nrow(meta) == 76)
ab  <- as.matrix(ab[wide$catalog_id, meta$sample])
grp <- factor(meta$group, levels = c("Control", "Case"))

cliff <- function(a, b) { m <- outer(a, b, "-"); (sum(m > 0) - sum(m < 0)) / length(m) }

# ---- intersection membership -------------------------------------------------
memb <- sapply(names(SETS), function(s) wide[[s]] == "positive")
key  <- apply(memb, 1, function(v) paste0(as.integer(v), collapse = ""))
wide$key <- key
wide$ancom <- factor(wide$ancom_direction, names(ACOL))
stopifnot(!any(is.na(wide$ancom)))

# ---- abundance carried by each intersection, per donor -----------------------
tot <- colSums(ab)
share <- t(rowsum(ab, key)) / tot * 100          # donors x intersections
per_donor <- as_tibble(share, rownames = "sample") %>%
  pivot_longer(-sample, names_to = "key", values_to = "share") %>%
  mutate(group = grp[match(sample, meta$sample)])

tests <- per_donor %>% group_by(key) %>%
  summarise(case_med = median(share[group == "Case"]),
            ctrl_med = median(share[group == "Control"]),
            case_lo = quantile(share[group == "Case"], .25),
            case_hi = quantile(share[group == "Case"], .75),
            ctrl_lo = quantile(share[group == "Control"], .25),
            ctrl_hi = quantile(share[group == "Control"], .75),
            delta = cliff(share[group == "Case"], share[group == "Control"]),
            p = wilcox.test(share ~ group, exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(q = p.adjust(p, "BH"))

counts <- wide %>% count(key, ancom, name = "n") %>% complete(key, ancom, fill = list(n = 0))
n_tot  <- wide %>% count(key, name = "n_genomes")
acid   <- wide %>% filter(strict_acid_uptake != "unknown") %>%
  group_by(key) %>% summarise(acid_pos = sum(strict_acid_uptake == "positive"),
                              acid_n = n(), .groups = "drop")

# Columns run from the largest control-side contributor to the smallest, so the
# reading order is "what healthy communities are mostly made of" first.
ord <- tests %>% arrange(desc(ctrl_med)) %>% pull(key)
lev <- function(d) mutate(d, key = factor(key, ord))
tests <- lev(tests); counts <- lev(counts); n_tot <- lev(n_tot)
per_donor <- lev(per_donor); acid <- lev(acid); wide <- lev(wide)

# ---- reporting ---------------------------------------------------------------
say("Trait-combination UpSet: abundance share, ANCOM-BC2 status, 584 genomes / 76 donors")
say(sprintf("%-6s %-42s %7s %8s %8s %7s %9s %s", "key", "combination", "genomes",
            "ctrl%", "case%", "delta", "q", "dep/ns/enr"))
for (k in ord) {
  bits <- strsplit(k, "")[[1]] == "1"
  lab  <- if (any(bits)) paste(names(SETS)[bits], collapse = " + ") else "none of the four"
  tt <- tests[tests$key == k, ]; cc <- counts[counts$key == k, ]
  a  <- acid[acid$key == k, ]
  say(sprintf("%-6s %-42s %7d %8.2f %8.2f %+7.2f %9.3g %d/%d/%d  acid %s",
              k, substr(lab, 1, 42), n_tot$n_genomes[n_tot$key == k],
              tt$ctrl_med, tt$case_med, tt$delta, tt$q,
              cc$n[cc$ancom == "Depleted"], cc$n[cc$ancom == "Not significant"],
              cc$n[cc$ancom == "Enriched"],
              if (nrow(a)) sprintf("%d/%d", a$acid_pos, a$acid_n) else "0/0"))
}

# ---- block A: abundance share, cases against controls ------------------------
barsA <- tests %>%
  select(key, Control = ctrl_med, Case = case_med) %>%
  pivot_longer(-key, names_to = "group", values_to = "med") %>%
  left_join(tests %>% select(key, ctrl_lo, ctrl_hi, case_lo, case_hi), by = "key") %>%
  mutate(group = factor(group, names(GCOL)),
         lo = ifelse(group == "Case", case_lo, ctrl_lo),
         hi = ifelse(group == "Case", case_hi, ctrl_hi))
starA <- tests %>% mutate(
  star = ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", ""))),
  lab  = ifelse(star == "", "", sprintf("%s  %+.2f", star, delta)),
  ytop = pmax(case_hi, ctrl_hi))
# Precomputed, not evaluated inside aes(): the saved plot object is reloaded by
# fig4_square_panels_n76.R, where YMAX would not exist.
starA$ystar <- (sqrt(starA$ytop) + sqrt(max(barsA$hi) * 1.55) * .045)^2
BRK  <- c(0, 0.5, 2, 5, 10, 20, 40, 60)
YMAX <- max(barsA$hi) * 1.55

pA <- ggplot(barsA, aes(key, med, fill = group)) +
  geom_col(position = position_dodge(width = .60), width = .56, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = .60),
                width = .15, linewidth = .18, colour = "grey30") +
  geom_text(data = starA, aes(key, ystar, label = lab),
            inherit.aes = FALSE, size = 5.2 / PT, colour = "grey20", vjust = 0) +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_y_sqrt(limits = c(0, max(barsA$hi) * 1.55), breaks = BRK, labels = as.character(BRK),
               expand = expansion(mult = c(0, .02))) +
  labs(y = "Share of community (%)\nsquare-root scale", x = NULL) +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 6), axis.text.y = element_text(size = 5.5),
        legend.background = element_blank(),
        legend.key.size = grid::unit(6, "pt"), legend.text = element_text(size = 5.5),
        plot.margin = margin(1, 2, 0, 2))

# ---- block B: genome counts split by ANCOM-BC2 direction ---------------------
pB <- ggplot(counts, aes(key, n, fill = ancom)) +
  geom_col(width = .56, colour = NA) +
  geom_text(data = n_tot, aes(key, n_genomes, label = n_genomes), inherit.aes = FALSE,
            vjust = -0.35, size = 5 / PT, colour = "grey25") +
  scale_fill_manual(values = ACOL, name = NULL) +
  scale_y_continuous(breaks = c(0, 100, 200), expand = expansion(mult = c(0, .18))) +
  labs(y = "Genomes", x = NULL) +
  theme_classic(base_size = 6.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 6), axis.text.y = element_text(size = 5.5),
        legend.background = element_blank(),
        legend.key.size = grid::unit(6, "pt"), legend.text = element_text(size = 5.5),
        plot.margin = margin(0, 2, 0, 2))

# ---- block C: intersection matrix -------------------------------------------
mat <- expand_grid(key = factor(ord, ord), set = names(SETS)) %>%
  mutate(on = mapply(function(k, s) strsplit(as.character(k), "")[[1]][match(s, names(SETS))] == "1",
                     key, set),
         set = factor(set, rev(names(SETS))))
seg <- mat %>% filter(on) %>% group_by(key) %>%
  summarise(ymin = min(as.integer(set)), ymax = max(as.integer(set)), .groups = "drop") %>%
  filter(ymax > ymin)
setlab <- sprintf("%s  (%d)", SETS[rev(names(SETS))],
                  colSums(memb)[rev(names(SETS))])

pC <- ggplot(mat, aes(key, set)) +
  geom_point(colour = "#DEDEDE", size = 1.35) +
  geom_segment(data = seg, aes(x = key, xend = key, y = ymin, yend = ymax),
               inherit.aes = FALSE, colour = "#2E3B40", linewidth = .32) +
  geom_point(data = filter(mat, on), colour = "#2E3B40", size = 1.5) +
  scale_y_discrete(labels = setlab) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.5, colour = "black", hjust = 1),
        plot.margin = margin(0, 2, 0, 2))

# ---- block D: acid-requirement strip over the modelled members ---------------
acid_strip <- acid %>% mutate(frac = acid_pos / acid_n,
                              lab = sprintf("%d/%d", acid_pos, acid_n))
pD <- ggplot(acid_strip, aes(key, 1, fill = frac)) +
  # Zero fractions remain explicitly labelled; omit their invisible white tiles.
  geom_tile(data = filter(acid_strip, frac > 0), colour = NA, width = .70, height = .8) +
  geom_text(aes(label = lab), size = 5 / PT, colour = "grey10") +
  scale_fill_gradient(low = "white", high = "#F2B705", limits = c(0, 1), guide = "none") +
  scale_y_continuous(breaks = 1, labels = "Model-derived\nacid requirement") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5, colour = "black", hjust = 1, lineheight = .95),
        plot.margin = margin(0, 2, 1, 2))

# All panels share the same column centres. A narrower page compacts the column
# spacing without changing their membership/order. No white background or bar
# outlines: keep individual vector objects accessible in Affinity.
editable_theme <- theme(
  plot.background = element_blank(), panel.background = element_blank(),
  legend.background = element_blank(), legend.box.background = element_blank(),
  legend.key = element_blank())
pA <- pA + coord_cartesian(clip = "off") + editable_theme
pB <- pB + coord_cartesian(clip = "off") + editable_theme
pC <- pC + coord_cartesian(clip = "off") + editable_theme
pD <- pD + coord_cartesian(clip = "off") + editable_theme

p4e <- (pA / pB / pC / pD) +
  plot_layout(heights = c(1.70, 1.00, .85, .38), guides = "collect") +
  plot_annotation(
    title = "Gut functional trait combinations, weighted by abundance\nand split by differential-abundance status",
    caption = paste("Bars: medians and interquartile ranges. Labels: Cliff's delta and Benjamini-Hochberg q",
                    sprintf("across %d combinations.", length(ord)),
                    "\n*** q < 0.001, ** q < 0.01, * q < 0.05.",
                    "\nEncoded capacity weighted by abundance, not expression or flux."),
    theme = theme(legend.position = "top", legend.justification = "right",
                  plot.background = element_blank(),
                  legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0),
                  legend.key.size = grid::unit(6, "pt"),
                  legend.text = element_text(size = 5.5),
                  plot.title = element_text(face = "bold", size = 7.5),
                  plot.caption = element_text(size = 5, colour = "grey35", hjust = 0),
                  plot.title.position = "plot",
                  # left inset keeps the title clear of the figure-level "e" tag
                  plot.margin = margin(2, 2, 2, 16)))

ggsave(file.path(OUTDIR, "Fig4e_TraitUpSet_n76.pdf"), p4e, width = 5.8, height = 4.0,
       units = "in", device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "Fig4e_TraitUpSet_n76.png"), p4e, width = 5.8, height = 4.0,
       units = "in", dpi = 600, bg = "transparent", type = "cairo-png")
saveRDS(p4e, file.path(OUTDIR, "Fig4e_TraitUpSet_n76.rds"))
write_csv(tests %>% left_join(n_tot, "key") %>% left_join(acid, "key"),
          file.path(OUTDIR, "Fig4e_trait_upset_tests_n76.csv"))
write_csv(counts, file.path(OUTDIR, "Fig4e_trait_upset_ancom_counts_n76.csv"))
write_csv(per_donor, file.path(OUTDIR, "Fig4e_trait_upset_per_donor_n76.csv"))
writeLines(STATS, file.path(OUTDIR, "FIG4E_TRAIT_UPSET_n76.txt"))
message("DONE -> Fig4e_TraitUpSet_n76.pdf")
