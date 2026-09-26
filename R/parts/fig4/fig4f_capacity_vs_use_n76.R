#!/usr/bin/env Rscript
# Figure 4f (combined): predicted capacity against predicted use.
#
# Three numbers exist for each fermentation product and they are not the same
# number. The published 4f showed only the first.
#
#   Capacity   per genome, the maximum formation flux of the cytosolic product
#              on the colonic medium with biomass held at 10% of optimum, then
#              averaged over the MAGs detected in a donor with abundance
#              weights. One genome at a time, no competition, no shared pool.
#              A ceiling.
#
#   Produced   the community solution. Every detected modelled MAG placed in one
#              shared medium, abundance weighted, solved together by cooperative
#              tradeoff. This is the secretion each member actually carries in
#              that solution, summed over members, so it still counts a molecule
#              that a neighbour immediately eats.
#
#   Exported   the same solution's net exchange with the medium. What is left
#              after the community has cross-fed on its own output, and the only
#              one of the three that a stool sample could in principle see.
#
# The gap between them is the result. Net export is 0.02% to 0.1% of gross
# production, so almost everything these communities make is consumed inside the
# community, and three products with real capacity, butyrate, propionate and
# D-lactate, are never produced at all in the community solution: the draft
# reconstructions carry the enzymes but not a transport route out of the cell,
# and formation is not coupled to biomass, so a growth objective has no reason
# to carry the flux. That is the same failure the external IBDMDB test found
# from the other end, where gene abundance predicted transcript abundance but
# not the measured faecal acid.
#
# Sub-panel (i) is the substrate side of the same solution: what the community
# eats. Sub-panel (ii) is the three-stage flow. They are lettered (i) and (ii)
# because the whole figure is one panel of Figure 4, namely 4f.
#
# Ribbon widths are mean shares within a stage, which is what makes a flow
# diagram add up; the medians and the case-control tests are in the stats file.
# Stages are separately normalised because their units differ, so the figure
# compares composition between stages, never magnitude.
#
# Outputs -> result2/n76/fig4/Fig4f_CapacityVsUse_n76.{pdf,png,rds}
# Leaves Fig4f_FluxCapacity_n76.* untouched.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(ggforce)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- Sys.getenv("FIG4F_OUTDIR", file.path(N76, "fig4"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
GF     <- file.path(N76, "metabolic_support/gapfill_v1")
STATS  <- character(); say <- function(s) { message(s); STATS <<- c(STATS, s) }
GCOL   <- c(Control = "#00BFC4", Case = "#F8766D")

NUT  <- c(uptake_fibre_derived = "Fibre-derived sugars",
          uptake_host_derived  = "Host-derived glycans",
          uptake_simple_acid   = "Simple acids",
          uptake_amino_acid    = "Amino acids")
# Stacking order is fixed across all three stages so ribbons never cross.
PROD <- c(acetate = "Acetate", formate = "Formate", ethanol = "Ethanol",
          succinate = "Succinate", lactate_L = "L-lactate",
          lactate_D = "D-lactate", propionate = "Propionate", butyrate = "Butyrate")
SEC  <- c(secrete_ac = "acetate", secrete_for = "formate", secrete_etoh = "ethanol",
          secrete_succ = "succinate", secrete_lac__L = "lactate_L")
PCOL <- c(Acetate = "#4C6E7A", Formate = "#7FA8B5", Ethanol = "#C9D8DC",
          Succinate = "#E3A34B", `L-lactate` = "#B5C99A", `D-lactate` = "#8FAE6E",
          Propionate = "#9B6A9E", Butyrate = "#C2453B")

raw <- read_tsv(file.path(GF, "taxon_substrate_use.tsv.gz"), show_col_types = FALSE) %>%
  mutate(across(starts_with(c("uptake_", "secrete_")), ~ replace_na(.x, 0)))
net <- read_csv(file.path(GF, "micom_net_exchange_71.csv"), show_col_types = FALSE)
cap <- read_csv(file.path(N76, "fig4/Fig4f_flux_capacity_per_donor_n76.csv"),
                show_col_types = FALSE)

donors <- sort(net$sample)
grp <- setNames(factor(net$group[match(donors, net$sample)],
                       levels = c("Control", "Case")), donors)
say(sprintf("donors with a community solution: %d (%d case, %d control)",
            length(donors), sum(grp == "Case"), sum(grp == "Control")))

# ---- the three stages, one matrix each, donors x products --------------------
mat_capacity <- cap %>% filter(sample %in% donors) %>%
  select(sample, product, capacity) %>% pivot_wider(names_from = product, values_from = capacity)
mat_capacity <- as.matrix(mat_capacity[match(donors, mat_capacity$sample), names(PROD)])

produced <- raw %>% group_by(sample_id) %>%
  summarise(across(all_of(names(SEC)), sum), .groups = "drop") %>%
  rename(all_of(setNames(names(SEC), SEC)))
mat_produced <- matrix(0, length(donors), length(PROD),
                       dimnames = list(donors, names(PROD)))
mat_produced[, unname(SEC)] <- as.matrix(
  produced[match(donors, produced$sample_id), unname(SEC)])

mat_exported <- as.matrix(net[match(donors, net$sample), names(PROD)])
mat_exported[mat_exported < 0] <- 0   # a net-consumed metabolite is not an export

stopifnot(all(rowSums(mat_capacity) > 0), all(rowSums(mat_produced) > 0))

say(sprintf("net export as a share of gross production, median over donors: %.3f%%",
            100 * median(rowSums(mat_exported) / rowSums(mat_produced))))
say(sprintf("donors with any net export: %d of %d", sum(rowSums(mat_exported) > 0), length(donors)))
say("")
say("Community totals per donor, mmol per gDW per hour (median)")
say(sprintf("  %-11s %12s %12s %12s", "product", "capacity*", "produced", "exported"))
for (p in names(PROD)) say(sprintf("  %-11s %12.4f %12.2f %12.4f", PROD[[p]],
  median(mat_capacity[, p]), median(mat_produced[, p]), median(mat_exported[, p])))
say("  * capacity is an abundance-weighted per-genome mean, not a community total")

# ---- shares, and the case-control test on each stage -------------------------
shares <- function(m) {
  s <- rowSums(m); out <- sweep(m, 1, ifelse(s > 0, s, NA), "/") * 100
  out[is.na(out)] <- 0; out
}
cliff <- function(a, b) { m <- outer(a, b, "-"); (sum(m > 0) - sum(m < 0)) / length(m) }

stage_tests <- lapply(c("Capacity", "Produced", "Exported"), function(st) {
  m <- switch(st, Capacity = mat_capacity, Produced = mat_produced, Exported = mat_exported)
  keep <- rowSums(m) > 0
  sh <- shares(m)[keep, , drop = FALSE]; g <- grp[keep]
  lapply(names(PROD), function(p) {
    v <- sh[, p]
    if (length(unique(v)) == 1) return(tibble(stage = st, product = PROD[[p]],
      ctrl = 0, case = 0, delta = 0, p = NA_real_, n = sum(keep)))
    tibble(stage = st, product = PROD[[p]], n = sum(keep),
           ctrl = median(v[g == "Control"]), case = median(v[g == "Case"]),
           delta = cliff(v[g == "Case"], v[g == "Control"]),
           p = wilcox.test(v ~ g, exact = FALSE)$p.value)
  }) %>% bind_rows()
}) %>% bind_rows() %>% group_by(stage) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()

say("")
say("Share of each stage held by each product, median percent per donor")
say(sprintf("  %-10s %-11s %4s %9s %9s %8s %10s", "stage", "product", "n",
            "control", "case", "delta", "q"))
for (i in seq_len(nrow(stage_tests))) with(stage_tests[i, ],
  say(sprintf("  %-10s %-11s %4d %8.2f%% %8.2f%% %+8.2f %10.3g",
              stage, product, n, ctrl, case, delta, ifelse(is.na(q), 1, q))))

# ---- 4f (i): substrate use --------------------------------------------------
don <- raw %>% group_by(sample_id, group) %>%
  summarise(across(all_of(names(NUT)), sum), .groups = "drop") %>%
  mutate(named_total = rowSums(across(all_of(names(NUT)))))
don$group <- factor(don$group, levels = c("Control", "Case"))
nut_t <- lapply(names(NUT), function(c) {
  v <- 100 * don[[c]] / don$named_total
  tibble(label = NUT[[c]], ctrl = median(v[don$group == "Control"]),
         case = median(v[don$group == "Case"]),
         delta = cliff(v[don$group == "Case"], v[don$group == "Control"]),
         p = wilcox.test(v ~ don$group, exact = FALSE)$p.value)
}) %>% bind_rows() %>% mutate(q = p.adjust(p, "BH")) %>% arrange(delta)
nut_l <- don %>% pivot_longer(all_of(names(NUT)), names_to = "feature", values_to = "flux") %>%
  mutate(share = 100 * flux / named_total, label = NUT[feature])

say("")
say("Nutrient use, share of the named nutrient pool taken up")
say(sprintf("  %-22s %9s %9s %8s %10s", "class", "control", "case", "delta", "q"))
for (i in seq_len(nrow(nut_t))) with(nut_t[i, ],
  say(sprintf("  %-22s %8.2f%% %8.2f%% %+8.2f %10.3g", label, ctrl, case, delta, q)))

nut_t <- nut_t %>% mutate(txt = sprintf("δ = %+.2f, %s", delta,
  ifelse(q < 1e-4, "q < 1e-4", sprintf("q = %.3g", q))))
ord <- nut_t$label
nut_l$label <- factor(nut_l$label, ord); nut_t$label <- factor(nut_t$label, ord)
xmax <- max(nut_l$share); nut_t$x <- xmax * 1.02

pa <- ggplot(nut_l, aes(share, label, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = .62, linewidth = .22,
               position = position_dodge(width = .72), alpha = .85) +
  geom_point(aes(colour = group), position = position_jitterdodge(
    jitter.width = .16, dodge.width = .72, seed = 42), size = .5, alpha = .55,
    show.legend = FALSE) +
  geom_text(data = nut_t, aes(x = x, y = label, label = txt), inherit.aes = FALSE,
            hjust = 0, size = 2.0, colour = "grey25") +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_colour_manual(values = GCOL) +
  scale_x_continuous(limits = c(0, xmax * 1.5), expand = expansion(mult = c(.01, 0))) +
  coord_cartesian(clip = "off") +
  labs(title = "(i)  What the community consumes",
       subtitle = "Share of the named nutrient pool taken up, per donor",
       x = "% of fibre + host-derived + amino acid + simple acid uptake", y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(size = 6.5),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        plot.title.position = "plot",
        legend.position = "bottom", legend.justification = "center",
        legend.background = element_blank(), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 7), legend.key.size = grid::unit(8, "pt"),
        plot.margin = margin(4, 6, 4, 4))

# ---- 4f (ii): capacity -> produced -> exported -------------------------------
# Mean share, computed only over donors that have any flux at that stage, so the
# five communities that export nothing do not silently shrink the export column.
stage_mean <- function(m, arm) {
  keep <- rowSums(m) > 0 & grp == arm
  colMeans(shares(m)[keep, , drop = FALSE])
}
nodes <- lapply(c("Control", "Case"), function(arm) {
  lapply(seq_along(c("Capacity", "Produced", "Exported")), function(i) {
    st <- c("Capacity", "Produced", "Exported")[i]
    m <- switch(st, Capacity = mat_capacity, Produced = mat_produced, Exported = mat_exported)
    v <- stage_mean(m, arm)
    tibble(arm = arm, stage = st, x = i, product = PROD[names(PROD)],
           w = as.numeric(v[names(PROD)])) %>%
      mutate(ymax = cumsum(w), ymin = ymax - w)
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  mutate(arm = factor(arm, c("Control", "Case")),
         product = factor(product, PROD))

# One curved band per product per consecutive stage pair.
bands <- bind_rows(lapply(c(1, 2), function(i) {
  a <- filter(nodes, x == i); b <- filter(nodes, x == i + 1)
  j <- inner_join(a, b, c("arm", "product"), suffix = c("1", "2"))
  bind_rows(
    transmute(j, arm, product, id = paste(arm, product, i), x = x1 + .18, y = ymin1),
    transmute(j, arm, product, id = paste(arm, product, i), x = x1 + .18, y = ymax1),
    transmute(j, arm, product, id = paste(arm, product, i), x = x2 - .18, y = ymax2),
    transmute(j, arm, product, id = paste(arm, product, i), x = x2 - .18, y = ymin2))
})) %>% mutate(arm = factor(arm, c("Control", "Case")), product = factor(product, PROD))

lab <- nodes %>% filter(w >= 6) %>%
  mutate(y = (ymin + ymax) / 2, txt = sprintf("%s %.0f%%", product, w))
# Butyrate is the product this section is about and its capacity band is too
# thin to hold a label, so it is called out beside the column instead of being
# the one unlabelled band in the figure.
lab_out <- nodes %>% filter(stage == "Capacity", w < 6) %>%
  mutate(y = (ymin + ymax) / 2, txt = sprintf("%s %.1f%%", product, w))

pb <- ggplot() +
  geom_diagonal_wide(data = bands, aes(x, y, group = id, fill = product),
                     alpha = .42, colour = NA, strength = .55) +
  geom_rect(data = nodes, aes(xmin = x - .18, xmax = x + .18,
                              ymin = ymin, ymax = ymax, fill = product),
            colour = "white", linewidth = .18) +
  geom_text(data = lab, aes(x = x, y = y, label = txt), size = 1.65, colour = "grey15") +
  geom_segment(data = lab_out, aes(x = x - .18, xend = x - .26, y = y, yend = y),
               linewidth = .2, colour = "grey35") +
  geom_text(data = lab_out, aes(x = x - .29, y = y, label = txt), size = 1.65,
            hjust = 1, colour = "grey15") +
  facet_wrap(~arm, nrow = 1) +
  scale_fill_manual(values = PCOL, name = NULL) +
  scale_x_continuous(breaks = 1:3, labels = c("Capacity", "Produced", "Exported"),
                     expand = expansion(mult = c(.42, .14))) +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = .02)) +
  coord_cartesian(clip = "off") +
  labs(title = "(ii)  Capacity against predicted use",
       subtitle = paste("Mean share of each stage, per arm. Stages are normalised separately:",
                        "capacity is a per-genome ceiling,\nproduced and exported are the",
                        "community solution before and after its own cross-feeding."),
       x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(size = 7, colour = "black"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7.5),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.2, colour = "grey35"),
        plot.title.position = "plot",
        legend.position = "bottom", legend.text = element_text(size = 6.2),
        legend.key.size = grid::unit(7, "pt"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(4, 6, 4, 4))

# Keep native text and vector objects directly editable in Affinity. Background
# rectangles and panel clipping are unnecessary here; geometry/data are unchanged.
editable_theme <- theme(
  plot.background = element_blank(), panel.background = element_blank(),
  legend.background = element_blank(), legend.box.background = element_blank(),
  legend.key = element_blank())
pa <- pa + editable_theme
pb <- pb + editable_theme

p <- pa / pb + plot_layout(heights = c(4, 6)) +
  plot_annotation(caption = paste(
    "MICOM cooperative tradeoff at fraction 0.5 on a western-diet colonic medium,",
    "71 of 76 donors solved.\nNet export is a median 0.02% of gross production, so",
    "almost everything produced is consumed within the community.\nPredicted fluxes",
    "under one stated diet, not rates and not measurements."),
    theme = theme(plot.background = element_blank(),
                  plot.caption = element_text(size = 5.5, colour = "grey35", hjust = 0)))

ggsave(file.path(OUTDIR, "Fig4f_CapacityVsUse_n76.pdf"), p, width = 6.4, height = 6.2,
       units = "in", device = cairo_pdf, bg = "transparent")
ggsave(file.path(OUTDIR, "Fig4f_CapacityVsUse_n76.png"), p, width = 6.4, height = 6.2,
       units = "in", dpi = 600, bg = "transparent", type = "cairo-png")
saveRDS(p, file.path(OUTDIR, "Fig4f_CapacityVsUse_n76.rds"))
write_csv(stage_tests, file.path(OUTDIR, "Fig4f_capacity_vs_use_stage_tests_n76.csv"))
write_csv(nodes, file.path(OUTDIR, "Fig4f_capacity_vs_use_flow_n76.csv"))
writeLines(STATS, file.path(OUTDIR, "FIG4F_CAPACITY_VS_USE_n76.txt"))
message("DONE -> Fig4f_CapacityVsUse_n76.pdf")
