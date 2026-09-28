#!/usr/bin/env Rscript
# Figure 4f (replacement): what each donor's community actually eats and excretes.
#
# The panel this replaces reported per-genome formation capacity, summed with
# abundance weights. That number is a ceiling: it asks what a genome could make
# if its whole network were pushed at one product, one genome at a time, with no
# competition for substrate and no shared pool. It cannot say whether the
# community would do it.
#
# This panel reports the community solution instead. Every modelled MAG detected
# in a donor is placed in one shared medium, weighted by its relative abundance,
# and solved together by cooperative tradeoff (MICOM 0.39.1, fraction 0.5). Each
# exchange flux is then signed and real: a negative flux is uptake from the
# shared pool, a positive one is excretion into it, competition is resolved, and
# cross-feeding between members is already netted out.
#
# Two quantities are shown.
#
#   a. Where the carbon and nitrogen come from. Uptake fluxes are binned into
#      fibre-derived sugars, host-derived glycans, amino acids and simple acids,
#      and expressed as a share of that named nutrient pool. The named pool is
#      only 12% (controls) to 17% (patients) of all exchange uptake; the rest is
#      water, protons, ions and trace cofactors, which carry no dietary meaning
#      and would swamp the shares if left in the denominator.
#
#   b. What comes back out, per unit of substrate consumed. Product excretion is
#      divided by that donor's total uptake flux, so it is a yield and not a
#      community size. Absolute excretion tracks how many members a community
#      has, and control communities are the larger ones (252 against 115 MAGs
#      detected at the median), so the raw fluxes would report community size.
#
# Neither quantity is a measurement. Both are the flux distribution of a
# constrained optimisation on one stated diet.
#
# Outputs -> result2/n76/fig4/Fig4f_SubstrateUse_n76.{pdf,png,rds}
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "fig4")
GF     <- file.path(N76, "metabolic_support/gapfill_v1")
STATS  <- character(); say <- function(s) { message(s); STATS <<- c(STATS, s) }
GCOL   <- c(Control = "#00BFC4", Case = "#F8766D")

NUT <- c(uptake_fibre_derived = "Fibre-derived sugars",
         uptake_host_derived  = "Host-derived glycans",
         uptake_simple_acid   = "Simple acids",
         uptake_amino_acid    = "Amino acids")
PROD <- c(secrete_ac = "Acetate", secrete_for = "Formate", secrete_lac__L = "L-lactate",
          secrete_succ = "Succinate", secrete_etoh = "Ethanol",
          secrete_co2 = "CO2", secrete_h2 = "H2")

# Axis labels are drawn as plotmath so the two gases carry real subscripts; the
# embedded Helvetica subset has no U+2082, so a literal CO₂ would drop a glyph.
AXIS_EXPR <- c(CO2 = "CO[2]", H2 = "H[2]")
axis_labeller <- function(x) parse(text = ifelse(
  x %in% names(AXIS_EXPR), AXIS_EXPR[x], sprintf("'%s'", x)))

raw <- read_tsv(file.path(GF, "taxon_substrate_use.tsv.gz"), show_col_types = FALSE) %>%
  mutate(across(starts_with(c("uptake_", "secrete_")), ~ replace_na(.x, 0)))

say(sprintf("taxon-donor rows %d | donors %d | families %d",
            nrow(raw), n_distinct(raw$sample_id), n_distinct(raw$family)))

don <- raw %>% group_by(sample_id, group) %>%
  summarise(across(starts_with(c("uptake_", "secrete_")), sum), .groups = "drop") %>%
  mutate(named_total = rowSums(across(all_of(names(NUT)))),
         product_total = rowSums(across(all_of(names(PROD)))))
stopifnot(nrow(don) == 76, all(don$named_total > 0), all(don$uptake_total > 0))
don$group <- factor(don$group, levels = c("Control", "Case"))

say(sprintf("donors solved: %d case, %d control",
            sum(don$group == "Case"), sum(don$group == "Control")))
say(sprintf("named nutrient uptake as a share of all exchange uptake: control median %.1f%%, case median %.1f%%",
            100 * median((don$named_total / don$uptake_total)[don$group == "Control"]),
            100 * median((don$named_total / don$uptake_total)[don$group == "Case"])))

cliff <- function(a, b) { m <- outer(a, b, "-"); (sum(m > 0) - sum(m < 0)) / length(m) }

# One test family per panel: the four nutrient classes, then the seven products.
# They answer different questions on different denominators, so pooling them
# under one BH correction would borrow significance across the two.
tabulate <- function(cols, denom, labels) {
  lapply(names(cols), function(c) {
    v <- 100 * don[[c]] / don[[denom]]
    tibble(feature = c, label = labels[[c]],
           ctrl = median(v[don$group == "Control"]), case = median(v[don$group == "Case"]),
           delta = cliff(v[don$group == "Case"], v[don$group == "Control"]),
           p = wilcox.test(v ~ don$group, exact = FALSE)$p.value)
  }) %>% bind_rows() %>% mutate(q = p.adjust(p, "BH")) %>% arrange(delta)
}

long <- function(cols, denom, labels) {
  don %>% select(sample_id, group, all_of(c(names(cols), denom))) %>%
    pivot_longer(all_of(names(cols)), names_to = "feature", values_to = "flux") %>%
    mutate(share = 100 * flux / .data[[denom]], label = labels[feature])
}

nut_t <- tabulate(NUT, "named_total", NUT)
pro_t <- tabulate(PROD, "uptake_total", PROD)
nut_l <- long(NUT, "named_total", NUT)
pro_l <- long(PROD, "uptake_total", PROD)

report <- function(title, t, unit) {
  say(""); say(title)
  say(sprintf("  %-22s %9s %9s %8s %10s", "feature", "control", "case", "delta", "q"))
  for (i in seq_len(nrow(t))) with(t[i, ],
    say(sprintf("  %-22s %8.2f%% %8.2f%% %+8.2f %10.3g", label, ctrl, case, delta, q)))
  say(sprintf("  unit: %s", unit))
}
report("Nutrient use, share of the named nutrient pool taken up", nut_t,
       "percent of fibre + host-derived + amino acid + simple acid uptake")
report("Product yield, excretion per unit of substrate consumed", pro_t,
       "percent of that donor's total exchange uptake flux")

# ---- who eats the amino acids ------------------------------------------------
# Panel a's largest shift is amino acid use. The community solution attributes
# every flux to a member, so the family carrying that shift is recoverable
# without a second model run.
fam <- raw %>% group_by(sample_id, group, family) %>%
  summarise(aa = sum(uptake_amino_acid), .groups = "drop") %>%
  group_by(sample_id) %>% mutate(share = 100 * aa / sum(aa)) %>% ungroup()
top <- fam %>% group_by(family) %>% summarise(total = sum(aa)) %>%
  slice_max(total, n = 10) %>% pull(family)
fam_t <- lapply(top, function(f) {
  v <- fam %>% filter(family == f) %>%
    right_join(distinct(don, sample_id, group), "sample_id") %>%
    mutate(share = replace_na(share, 0), group = coalesce(group.y, group.x))
  tibble(family = f,
         ctrl = median(v$share[v$group == "Control"]),
         case = median(v$share[v$group == "Case"]),
         delta = cliff(v$share[v$group == "Case"], v$share[v$group == "Control"]),
         p = wilcox.test(v$share ~ v$group, exact = FALSE)$p.value)
}) %>% bind_rows() %>% mutate(q = p.adjust(p, "BH")) %>% arrange(desc(delta))
say(""); say("Family share of community amino acid uptake, median percent per donor")
say(sprintf("  %-24s %9s %9s %8s %10s", "family", "control", "case", "delta", "q"))
for (i in seq_len(nrow(fam_t))) with(fam_t[i, ],
  say(sprintf("  %-24s %8.2f%% %8.2f%% %+8.2f %10.3g", family, ctrl, case, delta, q)))

# ---- panels ------------------------------------------------------------------
lab_of <- function(t) t %>% mutate(
  txt = sprintf("δ = %+.2f, %s", delta,
                ifelse(q < 1e-4, "q < 1e-4", sprintf("q = %.3g", q))))

panel <- function(dat, t, xlab, title, subtitle) {
  t <- lab_of(t)
  ord <- t$label
  dat$label <- factor(dat$label, ord); t$label <- factor(t$label, ord)
  xmax <- max(dat$share)
  t$x <- xmax * 1.02
  ggplot(dat, aes(share, label, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = .62, linewidth = .22,
                 position = position_dodge(width = .72), alpha = .85) +
    geom_point(aes(colour = group), position = position_jitterdodge(
      jitter.width = .16, dodge.width = .72, seed = 42), size = .5, alpha = .55,
      show.legend = FALSE) +
    geom_text(data = t, aes(x = x, y = label, label = txt), inherit.aes = FALSE,
              hjust = 0, size = 2.0, colour = "grey25") +
    scale_fill_manual(values = GCOL, name = NULL) +
    scale_colour_manual(values = GCOL) +
    scale_y_discrete(labels = axis_labeller) +
    scale_x_continuous(limits = c(0, xmax * 1.52),
                       expand = expansion(mult = c(.01, 0))) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text = element_text(size = 7, colour = "black"),
          axis.title.x = element_text(size = 7),
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 6.5, colour = "grey35"),
          plot.title.position = "plot",
          legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(size = 7), legend.key.size = grid::unit(8, "pt"),
          plot.margin = margin(4, 6, 4, 4))
}

# One shared legend at the foot: an in-panel key sat on the effect-size labels,
# which are drawn outside the panel on the right where the key had to go.
pa <- panel(nut_l, nut_t, "% of the named nutrient pool taken up",
            "What the community consumes",
            "Per donor, from the community flux solution")
pb <- panel(pro_l, pro_t, "% of total exchange uptake flux (yield)",
            "What it excretes per unit consumed",
            "Net community exchange, so cross-feeding is already subtracted")

p4f <- pa / pb + plot_layout(heights = c(4, 7), guides = "collect") +
  plot_annotation(caption = paste(
    "MICOM cooperative tradeoff at fraction 0.5 on a western-diet colonic medium,",
    sprintf("%d of 76 donors solved.\nPredicted community fluxes under one stated diet,", nrow(don)),
    "not rates and not measurements."),
    theme = theme(plot.caption = element_text(size = 5.5, colour = "grey35", hjust = 0),
                  legend.position = "bottom"))

ggsave(file.path(OUTDIR, "Fig4f_SubstrateUse_n76.pdf"), p4f, width = 5.2, height = 5.0,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(OUTDIR, "Fig4f_SubstrateUse_n76.png"), p4f, width = 5.2, height = 5.0,
       units = "in", dpi = 600, bg = "white", type = "cairo-png")
saveRDS(p4f, file.path(OUTDIR, "Fig4f_SubstrateUse_n76.rds"))
write_csv(bind_rows(mutate(nut_t, panel = "nutrient"), mutate(pro_t, panel = "product")),
          file.path(OUTDIR, "Fig4f_substrate_use_tests_n76.csv"))
write_csv(fam_t, file.path(OUTDIR, "Fig4f_substrate_amino_by_family_n76.csv"))
write_csv(don, file.path(OUTDIR, "Fig4f_substrate_use_per_donor_n76.csv"))
writeLines(STATS, file.path(OUTDIR, "FIG4F_SUBSTRATE_USE_n76.txt"))
message("DONE -> Fig4f_SubstrateUse_n76.pdf")
