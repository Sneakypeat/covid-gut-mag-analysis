#!/usr/bin/env Rscript
# Figure 4f: abundance-weighted fermentation-product formation capacity per donor.
#
# Panel e answers "who carries the genes, and how abundant are they". This one
# answers the question that caption could not: given the diet, can the network
# those genes sit in actually form the compound, and how much of that capacity
# does each donor's community carry.
#
# Per-genome capacity is the maximum formation flux of the cytosolic product on
# the western-diet colonic medium, from gapfilled CarveMe models. A demand
# reaction is used rather than the secretion reaction because the export step is
# the part of these reconstructions that is missing: all eight modelled genomes
# carrying ptb and buk have both enzymes and only one has a butyrate transporter.
# Methods frozen in gapfill_v1/FLUX_METHODS_LOCK.md.
#
# Per donor, capacity is the abundance-weighted mean over the modelled MAGs
# detected in that donor, so it is a property of the community, not of a genome.
#
# This is predicted capacity under one stated diet. Not a rate, not expression,
# not a measurement.
#
# Outputs -> result2/n76/fig4/Fig4f_FluxCapacity_n76.{pdf,png,rds}
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "fig4")
GF     <- file.path(N76, "metabolic_support/gapfill_v1")
FLOOR  <- 0.01      # percent, the detection floor used throughout this project
STATS  <- character(); say <- function(s) { message(s); STATS <<- c(STATS, s) }
PT <- 2.845276

LAB <- c(butyrate = "Butyrate", acetate = "Acetate", propionate = "Propionate",
         lactate_L = "L-lactate", lactate_D = "D-lactate",
         succinate = "Succinate", formate = "Formate", ethanol = "Ethanol")
GCOL <- c(Control = "#00BFC4", Case = "#F8766D")

# FVA table over the full 584-model catalogue. `<p>_max` is the growth-coupled
# maximum (FVA at fraction_of_optimum = 0.10), the same quantity the earlier
# 339-model scan called `<p>_growth_coupled`; `<p>_min` is the other end of the
# interval and is reported in the stats file.
cap  <- read_csv(file.path(GF, "fva_584.csv"), show_col_types = FALSE) %>%
  rename_with(~ sub("_max$", "_growth_coupled", .x), ends_with("_max"))
ab   <- read.delim(file.path(INDIR, "MAG_relative_abundance_percent.tsv"),
                   check.names = FALSE, row.names = 1)
meta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE)
stopifnot(nrow(meta) == 76, all(cap$catalog_id %in% rownames(ab)))

A   <- as.matrix(ab[cap$catalog_id, meta$sample])   # modelled MAGs x donors
det <- A > FLOOR
grp <- factor(meta$group, levels = c("Control", "Case"))

say(sprintf("modelled MAGs %d | donors %d | detected per donor: case median %.0f, control median %.0f",
            nrow(A), ncol(A), median(colSums(det)[grp == "Case"]),
            median(colSums(det)[grp == "Control"])))
say(sprintf("growth on the colonic medium: %d of %d models grow (median %.4f/h)",
            sum(cap$grows == 1), nrow(cap), median(cap$growth_gut, na.rm = TRUE)))
say("FVA interval: a product whose minimum exceeds zero must carry flux at 10% of optimal growth")
for (p in names(LAB)) {
  mn <- cap[[paste0(p, "_min")]]; st <- cap[[paste0(p, "_status")]]
  say(sprintf("  %-11s obligatory in %d of %d solved models", p,
              sum(st == "ok" & !is.na(mn) & mn > 1e-6), sum(st == "ok")))
}

# A blank capacity is three different things and they must not be averaged
# together. `no_metabolite` means the cytosolic product is absent from the
# network, which is a real zero and belongs in the denominator. `no_growth`,
# `infeasible` and `solver_error` mean the number is unknown, so those genomes
# leave both the numerator and the denominator for that product rather than
# being silently counted as zero.
USABLE <- c("ok", "no_metabolite")
weighted <- function(p, keep = rep(TRUE, nrow(cap))) {
  st <- cap[[paste0(p, "_status")]]
  v  <- cap[[paste0(p, "_growth_coupled")]]
  v[st == "no_metabolite"] <- 0
  use <- st %in% USABLE & !is.na(v) & keep
  W <- A * det
  W[!use, ] <- 0
  cs <- colSums(W)
  W <- sweep(W, 2, ifelse(cs > 0, cs, 1), "/")
  list(value = as.numeric(crossprod(ifelse(is.na(v), 0, v), W)),
       n_used = sum(use), n_dropped = sum(!(st %in% USABLE)),
       donors_empty = sum(cs <= 0))
}

per_donor <- lapply(names(LAB), function(p) {
  if (!paste0(p, "_status") %in% names(cap)) return(NULL)
  w <- weighted(p)
  st <- cap[[paste0(p, "_status")]]
  v  <- cap[[paste0(p, "_growth_coupled")]]
  tibble(product = p, sample = meta$sample, group = grp, capacity = w$value,
         n_used = w$n_used, n_dropped = w$n_dropped,
         n_capable = sum(st == "ok" & !is.na(v) & v > 1e-6))
}) %>% bind_rows()

say("")
say("Model coverage per product (genomes contributing to the weighting)")
say(sprintf("  %-11s %7s %9s %9s %9s", "product", "used", "true zero", "unknown", "can form"))
for (p in names(LAB)) {
  st <- cap[[paste0(p, "_status")]]
  d <- per_donor[per_donor$product == p, ][1, ]
  say(sprintf("  %-11s %7d %9d %9d %9d", p, d$n_used, sum(st == "no_metabolite"),
              d$n_dropped, d$n_capable))
}
stopifnot(all(per_donor$n_used > 0))

cliff <- function(a, b) { m <- outer(a, b, "-"); (sum(m > 0) - sum(m < 0)) / length(m) }
tests <- per_donor %>% group_by(product, n_capable, n_used, n_dropped) %>%
  summarise(ctrl_med = median(capacity[group == "Control"]),
            case_med = median(capacity[group == "Case"]),
            delta = cliff(capacity[group == "Case"], capacity[group == "Control"]),
            p = wilcox.test(capacity ~ group, exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(q = p.adjust(p, "BH")) %>% arrange(delta)

say("")
say("Abundance-weighted formation capacity, mmol per gDW per hour, growth-coupled")
say(sprintf("%-11s %8s %10s %10s %7s %10s", "product", "models", "control", "case", "delta", "q"))
for (i in seq_len(nrow(tests))) with(tests[i, ],
  say(sprintf("%-11s %8d %10.4f %10.4f %+7.2f %10.3g  (%d models used, %d unknown)",
              product, n_capable, ctrl_med, case_med, delta, q, n_used, n_dropped)))

# ---- sensitivity: is the rise simply Enterobacteriaceae? ---------------------
# Every product that rises in patients is a mixed-acid fermentation product, which
# is the Enterobacteriaceae signature, and that family is the expansion this
# manuscript is about. Dropping those genomes from the weighting asks whether any
# of the pattern survives without them.
# The 584-genome eligibility table, not the 339-model inventory: the latter has
# no row for the 245 later-built genomes, so a match against it silently scored
# 3 of the 12 Enterobacteriaceae as non-family.
inv <- read.delim(file.path(N76, "metabolic_support/MAG_eligibility_all584.tsv"),
                  check.names = FALSE)
fam <- inv$family[match(cap$catalog_id, inv$catalog_id)]
ent <- !is.na(fam) & fam == "Enterobacteriaceae"
say("")
say(sprintf("Enterobacteriaceae among the modelled MAGs: %d of %d", sum(ent), length(ent)))

sens <- lapply(names(LAB), function(p) {
  if (!paste0(p, "_status") %in% names(cap)) return(NULL)
  x <- weighted(p, keep = !ent)$value
  tibble(product = p,
         ctrl_med = median(x[grp == "Control"]), case_med = median(x[grp == "Case"]),
         delta = cliff(x[grp == "Case"], x[grp == "Control"]),
         p = wilcox.test(x ~ grp, exact = FALSE)$p.value)
}) %>% bind_rows() %>% mutate(q = p.adjust(p, "BH"))
say("Same test with Enterobacteriaceae genomes removed from the weighting")
say(sprintf("%-11s %10s %10s %7s %10s", "product", "control", "case", "delta", "q"))
for (i in seq_len(nrow(sens))) with(sens[i, ],
  say(sprintf("%-11s %10.4f %10.4f %+7.2f %10.3g", product, ctrl_med, case_med, delta, q)))
write_csv(sens, file.path(OUTDIR, "Fig4f_flux_capacity_noEntero_n76.csv"))

ord <- tests$product
per_donor$product <- factor(per_donor$product, ord)
tests$product <- factor(tests$product, ord)
tests <- tests %>% mutate(
  label = sprintf("%s, %s%+.2f", ifelse(q < 1e-4, "q < 1e-4", sprintf("q = %.3g", q)),
                  "δ = ", delta))

# Precomputed, not evaluated inside aes(): the saved plot object is reloaded by
# fig4_square_panels_n76.R, where these locals would not exist.
XMAX  <- max(per_donor$capacity)
tests$xlab <- (sqrt(XMAX) * 1.04)^2

p4f <- ggplot(per_donor, aes(capacity, product, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = .62, linewidth = .22,
               position = position_dodge(width = .72), alpha = .85) +
  geom_point(aes(colour = group), position = position_jitterdodge(
    jitter.width = .16, dodge.width = .72, seed = 42), size = .5, alpha = .55,
    show.legend = FALSE) +
  geom_text(data = tests, aes(x = xlab, y = product, label = label), inherit.aes = FALSE,
            hjust = 0, size = 2.0, colour = "grey25") +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_colour_manual(values = GCOL) +
  scale_y_discrete(labels = LAB) +
  scale_x_sqrt(limits = c(0, (sqrt(XMAX) * 1.62)^2),
               expand = expansion(mult = c(.01, 0))) +
  coord_cartesian(clip = "off") +
  labs(title = "Predicted fermentation-product formation capacity",
       subtitle = "Per donor: abundance-weighted maximum formation flux on a western-diet colonic medium",
       x = "mmol per gDW per hour (square-root scale)", y = NULL,
       caption = paste("Gapfilled genome-scale models, growth held at 10% of optimum.",
                       "Predicted capacity under one stated diet, not a rate or a measurement.")) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        plot.caption = element_text(size = 5.5, colour = "grey35", hjust = 0),
        plot.title.position = "plot",
        legend.position = c(.86, .12), legend.background = element_blank(),
        legend.text = element_text(size = 7), legend.key.size = grid::unit(8, "pt"),
        plot.margin = margin(4, 6, 4, 4))

ggsave(file.path(OUTDIR, "Fig4f_FluxCapacity_n76.pdf"), p4f, width = 5.2, height = 4.4,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(OUTDIR, "Fig4f_FluxCapacity_n76.png"), p4f, width = 5.2, height = 4.4,
       units = "in", dpi = 600, bg = "white", type = "cairo-png")
saveRDS(p4f, file.path(OUTDIR, "Fig4f_FluxCapacity_n76.rds"))
write_csv(tests, file.path(OUTDIR, "Fig4f_flux_capacity_tests_n76.csv"))
write_csv(per_donor, file.path(OUTDIR, "Fig4f_flux_capacity_per_donor_n76.csv"))
writeLines(STATS, file.path(OUTDIR, "FIG4F_FLUX_CAPACITY_n76.txt"))
message("DONE -> Fig4f_FluxCapacity_n76.pdf")
