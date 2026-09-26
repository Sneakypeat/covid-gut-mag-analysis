#!/usr/bin/env Rscript
# Supplementary figure: what the flux analysis rests on.
#
# Panel f of Figure 4 reports one number per donor. This figure shows the three
# things a reader has to accept before that number means anything: that the
# zero-flux result in the draft models was a transport gap rather than a missing
# pathway, that gapfilling was a light touch rather than a rebuild, and that the
# model-predicted capacity agrees with the independent gene call it is meant to
# supersede.
#
# Outputs -> result2/n76/supp_fig/flux/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
GF     <- file.path(N76, "metabolic_support/gapfill_v1")
OUTDIR <- file.path(N76, "supp_fig/sources/flux")   # working outputs; the final
# numbered PDF is copied to supp_fig/Supplementary_Figure_7_n76.pdf by hand
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS  <- character(); say <- function(s) { message(s); STATS <<- c(STATS, s) }
PT <- 2.845276
GCOL <- c(Control = "#00BFC4", Case = "#F8766D")

# Full 584-model FVA table; `_max` is the growth-coupled maximum.
cap   <- read_csv(file.path(GF, "fva_584.csv"), show_col_types = FALSE) %>%
  rename_with(~ sub("_max$", "_growth_coupled", .x), ends_with("_max"))
route <- read_csv(file.path(GF, "butyrate_route_audit.csv"), show_col_types = FALSE)
addn  <- read_csv(file.path(GF, "gapfill_footprint.csv"), show_col_types = FALSE)
wide  <- read_tsv(file.path(N76, "metabolic_support/gut_traits_v2/genome_traits_wide.tsv"),
                  show_col_types = FALSE)

# ---- a. the transport gap ----------------------------------------------------
pa <- ggplot(route, aes(reorder(reaction, n), n, fill = role)) +
  geom_col(width = .68) +
  geom_text(aes(label = n), hjust = -0.3, size = 5.5 / PT, colour = "grey20") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(enzyme = "#2E6E8E", transport = "#C2453B"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  # Rebuilt over the completed 584-model catalogue the panel no longer says
  # "the export step is simply absent": 18 of the 47 do carry an exchange.
  # What it does show is that the route is incomplete in most of them, and
  # that 14 lose ptb/buk between the gene call and the reconstruction.
  labs(title = "The butyrate route is incomplete in the models that carry the genes",
       subtitle = sprintf("Butyrate-route reactions among the %d modelled genomes DRAM calls ptb/buk-positive",
                          route$n_genomes[1]),
       x = NULL, y = "Genomes carrying the reaction") +
  theme_classic(base_size = 7)

# ---- b. gapfill footprint ----------------------------------------------------
pb <- ggplot(addn, aes(added)) +
  geom_histogram(binwidth = 1, fill = "#4C6E7A", colour = "white", linewidth = .15) +
  labs(title = "Gapfilling was a light touch",
       subtitle = sprintf("Reactions added per model to grow on the colonic medium (median %d, max %d)",
                          median(addn$added), max(addn$added)),
       x = "Reactions added", y = "Models") +
  theme_classic(base_size = 7)

# ---- c. predicted capacity against the independent gene call -----------------
gene <- wide %>% select(catalog_id, butyrate_terminal_pair, acetate_terminal_pair)
cc <- cap %>% inner_join(gene, "catalog_id") %>%
  transmute(catalog_id,
            butyrate = butyrate_growth_coupled, acetate = acetate_growth_coupled,
            but_gene = ifelse(butyrate_terminal_pair == "positive", "ptb/buk called", "not called"),
            ac_gene  = ifelse(acetate_terminal_pair == "positive", "pta/ackA called", "not called"))
# Only genomes whose optimisation actually solved enter the comparison against
# the gene call; a `no_growth` model would otherwise read as a confident zero.
cc$but_ok <- cap$butyrate_status[match(cc$catalog_id, cap$catalog_id)] %in% c("ok", "no_metabolite")
cc$ac_ok  <- cap$acetate_status[match(cc$catalog_id, cap$catalog_id)] %in% c("ok", "no_metabolite")
cl <- bind_rows(
  cc %>% filter(but_ok) %>% transmute(product = "Butyrate", capacity = butyrate, call = but_gene),
  cc %>% filter(ac_ok)  %>% transmute(product = "Acetate",  capacity = acetate,  call = ac_gene)) %>%
  filter(!is.na(capacity))
# `call` is relabelled to a shared two-level factor so each facet shows only its
# own two boxes: as a four-level factor both facets drew an empty category.
cl <- cl %>% mutate(called = factor(ifelse(grepl("called$", call) & call != "not called",
                                           "terminal pair called", "not called"),
                                    c("not called", "terminal pair called")))
say("Growth-coupled formation capacity by gene call, median mmol per gDW per hour")
for (pr in unique(cl$product)) for (lv in levels(cl$called)) {
  x <- cl$capacity[cl$product == pr & cl$called == lv]
  say(sprintf("  %-9s %-22s n = %3d  median %.3f", pr, lv, length(x), median(x)))
}
for (pr in unique(cl$product)) {
  d <- cl[cl$product == pr, ]
  w <- wilcox.test(capacity ~ called, data = d, exact = FALSE)
  say(sprintf("  %-9s Wilcoxon p = %.3g", pr, w$p.value))
}

pc <- ggplot(cl, aes(called, capacity, fill = called)) +
  geom_boxplot(outlier.shape = NA, width = .55, linewidth = .22, alpha = .85) +
  geom_jitter(width = .14, size = .5, alpha = .5, colour = "grey25") +
  facet_wrap(~product, scales = "free_y") +
  scale_fill_manual(values = c("not called" = "#BDBDBD",
                               "terminal pair called" = "#2E6E8E"), guide = "none") +
  labs(title = "Model-predicted capacity against the independent gene call",
       subtitle = "Growth-coupled maximum formation flux per genome, by whether DRAM called the terminal pair",
       x = NULL, y = "mmol per gDW per hour") +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7))

# ---- d. who carries the predicted capacity ----------------------------------
# The positive form of the Enterobacteriaceae test. Rather than removing the
# family and reporting that the difference goes away, this asks what share of the
# community's whole predicted capacity those genomes hold, which is the mechanism
# the manuscript argues for rather than a nuisance term.
ab   <- read.delim(file.path(INDIR, "MAG_relative_abundance_percent.tsv"),
                   check.names = FALSE, row.names = 1)
meta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE)
inv  <- read.delim(file.path(N76, "metabolic_support/MAG_eligibility_all584.tsv"), check.names = FALSE)
A    <- as.matrix(ab[cap$catalog_id, meta$sample])
det  <- A > 0.01
grp  <- factor(meta$group, levels = c("Control", "Case"))
ent  <- !is.na(inv$family[match(cap$catalog_id, inv$catalog_id)]) &
        inv$family[match(cap$catalog_id, inv$catalog_id)] == "Enterobacteriaceae"
PRODS <- c(succinate = "Succinate", propionate = "Propionate", formate = "Formate",
           ethanol = "Ethanol", acetate = "Acetate", lactate_D = "D-lactate",
           lactate_L = "L-lactate", butyrate = "Butyrate")
# Same status handling as panel f: `no_metabolite` is a real zero and stays in,
# `no_growth` and solver failures are unknown and leave the calculation.
attrib <- lapply(names(PRODS), function(pr) {
  st <- cap[[paste0(pr, "_status")]]
  v  <- cap[[paste0(pr, "_growth_coupled")]]
  v[st == "no_metabolite"] <- 0
  use <- st %in% c("ok", "no_metabolite") & !is.na(v)
  W <- A * det; W[!use, ] <- 0
  tot <- as.numeric(crossprod(ifelse(is.na(v), 0, v), W))
  We <- W; We[!ent, ] <- 0
  tibble(product = PRODS[[pr]], sample = meta$sample, group = grp,
         # A donor whose whole community capacity is numerically zero has no
         # meaningful share: one genome at 1e-37 would otherwise read as 100%.
         share = 100 * ifelse(tot > 1e-6,
                              as.numeric(crossprod(ifelse(is.na(v), 0, v), We)) / tot, 0))
}) %>% bind_rows() %>% mutate(product = factor(product, PRODS))
say("")
say("Enterobacteriaceae share of the community's predicted capacity, median % per donor")
for (pr in levels(attrib$product)) {
  d <- attrib[attrib$product == pr, ]
  say(sprintf("  %-11s control %5.1f%%  case %5.1f%%  case donors over 25%%: %d/38",
              pr, median(d$share[d$group == "Control"]), median(d$share[d$group == "Case"]),
              sum(d$share[d$group == "Case"] > 25)))
}
pd <- ggplot(attrib, aes(share, product, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = .62, linewidth = .22,
               position = position_dodge(width = .72), alpha = .85) +
  geom_point(aes(colour = group), position = position_jitterdodge(
    jitter.width = .15, dodge.width = .72, seed = 42), size = .5, alpha = .55,
    show.legend = FALSE) +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_colour_manual(values = GCOL) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(.01, .02))) +
  labs(title = sprintf("%d genomes carry about half of what patients are predicted to ferment", sum(ent)),
       subtitle = "Enterobacteriaceae share of each donor's total predicted formation capacity",
       x = "Share of the donor's predicted capacity (%)", y = NULL) +
  theme_classic(base_size = 7) +
  theme(legend.position = c(.9, .9), legend.background = element_blank(),
        legend.key.size = grid::unit(7, "pt"))
write_csv(attrib, file.path(OUTDIR, "entero_attribution_per_donor_n76.csv"))

# ---- e. the absolute scale --------------------------------------------------
# Figure 4f now reports each product's share of the capacity a community
# carries, because the three stages it compares have different units and only
# composition is comparable between them. That leaves the actual mmol values
# with nowhere to appear, so the panel they used to occupy is kept here.
# Its own caption repeated what this figure's legend already says and was being
# clipped at half width; the left margin leaves room for the "e" tag, which
# otherwise printed on top of the title.
pe <- wrap_elements(full = readRDS(file.path(N76, "fig4/Fig4f_FluxCapacity_n76.rds")) +
                      labs(caption = NULL) +
                      theme(plot.margin = margin(4, 8, 4, 24))) +
  theme(plot.tag.position = c(0.005, 0.995))

# Five panels in one column made a 7 x 15 in strip, which stretches every panel
# into a letterbox. Paired into rows instead: the two small diagnostics side by
# side, the two-facet comparison across the full width, then the two
# eight-product boxplots side by side.
fig <- ((pa | pb) / pc / (pd | pe)) +
  plot_layout(heights = c(1, 1.05, 1.45)) +
  plot_annotation(tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 10)))

ggsave(file.path(OUTDIR, "Supplementary_Figure_Flux_n76.pdf"), fig, width = 10.5, height = 11.5,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(OUTDIR, "Supplementary_Figure_Flux_n76.png"), fig, width = 10.5, height = 11.5,
       units = "in", dpi = 400, bg = "white", type = "cairo-png")
writeLines(STATS, file.path(OUTDIR, "SUPP_FLUX_n76.txt"))
message("DONE -> Supplementary_Figure_Flux_n76.pdf")
