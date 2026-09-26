# =============================================================================
# Figure 1, rebuilt on the 76-donor catalogue (584 MAGs, 38 case / 38 control)
#
# Supersedes both the published n = 33 figure and the donor-collapsed n = 16
# figure. The control arm is now 38 independent donors with technical runs
# merged before assembly, so there is no pseudoreplication and no collapse step.
#
# Panels and encodings follow the manuscript Figure 1 legend exactly. All
# plotting helpers (add_stats, get_tax_df, build_mag_df, plot_ancombc2_volcano,
# tidy_ancombc2_res) are adapted from ~/MAG_Analysis/Mags_taxanomy.R. The cached
# primary ANCOM-BC2 calls are used directly; no model is refitted here.
#
# Panel c (ZOE health index) is CONDITIONAL: it is drawn only when the rebuilt
# per-donor ZOE table exists. The 22 PREDICT controls and the new H5 donor had
# no published ZOE score, and the 15 legacy donors were re-QC'd with their runs
# merged, so the whole control arm is being re-profiled with MetaPhlAn. Until
# that lands the figure is built as a 7-panel a-g.
#
# Outputs -> result2/n76/
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(dplyr); library(tidyr)
  library(tibble); library(readr); library(ggplot2); library(ggrepel)
  library(patchwork); library(colorspace); library(scales)
  library(forcats); library(ggnewscale); library(aplot); library(ggstar)
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
source(file.path(BASE, "taxon_italics_n76.R"))   # expr_from_md(): taxon names in italics
RES    <- file.path(BASE, "results")
outdir <- file.path(BASE, "result2/n76")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# inputs stay in result2/n76; every Figure 1 output goes to its own folder
figdir <- file.path(outdir, "fig1")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

STATS <- c()
say <- function(...) { line <- paste0(...); message(line); STATS <<- c(STATS, line) }

# Prevalence = share of donors in which a MAG exceeds PREV_DETECTION relative
# abundance. The published rule (>= 1 read) saturates at this depth: control
# donors now carry a median 52 M mapped reads, so one stray read at 95% identity
# made E. coli, E. ruysiae, K. quasipneumoniae and C. portucalensis "present" in
# every control. A 0.01% floor makes presence comparable across 21-57 M reads.
# Sensitivity (prevalence-LFC Spearman, significant MAGs): >=1 read 0.31;
# rarefied 0.43; 1e-5 0.75; 1e-4 0.79; 1e-3 0.81.
PREV_DETECTION <- 1e-4

# ---- palette: MUST be phylum_colors_mags.rds (GTDB names) -------------------
pal <- readRDS(file.path(RES, "phylum_colors_mags.rds"))
phylum_pal <- pal

# Helvetica lacks the arrow/triangle glyphs used in the count labels; DejaVu Sans has them.
ARROW_FONT <- "DejaVu Sans"



ps_mags <- readRDS(file.path(outdir, "ps_mags_n76.rds"))
res76   <- read.csv(file.path(outdir, "ANCOMBC2_MAG_results_n76.csv"), stringsAsFactors = FALSE)

# Rebuilt per-donor ZOE table; absent until MetaPhlAn on the 38 controls finishes.
ZOE_NEW <- file.path(outdir, "ZOE_HealthIndex_PerDonor_n76.csv")
HAVE_ZOE <- file.exists(ZOE_NEW)

# ---- per-panel save ---------------------------------------------------------
# Panel letters depend only on whether the ZOE panel exists, and that is known
# above, so each panel can be written under its final name the moment it is
# built instead of waiting for the loop at the end of the script.
PANEL_KEYS <- c("PCoA", "DiversityInstability_marginal",
                if (HAVE_ZOE) "ZOE_HealthIndex",
                "Volcano", "Pie", "PrevalenceVsLFC", "Phylum", "Family")
PANEL_LETTER <- setNames(letters[seq_along(PANEL_KEYS)], PANEL_KEYS)
save_panel <- function(p, key, w, h) {
  f <- file.path(figdir, sprintf("Fig1%s_%s_n76.pdf", PANEL_LETTER[[key]], key))
  try(ggsave(f, p, width = w, height = h, device = cairo_pdf, bg = "transparent"), silent = TRUE)
  invisible(p)
}


say("=== FIGURE 1, n = 76 donors (38 case / 38 control) ===")
say(sprintf("samples: %d  (%s)", nsamples(ps_mags),
            paste(names(table(sample_data(ps_mags)$Group)),
                  table(sample_data(ps_mags)$Group), sep = "=", collapse = ", ")))

# =============================================================================
# ADAPTED HELPERS from Mags_taxanomy.R
# =============================================================================
`%||%` <- function(a,b) if (!is.null(a)) a else b

tidy_ancombc2_res <- function(res) {
  suff <- grep("^lfc_", names(res), value = TRUE) |> sub(pattern = "^lfc_", replacement = "")
  bind_rows(lapply(suff, function(s) {
    pick <- function(pref) { nm <- paste0(pref, s); if (nm %in% names(res)) res[[nm]] else NA_real_ }
    tibble(taxon = res$taxon %||% rownames(res),
           contrast = s,
           lfc = pick("lfc_"), se = pick("se_"), W = pick("W_"),
           p = pick("p_"), q = pick("q_"),
           diff = as.logical(res[[paste0("diff_", s)]]) %||% NA,
           passed_ss = as.logical(res[[paste0("passed_ss_", s)]]) %||% NA)
  })) %>% arrange(q)
}

add_stats <- function(df, ps) {
  if (is.null(ps)) return(df)
  ot  <- phyloseq::otu_table(ps); mat <- as.matrix(ot)
  if (!phyloseq::taxa_are_rows(ps)) mat <- t(mat)
  mat_rel <- sweep(mat, 2, pmax(colSums(mat), .Machine$double.eps), "/")
  tibble(taxon = rownames(mat),
         mean_abund = rowMeans(mat_rel, na.rm = TRUE),
         # verbatim helper, except detection uses PREV_DETECTION (see top of file)
         prevalence = rowMeans(mat_rel > PREV_DETECTION, na.rm = TRUE)) %>%
    right_join(df, by = "taxon")
}

get_tax_df <- function(ps) {
  tax_raw <- as.data.frame(phyloseq::tax_table(ps)) %>% tibble::rownames_to_column("taxon")
  if ("Species" %in% colnames(tax_raw)) {
    tax_raw$label <- gsub("^s__", "", tax_raw$Species)
  } else if ("Genus" %in% colnames(tax_raw)) {
    tax_raw$label <- paste(gsub("^g__", "", tax_raw$Genus), "sp.")
  } else tax_raw$label <- tax_raw$taxon
  tax_raw
}

build_mag_df <- function(ab2_long, ps, pal, q_thresh = 0.05) {
  tax_df <- get_tax_df(ps)
  ab2_long %>% tibble::as_tibble() %>%
    left_join(select(tax_df, taxon, Phylum, label), by = "taxon") %>%
    mutate(Phylum = if_else(is.na(Phylum) | !nzchar(Phylum), "Unknown", Phylum)) %>%
    add_stats(ps) %>%
    mutate(q_plot = pmax(q, .Machine$double.xmin),
           sig    = !is.na(diff) & diff,
           shape_key = case_when(sig & lfc > 0 ~ "Up", sig & lfc < 0 ~ "Down", TRUE ~ "NS"),
           Phylum = if_else(Phylum %in% names(pal), Phylum, "Unknown"))
}

plot_ancombc2_volcano <- function(df, outdir, pal, q_line = 0.05, label_y = 7.5,
                                  xlim_fixed = NULL, prev_floor = NULL,
                                  prev_roof = NULL, prev_gamma = 0.5) {
  # symmetric about zero so the enriched and depleted sides are read on the
  # same scale rather than the axis following whichever side reaches further
  if (is.null(xlim_fixed)) xlim_fixed <- c(-1, 1) * (max(abs(df$lfc), na.rm = TRUE) + 0.5)
  pf <- if (is.null(prev_floor)) min(df$prevalence, na.rm = TRUE) else prev_floor
  pr <- if (is.null(prev_roof))  max(df$prevalence, na.rm = TRUE) else prev_roof
  if (pf >= pr) { pf <- 0; pr <- 1 }
  p_scaled <- pmin(pmax((df$prevalence - pf) / (pr - pf), 0), 1)^prev_gamma
  shade_by_prev <- function(base_hex, p) {
    p <- pmin(pmax(p, 0), 1); rgb_m <- col2rgb(base_hex) / 255
    rgb((1-p)*1 + p*rgb_m[1,], (1-p)*1 + p*rgb_m[2,], (1-p)*1 + p*rgb_m[3,])
  }
  base_col <- pal[as.character(df$Phylum)]
  df <- df %>% mutate(
    fill_col   = shade_by_prev(base_col, p_scaled),
    stroke_col = if_else(sig, scales::alpha(base_col, 1) %>% colorspace::darken(0.3), "grey20"))
  abbreviate_binom <- function(x) sapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    parts <- strsplit(s, "\\s+")[[1]]
    if (length(parts) < 2) return(s)
    paste0(substr(parts[1],1,1), ". ", paste(parts[-1], collapse=" "))
  })
  lab_df <- df %>% filter(sig, -log10(q_plot) > label_y) %>%
    mutate(label_abbr = abbreviate_binom(label))
  ggplot(df, aes(x = lfc, y = -log10(q_plot))) +
    coord_cartesian(xlim = xlim_fixed) +
    scale_x_continuous(breaks = seq(-floor(xlim_fixed[2]/2)*2, floor(xlim_fixed[2]/2)*2, by = 2)) +
    annotate("rect", xmin=-1, xmax=1, ymin=-Inf, ymax=Inf, fill="grey80", alpha=0.12) +
    geom_hline(yintercept = -log10(q_line), linetype="dotted", colour="grey35") +
    geom_hline(yintercept = label_y, colour="red", size=0.3, linetype="dotted") +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50") +
    geom_point(data = subset(df, !sig), aes(size = mean_abund), shape=21,
               fill="#D0D0D0", colour="#7F7F7F", stroke=0.2, alpha=0.5, show.legend=FALSE) +
    geom_point(data = subset(df, sig),
               aes(size=mean_abund, fill=I(fill_col), shape=shape_key, colour=I(stroke_col)),
               stroke = 0.25) +
    ggrepel::geom_text_repel(data = lab_df, aes(label = label_abbr), size=3,
                             fontface="italic", show.legend=FALSE, max.overlaps=Inf) +
    scale_shape_manual(name="Direction", values=c(Down=25, Up=24, NS=21)) +
    scale_size_continuous(name="Mean abundance", range=c(1,4)) +
    theme_bw(base_size = 12) +
    theme(legend.position="none",
          panel.background=element_rect(fill="transparent", color=NA),
          plot.background =element_rect(fill="transparent", color=NA)) +
    labs(title="MAG Differential Abundance", x="Log Fold Change", y="-Log10(q)")
}

theme_mag <- theme_bw(base_size = 11) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_line(color="grey92"),
        panel.grid.major.y=element_line(color="grey92"),
        axis.text=element_text(color="black"),
        plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(color="grey40", size=8),
        plot.caption=element_text(color="grey50", size=7),
        legend.position="right",
        legend.key=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent", color=NA),
        plot.background =element_rect(fill="transparent", color=NA))

# =============================================================================
# CORE TABLES
# =============================================================================
mag_ancom_long <- tidy_ancombc2_res(res76) %>% filter(!grepl("Intercept", contrast))
mag_df <- build_mag_df(mag_ancom_long, ps_mags, pal, q_thresh = 0.05)

tax_df <- get_tax_df(ps_mags)
sig_mags_with_tax <- res76 %>% filter(diff_GroupCase) %>% left_join(tax_df, by = "taxon")

say(sprintf("MAGs tested: %d | significant: %d (enriched %d, depleted %d)",
            nrow(res76), nrow(sig_mags_with_tax),
            sum(sig_mags_with_tax$lfc_GroupCase > 0),
            sum(sig_mags_with_tax$lfc_GroupCase < 0)))

# =============================================================================
# PANEL A(i): Aitchison PCoA biplot + PERMANOVA
# =============================================================================
otu_raw <- as(otu_table(ps_mags), "matrix")
if (taxa_are_rows(ps_mags)) otu_raw <- t(otu_raw)
sample_data(ps_mags)$Shannon  <- vegan::diversity(otu_raw, index = "shannon")
sample_data(ps_mags)$Observed <- rowSums(otu_raw > 0)

SHANNON_LIM <- range(sample_data(ps_mags)$Shannon, na.rm = TRUE) + c(-0.05, 0.05)
SHANNON_BRK <- pretty(SHANNON_LIM, n = 4)
RICH_BRK    <- pretty(range(sample_data(ps_mags)$Observed, na.rm = TRUE), n = 4)
say(sprintf("Shannon range %.2f-%.2f ; observed richness range %d-%d",
            SHANNON_LIM[1], SHANNON_LIM[2],
            min(sample_data(ps_mags)$Observed), max(sample_data(ps_mags)$Observed)))

ps_clr <- transform_sample_counts(ps_mags, function(x) { xp <- x + 1; log(xp/exp(mean(log(xp)))) })
ord_ait <- ordinate(ps_clr, method = "PCoA", distance = "euclidean")
evals <- ord_ait$values$Eigenvalues; var_ait <- evals/sum(evals)
labs_ait <- c(sprintf("PCoA 1 (%.1f%%)", 100*var_ait[1]), sprintf("PCoA 2 (%.1f%%)", 100*var_ait[2]))
df_ait <- plot_ordination(ps_clr, ord_ait, justDF = TRUE)
df_ait$Group <- get_variable(ps_clr, "Group")

# Arrows come from phylum relative abundance, computed from the counts. The
# earlier version divided the CLR table by its own row sums; CLR values sum to
# zero by construction (measured: 1e-12), so every sample was divided by its own
# rounding error. That corrupted both which phyla were picked (Bacillota and
# Bacteroidota, three quarters of every sample, could never appear) and the
# correlations that set the arrow directions. arrow_scale below is unchanged.
ps_rel_ait  <- transform_sample_counts(ps_mags, function(x) x/sum(x))
ps_glom_ait <- tax_glom(ps_rel_ait, taxrank = "Phylum")
top5 <- names(sort(taxa_sums(ps_glom_ait), decreasing = TRUE)[1:5])
ps_top5 <- prune_taxa(top5, ps_glom_ait)
otu_top <- as(otu_table(ps_top5), "matrix"); if (taxa_are_rows(ps_top5)) otu_top <- t(otu_top)
common <- intersect(rownames(otu_top), rownames(ord_ait$vectors))
vec_ait <- cor(otu_top[common,], ord_ait$vectors[common, 1:2])
arrow_scale <- 2 * max(abs(df_ait$Axis.1), abs(df_ait$Axis.2))
arrow_df <- as.data.frame(vec_ait * arrow_scale)
arrow_df$Label <- as.character(tax_table(ps_top5)[rownames(arrow_df), "Phylum"])

otu_mat <- as(otu_table(ps_clr), "matrix"); if (taxa_are_rows(ps_clr)) otu_mat <- t(otu_mat)
meta_df <- as(sample_data(ps_clr), "data.frame")
set.seed(42)
perm_fit <- adonis2(dist(otu_mat, method="euclidean") ~ meta_df[["Group"]], permutations = 999)
r2_val <- round(perm_fit$R2[1], 3); p_val <- perm_fit$`Pr(>F)`[1]
p_lab <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", p_val))
stats_label <- paste0("PERMANOVA: R2 = ", r2_val, " | ", p_lab)
say(sprintf("PANEL A(i) PERMANOVA: R2 = %.4f, F = %.3f, p = %s   [manuscript n=33: R2 = 0.156, p = 0.001]",
            perm_fit$R2[1], perm_fit$F[1], format(p_val)))
say(sprintf("PANEL A(i) axes: %s ; %s", labs_ait[1], labs_ait[2]))

# ---- do the three control studies separate on the axes panel a displays? ----
# The control arm is three studies, and they do differ from one another overall.
# What matters for the case-control claim is whether that difference lies in the
# plotted plane. It does not: the studies are indistinguishable on both axes, so
# their overall difference sits in later components.
set.seed(42)
ctrl_i <- meta_df[["Group"]] == "Control"
perm_ctrl <- adonis2(dist(otu_mat[ctrl_i, ], method = "euclidean") ~ meta_df[["source"]][ctrl_i],
                     permutations = 999)
say(sprintf("PANEL A(i) PERMANOVA among the 3 control studies: R2 = %.4f, F = %.3f, p = %s",
            perm_ctrl$R2[1], perm_ctrl$F[1], format(perm_ctrl$`Pr(>F)`[1])))
say(sprintf("  R2 expected by chance for this design (df model / df total) = %.4f; adjusted R2 = %.4f",
            perm_ctrl$Df[1] / perm_ctrl$Df[nrow(perm_ctrl)],
            1 - (1 - perm_ctrl$R2[1]) * perm_ctrl$Df[nrow(perm_ctrl)] /
              perm_ctrl$Df[nrow(perm_ctrl) - 1]))

ax <- df_ait[match(rownames(meta_df), rownames(df_ait)), c("Axis.1", "Axis.2")]
stopifnot(nrow(ax) == nrow(meta_df), !anyNA(ax))
# PCoA axis signs are arbitrary; orient so controls sit on the positive side of axis 1.
if (median(ax$Axis.1[ctrl_i]) < median(ax$Axis.1[!ctrl_i])) ax$Axis.1 <- -ax$Axis.1
src_ctrl <- factor(meta_df[["source"]][ctrl_i])
for (a in c("Axis.1", "Axis.2")) {
  kw  <- kruskal.test(ax[[a]][ctrl_i] ~ src_ctrl)
  eta <- summary(aov(ax[[a]][ctrl_i] ~ src_ctrl))[[1]]
  say(sprintf("  %s among controls: Kruskal-Wallis chi2 = %.3f, p = %.3f (study explains %.1f%% of control variation on this axis)",
              sub("Axis.", "PCoA ", a, fixed = TRUE), kw$statistic, kw$p.value,
              100 * eta[1, 2] / sum(eta[, 2])))
}
a1 <- ax$Axis.1
say(sprintf("  PCoA 1 range: controls %.1f to %.1f ; patients %.1f to %.1f",
            min(a1[ctrl_i]), max(a1[ctrl_i]), min(a1[!ctrl_i]), max(a1[!ctrl_i])))
say(sprintf("  patients falling inside the control PCoA 1 range: %d of %d",
            sum(a1[!ctrl_i] >= min(a1[ctrl_i]) & a1[!ctrl_i] <= max(a1[ctrl_i])), sum(!ctrl_i)))
for (s in levels(src_ctrl)) {
  own <- a1[ctrl_i][src_ctrl == s]; rest <- a1[ctrl_i][src_ctrl != s]
  say(sprintf("  %-11s: %d of %d samples inside the PCoA 1 range of the other control studies",
              s, sum(own >= min(rest) & own <= max(rest)), length(own)))
}
say("  PCoA 1 is fitted to these data to maximise variance, so the range statement describes the ordination, not an independent test.")

# Case is one cohort; the control arm is three independent studies, so the
# shape channel carries cohort rather than just Group. ggplot2 only offers five
# fillable shapes (21-25), two of which are triangles, so ggstar is used for a
# fillable pentagon. Fill stays mapped to Shannon, so every shape must be
# fillable. Ellipses remain per Group, not per cohort.
COHORT_LAB <- c(
  "PRJNA1457742" = "Case (PRJNA1457742)",
  "PRJEB39223"   = "Control: PREDICT / ZOE (PRJEB39223)",
  "PRJEB7949"    = "Control: Cardiff (PRJEB7949)",
  "PRJEB7331"    = "Control: Sanger (PRJEB7331)")
# ggstar starshape codes: 15 circle, 13 square, 12 rhombus, 5 regular pentagon
COHORT_SHAPE <- c(
  "Case (PRJNA1457742)"                 = 15,
  "Control: PREDICT / ZOE (PRJEB39223)" = 13,
  "Control: Cardiff (PRJEB7949)"        = 12,
  "Control: Sanger (PRJEB7331)"         = 5)

df_ait$Cohort <- factor(unname(COHORT_LAB[as.character(get_variable(ps_clr, "source"))]),
                        levels = names(COHORT_SHAPE))
stopifnot(!anyNA(df_ait$Cohort))

# ggplot draws points in row order, so the 22 PREDICT/ZOE squares were landing
# on top of the smaller cohorts. Draw squares first, then pentagons, then
# diamonds last. Factor levels are untouched, so the legend order is unchanged.
DRAW_ORDER <- c("Case (PRJNA1457742)",
                "Control: PREDICT / ZOE (PRJEB39223)",
                "Control: Sanger (PRJEB7331)",
                "Control: Cardiff (PRJEB7949)")
df_ait <- df_ait[order(match(as.character(df_ait$Cohort), DRAW_ORDER)), ]
say(sprintf("PANEL A cohorts: %s",
            paste(names(table(df_ait$Cohort)), table(df_ait$Cohort), sep = " n=", collapse = ", ")))

my_line_colors <- c("Control" = "#00BFC4", "Case" = "#F8766D")

pA1 <- ggplot(df_ait, aes(x = Axis.1, y = Axis.2)) +
  stat_ellipse(aes(color = Group), type="norm", level=0.95, linetype=2, linewidth=0.8, show.legend=FALSE) +
  geom_star(aes(fill = Shannon, size = Observed, starshape = Cohort), color="black", starstroke=0.3) +
  scale_color_manual(values = my_line_colors) +
  scale_starshape_manual(values = COHORT_SHAPE, name = "Cohort") +
  scale_fill_distiller(palette="Blues", name="Shannon\nAlpha Diversity",
                       limits=SHANNON_LIM, breaks=SHANNON_BRK, direction=1) +
  scale_size_continuous(range=c(2,6), name="Observed\nRichness", breaks=RICH_BRK) +
  geom_segment(data=arrow_df, aes(x=0,y=0,xend=Axis.1,yend=Axis.2),
               arrow=arrow(length=unit(0.15,"cm")), color="black", linewidth=0.6) +
  geom_text_repel(data=arrow_df, aes(x=Axis.1, y=Axis.2, label=Label), size=3.5,
                  fontface="bold.italic", color="black", bg.color="white", bg.r=0.15,
                  box.padding=0.6, point.padding=0.3, force=5, max.overlaps=Inf) +
  # without an explicit override the size legend inherits ggstar's star glyph
  guides(starshape = guide_legend(override.aes = list(size = 3.2, fill = "grey40"),
                                  order = 1, nrow = 2),
         fill = guide_colourbar(order = 2, barwidth = unit(3.2, "cm"), barheight = unit(0.3, "cm")),
         size = guide_legend(override.aes = list(starshape = 15, fill = "grey60"),
                             order = 3, nrow = 1)) +
  labs(title="Aitchison Biplot: Compositional Difference",
       subtitle=paste0(stats_label, "\nShape: cohort | Fill: Shannon | Size: Richness"),
       x=labs_ait[1], y=labs_ait[2]) +
  theme_bw(base_size=12) +
  # legends moved below so they stop compressing the plotting area horizontally;
  # aspect.ratio pins the panel square regardless of the canvas
  theme(aspect.ratio=1, legend.position="bottom", legend.box="vertical",
        legend.box.just="left", panel.grid.minor=element_blank(),
        legend.text=element_text(size=7.5), legend.title=element_text(size=8.5),
        legend.key.size=unit(0.4, "cm"), legend.margin=margin(1, 4, 1, 4),
        legend.box.spacing=unit(0.2, "cm")) +
  coord_cartesian(clip="off")

save_panel(pA1, "PCoA", 6, 6.7)

# =============================================================================
# PANEL A(ii) + B: Shannon, beta dispersion
# =============================================================================
otu2 <- as(otu_table(ps_mags), "matrix"); if (taxa_are_rows(ps_mags)) otu2 <- t(otu2)
clr_mat <- apply(otu2 + 1, 1, function(x) log(x) - mean(log(x)))
ps_clr2 <- ps_mags; otu_table(ps_clr2) <- otu_table(t(clr_mat), taxa_are_rows = FALSE)
dist_ait <- phyloseq::distance(ps_clr2, method = "euclidean")
sample_data(ps_clr2)$Group <- as.factor(sample_data(ps_clr2)$Group)
meta <- as(sample_data(ps_clr2), "data.frame")
dispersion <- betadisper(dist_ait, meta$Group)

ps_counts <- ps_mags; otu_table(ps_counts) <- otu_table(round(otu_table(ps_counts), 0),
                                                        taxa_are_rows = taxa_are_rows(ps_mags))
alpha_df <- estimate_richness(ps_counts, measures="Shannon") %>% rownames_to_column("Sample")
alpha_df$Group <- meta$Group[match(alpha_df$Sample, rownames(meta))]

df_plot <- data.frame(Sample = names(dispersion$distances),
                      Group = dispersion$group,
                      DistCent = dispersion$distances) %>%
  left_join(alpha_df, by = c("Sample","Group"))

set.seed(42)
perm_disp <- permutest(dispersion, permutations = 999)
w_disp <- with(df_plot, wilcox.test(DistCent[Group == "Case"],
                                    DistCent[Group == "Control"], exact = FALSE))
w_shan <- with(alpha_df, wilcox.test(Shannon[Group == "Case"],
                                    Shannon[Group == "Control"], exact = FALSE))
say(sprintf("PANEL A(ii) Shannon Wilcoxon: W = %.1f, p = %.4g", w_shan$statistic, w_shan$p.value))
say(sprintf("  Shannon median Control = %.3f, Case = %.3f",
            median(alpha_df$Shannon[alpha_df$Group=="Control"]),
            median(alpha_df$Shannon[alpha_df$Group=="Case"])))
say(sprintf("PANEL B(i) beta dispersion permutest: F = %.3f, p = %.4g",
            perm_disp$tab$F[1], perm_disp$tab$`Pr(>F)`[1]))
say(sprintf("  dispersion Wilcoxon: W = %.1f, p = %.4g", w_disp$statistic, w_disp$p.value))
m <- lm(DistCent ~ Shannon * Group, data = df_plot); m_anova <- anova(m)
say(sprintf("PANEL B(ii) interaction model adj R2 = %.3f; sequential ANOVA: Shannon p = %.4g; Group p = %.4g; interaction p = %.4g",
            summary(m)$adj.r.squared, m_anova["Shannon","Pr(>F)"],
            m_anova["Group","Pr(>F)"], m_anova["Shannon:Group","Pr(>F)"]))

my_colors <- c("Control" = "#00BFC4", "Case" = "#F8766D")

# PANEL B: one composite, main scatter with marginal BOXPLOTS.
# Layout skeleton reused verbatim from Mags_taxanomy.R lines 640-655
# (Figure_Abundance_Prevalence_Smooth.pdf), density margins swapped for
# boxplots. The right margin absorbs what used to be a separate b(i).
b_main <- ggplot(df_plot, aes(x=Shannon, y=DistCent, color=Group)) +
  geom_point(size=1.8, alpha=0.75) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.6) +
  scale_color_manual(values=my_colors) +
  labs(title="Diversity vs dispersion",
       subtitle=sprintf("adj R² = %.3f; interaction p = %.3g",
                        summary(m)$adj.r.squared, m_anova["Shannon:Group","Pr(>F)"]),
       x="Shannon Diversity (Alpha)", y="Distance to Centroid (Aitchison)") +
  # square plotting area; the marginal tracks ride on top of it
  theme_mag + theme(legend.position="none", aspect.ratio = 1)

# Marginal boxplots carry jittered points, matching the project's own style
# reference results/Figure-1b_Diversity_vs_Instability.pdf. Group on the
# cross-axis so the two boxes sit on separate tracks rather than overplotting,
# and the jitter is seeded through position_jitter (geom_jitter(seed=) is
# silently ignored by current ggplot2).
b_top <- ggplot(df_plot, aes(x = Shannon, y = Group, fill = Group)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA, linewidth=0.35, width=0.62) +
  geom_point(position = position_jitter(height = 0.17, width = 0, seed = 42),
             size=0.75, alpha=0.8, shape=18, colour="grey15") +
  scale_fill_manual(values=my_colors) +
  theme_void() + theme(legend.position="none")

b_right <- ggplot(df_plot, aes(x = Group, y = DistCent, fill = Group)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA, linewidth=0.35, width=0.62) +
  geom_point(position = position_jitter(width = 0.17, height = 0, seed = 42),
             size=0.75, alpha=0.8, shape=18, colour="grey15") +
  scale_fill_manual(values=my_colors) +
  theme_void() + theme(legend.position="none")

# The patchwork `design` skeleton from Mags_taxanomy.R:640-655 cannot align the
# margins to the panel region (it recycles a 2-element width vector over a
# 3-column design, and patchwork does not share axes across design cells). The
# same script's own Figure-1b (line ~7361) uses aplot::insert_top/insert_right,
# which DOES align marginal plots to the main axes. Using that instead, since
# "sharing axes with the scatter" is the requirement.
pB_aplot <- b_main |>
  aplot::insert_top(b_top,   height = 0.18) |>
  aplot::insert_right(b_right, width  = 0.18)
# aplot objects are not ggplots; convert so ggsave/patchwork can consume them
pB <- ggplotify::as.ggplot(pB_aplot)

save_panel(pB, "DiversityInstability_marginal", 6.5, 6.3)

# =============================================================================
# PANEL C: ZOE microbiome health index  (CONDITIONAL)
#
# Drawn only once ZOE_HealthIndex_PerDonor_n76.csv exists. The published ZOE
# table covers 38 cases and the 33 legacy control RUNS; it has no score for the
# 22 PREDICT controls or the new H5 donor, and its legacy control scores were
# computed on per-run reads that no longer correspond to the donor-level reads
# used here. Rather than mix the two, the whole control arm is re-profiled.
# =============================================================================
if (HAVE_ZOE) {
  cliffs_delta <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))

  zoe_col <- read.csv(ZOE_NEW, stringsAsFactors = FALSE)
  stopifnot(all(c("donor", "health_index", "cohort") %in% names(zoe_col)))

  ci_ctrl  <- zoe_col$health_index[zoe_col$cohort == "Control"]
  ci_covid <- zoe_col$health_index[zoe_col$cohort == "COVID"]
  wt_new <- wilcox.test(ci_covid, ci_ctrl, exact = FALSE)
  cd_new <- cliffs_delta(ci_covid, ci_ctrl)
  say(sprintf("PANEL C ZOE (n=%d ctrl vs %d case): Wilcoxon p = %.3g, Cliff's delta = %.2f",
              length(ci_ctrl), length(ci_covid), wt_new$p.value, cd_new))
  say(sprintf("  ZOE medians: Control %.3f (range %.3f to %.3f) ; COVID %.3f",
              median(ci_ctrl), min(ci_ctrl), max(ci_ctrl), median(ci_covid)))

  # Canonical Figure 1 colours (Mags_taxanomy.R:579); ZOE "COVID" maps to Case.
  zoe_col$GroupF <- factor(ifelse(zoe_col$cohort == "COVID", "Case", "Control"),
                           levels = c("Control", "Case"))
  p_lab_zoe <- if (wt_new$p.value < 0.001) "p < 0.001" else sprintf("p = %.3g", wt_new$p.value)
  y_top_z <- max(zoe_col$health_index)

  pZ <- ggplot(zoe_col, aes(x = GroupF, y = health_index, fill = GroupF, color = GroupF)) +
    geom_violin(alpha = 0.22, width = 0.55, linewidth = 0.3) +
    geom_boxplot(width = 0.16, alpha = 0.55, outlier.shape = NA, linewidth = 0.4) +
    # position_jitter(seed=) not geom_jitter(seed=): the latter is silently ignored
    geom_point(position = position_jitter(width = 0.07, seed = 42), size = 1.5, alpha = 0.8) +
    scale_fill_manual(values = my_colors, guide = "none") +
    scale_color_manual(values = my_colors, guide = "none") +
    scale_x_discrete(expand = expansion(add = 0.55)) +
    coord_cartesian(clip = "off") +
    annotate("segment", x = 1, xend = 2, y = y_top_z + 0.06, yend = y_top_z + 0.06,
             color = "grey40") +
    annotate("text", x = 1.5, y = y_top_z + 0.11,
             label = sprintf("%s,  Cliff's delta = %.2f", p_lab_zoe, cd_new), size = 3.2) +
    labs(title = "Microbiome health index (ZOE)",
         subtitle = sprintf("%d control / %d case donors", length(ci_ctrl), length(ci_covid)),
         x = NULL, y = "Microbiome health index (ZOE)") +
    theme_mag + theme(legend.position = "none", plot.margin = margin(18, 6, 6, 6))
} else {
  pZ <- NULL
  say("PANEL C ZOE: SKIPPED - ZOE_HealthIndex_PerDonor_n76.csv not present yet.")
  say("  Figure built as 7 panels (a-g). Re-run once MetaPhlAn on the 38 controls completes.")
}

# =============================================================================
# PANEL C: volcano
# =============================================================================
if (!is.null(pZ)) save_panel(pZ, "ZOE_HealthIndex", 4.6, 5.2)

pC <- plot_ancombc2_volcano(df = mag_df, outdir = outdir, pal = pal)

save_panel(pC, "Volcano", 6, 6)

# =============================================================================
# PANEL D: pie
# =============================================================================
df_summary <- mag_df %>%
  mutate(Status = case_when(shape_key=="Up" ~ "Enriched", shape_key=="Down" ~ "Depleted", TRUE ~ "NS"),
         Status = factor(Status, levels=c("Enriched","Depleted","NS")))
pie_data <- df_summary %>% group_by(Status) %>% tally() %>% mutate(perc = n/sum(n)*100)
say(sprintf("PANEL D pie: %s", paste(pie_data$Status, pie_data$n,
                                     sprintf("(%.1f%%)", pie_data$perc), collapse="; ")))

pD <- ggplot(pie_data, aes(x="", y=n, fill=Status)) +
  geom_bar(stat="identity", width=1, color="white") + coord_polar("y", start=0) +
  scale_fill_manual(values=c(Enriched="grey25", Depleted="grey60", NS="grey92")) +
  theme_void() +
  labs(title="MAG Significance Distribution", subtitle=paste0("Total MAGs: ", nrow(df_summary))) +
  geom_text(aes(label=paste0(Status,"\n",n), colour=Status), position=position_stack(vjust=0.5),
            size=4, fontface="bold", show.legend=FALSE) +
  scale_colour_manual(values=c(Enriched="white", Depleted="black", NS="black")) +
  theme(plot.title=element_text(hjust=0.5, face="bold"), plot.subtitle=element_text(hjust=0.5))

save_panel(pD, "Pie", 5, 5)

# =============================================================================
# PANEL E (lettered f): prevalence vs effect size, significant MAGs
#
# Transcribed from Mags_taxanomy.R:2935-3082 ("Figure1F_Prevalence_Effect").
# Axes follow the ORIGINAL: x = log fold change, y = prevalence difference,
# with the red "core" and orange "opportunistic" zones drawn on the y axis.
# An earlier adaptation had these axes swapped, which stacked every MAG with a
# prevalence difference of exactly zero into one vertical column at x = 0.
# Prevalence uses PREV_DETECTION (see top of file), not the original > 0.
# =============================================================================
ps_rel_e  <- transform_sample_counts(ps_mags, function(x) x / sum(x))
otu_rel_e <- as.matrix(otu_table(ps_rel_e))
if (!taxa_are_rows(ps_rel_e)) otu_rel_e <- t(otu_rel_e)

meta_e <- sample_data(ps_rel_e)
case_samples <- rownames(meta_e)[meta_e$Group == "Case"]
ctrl_samples <- rownames(meta_e)[meta_e$Group == "Control"]

prev_df_e <- data.frame(
  taxon      = rownames(otu_rel_e),
  Prev_Case  = rowMeans(otu_rel_e[, case_samples] > PREV_DETECTION),
  Prev_Ctrl  = rowMeans(otu_rel_e[, ctrl_samples] > PREV_DETECTION),
  Mean_Abund = rowMeans(otu_rel_e) * 100
)

prev_plot_df <- sig_mags_with_tax %>%
  left_join(prev_df_e, by = "taxon") %>%
  mutate(Prev_Diff = Prev_Case - Prev_Ctrl)

n_points <- nrow(prev_plot_df)
ct <- suppressWarnings(cor.test(prev_plot_df$lfc_GroupCase, prev_plot_df$Prev_Diff,
                                method = "spearman", use = "pairwise.complete.obs"))
cor_label <- paste0("N = ", n_points,
                    "\nSpearman R = ", sprintf("%.2f", ct$estimate),
                    "\np = ", formatC(ct$p.value, format = "e", digits = 1))
say(sprintf("PANEL F Spearman: rho = %.3f, p = %.3g, n = %d   [manuscript n=33: rho = 0.65, p = 2.8e-29, n = 229]",
            ct$estimate, ct$p.value, n_points))
say(sprintf("  prevalence difference exactly 0: %d of %d (present in every donor of both groups: %d)",
            sum(prev_plot_df$Prev_Diff == 0), n_points,
            sum(prev_plot_df$Prev_Case == 1 & prev_plot_df$Prev_Ctrl == 1)))
n_core <- sum(prev_plot_df$lfc_GroupCase > 0 & prev_plot_df$Prev_Diff > 0.2)
n_opp  <- sum(prev_plot_df$lfc_GroupCase > 0 & prev_plot_df$Prev_Diff > -0.1 & prev_plot_df$Prev_Diff <= 0.2)
say(sprintf("  red zone (LFC>0, prev diff>0.2): %d ; orange zone (-0.1 to 0.2): %d", n_core, n_opp))

x_min <- min(prev_plot_df$lfc_GroupCase, na.rm = TRUE)
x_max <- max(prev_plot_df$lfc_GroupCase, na.rm = TRUE)
y_min <- min(prev_plot_df$Prev_Diff,     na.rm = TRUE)
y_max <- max(prev_plot_df$Prev_Diff,     na.rm = TRUE)
x_pos <- x_min + 0.05 * (x_max - x_min)
y_pos <- y_max - 0.05 * (y_max - y_min)

SHOW_LABELS <- FALSE

pE <- ggplot(prev_plot_df, aes(x = lfc_GroupCase, y = Prev_Diff)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0.2,  ymax = 1,
           fill = "#E41A1C", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -0.1, ymax = 0.2,
           fill = "#FF7F00", alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotdash",
              linewidth = 0.8, alpha = 1) +
  geom_point(aes(size = Mean_Abund, fill = Phylum),
             shape = 21, colour = "grey30", stroke = 0.35, alpha = 0.9) +
  annotate("text", x = x_pos, y = y_pos, label = cor_label,
           hjust = 0, vjust = 1, size = 3, fontface = "bold") +
  scale_fill_manual(values = pal, name = "Phylum (GTDB)",
                    labels = function(x) expr_from_md(md_taxon(x))) +
  scale_size_continuous(range = c(2, 6), name = "Mean \nrelative abundance (%)",
                        labels = scales::label_number(accuracy = 0.01)) +
  labs(title = "Prevalence vs. Effect Size",
       subtitle = "Red: prevalence difference >0.2 | Orange: prevalence difference -0.1 to 0.2",
       x = "Log Fold Change (Effect Size)",
       y = "Prevalence Difference (Case - Control)") +
  theme_bw(base_size = 9) +
  theme(aspect.ratio = 1, legend.position  = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA))

if (SHOW_LABELS) {
  pE <- pE + geom_text_repel(
    data = head(arrange(prev_plot_df, desc(abs(lfc_GroupCase))), 30),
    aes(label = sub("^(\\S)\\S+ ", "\\1. ", gsub("_", " ", Species))),
    size = 2, fontface = "italic", box.padding = 0.5,
    max.overlaps = Inf, min.segment.length = 0)
}

save_panel(pE, "PrevalenceVsLFC", 6.5, 5.5)

# =============================================================================
# PANELS F-G: clade direction summaries; Fisher conditional on primary DA MAGs
# =============================================================================
source(file.path(BASE, "clade_fisher_figure_n76.R"), local = TRUE)
clade_figs <- build_clade_fisher_panels(BASE, save = TRUE)
pF <- clade_figs$p_phylum
pG <- clade_figs$p_family
pG_inset <- clade_figs$p_inset
for (line in clade_figs$stats) say(line)

# =============================================================================
# SAVE PANELS + COMPOSED FIGURE
#
# Panel letters are assigned from the list actually built, so dropping the
# conditional ZOE panel re-letters d-h down to c-g instead of leaving a gap.
# =============================================================================
panels <- list(Fig1_PCoA = pA1, Fig1_DiversityInstability_marginal = pB)
dims   <- list(c(6, 6.7), c(6.5, 6.3))   # PCoA: square panel plus its bottom legend; b: square scatter
if (HAVE_ZOE) { panels$Fig1_ZOE_HealthIndex <- pZ; dims <- c(dims, list(c(4.6, 5.2))) }
panels <- c(panels, list(Fig1_Volcano = pC, Fig1_Pie = pD,
                         Fig1_PrevalenceVsLFC = pE, Fig1_Phylum = pF, Fig1_Family = pG))
dims   <- c(dims, list(c(6, 6), c(5, 5), c(6.5, 5.5), c(5.0, 4.2), c(11.5, 6)))   # g: width carries the stats column; height tightened so the 13 phylum rows sit close

letters_used <- letters[seq_along(panels)]
names(panels) <- paste0("Fig1", letters_used, sub("^Fig1", "", names(panels)))

for (i in seq_along(panels)) {
  ggsave(file.path(figdir, paste0(names(panels)[i], "_n76.pdf")), panels[[i]],
         width = dims[[i]][1], height = dims[[i]][2], device = cairo_pdf, bg = "transparent")
}
say(sprintf("panels written: %s", paste(names(panels), collapse = ", ")))

source(file.path(BASE, "figure1_compose_n76.R"), local = TRUE)
fig1 <- compose_figure1_n76()
ggsave(file.path(figdir, "Figure_1_n76.pdf"), fig1, width = 148.1667,
       height = 222, units = "mm", device = cairo_pdf, bg = "transparent")
say("Combined figure: 148.17 x 222 mm; designed at the DOCX insertion width; text 5-7 pt.")
say("Family panel h uses two consecutive ordered strips; all 584 MAGs retained. Full test descriptions are in the figure legend and standalone panels.")

writeLines(STATS, file.path(figdir, "Figure_1_n76_STATS.txt"))
message("\nDONE -> ", file.path(figdir, "Figure_1_n76.pdf"))
