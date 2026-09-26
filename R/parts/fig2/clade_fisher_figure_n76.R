# Rebuild only Figure 1g, 1h and the family-count inset. No ANCOM-BC2 refit.
# Run: Rscript clade_fisher_figure_n76.R
# Also sourced by both Figure 1 entry points; sourcing does not execute it.
#
# Estimand: enriched/depleted composition among the 290 primary DA MAGs.
# Each two-sided Fisher table is [clade up, clade down; other up, other down].
# Eligibility: >=4 catalogue MAGs and >=1 DA MAG, fixed for both ranks.
# BH is applied within rank over eligible tests only. All 584 MAGs are plotted.
# These are exploratory catalogue-level summaries, not clade abundance tests.


build_clade_fisher_panels <- function(
    base = Sys.getenv("COVID_MAG_BASE", unset = getwd()),
    save = TRUE) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
  })
  n76 <- file.path(base, "result2/n76")
  inputs <- file.path(base, "result2/inputfiles/new")
  figdir <- file.path(n76, "fig1")
  tabledir <- file.path(n76, "taxon_level")
  anc_file <- file.path(n76, "ANCOMBC2_MAG_results_n76.csv")
  tax_file <- file.path(inputs, "MAG_quality_taxonomy.tsv")
  count_file <- file.path(inputs, "MAG_read_count.tsv")
  palette_file <- file.path(base, "results/phylum_colors_mags.rds")
  anc <- read_csv(anc_file, show_col_types = FALSE)
  tax <- read.delim(tax_file, check.names = FALSE, quote = "", comment.char = "")
  counts <- read.delim(count_file, check.names = FALSE, quote = "", comment.char = "")
  stopifnot(!anyDuplicated(anc$taxon), !anyDuplicated(tax$catalog_id),
            !anyDuplicated(counts[[1]]), setequal(anc$taxon, tax$catalog_id),
            setequal(anc$taxon, counts[[1]]), nrow(anc) == 584L)
  rank_part <- function(x, prefix) vapply(strsplit(x, ";", fixed = TRUE), function(z) {
    hit <- z[startsWith(z, prefix)]
    if (length(hit) != 1L) return(NA_character_)
    sub(prefix, "", hit, fixed = TRUE)
  }, character(1))
  taxonomy <- tibble(taxon = tax$catalog_id,
                    Phylum = rank_part(tax$gtdbtk_classification, "p__"),
                    Family = rank_part(tax$gtdbtk_classification, "f__"))
  count_matrix <- as.matrix(counts[, -1, drop = FALSE])
  storage.mode(count_matrix) <- "numeric"
  stopifnot(ncol(count_matrix) == 76L, all(is.finite(count_matrix)),
            all(count_matrix >= 0))
  prevalence <- tibble(taxon = counts[[1]],
                       Prevalence = rowMeans(count_matrix > 0))
  dat <- anc %>%
    transmute(taxon, lfc = lfc_GroupCase, mag_q = q_GroupCase,
              sig = diff_GroupCase %in% TRUE) %>%
    left_join(taxonomy, by = "taxon") %>%
    left_join(prevalence, by = "taxon") %>%
    mutate(Family = ifelse(is.na(Family) | Family == "", "Unassigned", Family),
           Status = case_when(sig & lfc > 0 ~ "Up", sig & lfc < 0 ~ "Down",
                              TRUE ~ "NS"))
  stopifnot(nrow(dat) == 584L, all(is.finite(dat$lfc)),
            !anyNA(dat$Phylum), all(nzchar(dat$Phylum)),
            !anyNA(dat$Prevalence), !any(dat$sig & dat$lfc == 0),
            sum(dat$Status == "Up") == 87L,
            sum(dat$Status == "Down") == 203L,
            sum(dat$Status == "NS") == 294L)
  # Keep the published GTDB palette; missing taxa fail instead of becoming Other.
  pal <- readRDS(palette_file)
  names(pal) <- trimws(names(pal))
  if (!all(unique(dat$Phylum) %in% names(pal)))
    stop("Unmapped phylum: ", paste(setdiff(unique(dat$Phylum), names(pal)), collapse = ", "))
  up_total <- sum(dat$Status == "Up")
  down_total <- sum(dat$Status == "Down")

  summarise_rank <- function(rank) {
    s <- dat %>% group_by(clade = .data[[rank]]) %>%
      summarise(n_total = n(), n_up = sum(Status == "Up"),
                n_down = sum(Status == "Down"), n_ns = sum(Status == "NS"),
                mean_lfc = mean(lfc), min_lfc = min(lfc), max_lfc = max(lfc),
                Phylum = first(Phylum), .groups = "drop") %>%
      mutate(n_da = n_up + n_down, other_up = up_total - n_up,
             other_down = down_total - n_down,
             eligible = n_total >= 4 & n_da > 0 & clade != "Unassigned",
             odds_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
             p = NA_real_, q = NA_real_)
    for (i in which(s$eligible)) {
      tab <- matrix(c(s$n_up[i], s$n_down[i], s$other_up[i], s$other_down[i]),
                    nrow = 2, byrow = TRUE,
                    dimnames = list(c("clade", "other"), c("up", "down")))
      stopifnot(sum(tab) == 290L, all(tab >= 0), all(rowSums(tab) > 0))
      f <- stats::fisher.test(tab, alternative = "two.sided", conf.level = 0.95)
      s$odds_ratio[i] <- unname(f$estimate)
      s$ci_lower[i] <- f$conf.int[1]
      s$ci_upper[i] <- f$conf.int[2]
      s$p[i] <- f$p.value
    }
    s$q[s$eligible] <- p.adjust(s$p[s$eligible], method = "BH")
    s %>% mutate(log2_odds_ratio = log2(odds_ratio),
                 bias = case_when(!eligible ~ "Not tested",
                                  q < 0.05 & odds_ratio > 1 ~ "Enrichment bias",
                                  q < 0.05 & odds_ratio < 1 ~ "Depletion bias",
                                  TRUE ~ "No detected bias")) %>% arrange(mean_lfc, clade)
  }
  phy <- summarise_rank("Phylum")
  fam <- summarise_rank("Family")
  stopifnot(sum(phy$n_total) == 584L, sum(fam$n_total) == 584L,
            sum(phy$eligible) == 8L, sum(fam$eligible) == 33L,
            all(phy$n_total == phy$n_up + phy$n_down + phy$n_ns),
            all(fam$n_total == fam$n_up + fam$n_down + fam$n_ns))
  fmt_q <- function(x) ifelse(is.na(x), "NT", trimws(formatC(x, format = "g", digits = 2)))
  italic_taxa <- function(x) as.expression(lapply(x, function(z) {
    if (z == "Unassigned") z else bquote(italic(.(z)))
  }))
  font <- "Helvetica"
  pt_mm <- 72.27 / 25.4
  theme_clade <- function() {
    theme_bw(base_size = 6, base_family = font) +
      theme(text = element_text(size = 6, colour = "black"),
            axis.text = element_text(size = 5.7, colour = "black"),
            axis.title = element_text(size = 6),
            plot.title = element_text(size = 7, face = "bold"),
            plot.subtitle = element_text(size = 5.5, colour = "grey25"),
            plot.caption = element_text(size = 5, colour = "grey30", hjust = 0,
                                        family = "DejaVu Sans"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
            panel.border = element_rect(colour = "grey35", linewidth = 0.3),
            axis.ticks = element_line(linewidth = 0.25),
            legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
            legend.key.size = grid::unit(2.8, "mm"),
            plot.margin = margin(5, 5, 5, 5),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.key = element_rect(fill = "transparent", colour = NA))
  }
  is_bias <- function(s) s$bias %in% c("Enrichment bias", "Depletion bias")

  # Phylum: individual LFCs, with one separate count/test line below each row.
  # Count text is not located on the LFC scale, and no Fisher OR is drawn as LFC.
  phy <- phy %>% mutate(row = row_number(),
                        detail = sprintf("up %d / down %d  |  q=%s",
                                         n_up, n_down, fmt_q(q)))
  pd <- dat %>% left_join(phy %>% select(clade, row), by = c("Phylum" = "clade"))
  phy_marks <- phy[is_bias(phy), ] %>%
    mutate(shape = ifelse(bias == "Enrichment bias", 24L, 25L),
           x = ifelse(bias == "Enrichment bias", max_lfc + 0.65, min_lfc - 0.65))
  x_limits <- range(dat$lfc) + c(-1.0, 1.0)
  # right-hand stats column, drawn outside the panel: the panel stops at the
  # data range and the text sits in the plot margin, so no empty axis is carried
  x_stats <- x_limits[2] + 0.03 * diff(x_limits)
  p_phylum <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_point(data = filter(pd, !sig), aes(lfc, row + 0.07, fill = Phylum),
               shape = 21, size = 1.15, stroke = 0.1, colour = "grey60", alpha = 0.22,
               position = position_jitter(height = 0.13, width = 0, seed = 42)) +
    geom_point(data = filter(pd, sig), aes(lfc, row + 0.07, fill = Phylum),
               shape = 21, size = 1.5, stroke = 0.15, colour = "grey25", alpha = 0.9,
               position = position_jitter(height = 0.13, width = 0, seed = 42)) +
    geom_point(data = phy_marks, aes(x, row + 0.07, shape = shape),
               size = 1.7, fill = "black", colour = "black", stroke = 0.15) +
    geom_text(data = phy, aes(x = x_stats, y = row + 0.07, label = detail),
              hjust = 0, family = font, size = 5 / pt_mm, colour = "grey25") +
    scale_fill_manual(values = pal, guide = "none") + scale_shape_identity() +
    scale_y_continuous(breaks = phy$row, labels = italic_taxa(phy$clade),
                       limits = c(0.45, nrow(phy) + 0.5), expand = c(0, 0)) +
    scale_x_continuous(breaks = c(-5, 0, 5), expand = c(0, 0)) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(title = "Distribution of MAG responses by phylum",
         subtitle = "584 MAGs: 87 up, 203 down, 294 nonsignificant",
         x = "MAG log fold change (case vs control)", y = NULL,
         caption = paste("Filled = DA; faint = nonsignificant MAGs.",
                         "▲ enrichment bias; ▼ depletion bias, relative to other DA MAGs.",
                         "Two-sided Fisher; BH across 8 eligible phyla (≥4 MAGs).",
                         "NT = clade not tested.", sep = "\n")) +
    theme_clade() + theme(panel.grid.major.y = element_blank(),
                          plot.margin = margin(5.5, 74, 5.5, 5.5, "pt"))

  # A cubic size ramp: prevalence is 0.118-1.00 with a median of 0.934, so a
  # linear map leaves nearly every dot the same size.
  pow3_trans <- scales::trans_new("pow3", transform = function(x) x^3,
                                  inverse = function(x) x^(1/3))

  # Family: retain the zero-centred MAG/stem design. Include all nonsignificant
  # MAGs faintly and add per-clade counts without hiding small/no-DA families.
  family_order <- fam$clade
  fd <- dat %>% mutate(Family = factor(Family, family_order),
                       PrevSize = pmin(pmax(Prevalence, 0.25), 1))
  fam <- fam %>% mutate(Family = factor(clade, family_order))
  family_marks <- fam[is_bias(fam), ] %>%
    mutate(shape = ifelse(bias == "Enrichment bias", 24L, 25L),
           y = ifelse(bias == "Enrichment bias", pmax(0, max_lfc) + 0.45,
                      pmin(0, min_lfc) - 0.45))
  summary_y <- max(dat$lfc) - 0.6
  family_up <- sum(fam$bias == "Enrichment bias")
  family_down <- sum(fam$bias == "Depletion bias")
  family_total <- sum(fam$clade != "Unassigned")
  family_tested <- sum(fam$eligible)
  p_family <- ggplot(fd, aes(Family, lfc)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_segment(data = filter(fd, sig), aes(xend = Family, y = 0, yend = lfc, colour = Phylum),
                 linewidth = 0.22, alpha = 0.45) +
    geom_point(data = filter(fd, !sig), aes(fill = Phylum, size = PrevSize),
               shape = 21, colour = "grey60", stroke = 0.1, alpha = 0.22,
               position = position_jitter(width = 0.12, height = 0, seed = 42)) +
    geom_point(data = filter(fd, sig), aes(fill = Phylum, size = PrevSize),
               shape = 21, colour = "grey25", stroke = 0.18, alpha = 0.9) +
    geom_point(data = family_marks, aes(Family, y, shape = shape),
               inherit.aes = FALSE, size = 1.8, fill = "black", colour = "black") +
    annotate("text", x = 1, y = summary_y, label = paste0(family_up, "↑"),
             hjust = 0, size = 6 / pt_mm, family = "DejaVu Sans", fontface = "bold",
             colour = "#FF6347") +
    annotate("text", x = 2.8, y = summary_y, label = paste0(family_down, "↓"),
             hjust = 0, size = 6 / pt_mm, family = "DejaVu Sans", fontface = "bold",
             colour = "#377EB8") +
    annotate("text", x = 1, y = summary_y - 0.8,
             label = sprintf("n = %d families (%d tested)", family_total, family_tested),
             hjust = 0, size = 6 / pt_mm, family = font, fontface = "bold") +
    scale_fill_manual(values = pal, name = "Phylum", labels = italic_taxa) +
    scale_colour_manual(values = pal, guide = "none") + scale_shape_identity() +
    scale_size_continuous(range = c(0.7, 3.6), limits = c(0.25, 1),
                          transform = pow3_trans,
                          breaks = c(0.25, 0.5, 0.75, 1),
                          labels = c("\u2264 25%", "50%", "75%", "100%"),
                          name = "Detection prevalence") +
    scale_x_discrete(limits = family_order, labels = italic_taxa,
                     expand = expansion(add = 0.55)) +
    scale_y_continuous(limits = c(min(dat$lfc) - 1.1, max(dat$lfc) + 1.1),
                       breaks = c(-5, 0, 5), expand = c(0, 0)) +
    labs(title = "Distribution of MAG-level responses across families",
         subtitle = "584 MAGs in 78 named families + 1 unassigned group",
         x = NULL, y = "MAG log fold change (case vs control)",
         caption = paste("Filled = DA; faint = nonsignificant. ▲ / ▼ = enrichment / depletion bias relative to other DA MAGs (not aggregate family abundance).",
                         "Two-sided Fisher among 87 enriched and 203 depleted MAGs; BH across 33 eligible families (≥4 MAGs, ≥1 DA). Other families are shown but untested.",
                         "Point size: donors with >0 mapped reads / 76, scaled from 25% on a cubic size ramp so the 0.7-1.0 band where most MAGs sit stays separable. Exact odds ratios, 95% CIs and q-values are in the accompanying tables.", sep = "\n")) +
    theme_clade() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.7,
                                     margin = margin(t = 0.5)),
          axis.ticks.length.x = grid::unit(1.2, "pt"),
          panel.grid.major.x = element_blank(), legend.position = "right",
          legend.spacing.y = grid::unit(1, "mm")) +
    guides(fill = guide_legend(override.aes = list(size = 1.8, alpha = 1), order = 1),
           size = guide_legend(override.aes = list(fill = "grey65", alpha = 1), order = 2))

  # The inset now counts Fisher classifications, not signs of fitted BLUPs or
  # signs of the within-family mean. Exclude the unassigned group from families.
  bias_levels <- c("Not tested", "No detected bias", "Depletion bias", "Enrichment bias")
  inset_data <- fam %>% filter(clade != "Unassigned") %>%
    count(bias, Phylum, name = "n_families") %>%
    mutate(bias = factor(bias, bias_levels))
  inset_totals <- inset_data %>% group_by(bias) %>%
    summarise(total = sum(n_families), .groups = "drop")
  stopifnot(sum(inset_totals$total) == 78L)
  p_inset <- ggplot(inset_data, aes(n_families, bias, fill = Phylum)) +
    geom_col(width = 0.58, colour = "white", linewidth = 0.2) +
    geom_text(data = inset_totals, aes(x = total + 0.8, y = bias, label = total),
              inherit.aes = FALSE, hjust = 0, family = font, size = 6 / pt_mm) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(expand = expansion(add = c(0, 5)), breaks = seq(0, 50, 10)) +
    labs(title = "Families by directional Fisher result", x = "Number of families", y = NULL,
         subtitle = "78 named families; 33 eligible tests",
         caption = paste("Enrichment/depletion bias: BH q < 0.05 relative to other DA MAGs.",
                         "No detected bias: q ≥ 0.05. Not tested: <4 MAGs or no DA MAGs.",
                         "Stack colours use the family-panel phylum key; unassigned group excluded.", sep = "\n")) +
    theme_clade() + theme(panel.grid.major.y = element_blank())

  stats <- c("Figure 1g/h: directional Fisher clade summaries",
             "584 MAGs; 76 donors; primary DA = 87 enriched + 203 depleted; NS = 294.",
             "Fisher matrix rows: clade / remaining catalogue; columns: enriched / depleted.",
             "Background = all 290 primary DA MAGs, including members of untested clades.",
             "Two-sided exact Fisher; conditional odds ratio and unadjusted 95% exact CI.",
             "Eligibility: >=4 total MAGs, >=1 DA MAG and assigned rank; BH within rank.",
             "No continuity correction; zero/infinite odds ratios are retained in tables.",
             "OR>1 is enrichment bias versus the remaining DA MAGs, not necessarily an up majority.",
             "Exploratory catalogue-level comparison; MAG calls share donors and phylogenetic dependence.",
             "Not a test of summed clade abundance and not independent patient-level replication.",
             "All MAGs plotted; only annotated eligible clades receive Fisher triangles.",
             "Family point size retains detection prevalence at >0 mapped reads / 76 donors.",
             sprintf("Family corner summary: %d enrichment-biased (orange), %d depletion-biased (blue); %d named families, %d tested.",
                     family_up, family_down, family_total, family_tested),
             "Phylum canvas: 4 x 6 inches; family: 14 x 6; family inset: 5.2 x 3.2.",
             "Visible text: 5-7 pt; canonical phylum_colors_mags.rds palette; transparent canvas.")
  for (rank in c("Phylum", "Family")) {
    s <- if (rank == "Phylum") phy else fam
    stats <- c(stats, sprintf("%s: %d eligible; %d enrichment-biased; %d depletion-biased.",
                             rank, sum(s$eligible), sum(s$bias == "Enrichment bias"),
                             sum(s$bias == "Depletion bias")))
    for (i in which(is_bias(s))) stats <- c(stats, sprintf(
      "  %s: up=%d down=%d NS=%d total=%d; OR=%g CI=[%g,%g] p=%.8g q=%.8g",
      s$clade[i], s$n_up[i], s$n_down[i], s$n_ns[i], s$n_total[i], s$odds_ratio[i],
      s$ci_lower[i], s$ci_upper[i], s$p[i], s$q[i]))
  }
  input_hashes <- tools::md5sum(c(anc_file, tax_file, count_file, palette_file))
  stats <- c(stats, "Input MD5:", paste(names(input_hashes), unname(input_hashes)))
  if (save) {
    dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(tabledir, recursive = TRUE, showWarnings = FALSE)
    write_csv(phy, file.path(tabledir, "directional_fisher_phylum_n76.csv"))
    write_csv(fam %>% select(-Family), file.path(tabledir, "directional_fisher_family_n76.csv"))
    write_csv(dat, file.path(tabledir, "directional_fisher_MAG_plot_data_n76.csv"))
    write_csv(inset_totals, file.path(tabledir, "directional_fisher_family_counts_n76.csv"))
    ggplot2::ggsave(file.path(figdir, "Fig1g_Phylum_n76.pdf"), p_phylum,
                    width = 4, height = 6, units = "in", device = cairo_pdf, bg = "transparent")
    ggplot2::ggsave(file.path(figdir, "Fig1h_Family_n76.pdf"), p_family,
                    width = 14, height = 6, units = "in", device = cairo_pdf, bg = "transparent")
    ggplot2::ggsave(file.path(figdir, "Fig1h_Family_inset_counts_n76.pdf"), p_inset,
                    width = 5.2, height = 3.2, units = "in", device = cairo_pdf, bg = "transparent")
    writeLines(stats, file.path(figdir, "Fig1gh_DIRECTIONAL_FISHER_STATS_n76.txt"))
    message(paste(stats[1:16], collapse = "\n"))
  }
  invisible(list(p_phylum = p_phylum, p_family = p_family, p_inset = p_inset,
                 phylum = phy, family = fam, data = dat, stats = stats))
}

if (sys.nframe() == 0L) build_clade_fisher_panels()
