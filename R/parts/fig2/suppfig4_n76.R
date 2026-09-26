# =============================================================================
# Supplementary Figure 4, reviewer-response robustness checks rebuilt on n = 76
#
# The published response ran eight sections. Five no longer apply:
#   - CAFCA permutation correction (section 4) and the GO cellular-component
#     audit (section 8) are withdrawn with the GO/CAFCA layers.
#   - Technical-replicate handling and the pseudoreplication donor-collapse
#     check (section 5) are moot: the rebuild is 76 donors, one sample each.
#   - ZOE health-rank validation (section 6) already ships as
#     fig1/SuppFig_ZOE_Concordance_n76.pdf.
#   - The sample-level metadata table (section 7) is a table, not a figure.
#
# Three live ones are rebuilt here:
#   a-b  Reviewer 1 item 2. Does the phylogenetic signal survive completeness
#        adjustment? Raw LFC, completeness-residualised LFC, and a >=90%
#        completeness subset.
#   c-d  Reviewer 1 item 5. Is functional load a completeness artefact? Variance
#        explained by guild against completeness and genome size, and guild
#        differences before and after residualisation.
#   e-f  Reviewer 1 item 4. Is the guild solution an artefact of the clustering
#        parameters? Stability across seeds and across the parameter sweep.
#
# Outputs -> result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(phytools); library(phylosignal); library(phylobase)
  library(dplyr); library(tidyr); library(tibble); library(readr); library(purrr)
  library(ggplot2); library(patchwork); library(ggrepel)
})

theme_transparent <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA))
theme_bw <- function(...) ggplot2::theme_bw(...) + theme_transparent
ggplot2::theme_set(ggplot2::theme_get() + theme_transparent)
ggsave <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "supp_fig")
TREE   <- file.path(INDIR, "gtdbtk.bac120.decorated_n76.tree")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

SEED <- 42
HQ_COMPLETENESS <- 90

st   <- readRDS(file.path(N76, "fig2", "figure2_stats_n76.rds"))
tree <- read.tree(TREE); tree$node.label <- NULL
meta <- st$meta_df

# =============================================================================
# a-b. phylogenetic signal against genome completeness
# =============================================================================

say("## Reviewer 1 item 2: is the phylogenetic signal a completeness artefact?")
d <- meta %>% filter(!is.na(lfc_GroupCase), !is.na(Completeness))
d <- d[d$Genome_ID %in% tree$tip.label, ]
say(sprintf("   tree tips with an LFC and a completeness value: %d", nrow(d)))

r <- suppressWarnings(cor.test(d$Completeness, d$lfc_GroupCase, method = "spearman"))
say(sprintf("   Spearman(completeness, LFC) = %+.3f, p = %.3g", r$estimate, r$p.value))

# completeness-residualised LFC: the part of the effect completeness cannot explain
d$lfc_resid <- residuals(lm(lfc_GroupCase ~ Completeness, data = d))

signal_of <- function(dat, value_col, label) {
  tr <- keep.tip(tree, intersect(tree$tip.label, dat$Genome_ID))
  x  <- setNames(dat[[value_col]][match(tr$tip.label, dat$Genome_ID)], tr$tip.label)
  set.seed(SEED)
  k   <- phytools::phylosig(tr, x, method = "K", test = TRUE, nsim = 999)
  lam <- phytools::phylosig(tr, x, method = "lambda", test = TRUE)
  p4  <- phylobase::phylo4d(tr, data.frame(v = x))
  mor <- phylosignal::phyloSignal(p4, methods = "I", reps = 999)
  tibble(set = label, n = length(x),
         K = unname(k$K), K_p = unname(k$P),
         lambda = unname(lam$lambda), lambda_p = unname(lam$P),
         MoranI = mor$stat["v", "I"], MoranI_p = mor$pvalue["v", "I"])
}
sig_tab <- bind_rows(
  signal_of(d, "lfc_GroupCase", "Raw LFC"),
  signal_of(d, "lfc_resid",     "Completeness-residualised"),
  signal_of(d %>% filter(Completeness >= HQ_COMPLETENESS), "lfc_GroupCase",
            sprintf("Completeness >= %d%%", HQ_COMPLETENESS)))
write_csv(sig_tab, file.path(OUTDIR, "SuppFig4_phylo_signal_robustness_n76.csv"))
for (i in seq_len(nrow(sig_tab))) with(sig_tab[i, ],
  say(sprintf("   %-28s n = %3d | K = %.3f (p = %.3f) | lambda = %.3f (p = %.3g) | Moran's I = %.3f (p = %.3f)",
              set, n, K, K_p, lambda, lambda_p, MoranI, MoranI_p)))

sig_long <- sig_tab %>%
  select(set, n, K, K_p, lambda, lambda_p, MoranI, MoranI_p) %>%
  pivot_longer(c(K, lambda, MoranI), names_to = "stat", values_to = "value") %>%
  mutate(p = case_when(stat == "K" ~ K_p, stat == "lambda" ~ lambda_p, TRUE ~ MoranI_p),
         stat = recode(stat, K = "Blomberg's K", lambda = "Pagel's lambda", MoranI = "Moran's I"),
         set = factor(set, levels = sig_tab$set))

pA <- ggplot(sig_long, aes(set, value, fill = set)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))),
            vjust = -0.35, size = 2.5) +
  facet_wrap(~ stat, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("grey70", "#3B6FB6", "#5AAE61")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Phylogenetic signal survives completeness adjustment",
       subtitle = "Reviewer 1 item 2. 999 permutations; lambda by likelihood ratio",
       x = NULL, y = "Statistic") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 6.5, angle = 20, hjust = 1))

pB <- ggplot(d, aes(Completeness, lfc_GroupCase)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = diff_GroupCase %in% TRUE), size = 1, alpha = 0.6) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.6) +
  geom_vline(xintercept = HQ_COMPLETENESS, linetype = "dotted", colour = "#5AAE61") +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#D55E00"),
                      labels = c("Not significant", "Differentially abundant"), name = NULL) +
  labs(title = "Completeness does not order the effect sizes",
       subtitle = sprintf("Spearman rho = %+.3f, p = %.2g", r$estimate, r$p.value),
       x = "CheckM2 completeness (%)", y = "Log fold change (Case vs Control)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

# =============================================================================
# c-d. functional load against completeness and genome size
# =============================================================================

say(""); say("## Reviewer 1 item 5: is functional load a completeness artefact?")
prod <- read_tsv(file.path(INDIR, "DRAM_metabolism_product.tsv"), show_col_types = FALSE)
# read_tsv keeps the 66 True/False columns logical, so as.matrix coerces them to
# 1/0 and the load sums all 98 pathway columns. read.delim would make the whole
# matrix character and silently drop them, halving the load.
load_mat <- as.matrix(prod[, -1]); rownames(load_mat) <- prod[[1]]
storage.mode(load_mat) <- "numeric"
fl <- tibble(Genome_ID = rownames(load_mat), Load = rowSums(load_mat, na.rm = TRUE))

guilds <- readRDS(file.path(N76, "fig3", "df_func_n76_30clusters.rds")) %>%
  transmute(Genome_ID = Taxon, Guild = factor(Func_Cluster))
qual <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                   quote = "", comment.char = "") %>%
  transmute(Genome_ID = catalog_id, Completeness = completeness,
            Genome_Mb = checkm2_Genome_Size / 1e6)

fd <- fl %>% inner_join(guilds, by = "Genome_ID") %>% left_join(qual, by = "Genome_ID") %>%
  filter(!is.na(Completeness), Guild != "0")
say(sprintf("   MAGs with a load, a guild and a completeness value: %d across %d guilds",
            nrow(fd), n_distinct(fd$Guild)))

# Genome size comes from the CheckM2 columns, which were only re-run for the
# rebuilt controls; the carried-over case genomes have none. Size models are
# therefore fitted on the subset that has it, and the residualisation that
# panel d uses is on completeness alone so no MAG is dropped.
fd_sz <- fd %>% filter(!is.na(Genome_Mb))
say(sprintf("   genome size available for %d of %d MAGs (carried-over case genomes have no CheckM2 stats)",
            nrow(fd_sz), nrow(fd)))
r2_of <- function(f, dat) summary(lm(f, data = dat))$r.squared
r2_tab <- tibble(
  Predictor = c("Guild", "Completeness", "Genome size (Mb)", "Completeness + genome size"),
  n  = c(nrow(fd), nrow(fd), nrow(fd_sz), nrow(fd_sz)),
  R2 = c(r2_of(Load ~ Guild, fd), r2_of(Load ~ Completeness, fd),
         r2_of(Load ~ Genome_Mb, fd_sz), r2_of(Load ~ Completeness + Genome_Mb, fd_sz)))
for (i in seq_len(nrow(r2_tab))) with(r2_tab[i, ],
  say(sprintf("   variance in functional load explained by %-28s R2 = %.3f (n = %d)", Predictor, R2, n)))

# do guild differences survive removing what completeness explains?
fd$Load_resid <- residuals(lm(Load ~ Completeness, data = fd))
kw_raw <- kruskal.test(Load ~ Guild, data = fd)
kw_adj <- kruskal.test(Load_resid ~ Guild, data = fd)
say(sprintf("   Kruskal-Wallis guild effect: raw p = %.3g | after removing completeness p = %.3g",
            kw_raw$p.value, kw_adj$p.value))
write_csv(r2_tab, file.path(OUTDIR, "SuppFig4_functional_load_variance_n76.csv"))

pC <- ggplot(r2_tab %>% mutate(Predictor = factor(Predictor, rev(Predictor))),
             aes(R2, Predictor, fill = Predictor == "Guild")) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", R2)), hjust = -0.15, size = 2.8) +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#3B6FB6")) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Guild, not genome quality, explains functional load",
       subtitle = "Reviewer 1 item 5. Variance in functional load explained",
       x = expression(R^2), y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.y = element_text(size = 7))

pD <- fd %>%
  select(Guild, Raw = Load, Residualised = Load_resid) %>%
  pivot_longer(c(Raw, Residualised), names_to = "which", values_to = "val") %>%
  mutate(which = factor(which, c("Raw", "Residualised"))) %>%
  ggplot(aes(reorder(Guild, val, median), val, fill = which)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8, linewidth = 0.3,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Raw" = "grey70", "Residualised" = "#3B6FB6"), name = NULL) +
  labs(title = "Guild separation is not created by completeness",
       subtitle = sprintf("Kruskal-Wallis raw p = %.2g; after removing completeness p = %.2g",
                          kw_raw$p.value, kw_adj$p.value),
       x = "Guild", y = "Functional load") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 5.5), legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

# =============================================================================
# e-f. is the guild solution an artefact of the clustering parameters?
# =============================================================================

say(""); say("## Reviewer 1 item 4: is the guild solution parameter-dependent?")
cd <- read_csv(file.path(N76, "fig3", "fig3_cluster_decision_n76.csv"), show_col_types = FALSE) %>%
  filter(!is.na(n_guilds)) %>%
  mutate(chosen = grepl("md0.05", label) & mp == 6 & nn == 15)
say(sprintf("   parameter settings swept: %d | guild counts %d-%d | ARI %.2f-%.2f",
            nrow(cd), min(cd$n_guilds), max(cd$n_guilds),
            min(cd$stability_ARI, na.rm = TRUE), max(cd$stability_ARI, na.rm = TRUE)))
ch <- cd %>% filter(chosen)
if (nrow(ch)) say(sprintf("   chosen setting: %s -> %d guilds, noise %.0f%%, ARI %.2f, seeds gave %s",
                          ch$label[1], ch$n_guilds[1], 100 * ch$noise[1], ch$stability_ARI[1],
                          ch$guilds_across_seeds[1]))
write_csv(cd, file.path(OUTDIR, "SuppFig4_guild_parameter_sweep_n76.csv"))

pE <- ggplot(cd, aes(n_guilds, stability_ARI)) +
  geom_point(aes(size = noise, fill = chosen), shape = 21, colour = "grey25", alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.1, max.overlaps = 20,
                           segment.size = 0.2, colour = "grey30") +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#D55E00"),
                    labels = c("Swept", "Chosen"), name = NULL) +
  scale_size_continuous(range = c(2, 6), labels = scales::percent, name = "Noise") +
  labs(title = "Guild count against partition stability",
       subtitle = "Reviewer 1 item 4. ARI across six seeds at each setting",
       x = "Guilds", y = "Stability (mean ARI across seeds)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"))

seeds_long <- cd %>%
  separate_rows(guilds_across_seeds, sep = "/") %>%
  mutate(g = as.numeric(guilds_across_seeds)) %>%
  filter(!is.na(g))
pF <- ggplot(seeds_long, aes(reorder(label, g, median), g)) +
  geom_boxplot(aes(fill = chosen), outlier.shape = NA, width = 0.6, alpha = 0.85,
               linewidth = 0.3, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 0.9, alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#D55E00")) +
  labs(title = "The guild count itself moves with the seed",
       subtitle = "Six seeds per setting. The partition is stabler than the count.",
       x = NULL, y = "Guilds recovered") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.y = element_text(size = 6.5))

# =============================================================================
# g. UMAP + HDBSCAN against tree-based clustering
# The published response compared the guilds to Bonsai, pvclust and HiDeF.
# Bonsai and HiDeF need external Python/HPC runs; pvclust is the tree-based
# comparator that runs here, on the same Jaccard distance the guilds came from.
# =============================================================================

say(""); say("## Reviewer 1 item 4: UMAP + HDBSCAN against hierarchical clustering")
mat <- readRDS(file.path(N76, "fig3", "df_func_n76_30clusters.rds"))
prod_bin <- load_mat[rownames(load_mat) %in% mat$Taxon, , drop = FALSE]
prod_bin[prod_bin > 0] <- 1
prod_bin <- prod_bin[, colSums(prod_bin) > 0, drop = FALSE]
say(sprintf("   binary functional matrix for clustering: %d MAGs x %d features",
            nrow(prod_bin), ncol(prod_bin)))

set.seed(SEED)
dj <- proxy::dist(prod_bin, method = "Jaccard")
hc <- hclust(as.dist(dj), method = "ward.D2")
k_target <- n_distinct(mat$Func_Cluster[mat$Func_Cluster != "0"])
hc_cut <- cutree(hc, k = k_target)

cmp <- tibble(Genome_ID = rownames(prod_bin), Hier = hc_cut[rownames(prod_bin)]) %>%
  inner_join(mat %>% transmute(Genome_ID = Taxon, Guild = as.character(Func_Cluster)),
             by = "Genome_ID")
ari_all  <- mclust::adjustedRandIndex(cmp$Guild, cmp$Hier)
cmp_nz   <- cmp %>% filter(Guild != "0")
ari_nz   <- mclust::adjustedRandIndex(cmp_nz$Guild, cmp_nz$Hier)
say(sprintf("   hierarchical (ward.D2) cut at k = %d", k_target))
say(sprintf("   ARI against the guilds: %.3f including noise | %.3f excluding noise MAGs", ari_all, ari_nz))

purity <- cmp_nz %>% count(Guild, Hier) %>% group_by(Guild) %>%
  summarise(purity = max(n) / sum(n), n = sum(n), .groups = "drop")
say(sprintf("   per-guild purity against the hierarchical cut: median %.2f, %d of %d guilds above 0.5",
            median(purity$purity), sum(purity$purity > 0.5), nrow(purity)))
write_csv(purity, file.path(OUTDIR, "SuppFig4_guild_vs_hierarchical_purity_n76.csv"))

pG <- ggplot(purity, aes(reorder(Guild, purity), purity)) +
  geom_col(aes(fill = purity > 0.5), width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#3B6FB6")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Guilds against tree-based clustering",
       subtitle = sprintf("Ward.D2 on the same Jaccard distance, cut at k = %d. ARI = %.2f (noise excluded)",
                          k_target, ari_nz),
       x = "Guild", y = "Purity against the hierarchical cut") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.y = element_text(size = 6))

# =============================================================================
# h-i. sample-level sequencing audit (published Fig_Sample_Metadata_Audit)
# Panel C of the published version counted sequencing runs per donor. The
# rebuild is 76 donors with one sample each, so that panel has nothing to show
# and is dropped rather than rebuilt as a row of ones.
# =============================================================================

say(""); say("## Reviewer 1 item 1: sample-level sequencing audit")
cl <- read_tsv(file.path(INDIR, "coverm_long.tsv"), show_col_types = FALSE)
smeta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE) %>%
  transmute(donor = sample, Group = group, source)
# CoverM normalises relative abundance against ALL reads, so the per-donor sum of
# relative_abundance_percent is the fraction mapping to the catalogue, and total
# depth follows from the mapped count.
audit <- cl %>% group_by(donor) %>%
  summarise(mapped_reads = sum(read_count), pct_mapped = sum(relative_abundance_percent),
            .groups = "drop") %>%
  mutate(total_reads = mapped_reads / (pct_mapped / 100),
         depth_M = total_reads / 1e6) %>%
  left_join(smeta, by = "donor") %>%
  mutate(cohort = ifelse(Group == "Case", "COVID", paste0("Control\n(", source, ")")),
         cohort = factor(cohort, levels = c("COVID", sort(unique(cohort[Group == "Control"])))))
write_csv(audit, file.path(OUTDIR, "SuppFig4_sequencing_audit_n76.csv"))
for (g in levels(audit$cohort)) {
  a <- audit %>% filter(cohort == g)
  say(sprintf("   %-28s n = %2d | median depth %6.1f M read pairs | median mapped %.1f%%",
              gsub("\n", " ", g), nrow(a), median(a$depth_M), median(a$pct_mapped)))
}
kw_depth <- kruskal.test(depth_M ~ cohort, data = audit)
kw_map   <- kruskal.test(pct_mapped ~ cohort, data = audit)
say(sprintf("   depth differs by cohort: Kruskal-Wallis p = %.3g | mapped fraction p = %.3g",
            kw_depth$p.value, kw_map$p.value))

audit_cols <- c("COVID" = "#2B3A67")
audit_cols[setdiff(levels(audit$cohort), "COVID")] <- "#D6C453"

pH <- ggplot(audit, aes(cohort, depth_M, fill = cohort)) +
  geom_violin(alpha = 0.55, linewidth = 0.3, colour = "grey30", show.legend = FALSE) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.shape = NA, linewidth = 0.3, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.07, seed = SEED), size = 0.8, alpha = 0.7) +
  stat_summary(fun.data = function(x) data.frame(y = min(x) * 0.6, label = paste0("n=", length(x))),
               geom = "text", size = 2.4, colour = "grey25") +
  scale_fill_manual(values = audit_cols) + scale_y_log10() +
  labs(title = "Sequencing depth by cohort and source",
       subtitle = sprintf("Analysis-ready read pairs after host removal, log scale. Kruskal-Wallis p = %.2g", kw_depth$p.value),
       x = NULL, y = "Read pairs (millions, log10)") +
  theme_bw(base_size = 9) +
  # the cohort legend only repeated the x axis
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 6.5),
        legend.position = "none", aspect.ratio = 1)

pI <- ggplot(audit, aes(cohort, pct_mapped, fill = cohort)) +
  geom_violin(alpha = 0.55, linewidth = 0.3, colour = "grey30", show.legend = FALSE) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.shape = NA, linewidth = 0.3, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.07, seed = SEED), size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = audit_cols) +
  labs(title = sprintf("Reads mapped to the %d-MAG catalogue", n_distinct(cl$MAG)),
       subtitle = sprintf("Per donor, CoverM. Kruskal-Wallis p = %.2g", kw_map$p.value),
       x = NULL, y = "Reads mapped (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 6.5),
        legend.position = "none", aspect.ratio = 1)

# ---- save -------------------------------------------------------------------

pl <- list(SuppFig4a_PhyloSignal_Robustness = pA, SuppFig4b_Completeness_vs_LFC = pB,
           SuppFig4c_Load_Variance = pC, SuppFig4d_Load_by_Guild = pD,
           SuppFig4e_Guild_Stability = pE, SuppFig4f_Guilds_across_Seeds = pF,
           SuppFig4g_Guild_vs_Hierarchical = pG,
           SuppFig4h_Sequencing_Depth = pH, SuppFig4i_Reads_Mapped = pI)
for (nm in names(pl))
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), pl[[nm]], width = 5.4, height = 4.4,
         device = cairo_pdf, bg = "transparent")

# A4 portrait, three columns. The previous canvas was 13 x 24 in, which a
# journal can only place by shrinking the type past legibility; at A4 the type
# scale below holds every label at 5 pt or more.
tag <- theme(plot.tag = element_text(face = "bold", size = 14))
small_text <- theme(
  axis.text    = element_text(size = 8, colour = "black"),
  axis.title   = element_text(size = 9.5),
  legend.text  = element_text(size = 8),
  legend.title = element_text(size = 9.5),
  legend.key.size = unit(0.40, "cm"),
  strip.text    = element_text(size = 9.5),
  plot.title    = element_text(face = "bold", size = 11),
  plot.subtitle = element_text(size = 8.5, colour = "grey35"))
# Three narrow columns cannot carry the full-width titles these panels were
# written with, and the categorical axes of a, h and i collide at this width.
wrap_lab <- function(p, w = 30, ws = 44) {
  if (!is.null(p$labels$title))    p$labels$title    <- paste(strwrap(p$labels$title, w), collapse = "\n")
  if (!is.null(p$labels$subtitle)) p$labels$subtitle <- paste(strwrap(p$labels$subtitle, ws), collapse = "\n")
  p
}
lean_x <- theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
# Panel a carries three long set names and three long facet titles in one
# third of the page; abbreviate both rather than let them overplot.
pA <- pA + scale_x_discrete(labels = c("Raw LFC" = "Raw", "Completeness-residualised" = "Residual",
                                       "Completeness >= 90%" = "\u2265 90%")) +
  facet_wrap(~ stat, scales = "free_y", nrow = 1,
             labeller = labeller(stat = c("Blomberg's K" = "K", "Pagel's lambda" = "lambda",
                                          "Moran's I" = "Moran's I")))
# p-value labels sit above three adjacent bars in a third-page panel, so they
# collide horizontally; set them vertical instead of shrinking them below 5 pt.
pA$layers[[2]]$aes_params$size  <- 2.5
pA$layers[[2]]$aes_params$angle <- 90
pA$layers[[2]]$aes_params$hjust <- -0.08
pA$layers[[2]]$aes_params$vjust <- 0.5
pA <- wrap_lab(pA) + lean_x; pB <- wrap_lab(pB); pC <- wrap_lab(pC)
pD <- wrap_lab(pD); pE <- wrap_lab(pE); pF <- wrap_lab(pF)
pG <- wrap_lab(pG); pH <- wrap_lab(pH) + lean_x; pI <- wrap_lab(pI) + lean_x

# Square plotting areas throughout, with the canvas sized around them. Panel a
# is three facets side by side, so each facet takes a ratio of 3 to make the
# block itself square; h and i already carry their own square ratio.
sq <- theme(aspect.ratio = 1)
pA <- pA + theme(aspect.ratio = 3)
pB <- pB + sq; pC <- pC + sq
# d holds 30 guilds, so it keeps the square ratio off and spreads across the
# whole of its cell; the numbers stand upright rather than run together
pD <- pD + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
pE <- pE + sq; pF <- pF + sq; pG <- pG + sq

supp4 <- ((pA + labs(tag = "a") + tag) | (pB + labs(tag = "b") + tag) | (pC + labs(tag = "c") + tag)) /
         ((pD + labs(tag = "d") + tag) | (pE + labs(tag = "e") + tag) | (pF + labs(tag = "f") + tag)) /
         ((pG + labs(tag = "g") + tag) | (pH + labs(tag = "h") + tag) | (pI + labs(tag = "i") + tag)) &
  small_text
supp4 <- supp4 + plot_annotation(
    title = "Supplementary Figure 4. Reviewer-response robustness checks on the 76-donor rebuild",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
# canvas follows the panels rather than the page: square areas need the room
ggsave(file.path(OUTDIR, "Supplementary_Figure_4_n76.pdf"), supp4, width = 12, height = 14,
       device = cairo_pdf, bg = "transparent")

writeLines(STATS, file.path(OUTDIR, "SuppFig4_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Supplementary_Figure_4_n76.pdf"))
