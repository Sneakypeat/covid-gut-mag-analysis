# =============================================================================
# Supplementary Figure 2, 76-donor catalogue
#
#   a      heatmap of species-level functional drivers, sub-clustered by status
#   b-d    t-SNE overlays of cytochrome completeness (Complex III, high-, low-affinity)
#   e      t-SNE coloured by log fold change
#   f,g,i  binary t-SNE overlays (nitrate/TMAO, sulfur oxyanion redox, butyrate)
#   h      phylum overlay
#   j      raincloud of LFC for module-carrying MAGs, Wilcoxon against zero
#
# The t-SNE embedding is the one Figure 3 drew (fig3/df_tsne_plot_n76.rds), so
# every panel sits on identical coordinates rather than a re-run of Rtsne.
#
# Outputs -> result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(purrr)
  library(ggplot2); library(patchwork); library(pheatmap); library(reshape2)
  library(viridis); library(ggdist); library(grid)
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
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "supp_fig")
source(file.path(BASE, "taxon_italics_n76.R"))   # md_taxon(): taxon names in italics
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

phylum_colors <- readRDS(file.path(BASE, "results/phylum_colors_mags.rds"))
df_tsne_plot <- readRDS(file.path(N76, "fig3", "df_tsne_plot_n76.rds"))
res_shift <- read_csv(file.path(N76, "fig3", "fig3_functional_drivers_n76.csv"), show_col_types = FALSE)
tax <- read.delim(file.path(INDIR, "MAG_quality_taxonomy.tsv"), check.names = FALSE,
                  quote = "", comment.char = "")

dram <- read.delim(file.path(INDIR, "DRAM_metabolism_product.tsv"), check.names = FALSE,
                   quote = "", comment.char = "")
rownames(dram) <- dram[[1]]; dram[[1]] <- NULL
mat_completeness <- as.matrix(data.frame(lapply(dram, function(x) {
  x <- as.character(x); x[x == "True"] <- "1"; x[x == "False"] <- "0"
  suppressWarnings(as.numeric(x))
}), check.names = FALSE, row.names = rownames(dram)))
mat_completeness[is.na(mat_completeness)] <- 0
cn <- colnames(mat_completeness)
score_module <- function(mat, cols, mode = "max") {
  if (length(cols) == 0) return(rep(0, nrow(mat)))
  sub <- mat[, cols, drop = FALSE]
  if (mode == "max") apply(sub, 1, max, na.rm = TRUE) else rowMeans(sub, na.rm = TRUE)
}
say(sprintf("MAGs on the t-SNE: %d | DRAM features: %d", nrow(df_tsne_plot), ncol(mat_completeness)))

# ---- a. species-level driver heatmap ---------------------------------------

blank <- function(x) is.na(x) | x == ""
df_driver_meta <- tax %>%
  transmute(Taxon = catalog_id, Genus = genus, Species_raw = species) %>%
  mutate(Display_Name = case_when(!blank(Species_raw) ~ gsub("_", " ", Species_raw),
                                  !blank(Genus) ~ paste(gsub("_", " ", Genus), "sp."),
                                  TRUE ~ Taxon)) %>%
  select(Taxon, Species = Display_Name) %>%
  left_join(df_tsne_plot %>% select(Taxon, LFC_Status, Phylum_Plot, Func_Cluster), by = "Taxon") %>%
  mutate(Species = make.unique(Species))

meta_lfc <- df_tsne_plot %>% select(Taxon, LFC_capped)
mat_sub <- mat_completeness[meta_lfc$Taxon, ]
mat_filtered <- mat_sub[, colSums(mat_sub > 0) >= 5]
target_functions <- res_shift %>% filter(P_Val < 0.05) %>% arrange(desc(Diff))   # same rule as Fig. 3c

driver_data_long <- mat_sub[, colnames(mat_sub) %in% target_functions$Function] %>%   # all 35, as in Fig. 3c
  as.data.frame() %>% rownames_to_column("Taxon") %>%
  inner_join(df_driver_meta, by = "Taxon") %>%
  pivot_longer(cols = any_of(target_functions$Function), names_to = "Function",
               values_to = "Completeness")

# Keep species that are significant AND actually drive a pathway (a MAG >80% complete)
active_species_summary <- driver_data_long %>%
  group_by(Species, Function, LFC_Status, Phylum_Plot) %>%
  summarise(MAG_Count = sum(Completeness > 0.8), .groups = "drop") %>%
  filter(LFC_Status %in% c("Enriched", "Depleted"), MAG_Count > 0)
valid_species <- unique(active_species_summary$Species)
say(sprintf("PANEL a: %d significant driver functions x %d active species",
            nrow(target_functions), length(valid_species)))

plot_data_long <- driver_data_long %>% filter(Species %in% valid_species) %>% mutate(Label = Function)
plot_mat <- acast(plot_data_long, Species ~ Label, value.var = "Completeness", fun.aggregate = mean)
plot_mat[plot_mat == 0] <- NA

species_status_map <- df_driver_meta %>% filter(Species %in% rownames(plot_mat)) %>%
  select(Species, Status = LFC_Status) %>% distinct(Species, .keep_all = TRUE)
get_clustered_order <- function(mat) {
  if (nrow(mat) < 2) return(mat)
  mz <- mat; mz[is.na(mz)] <- 0
  mat[hclust(dist(mz, method = "euclidean"), method = "ward.D2")$order, , drop = FALSE]
}
mat_enr <- get_clustered_order(plot_mat[rownames(plot_mat) %in%
                                          species_status_map$Species[species_status_map$Status == "Enriched"], , drop = FALSE])
mat_dep <- get_clustered_order(plot_mat[rownames(plot_mat) %in%
                                          species_status_map$Species[species_status_map$Status == "Depleted"], , drop = FALSE])
plot_mat_ordered <- rbind(mat_enr, mat_dep)
gap_row_pos <- nrow(mat_enr)
say(sprintf("   enriched species %d | depleted species %d", nrow(mat_enr), nrow(mat_dep)))

ann_row <- df_driver_meta %>% filter(Species %in% rownames(plot_mat_ordered)) %>%
  select(Species, Phylum = Phylum_Plot, Status = LFC_Status) %>%
  distinct(Species, .keep_all = TRUE) %>% column_to_rownames("Species")
ann_row <- ann_row[rownames(plot_mat_ordered), , drop = FALSE]
ann_row$Status <- droplevels(factor(ann_row$Status))
ann_colors <- list(Status = c("Enriched" = "#D55E00", "Depleted" = "#0072B2"),
                   Phylum = phylum_colors[names(phylum_colors) %in% unique(ann_row$Phylum)])

valid_cols <- intersect(target_functions$Function, colnames(plot_mat_ordered))
plot_mat_ordered <- plot_mat_ordered[, valid_cols, drop = FALSE]
magma_cols <- viridis(100, option = "magma", begin = 0.2, end = 0.9, direction = -1)

# ---- a. species x function heatmap, square tiles ---------------------------
# Drawn as one tile per cell at a fixed physical size (ggh4x), not stretched to
# fill a panel: 232 species over one column would need 16.7 in of height to
# label at 5 pt, so the species are split into two columns of 116. Status and
# phylum ride as two annotation columns on the left of each block, each with its
# own fill scale (ggnewscale).
suppressPackageStartupMessages({ library(ggh4x); library(ggnewscale) })
TILE_A  <- 0.072                      # inches; 5 pt type needs ~0.069 in of row
SP_PER  <- ceiling(nrow(plot_mat_ordered) / 2)
clip_at <- function(x, n) ifelse(nchar(x) <= n, x,
                                 paste0(sub("\\s+\\S*$", "", substr(x, 1, n)), "\u2026"))

sp_order <- rownames(plot_mat_ordered)
hm_long <- as.data.frame(plot_mat_ordered) %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Function", values_to = "Completeness") %>%
  mutate(rank  = match(Species, sp_order),
         block = ifelse(rank <= SP_PER, sprintf("Species 1-%d", SP_PER),
                        sprintf("Species %d-%d", SP_PER + 1, length(sp_order))),
         Function = factor(clip_at(Function, 46),
                           levels = clip_at(colnames(plot_mat_ordered), 46)),
         Species  = factor(Species, levels = rev(sp_order)),
         lab      = clip_at(as.character(Species), 30))
sp_labels <- setNames(clip_at(rev(sp_order), 30), rev(sp_order))

ann_long <- ann_row %>% rownames_to_column("Species") %>%
  mutate(rank  = match(Species, sp_order),
         block = ifelse(rank <= SP_PER, sprintf("Species 1-%d", SP_PER),
                        sprintf("Species %d-%d", SP_PER + 1, length(sp_order))),
         Species = factor(Species, levels = rev(sp_order)))

p_hm <- ggplot(hm_long, aes(Function, Species)) +
  geom_tile(aes(fill = Completeness), colour = NA) +
  scale_fill_gradientn(colours = magma_cols, na.value = "grey97",
                       name = "Pathway\ncompleteness", limits = c(0, 1)) +
  new_scale_fill() +
  geom_tile(data = ann_long, aes(x = "Status", y = Species, fill = Status), inherit.aes = FALSE) +
  scale_fill_manual(values = c("Enriched" = "#D55E00", "Depleted" = "#0072B2"), name = "Status") +
  new_scale_fill() +
  geom_tile(data = ann_long, aes(x = "Phylum", y = Species, fill = Phylum), inherit.aes = FALSE) +
  scale_fill_manual(values = phylum_colors, name = "Phylum (GTDB)",
                    labels = function(x) expr_from_md(md_taxon(x))) +
  scale_x_discrete(limits = c("Status", "Phylum", levels(hm_long$Function)), expand = c(0, 0)) +
  scale_y_discrete(labels = sp_labels, expand = c(0, 0)) +
  facet_wrap(~ block, nrow = 1, scales = "free_y") +
  force_panelsizes(rows = unit(SP_PER * TILE_A, "in"),
                   cols = unit((ncol(plot_mat_ordered) + 2) * TILE_A, "in")) +
  labs(title = "Species-level functional drivers",
       subtitle = sprintf("%d active species x %d significant driver functions; species ordered by response status, then clustered",
                          nrow(plot_mat_ordered), ncol(plot_mat_ordered)),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(size = 5, colour = "black", face = "italic"),
        panel.grid = element_blank(), strip.text = element_text(size = 6, face = "bold"),
        legend.title = element_text(size = 6, face = "bold"), legend.text = element_text(size = 5),
        legend.key.size = unit(0.22, "cm"), legend.position = "right",
        plot.title = element_text(face = "bold", size = 8),
        plot.subtitle = element_text(size = 5.5, colour = "grey35"))

hm_file <- file.path(OUTDIR, "SuppFig2a_Driver_Heatmap_n76.pdf")
ggsave(hm_file, p_hm, width = 8.27, height = 11.69, device = cairo_pdf,
       bg = "transparent", limitsize = FALSE)

# ---- b-d. cytochrome overlays ----------------------------------------------

cyto_cols_all <- cn[grepl("Complex III|High affinity|Low affinity", cn)]
df_cyto_final <- mat_completeness %>% as.data.frame() %>% select(all_of(cyto_cols_all)) %>%
  rownames_to_column("Taxon") %>%
  pivot_longer(-Taxon, names_to = "Raw_Function", values_to = "Completeness") %>%
  filter(Completeness > 0) %>%
  mutate(Guild = case_when(grepl("High affinity", Raw_Function) ~ "High-Affinity Cytochrome",
                           grepl("Low affinity", Raw_Function)  ~ "Low-Affinity Cytochrome",
                           grepl("Complex III", Raw_Function)   ~ "Complex III")) %>%
  group_by(Taxon, Guild) %>% summarise(Completeness = max(Completeness), .groups = "drop") %>%
  left_join(df_tsne_plot, by = "Taxon")

base_theme <- theme_bw(base_size = 9) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9))
create_cyto_subplot <- function(guild_name) {
  ggplot() +
    geom_point(data = df_tsne_plot, aes(tSNE1, tSNE2), color = "grey50", size = 1.2, alpha = 0.5) +
    geom_point(data = filter(df_cyto_final, Guild == guild_name),
               aes(tSNE1, tSNE2, color = Completeness), size = 2.3, alpha = 0.9) +
    scale_color_viridis_c(option = "mako", direction = -1, limits = c(0, 1)) +
    labs(title = guild_name) + base_theme
}
p_c3   <- create_cyto_subplot("Complex III")
p_high <- create_cyto_subplot("High-Affinity Cytochrome")
p_low  <- create_cyto_subplot("Low-Affinity Cytochrome")
for (g in c("Complex III", "High-Affinity Cytochrome", "Low-Affinity Cytochrome"))
  say(sprintf("   %s: carried by %d MAGs", g, sum(df_cyto_final$Guild == g)))

# ---- e. LFC overlay ---------------------------------------------------------

# ggplot draws in row order, so the near-zero greys were landing on top of the
# strongly shifted genomes and hiding them. Draw weakest first, strongest last.
p_lfc <- ggplot(dplyr::arrange(df_tsne_plot, abs(LFC_capped)),
                aes(tSNE1, tSNE2, color = LFC_capped)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_gradient2(low = "#00BFC4", mid = "grey90", high = "#F8766D", midpoint = 0,
                        limits = c(-2, 2), name = "Log Fold\nChange") +
  labs(title = "Differential abundance (LFC)") + base_theme

# ---- f, g, i. binary functional axes ---------------------------------------

axis_defs <- list(
  "Nitrate/TMAO respiration" = c("Nitrogen metabolism: nitrate => nitrite",
                                 "Nitrogen metabolism: nitrite => nitric oxide",
                                 "Other Reductases: TMAO reductase"),
  "Sulfur oxyanion redox" = c("Sulfur metabolism: tetrathionate => thiosulfate",
                              "Sulfur metabolism: Thiosulfate oxidation by SOX complex, thiosulfate => sulfate",
                              "Sulfur metabolism: thiosulfate => sulfite",
                              "Sulfur metabolism: dissimilatory sulfate reduction (and oxidation) sulfate => sulfide"),
  "Butyrate pathway modules" = c("SCFA and alcohol conversions: Butyrate, pt 1",
                                 "SCFA and alcohol conversions: Butyrate, pt 2"))
axis_defs <- discard(map(axis_defs, ~ intersect(.x, cn)), ~ length(.x) == 0)
df_axes_scores <- imap_dfr(axis_defs, function(cols, axis_name)
  tibble(Taxon = rownames(mat_completeness), Axis = axis_name,
         Completeness = score_module(mat_completeness, cols, "max"))) %>%
  left_join(df_tsne_plot, by = "Taxon")

plot_functional_axis <- function(target_axis, plot_title, hex_color = "#00BFA5") {
  sub <- df_axes_scores %>% filter(Axis == target_axis, Completeness > 0) %>% mutate(Status = "Present")
  say(sprintf("   %s: present in %d MAGs", target_axis, nrow(sub)))
  ggplot() +
    geom_point(data = df_tsne_plot, aes(tSNE1, tSNE2), color = "grey50", size = 1.2, alpha = 0.5) +
    geom_point(data = sub, aes(tSNE1, tSNE2, color = Status), size = 2.2, alpha = 0.9) +
    scale_color_manual(values = c("Present" = hex_color), name = "") +
    labs(title = plot_title) + base_theme + theme(legend.position = "none")
}
p_nitrate  <- plot_functional_axis("Nitrate/TMAO respiration", "Nitrate / TMAO")
p_sulfur   <- plot_functional_axis("Sulfur oxyanion redox", "Sulfur (SOX/DSR)")
p_butyrate <- plot_functional_axis("Butyrate pathway modules", "Butyrate")

# ---- h. phylum overlay ------------------------------------------------------

p_phylum <- ggplot(df_tsne_plot, aes(tSNE1, tSNE2, fill = Phylum_Plot)) +
  geom_point(shape = 21, size = 2.2, stroke = 0.2, colour = "grey30", alpha = 0.9) +
  scale_fill_manual(values = phylum_colors, name = "Phylum (GTDB)",
                    labels = function(x) expr_from_md(md_taxon(x))) +
  labs(title = "Phylum") + base_theme +
  theme(legend.position = "right", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))

# ---- j. guild raincloud -----------------------------------------------------

guild_defs <- list(
  "Low-affinity terminal oxidases" = c("Complex IV Low affinity: Cytochrome aa3-600 menaquinol oxidase",
                                       "Complex IV Low affinity: Cytochrome c oxidase",
                                       "Complex IV Low affinity: Cytochrome c oxidase, prokaryotes",
                                       "Complex IV Low affinity: Cytochrome o ubiquinol oxidase"),
  "Cytochrome bd-associated respiration" = c("Complex III: Cytochrome bd ubiquinol oxidase",
                                             "Complex IV High affinity: Cytochrome bd ubiquinol oxidase"),
  "Nitrate/TMAO respiration" = axis_defs[["Nitrate/TMAO respiration"]],
  "Sulfur oxyanion redox"    = axis_defs[["Sulfur oxyanion redox"]],
  "Butyrate pathway modules" = axis_defs[["Butyrate pathway modules"]],
  "Propionate pathway modules" = c("SCFA and alcohol conversions: Propionate, pt 1",
                                   "SCFA and alcohol conversions: Propionate, pt 2"))
guild_defs <- discard(map(guild_defs, ~ intersect(.x, cn)), ~ length(.x) == 0)

guild_long <- imap_dfr(guild_defs, function(cols, g)
  tibble(Taxon = rownames(mat_completeness), Guild = g,
         Score = score_module(mat_completeness, cols, "mean"))) %>%
  left_join(df_tsne_plot %>% select(Taxon, LFC_plot = LFC_capped), by = "Taxon") %>%
  filter(!is.na(LFC_plot)) %>% mutate(In_guild = Score > 0)

stats_res <- guild_long %>% filter(In_guild) %>% group_by(Guild) %>%
  summarise(n_in_guild = n(), n_case = sum(LFC_plot > 0), n_ctrl = sum(LFC_plot < 0),
            median_LFC = median(LFC_plot),
            p_vs_zero = tryCatch(wilcox.test(LFC_plot, mu = 0, conf.int = TRUE)$p.value,
                                 error = function(e) NA_real_), .groups = "drop") %>%
  mutate(q_vs_zero = p.adjust(p_vs_zero, "BH"),
         plot_significance = case_when(q_vs_zero < 0.001 ~ "***", q_vs_zero < 0.01 ~ "**",
                                       q_vs_zero < 0.05 ~ "*", TRUE ~ "ns"),
         Label_Text = paste0(Guild, "\n(n=", n_in_guild, "; +LFC=", n_case, " | -LFC=", n_ctrl, ")"))
write_csv(stats_res, file.path(OUTDIR, "SuppFig2j_Guild_LFC_Stats_n76.csv"))
say("PANEL j: guild-level LFC, one-sample Wilcoxon against zero (BH)")
for (i in seq_len(nrow(stats_res))) with(stats_res[i, ],
  say(sprintf("   %-38s n=%3d (+%3d/-%3d) median LFC %+.2f q=%.3g %s",
              Guild, n_in_guild, n_case, n_ctrl, median_LFC, q_vs_zero, plot_significance)))

plot_guilds <- intersect(c("Low-affinity terminal oxidases", "Cytochrome bd-associated respiration",
                           "Nitrate/TMAO respiration", "Sulfur oxyanion redox",
                           "Butyrate pathway modules", "Propionate pathway modules"),
                         names(guild_defs))
lev <- stats_res %>% filter(Guild %in% plot_guilds) %>%
  mutate(Guild = factor(Guild, levels = plot_guilds)) %>% arrange(Guild) %>% pull(Label_Text)
df_rain <- guild_long %>% filter(In_guild, Guild %in% plot_guilds) %>%
  left_join(stats_res %>% select(Guild, Label_Text, plot_significance), by = "Guild") %>%
  mutate(Label_Text = factor(Label_Text, levels = lev))
y_pad <- diff(range(df_rain$LFC_plot, na.rm = TRUE)) * 0.08
stats_plot <- df_rain %>% group_by(Label_Text, plot_significance) %>%
  summarise(y_pos = max(LFC_plot, na.rm = TRUE) + y_pad, .groups = "drop")

p_rain <- ggplot(df_rain, aes(x = Label_Text, y = LFC_plot)) +
  ggdist::stat_halfeye(adjust = 0.6, width = 0.65, .width = 0, justification = -0.25,
                       fill = "grey78", color = "grey45", alpha = 0.9, point_colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "grey96", color = "grey25",
               linewidth = 0.35, alpha = 0.95) +
  geom_point(position = position_jitter(width = 0.055, height = 0, seed = 42),
             size = 1.1, alpha = 0.6, color = "grey15") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_text(data = stats_plot, aes(x = Label_Text, y = y_pos, label = plot_significance),
            inherit.aes = FALSE, size = 3.2, fontface = "bold", color = "grey10") +
  coord_flip(clip = "off") +
  labs(title = "Guild-level LFC distributions",
       subtitle = "Stars: BH-adjusted one-sample Wilcoxon against zero LFC",
       x = NULL, y = "Log fold change (case vs control)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, color = "grey35"),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

# ---- save panels + composite ------------------------------------------------

pl <- list(SuppFig2b_ComplexIII = p_c3, SuppFig2c_HighAffinity = p_high,
           SuppFig2d_LowAffinity = p_low, SuppFig2e_LFC = p_lfc,
           SuppFig2f_NitrateTMAO = p_nitrate, SuppFig2g_Sulfur = p_sulfur,
           SuppFig2h_Phylum = p_phylum, SuppFig2i_Butyrate = p_butyrate,
           SuppFig2j_Raincloud = p_rain)
for (nm in names(pl))
  ggsave(file.path(OUTDIR, paste0(nm, "_n76.pdf")), pl[[nm]],
         width = if (grepl("Raincloud", nm)) 6.5 else 4.6,
         height = if (grepl("Raincloud", nm)) 6 else 4, device = cairo_pdf, bg = "transparent")

# Two A4 landscape pages. One page cannot hold a 232-column heatmap and nine
# scatter panels without shrinking something past legibility, and stretching
# the heatmap to fill a shared row is what distorted it before.
tag <- theme(plot.tag = element_text(face = "bold", size = 11))
# The t-SNE panels are an embedding: x and y are the same kind of quantity, so
# they have to be drawn square or the cloud shape misleads.
square_tsne <- theme(aspect.ratio = 1)
for (nm in c("p_c3", "p_high", "p_low", "p_lfc", "p_nitrate", "p_sulfur", "p_phylum", "p_butyrate"))
  assign(nm, get(nm) + square_tsne)

small_text <- theme(
  axis.text    = element_text(size = 5.5, colour = "black"),
  axis.title   = element_text(size = 6.5),
  legend.text  = element_text(size = 5.5),
  legend.title = element_text(size = 6.5),
  legend.key.size = unit(0.28, "cm"),
  plot.title    = element_text(face = "bold", size = 7.5),
  plot.subtitle = element_text(size = 5.8, colour = "grey35"))

page1 <- (wrap_elements(full = p_hm) + labs(tag = "a") + tag) +
  plot_annotation(
    title = "Functional and metabolic landscape of disease-associated and healthy commensal guilds",
    theme = theme(plot.title = element_text(face = "bold", size = 9)))

page2 <- (((p_c3 + labs(tag = "b") + tag) | (p_high + labs(tag = "c") + tag) |
           (p_low + labs(tag = "d") + tag)) /
          ((p_lfc + labs(tag = "e") + tag) | (p_nitrate + labs(tag = "f") + tag) |
           (p_sulfur + labs(tag = "g") + tag)) /
          ((p_phylum + labs(tag = "h") + tag) | (p_butyrate + labs(tag = "i") + tag) |
           (p_rain + labs(tag = "j") + tag))) & small_text
# No title on page 2: page 1 carries it, and repeating it reads as a duplicate.
page2 <- page2 + plot_annotation(
    title = NULL, theme = theme(plot.margin = margin(6, 6, 6, 6)))

f1 <- file.path(tempdir(), "supp2_page1.pdf"); f2 <- file.path(tempdir(), "supp2_page2.pdf")
ggsave(f1, page1, width = 8.27, height = 11.69, device = cairo_pdf, bg = "transparent")   # portrait
ggsave(f2, page2, width = 11.69, height = 8.27, device = cairo_pdf, bg = "transparent")   # landscape
out2 <- file.path(OUTDIR, "Supplementary_Figure_4_n76.pdf")
if (nzchar(Sys.which("pdfunite"))) {
  system2("pdfunite", c(shQuote(f1), shQuote(f2), shQuote(out2)))
} else {
  file.copy(f1, out2, overwrite = TRUE)
  warning("pdfunite not found: only page 1 written")
}

writeLines(STATS, file.path(OUTDIR, "SuppFig2_STATS_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "Supplementary_Figure_4_n76.pdf"))
