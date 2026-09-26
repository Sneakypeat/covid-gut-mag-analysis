# =============================================================================
# KEGG category and module ORA (module result drawn as Figure 4d by
# fig4_square_panels_n76.R; the category result is a table only)
#
#   c  KEGG pathway ORA summarised into functional categories (KEGG BRITE
#      br08901 level 2). Each category is tested as one gene set (the union of
#      the KOs of its pathways) with the SAME hypergeometric test as the pathway
#      ORA, so the summary is a test, not an average of p-values.
#   d  KEGG module ORA (KEGG modules, M-numbers), same test, grouped by the
#      module hierarchy (br:ko00002).
#   Both drawn as single-column heatmaps in the style of Figure 3c.
#
#   The raw-p < 0.05 butterfly (80 pathways) that sat in Supplementary Fig 6
#   was dropped: Figure 4c already shows every pathway at BH q < 0.05.
#
# Test, reproduced exactly from KEGG_KO_hypergeometric_ORA_n76.csv (checked on
# "Metabolic pathways", both directions): universe N = tested KOs carrying at
# least one annotation of the kind being tested; K = KOs higher in that group
# (LinDA) within the universe; M = set size; k = overlap;
# p = phyper(k - 1, M, N - M, K, lower.tail = FALSE); BH within direction.
#
# Categories removed from the pathway summary, stated rather than hidden:
#   - "Global and overview maps" (Metabolic pathways, Microbial metabolism in
#     diverse environments, ...): catch-alls that overlap every other category.
#   - "Organismal Systems" and "Human Diseases", except the bacterial
#     "Drug resistance: antimicrobial" subclass: eukaryotic physiology, not
#     interpretable for bacterial genes.
#
# Outputs -> result2/n76/fig4/ and result2/n76/supp_fig/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(tibble)
  library(jsonlite); library(ggplot2); library(forcats); library(patchwork)
})

theme_transparent <- ggplot2::theme(
  plot.background       = element_rect(fill = "transparent", color = NA),
  panel.background      = element_rect(fill = "transparent", color = NA),
  legend.background     = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA),
  legend.key            = element_rect(fill = "transparent", color = NA))
theme_bw      <- function(...) ggplot2::theme_bw(...)      + theme_transparent
theme_minimal <- function(...) ggplot2::theme_minimal(...) + theme_transparent
ggplot2::theme_set(ggplot2::theme_get() + theme_transparent)
ggsave <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

BASE  <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
KF    <- file.path(BASE, "result2/n76/fig4/kegg_fast")
REF   <- file.path(BASE, "result2/n76/fig4/gene_catalogue/kegg_reference_20260919")
FIG4  <- file.path(BASE, "result2/n76/fig4")
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

# ---- inputs -----------------------------------------------------------------

linda <- read_csv(file.path(KF, "linda_results/LinDA_KO_all_results_n76.csv"), show_col_types = FALSE) %>%
  transmute(KO = Feature, Status)
tested  <- unique(linda$KO)
up_all  <- linda$KO[linda$Status == "Higher in Case"]
dn_all  <- linda$KO[linda$Status == "Higher in Control"]
say(sprintf("LinDA KOs tested %d | higher in Case %d | higher in Control %d",
            length(tested), length(up_all), length(dn_all)))

ko_pw <- read_tsv(file.path(BASE, "result2/n76/fig4/gene_catalogue/KEGG_KO_pathway_map_20260917.tsv"),
                  show_col_types = FALSE) %>%
  transmute(KO = KO_ID, pid = sub("^ko", "", pathway_id), pathway_name)

ko_md <- read_tsv(file.path(REF, "kegg_link_ko_module.tsv"), col_names = c("KO", "module"),
                  show_col_types = FALSE) %>%
  mutate(KO = sub("^ko:", "", KO), module = sub("^md:", "", module))
md_names <- read_tsv(file.path(REF, "kegg_list_module.tsv"), col_names = c("module", "module_name"),
                     show_col_types = FALSE) %>% mutate(module = sub("^md:", "", module))

# pathway -> level1 / level2 from BRITE br08901
br <- fromJSON(file.path(REF, "br08901_pathway_hierarchy.json"), simplifyVector = FALSE)
pw_class <- map_dfr(br$children, function(l1)
  map_dfr(l1$children, function(l2)
    map_dfr(l2$children, function(p)
      tibble(level1 = l1$name, level2 = l2$name,
             pid = sub("^(\\d{5}).*", "\\1", p$name)))))

# module -> category from br:ko00002 (level 2 = e.g. "Carbohydrate metabolism")
mh <- fromJSON(file.path(REF, "ko00002_module_hierarchy.json"), simplifyVector = FALSE)
md_class <- map_dfr(mh$children, function(l1)
  map_dfr(l1$children, function(l2)
    map_dfr(l2$children, function(l3)
      map_dfr(l3$children, function(m)
        tibble(module = sub("^(M\\d{5}).*", "\\1", m$name),
               mclass = l2$name, msub = l3$name)))))
say(sprintf("KEGG reference: %d pathways classified, %d modules classified",
            nrow(pw_class), n_distinct(md_class$module)))

# ---- the test, identical to the pathway ORA ---------------------------------

ora <- function(set_df, set_col) {
  univ <- intersect(tested, set_df$KO)
  up <- intersect(up_all, univ); dn <- intersect(dn_all, univ)
  N <- length(univ)
  set_df %>% filter(KO %in% univ) %>% group_by(.data[[set_col]]) %>%
    summarise(M = n_distinct(KO), k_case = sum(unique(KO) %in% up),
              k_control = sum(unique(KO) %in% dn), .groups = "drop") %>%
    mutate(N = N, K_case = length(up), K_control = length(dn),
           p_case    = phyper(k_case - 1,    M, N - M, K_case,    lower.tail = FALSE),
           p_control = phyper(k_control - 1, M, N - M, K_control, lower.tail = FALSE),
           q_case = p.adjust(p_case, "BH"), q_control = p.adjust(p_control, "BH"),
           dir = ifelse(q_case <= q_control, "Case", "Control"),
           q = pmin(q_case, q_control),
           signed = ifelse(dir == "Case", 1, -1) * -log10(q))
}

# sanity: the pathway-level test reproduces the table already on disk
pw_chk <- ora(ko_pw, "pid")
ref_tab <- read_csv(file.path(KF, "linda_results/KEGG_KO_hypergeometric_ORA_n76.csv"),
                    show_col_types = FALSE) %>% mutate(pid = sub("^ko", "", pathway_id))
chk <- pw_chk %>% inner_join(ref_tab %>% select(pid, p_case_ref = p_case), by = "pid")
say(sprintf("pathway ORA reproduction: %d pathways matched, max |log10 p difference| = %.2g",
            nrow(chk), max(abs(log10(chk$p_case) - log10(chk$p_case_ref)), na.rm = TRUE)))

# ---- c. pathway ORA summarised into functional categories -------------------

keep_cls <- pw_class %>%
  filter(level2 != "Global and overview maps",
         !level1 %in% c("Organismal Systems", "Human Diseases", "Drug Development") |
           level2 == "Drug resistance: antimicrobial")
ko_cls <- ko_pw %>% inner_join(keep_cls, by = "pid") %>% distinct(KO, level1, level2)
cls <- ora(ko_cls, "level2") %>% left_join(distinct(keep_cls, level1, level2), by = "level2")

# how many of each category's pathways were individually significant
pw_sig <- ref_tab %>% mutate(q = pmin(q_case, q_control),
                             dir = ifelse(q_case <= q_control, "Case", "Control")) %>%
  filter(q < 0.05) %>% inner_join(keep_cls, by = "pid") %>%
  count(level2, dir) %>% pivot_wider(names_from = dir, values_from = n, values_fill = 0)
cls <- cls %>% left_join(pw_sig, by = "level2") %>%
  mutate(across(c(Case, Control), ~ replace_na(.x, 0)),
         lab = sprintf("%d up / %d down", Case, Control))
write_csv(cls, file.path(FIG4, "KEGG_category_ORA_n76.csv"))
cls_sig <- cls %>% filter(q < 0.05)
say(""); say(sprintf("## c. functional categories: %d tested, %d at BH q < 0.05 (%d Case, %d Control)",
                     nrow(cls), nrow(cls_sig), sum(cls_sig$dir == "Case"), sum(cls_sig$dir == "Control")))
for (i in seq_len(nrow(cls_sig))) with(arrange(cls_sig, desc(signed))[i, ],
  say(sprintf("   %-48s %-7s q = %.2g  (%s pathways)", level2, dir, q, lab)))

# Single tile column sorted by signed -log10(FDR), exactly the Fig 3c layout.
# An earlier version faceted by KEGG level; the long strip labels consumed the
# whole panel width and the tiles did not render. The category of every row is
# kept in the CSV rather than drawn.
one_d <- function(df, row_col, title, subtitle, cap = NULL) {
  lim <- if (is.null(cap)) max(abs(df$signed)) else cap
  df %>% mutate(row = fct_reorder(.data[[row_col]], signed)) %>%
    ggplot(aes(x = "Shift", y = row)) +
    geom_tile(aes(fill = signed), colour = "white", linewidth = 0.4, width = 0.55) +
    scale_x_discrete(expand = expansion(mult = 0.6)) +
    scale_fill_gradient2(low = "#0f4480", mid = "grey98", high = "#b35806", midpoint = 0,
                         limits = c(-lim, lim), oob = scales::squish,
                         name = if (is.null(cap)) "Signed\n-log10(FDR)" else
                           sprintf("Signed\n-log10(FDR)\ncapped at %g", cap)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6.5, colour = "black"), panel.grid = element_blank(),
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6),
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"),
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 6.5, colour = "grey35"))
}

pC <- one_d(cls_sig, "level2",
            "KEGG functional categories",
            "Pathway ORA pooled by category, BH q < 0.05\nOrange: Case, blue: Control",
            cap = 10)

# ---- d. module ORA -----------------------------------------------------------

# Cut at the last space before the limit. Cutting names at their first comma, as
# an earlier version did, made M00651 and M00652 both read "Vancomycin
# resistance" (they are the D-Ala-D-Lac and D-Ala-D-Ser types).
shorten <- function(x, n) ifelse(nchar(x) <= n, x,
  paste0(sub("\\s+\\S*$", "", substr(x, 1, n)), "\u2026"))

md <- ora(ko_md, "module") %>% left_join(md_names, by = "module") %>%
  left_join(distinct(md_class, module, .keep_all = TRUE), by = "module") %>%
  mutate(module_name = gsub("\\s*\\[[^]]*\\]", "", module_name),
         short = sprintf("%s %s", module, shorten(module_name, 50)),
         mclass = replace_na(mclass, "Unclassified"))
write_csv(md, file.path(FIG4, "Fig4d_KEGG_module_ORA_n76.csv"))
md_sig <- md %>% filter(q < 0.05)
say(""); say(sprintf("## d. KEGG modules: %d tested, %d at BH q < 0.05 (%d Case, %d Control)",
                     nrow(md), nrow(md_sig), sum(md_sig$dir == "Case"), sum(md_sig$dir == "Control")))
for (i in seq_len(min(20, nrow(md_sig)))) with(arrange(md_sig, desc(abs(signed)))[i, ],
  say(sprintf("   %-50s %-7s q = %.2g  [%s]", substr(short, 1, 50), dir, q, mclass)))

pD <- one_d(md_sig, "short",
            "KEGG modules",
            "Module ORA, same test, BH q < 0.05\nOrange: Case, blue: Control")

# Panel PDFs for c and d are drawn by fig4_square_panels_n76.R (square tiles).
# Writing them here as well used to overwrite those with stretched slabs.
saveRDS(list(panel_c = pC, panel_d = pD, n_c = nrow(cls_sig), n_d = nrow(md_sig)),
        file.path(FIG4, "Fig4cd_KEGG_summary_panels_n76.rds"))

writeLines(STATS, file.path(FIG4, "KEGG_SUMMARY_STATS_n76.txt"))
message("DONE")
