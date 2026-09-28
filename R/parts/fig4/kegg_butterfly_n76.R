# =============================================================================
# KEGG pathway butterfly plot, n = 76
# Pathways are selected at BH q < 0.05 and encoded by signed -log10(FDR).
# Outputs -> result2/n76/fig4/
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(ggnewscale); library(forcats)
})
theme_transparent <- ggplot2::theme(
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA),
  legend.key = element_rect(fill = "transparent", color = NA))

BASE <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
OUT  <- file.path(BASE, "result2/n76/fig4")
STATS <- c(); say <- function(...) { l <- paste0(...); message(l); STATS <<- c(STATS, l) }

ora <- read_csv(file.path(BASE, "result2/n76/fig4/kegg_fast/linda_results/KEGG_KO_hypergeometric_ORA_n76.csv"),
                show_col_types = FALSE)
say(sprintf("pathways tested: %d", nrow(ora)))

# one row per pathway, taking whichever direction is the stronger signal
bf <- ora %>%
  mutate(dir = ifelse(q_case <= q_control, "Case", "Control"),
         q   = pmin(q_case, q_control),
         k   = ifelse(dir == "Case", k_case, k_control),
         K   = ifelse(dir == "Case", K_case, K_control),
         rich = k / K,
         signed = ifelse(dir == "Case", 1, -1) * -log10(q)) %>%
  filter(q < 0.05, is.finite(signed))
say(sprintf("FDR-significant (BH q < 0.05): %d pathways (%d case, %d control)",
            nrow(bf), sum(bf$dir == "Case"), sum(bf$dir == "Control")))
say(sprintf("rich factor range: %.2f-%.2f | |signed -log10(q)| max %.1f",
            min(bf$rich), max(bf$rich), max(abs(bf$signed))))

bf <- bf %>% mutate(pathway_name = fct_reorder(pathway_name, signed))
write_csv(bf, file.path(OUT, "Fig4d_KEGG_butterfly_data_n76.csv"))

up <- bf %>% filter(dir == "Case"); dn <- bf %>% filter(dir == "Control")
p_bf <- ggplot(bf, aes(signed, pathway_name)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = c(-log10(0.05), log10(0.05)), linetype = "dotted", colour = "grey60") +
  geom_segment(aes(x = 0, xend = signed, yend = pathway_name), colour = "grey75", linewidth = 0.3) +
  geom_point(data = up, aes(size = rich, colour = q)) +
  scale_colour_gradient(low = "#B2182B", high = "#FDBF6F", trans = "log10", name = "FDR (Case)") +
  new_scale_colour() +
  geom_point(data = dn, aes(size = rich, colour = q)) +
  scale_colour_gradient(low = "#2166AC", high = "#C6DBEF", trans = "log10", name = "FDR (Control)") +
  scale_size_continuous(range = c(1.4, 5), name = "Rich factor") +
  labs(title = "Pathway Enrichment (Butterfly Plot)",
       subtitle = sprintf("Left: enriched in Control | Right: enriched in Case. %d pathways at BH q < 0.05.", nrow(bf)),
       x = expression(Signed~-log[10](FDR)), y = NULL) +
  theme_bw(base_size = 9) + theme_transparent +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        axis.text.y = element_text(size = 6), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5))

h <- max(5, 0.16 * nrow(bf))
ggsave(file.path(OUT, "Fig4d_KEGG_Butterfly_n76.pdf"), p_bf, width = 8.5, height = h,
       device = cairo_pdf, bg = "transparent", limitsize = FALSE)
saveRDS(p_bf, file.path(OUT, "Fig4d_KEGG_butterfly_plot_n76.rds"))
writeLines(STATS, file.path(OUT, "KEGG_BUTTERFLY_STATS_n76.txt"))
message("DONE -> ", file.path(OUT, "Fig4d_KEGG_Butterfly_n76.pdf"))
