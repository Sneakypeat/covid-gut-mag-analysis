# Shared plotting setup.
#
# Figures are exported on a transparent canvas so they can be composed in a
# vector editor without a white rectangle behind every panel. theme_bw,
# theme_minimal and theme_void are wrapped rather than replaced so existing
# theme calls inherit it without being edited one by one.
suppressPackageStartupMessages({
  library(ggplot2)
})

theme_transparent <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA))

theme_bw      <- function(...) ggplot2::theme_bw(...)      + theme_transparent
theme_minimal <- function(...) ggplot2::theme_minimal(...) + theme_transparent
theme_void    <- function(...) ggplot2::theme_void(...)    + theme_transparent
ggsave        <- function(..., bg = "transparent") ggplot2::ggsave(..., bg = bg)

# Taxon names are italicised through plotmath rather than ggtext, because
# ggtext elements are dropped by ggh4x::force_panelsizes and cannot be merged
# by patchwork's & operator.
if (file.exists(file.path("R", "common_taxon_italics.R"))) {
  source(file.path("R", "common_taxon_italics.R"), local = FALSE)
}
