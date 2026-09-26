# Figure 4 and its supplementary figures: polyphenol capacity, KEGG shifts, traits and predicted flux
#
# Runs every step for this figure in order. Each part is executed in its own
# environment, because the parts were written as standalone scripts and reuse
# short names such as pA and OUTDIR; isolating them keeps that harmless.
#
#   export COVID_MAG_BASE=/path/to/unpacked/zenodo/deposit
#   Rscript R/figure4.R

source("R/config.R")
source("R/common.R")

PARTS <- c(
  "R/parts/fig4/camper_n76.R",
  "R/parts/fig4/camper_n76_checks.R",
  "R/parts/fig4/camper_n76_followup.R",
  "R/parts/fig4/camper_n76_taxonomy_adjustment.R",
  "R/parts/fig4/kegg_summary_n76.R",
  "R/parts/fig4/kegg_butterfly_n76.R",
  "R/parts/fig4/kegg_entero_adjust_n76.R",
  "R/parts/fig4/mge_n76.R",
  "R/parts/fig4/fig4e_trait_capacity_n76.R",
  "R/parts/fig4/fig4e_trait_upset_n76.R",
  "R/parts/fig4/fig4f_flux_capacity_n76.R",
  "R/parts/fig4/fig4f_micom_substrate_n76.R",
  "R/parts/fig4/fig4f_capacity_vs_use_n76.R",
  "R/parts/fig4/suppfig6_n76.R",
  "R/parts/fig4/suppfig_flux_n76.R",
  "R/parts/fig4/fig4_square_panels_n76.R"
)

for (p in PARTS) {
  if (!file.exists(p)) { message("skipping, not present: ", p); next }
  message("== ", p)
  sys.source(p, envir = new.env(parent = globalenv()))
}
