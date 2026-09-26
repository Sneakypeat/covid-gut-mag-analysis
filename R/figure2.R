# Figure 2 and its supplementary figures: MAG catalogue, novelty, phylogenetic structure, strain diversity
#
# Runs every step for this figure in order. Each part is executed in its own
# environment, because the parts were written as standalone scripts and reuse
# short names such as pA and OUTDIR; isolating them keeps that harmless.
#
#   export COVID_MAG_BASE=/path/to/unpacked/zenodo/deposit
#   Rscript R/figure2.R

source("R/config.R")
source("R/common.R")

PARTS <- c(
  "R/parts/fig2/novel_n76.R",
  "R/parts/fig2/gunc_n76.R",
  "R/parts/fig2/pg4_novelty_n76.R",
  "R/parts/fig2/taxon_level_da_n76.R",
  "R/parts/fig2/clade_summary_tests_n76.R",
  "R/parts/fig2/clade_fisher_figure_n76.R",
  "R/parts/fig2/strain_model_refit_n76.R",
  "R/parts/fig2/figure2_n76_stats.R",
  "R/parts/fig2/figure2_n76.R",
  "R/parts/fig2/suppfig3_n76.R",
  "R/parts/fig2/suppfig4_n76.R"
)

for (p in PARTS) {
  if (!file.exists(p)) { message("skipping, not present: ", p); next }
  message("== ", p)
  sys.source(p, envir = new.env(parent = globalenv()))
}
