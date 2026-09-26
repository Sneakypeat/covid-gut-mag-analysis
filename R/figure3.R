# Figure 3 and its supplementary figures: functional guilds and clustering robustness
#
# Runs every step for this figure in order. Each part is executed in its own
# environment, because the parts were written as standalone scripts and reuse
# short names such as pA and OUTDIR; isolating them keeps that harmless.
#
#   export COVID_MAG_BASE=/path/to/unpacked/zenodo/deposit
#   Rscript R/figure3.R

source("R/config.R")
source("R/common.R")

PARTS <- c(
  "R/parts/fig3/figure3_n76_tuning.R",
  "R/parts/fig3/figure3_n76_cluster_decision.R",
  "R/parts/fig3/figure3_n76.R",
  "R/parts/fig3/bonsai_n76.R",
  "R/parts/fig3/suppfig2_n76.R",
  "R/parts/fig3/suppfig5_n76.R"
)

for (p in PARTS) {
  if (!file.exists(p)) { message("skipping, not present: ", p); next }
  message("== ", p)
  sys.source(p, envir = new.env(parent = globalenv()))
}
