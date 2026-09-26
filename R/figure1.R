# Figure 1 and its supplementary figures: cohort, diversity, health index, external replication
#
# Runs every step for this figure in order. Each part is executed in its own
# environment, because the parts were written as standalone scripts and reuse
# short names such as pA and OUTDIR; isolating them keeps that harmless.
#
#   export COVID_MAG_BASE=/path/to/unpacked/zenodo/deposit
#   Rscript R/figure1.R

source("R/config.R")
source("R/common.R")

PARTS <- c(
  "R/parts/fig1/figure1_n76.R",
  "R/parts/fig1/zoe_n76.R",
  "R/parts/fig1/prevalence_legend_n76.R",
  "R/parts/fig1/alpha_controls_rarefied_n76.R",
  "R/parts/fig1/depth_rarefaction_n76.R",
  "R/parts/fig1/depth_rarefaction_seeds_n76.R",
  "R/parts/fig1/figure1_compose_n76.R",
  "R/parts/fig1/suppfig1_n76.R",
  "R/parts/fig1/suppfig7_n76.R"
)

for (p in PARTS) {
  if (!file.exists(p)) { message("skipping, not present: ", p); next }
  message("== ", p)
  sys.source(p, envir = new.env(parent = globalenv()))
}
