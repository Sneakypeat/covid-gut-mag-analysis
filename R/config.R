# Paths for every script in this repository.
#
# Nothing here is machine specific. Set COVID_MAG_BASE to wherever you have
# unpacked the Zenodo deposit (10.5281/zenodo.19859581); everything else is
# derived from it. The default is the working directory, so running the figure
# scripts from the unpacked deposit needs no configuration at all.
#
#   export COVID_MAG_BASE=/path/to/unpacked/deposit
#   Rscript R/figure1.R

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "01_inputs")        # analysis inputs from the deposit
OUTDIR <- Sys.getenv("COVID_MAG_OUT", unset = file.path(BASE, "results"))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(INDIR)) {
  stop("Input directory not found: ", INDIR, "\n",
       "Set COVID_MAG_BASE to the directory holding 01_inputs/ from the ",
       "Zenodo deposit 10.5281/zenodo.19859581")
}
