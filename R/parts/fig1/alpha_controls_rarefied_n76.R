# =============================================================================
# Do the three control studies still differ in alpha diversity once sequencing
# depth is equalised?
#
# On raw counts they do (Shannon, Kruskal-Wallis p = 5e-04). Control depth also
# differs sharply (p = 5e-07; medians 57.4M, 40.6M and 21.0M reads) and Shannon
# tracks depth at rho = +0.43, so depth is a candidate explanation. It is not a
# clean one: PRJEB7331 has the second highest depth and the lowest Shannon.
#
# This rarefies all 76 donors to the shallowest library and repeats the test
# across 25 seeds, so the answer does not rest on one random draw. Both the
# among-control test and the case-control test are reported, because an index
# that erases control heterogeneity by also erasing the case signal is useless.
#
# Outputs -> result2/n76/supp_fig/depth/
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(dplyr); library(tidyr); library(readr)
})

BASE   <- Sys.getenv("COVID_MAG_BASE", unset = getwd())
INDIR  <- file.path(BASE, "result2/inputfiles/new")
N76    <- file.path(BASE, "result2/n76")
OUTDIR <- file.path(N76, "supp_fig/depth")
N_SEED <- 25
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
STATS <- character(); say <- function(s) { message(s); STATS <<- c(STATS, s) }

ps   <- readRDS(file.path(N76, "ps_mags_n76.rds"))
meta <- read.delim(file.path(INDIR, "sample_metadata_76.tsv"), check.names = FALSE)
otu  <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
M <- t(otu)
rownames(meta) <- meta$sample; meta <- meta[rownames(M), ]
stopifnot(nrow(M) == 76, all(rownames(M) == meta$sample))

depth <- rowSums(M)
ctrl  <- meta$group == "Control"
src   <- factor(meta$source[ctrl])
grp   <- factor(meta$group, levels = c("Control", "Case"))

cliff <- function(a, b) { m <- outer(a, b, "-"); (sum(m > 0) - sum(m < 0)) / length(m) }
indices <- function(X) {
  S  <- specnumber(X); H <- diversity(X, "shannon")
  list(Observed = S, Chao1 = estimateR(round(X))["S.chao1", ], Shannon = H,
       Pielou = H / log(S), Simpson = diversity(X, "simpson"),
       InvSimpson = diversity(X, "invsimpson"))
}

say(sprintf("donors %d (case %d, control %d) | MAGs %d", nrow(M), sum(!ctrl), sum(ctrl), ncol(M)))
say(sprintf("depth median %.1fM, range %.1fM to %.1fM; rarefying every donor to %s reads",
            median(depth) / 1e6, min(depth) / 1e6, max(depth) / 1e6,
            format(min(depth), big.mark = ",")))
say(sprintf("control depth by study (median reads): %s",
            paste(sprintf("%s %.1fM", levels(src),
                          tapply(depth[ctrl], src, median) / 1e6), collapse = " | ")))
say(sprintf("among-control depth: Kruskal-Wallis p = %.3g", kruskal.test(depth[ctrl] ~ src)$p.value))

# ---- raw counts, for the comparison the rarefied result has to beat ----------
raw <- indices(M)
raw_rows <- lapply(names(raw), function(nm) {
  x <- raw[[nm]]
  tibble(index = nm, data = "raw", seed = NA_integer_,
         among_control_p = kruskal.test(x[ctrl] ~ src)$p.value,
         case_control_p  = wilcox.test(x ~ grp, exact = FALSE)$p.value,
         case_control_delta = cliff(x[!ctrl], x[ctrl]),
         depth_rho = suppressWarnings(cor(x, depth, method = "spearman")))
}) %>% bind_rows()

# ---- rarefied counts, 25 independent draws ----------------------------------
scale_to <- min(depth)
rar_rows <- lapply(seq_len(N_SEED), function(s) {
  set.seed(s)
  R <- rrarefy(M, sample = scale_to)
  stopifnot(all(abs(rowSums(R) - scale_to) < 1e-6))
  ix <- indices(R)
  lapply(names(ix), function(nm) {
    x <- ix[[nm]]
    tibble(index = nm, data = "rarefied", seed = s,
           among_control_p = kruskal.test(x[ctrl] ~ src)$p.value,
           case_control_p  = wilcox.test(x ~ grp, exact = FALSE)$p.value,
           case_control_delta = cliff(x[!ctrl], x[ctrl]),
           depth_rho = suppressWarnings(cor(x, depth, method = "spearman")))
  }) %>% bind_rows()
}) %>% bind_rows()

all_rows <- bind_rows(raw_rows, rar_rows)
summary_tbl <- all_rows %>%
  group_by(index, data) %>%
  summarise(seeds = sum(!is.na(seed)) + as.integer(all(is.na(seed))),
            among_p_median = median(among_control_p),
            among_p_min = min(among_control_p), among_p_max = max(among_control_p),
            among_sig_frac = mean(among_control_p < 0.05),
            case_p_median = median(case_control_p),
            case_delta_median = median(case_control_delta), .groups = "drop") %>%
  arrange(index, desc(data))

say("")
say("## among-control alpha diversity, raw versus depth-equalised")
say(sprintf("   %-11s %-9s %-26s %-24s %s", "index", "data",
            "among-control p", "case vs control", "seeds sig / n"))
for (i in seq_len(nrow(summary_tbl))) with(summary_tbl[i, ], {
  rng <- if (data == "rarefied") sprintf("median %.4f [%.4f, %.4f]", among_p_median, among_p_min, among_p_max)
         else sprintf("%.4f", among_p_median)
  say(sprintf("   %-11s %-9s %-26s delta %+.2f, p = %-8.1e %s",
              index, data, rng, case_delta_median, case_p_median,
              ifelse(data == "rarefied", sprintf("%d/%d", round(among_sig_frac * N_SEED), N_SEED), "")))
})

# Per-study medians on one representative draw, for the manuscript sentence.
set.seed(1); R1 <- rrarefy(M, sample = scale_to)
sh <- diversity(R1, "shannon")
per_study <- tibble(study = c(levels(src), "Cases"),
                    n = c(as.integer(table(src)), sum(!ctrl)),
                    median_shannon = c(tapply(sh[ctrl], src, median), median(sh[!ctrl])))
say("")
say("## rarefied Shannon by cohort (seed 1)")
for (i in seq_len(nrow(per_study))) with(per_study[i, ],
  say(sprintf("   %-11s n = %2d  median Shannon %.2f", study, n, median_shannon)))
say(sprintf("   largest gap between control studies %.2f | lowest control study minus cases %.2f",
            diff(range(per_study$median_shannon[per_study$study != "Cases"])),
            min(per_study$median_shannon[per_study$study != "Cases"]) -
              per_study$median_shannon[per_study$study == "Cases"]))

write_csv(all_rows, file.path(OUTDIR, "alpha_controls_rarefied_all_n76.csv"))
write_csv(summary_tbl, file.path(OUTDIR, "alpha_controls_rarefied_summary_n76.csv"))
write_csv(per_study, file.path(OUTDIR, "alpha_controls_rarefied_by_study_n76.csv"))
writeLines(STATS, file.path(OUTDIR, "ALPHA_CONTROLS_RAREFIED_n76.txt"))
message("DONE -> ", file.path(OUTDIR, "ALPHA_CONTROLS_RAREFIED_n76.txt"))
