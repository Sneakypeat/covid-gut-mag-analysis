# Genome-resolved metagenomics of the gut microbiome in SARS-CoV-2 infection

Analysis code for a 76-donor, 584-MAG study: 38 SARS-CoV-2 patients against 38
pre-pandemic controls drawn from three independent European cohorts.

**Code only.** No data is tracked here. Inputs, serialised R objects and result
tables are deposited on Zenodo (`10.5281/zenodo.19859581`); raw reads and the
MAG catalogue are in SRA under `PRJNA1457742`.

## Running it

Four entry points organize each main figure and its corresponding supplementary
analyses. Run commands from the repository root.

    export COVID_MAG_BASE=/path/to/unpacked/zenodo/deposit
    Rscript R/figure1.R     # Fig 1 + Supplementary Figures 1 and 2
    Rscript R/figure2.R     # Fig 2 + Supplementary Figures 3 and 8
    Rscript R/figure3.R     # Fig 3 + Supplementary Figures 4 and 5
    Rscript R/figure4.R     # Fig 4 + Supplementary Figures 6 and 7

`R/config.R` defines the shared paths. `COVID_MAG_BASE` points at the unpacked
deposit and `COVID_MAG_OUT` optionally redirects outputs. The scripts
are checked for syntax; they have not been run end to end on a fresh machine,
and some expect the directory layout used during the analysis.

The ZOE health-rank step additionally needs the external reference workbook
`41586_2025_9854_MOESM3_ESM.xlsx`. Place it under `references/` within
`COVID_MAG_BASE`, or set `COVID_MAG_ZOE_XLSX` to its full path.

Each entry point runs its parts in isolated environments, because the parts
were written as standalone scripts and reuse short names such as `pA` and
`OUTDIR`.

## Layout

| Path | Contents |
|---|---|
| `R/figure1.R` … `figure4.R` | the four entry points |
| `R/config.R`, `R/common.R` | paths, and the shared transparent-canvas theme |
| `R/parts/fig1` … `fig4` | the analysis steps each figure runs |
| `python/` | supplementary-table assembly, the KEGG module builder, the genome trait builder, the MICOM collector, and the two HMP2 external-validation tests |
| `metabolic_support/` | the support graph behind the acid-requirement row of Figure 4e, the genome trait calls, and their tests |


## Citation

Manuscript in submission. Cite the Zenodo deposit for the data and this
repository for the code.
