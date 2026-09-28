#!/usr/bin/env python3
"""Assemble the supplementary tables for the 76-donor manuscript.

The submitted version shipped six tables. Three of them describe analyses that
survive into this version and are rebuilt here on the 584-MAG catalogue; the
other three covered the GO network, the CAZyme/Pfam layers and CAFCA, all of
which were withdrawn, so their slots are reused rather than renumbered around.
Everything the rebuilt manuscript reports and does not already show in a figure
gets a home here.

Each workbook opens with a Cover sheet naming, for every data sheet, the source
file it was built from and its row count, so a reader can trace any number back
to the file that produced it without reading this script.

Tables that are matrices rather than result tables (per-donor module
completeness, per-taxon exchange fluxes) are written as gzipped CSV alongside
the workbooks instead of being forced into sheets: they are meant to be read by
software, and Excel's row ceiling is not the reason to subset a result.

Usage: build_supp_tables_n76.py [--outdir <dir>]
"""
import argparse
import gzip
import os
import re
import shutil
import sys

import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
N76 = os.path.join(BASE, "result2/n76")

# (sheet name, path relative to result2/n76, one-line description)
TABLES = {
    "S1_Cohort_sequencing_and_robustness": [
        ("Donor_depth", "supp_fig/sources/depth/depth_per_donor_n76.csv",
         "Sequencing depth and mapped reads per donor"),
        ("Sequencing_audit", "supp_fig/tables/SuppFig4_sequencing_audit_n76.csv",
         "Per-donor read audit by arm and source study"),
        ("ZOE_health_index", "ZOE_HealthIndex_PerDonor_n76.csv",
         "Per-donor ZOE gut-microbiome health index"),
        ("ZOE_concordance_per_SGB", "ZOE_Concordance_PerSGB_n76.csv",
         "Per-SGB ZOE health rank against correlation with case status"),
    ],
    "S2_MAG_catalogue_and_novelty": [
        ("MAG_eligibility", "metabolic_support/MAG_eligibility_all584.tsv",
         "All 584 catalogue MAGs: quality, taxonomy, modelling eligibility"),
        ("Novelty", "novel/all_MAGs_novelty_n76.csv",
         "Species-level assignment status for every MAG"),
        ("Novel_candidates", "novel/novel_MAGs_n76.csv",
         "The 20 putative novel species"),
        ("GUNC", "novel/gunc_584_annotated.csv",
         "GUNC chimerism assessment per MAG"),
    ],
    "S3_Differential_abundance": [
        ("ANCOMBC2_all_MAGs", "ANCOMBC2_MAG_results_n76.csv",
         "ANCOM-BC2 at MAG level, all 584 genomes"),
        ("ANCOMBC2_significant", "ANCOMBC2_SIGNIFICANT_with_Taxonomy_n76.csv",
         "The 290 differentially abundant MAGs with taxonomy"),
        ("By_catalogue_source", "significant_by_catalogue_source_n76.csv",
         "Significant MAGs split by the assembly they derive from"),
        ("Raw_vs_rarefied", "supp_fig/sources/depth/ancombc2_raw_vs_rarefied_n76.csv",
         "Same test on rarefied counts, as a depth sensitivity analysis"),
    ],
    "S3_Differential_abundance_PART2": [
        ("ANCOMBC2_phylum", "taxon_level/ancombc2_phylum_level_n76.csv",
         "Summed clade abundance tested at phylum level"),
        ("ANCOMBC2_family", "taxon_level/ancombc2_family_level_n76.csv",
         "Summed clade abundance tested at family level"),
        ("Fisher_phylum", "taxon_level/directional_fisher_phylum_n76.csv",
         "Conditional directional Fisher test, phylum"),
        ("Fisher_family", "taxon_level/directional_fisher_family_n76.csv",
         "Conditional directional Fisher test, family"),
        ("Clade_ORA_phylum", "taxon_level/clade_overrepresentation_phylum_n76.csv",
         "Clade over-representation among differentially abundant MAGs, phylum"),
        ("Clade_ORA_family", "taxon_level/clade_overrepresentation_family_n76.csv",
         "Clade over-representation among differentially abundant MAGs, family"),
        ("Phylum_GLMM_BLUPs", "phylum_glmm_blups_n76.tsv",
         "Phylum random effects from the mixed model"),
        ("Family_GLMM_BLUPs", "family_glmm_blups_n76.tsv",
         "Family random effects from the mixed model"),
    ],
    "S4_Phylogeny_and_functional_guilds": [
        ("Global_signal", "fig2/Table_S3_Global_Phylo_Stats_n76.csv",
         "Moran's I, Abouheif's Cmean, Pagel's lambda, Blomberg's K"),
        ("Local_categorisation", "fig2/Table_S2_Phylogenetic_Categorization_n76.csv",
         "Local Moran's I class per MAG"),
        ("Node_enrichment", "fig2/Table_S5_Node_Enrichment_n76.csv",
         "One-sided Fisher test per internal node"),
        ("Picante_clustering", "fig2/Table_S4_Picante_Clustering_n76.csv",
         "Mean pairwise and nearest-taxon distance SES"),
        ("Quality_phyloglm", "supp_fig/sources/gunc/gunc_phyloglm_models_n76.csv",
         "Phylogenetic logistic regression of enrichment on genome quality"),
        ("Contamination_models", "supp_fig/tables/SuppFig3_contamination_models_n76.csv",
         "Contamination sensitivity variants of that regression"),
    ],
    "S4_Phylogeny_and_functional_guilds_PART2": [
        ("Guild_labels", "fig3/fig3_guild_labels_n76.csv",
         "Guild assignment and size"),
        ("Functional_drivers", "fig3/fig3_functional_drivers_n76.csv",
         "Pathway completeness between case-enriched and case-depleted MAGs"),
        ("Drivers_carriage3", "fig3/fig3_functional_drivers_carriage3_n76.csv",
         "Same test restricted to features carried by at least three MAGs"),
        ("Niche_map", "fig3/fig3_niche_map_n76.csv",
         "Phylum by guild hypergeometric enrichment"),
        ("Cluster_sweep", "fig3/fig3_cluster_sweep_n76.csv",
         "HDBSCAN and UMAP parameter sweep"),
        ("Method_ARI", "bonsai/bonsai_method_ari_n76.csv",
         "Adjusted Rand index between partitioning methods"),
        ("Guild_vs_hierarchical", "supp_fig/tables/SuppFig4_guild_vs_hierarchical_purity_n76.csv",
         "Guild purity against hierarchical clustering"),
        ("Guild_LFC_stats", "supp_fig/tables/SuppFig2j_Guild_LFC_Stats_n76.csv",
         "Log fold change distribution by guild"),
    ],
    "S5_Functional_capacity": [
        ("Module_tests", "fig4/camper_module_tests_n76.csv",
         "Abundance-weighted CAMPER module completeness, per module"),
        ("Class_tests", "fig4/camper_class_tests_n76.csv",
         "Same, pooled by compound class"),
        ("Oxygen_tests", "fig4/camper_oxygen_tests_n76.csv",
         "Same, pooled by oxygen requirement"),
        ("Entero_contribution", "fig4/camper_entero_contribution_n76.csv",
         "Enterobacteriaceae share of the polyphenol capacity increase"),
        ("Carriage_by_direction", "fig4/camper_carriage_by_direction_n76.csv",
         "Module carriage split by ANCOM-BC2 direction"),
        ("Sensitivity_highquality", "fig4/camper_sensitivity_highquality_n76.csv",
         "Restricted to high-quality MAGs"),
        ("Sensitivity_unnormalised", "fig4/camper_sensitivity_unnormalised_n76.csv",
         "Without within-donor renormalisation"),
    ],
    "S5_Functional_capacity_PART2": [
        ("Category_ORA", "fig4/KEGG_category_ORA_n76.csv",
         "Over-representation among differentially abundant KOs, category level"),
        ("Pathway_ORA", "fig4/Fig4c_KEGG_pathway_ORA_q05_n76.csv",
         "The 48 pathways at BH q < 0.05"),
        ("Module_ORA", "fig4/Fig4d_KEGG_module_ORA_n76.csv",
         "Module-level over-representation"),
        ("Butterfly", "fig4/Fig4d_KEGG_butterfly_data_n76.csv",
         "Signed enrichment used for the butterfly panel"),
    ],
    "S6_Modelling_strain_and_replication": [
        ("Genome_traits", "metabolic_support/gut_traits_v2/genome_traits_wide.tsv",
         "Per-genome carbohydrate and fermentation marker calls"),
        ("Trait_combinations", "fig4/Fig4e_trait_upset_tests_n76.csv",
         "Donor share of each observed marker combination"),
        ("Trait_ANCOM_counts", "fig4/Fig4e_trait_upset_ancom_counts_n76.csv",
         "MAG counts per combination by differential-abundance direction"),
        ("Flux_capacity_tests", "fig4/Fig4f_flux_capacity_tests_n76.csv",
         "Abundance-weighted formation capacity per product"),
        ("Flux_capacity_noEntero", "fig4/Fig4f_flux_capacity_noEntero_n76.csv",
         "Same with Enterobacteriaceae removed from the weighting"),
        ("Entero_attribution", "fig4/Fig4f_entero_attribution_n76.csv",
         "Enterobacteriaceae share of predicted capacity per product"),
        ("Community_growth_exchange", "metabolic_support/gapfill_v1/micom_net_exchange_76.csv",
         "MICOM community growth rate and net exchange, 71 donors"),
        ("Substrate_use_tests", "fig4/Fig4f_substrate_use_tests_n76.csv",
         "Community nutrient use and product yield"),
        ("Substrate_by_family", "fig4/Fig4f_substrate_amino_by_family_n76.csv",
         "Family share of community amino acid uptake"),
        ("Capacity_vs_use_stages", "fig4/Fig4f_capacity_vs_use_stage_tests_n76.csv",
         "Product share at the capacity, produced and exported stages"),
        ("Butyrate_route_audit", "metabolic_support/gapfill_v1/butyrate_route_audit.csv",
         "Butyrate-route reactions among the ptb/buk-positive models"),
        ("Gapfill_footprint", "metabolic_support/gapfill_v1/gapfill_footprint.csv",
         "Reactions added per model by gapfilling"),
    ],
    "S6_Modelling_strain_and_replication_PART2": [
        ("Strain_per_genome", "fig2/Fig2b_strain_per_genome_n76.csv",
         "Within-species nucleotide diversity per genome and arm"),
        ("External_effects", "supp_fig/sources/external_validation/analysis/"
                             "external_validation_effects.tsv",
         "Per-cohort and pooled effect sizes for the replicated endpoints"),
        ("External_sample_metrics", "supp_fig/sources/external_validation/analysis/"
                                    "external_validation_sample_metrics.tsv",
         "Per-sample endpoint values in the three external cohorts"),
        ("Severity_trend", "supp_fig/sources/external_validation/analysis/"
                           "PRJEB43555_severity_trend.tsv",
         "Endpoint correlation with reported severity in PRJEB43555"),
        ("IBDMDB_expression", "metabolic_support/expression_validation_v1/ibdmdb_ec_dna_rna.csv",
         "External test: gene abundance against transcript abundance, HMP2"),
        ("IBDMDB_metabolite", "metabolic_support/expression_validation_v1/"
                              "ibdmdb_metabolite_tests.csv",
         "External test: pathway abundance against the measured faecal acid"),
    ],
    "S1_Cohort_sequencing_and_robustness_PART2": [
        ("Rarefaction_stability", "supp_fig/sources/depth/rarefaction_seed_stability_n76.csv",
         "Alpha diversity across rarefaction seeds"),
        ("Alpha_rarefied_by_study", "supp_fig/sources/depth/alpha_controls_rarefied_by_study_n76.csv",
         "Control alpha diversity by source study after rarefaction"),
        ("Imbalance_provenance", "supp_fig/tables/SuppFig3_imbalance_provenance_n76.csv",
         "Catalogue imbalance between case and control assemblies"),
        ("Guild_parameter_sweep", "supp_fig/tables/SuppFig4_guild_parameter_sweep_n76.csv",
         "Guild count across parameter settings and seeds"),
        ("ARI_matrix", "supp_fig/tables/SuppFig5_ari_matrix_n76.csv",
         "Pairwise ARI between all partitioning methods"),
        ("MGE_per_MAG", "mge/mge_per_mag_annotated_n76.tsv",
         "Mobile genetic element content per MAG"),
        ("MGE_phyloglm", "mge/mge_phyloglm_models_n76.csv",
         "Phylogenetic regression of enrichment on MGE content"),
        ("COPTR_replication", "coptr/coptr_mag_tests.tsv",
         "Peak-to-trough replication rate per MAG"),
    ],
}

# Matrices too large to be worth a worksheet. These are deposited in the Zenodo
# archive rather than shipped beside the workbooks, so the builder only records
# where they live; it no longer copies them.
DATA_FILES = [
    ("Supplementary_Data_1_camper_module_completeness.csv.gz",
     "fig4/camper_module_completeness_n76.tsv",
     "Per-genome CAMPER module completeness, 584 MAGs"),
    ("Supplementary_Data_2_camper_donor_capacity.csv.gz",
     "fig4/camper_donor_module_capacity_n76.tsv",
     "Abundance-weighted CAMPER module capacity per donor"),
    ("Supplementary_Data_3_taxon_substrate_use.csv.gz",
     "metabolic_support/gapfill_v1/taxon_substrate_use.tsv.gz",
     "Per-taxon uptake and secretion fluxes from the community solution"),
    ("Supplementary_Data_4_coptr_ptr_long.csv.gz",
     "coptr/coptr_ptr_long.tsv",
     "Peak-to-trough ratios, long format"),
    ("Supplementary_Data_5_micom_members.csv.gz",
     "metabolic_support/gapfill_v1/micom_members_n76.csv",
     "Per-member growth rates in each donor community"),
]


def merge_parts(tables):
    """Fold the _PART2 groups into their parent so eleven groups ship as six.

    The groups are declared separately because each one is a coherent block of
    related sheets and keeping them apart makes the mapping to source files
    readable; they are merged here because a reviewer would rather open six
    workbooks than eleven.
    """
    out = {}
    for name, sheets in tables.items():
        key = name[:-len("_PART2")] if name.endswith("_PART2") else name
        out.setdefault(key, []).extend(sheets)
    return out


def read_any(path):
    sep = "\t" if path.endswith((".tsv", ".tsv.gz")) else ","
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        df = pd.read_csv(fh, sep=sep, low_memory=False)
    return strip_paths(df)


ABS_PATH = re.compile(r"^(?:/|[A-Za-z]:\\)\S*[/\\]")


def strip_paths(df):
    """Reduce absolute paths to basenames.

    Some inputs carry the cluster location of each read file, assembly or bin
    (`read1`, `assembly_dir`, `source_fasta`), which names the account the
    pipeline ran under. The basename keeps the sample or bin identifier, so
    joins still work; the directory is private and is dropped.
    """
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(
                lambda v: os.path.basename(v.rstrip("/\\"))
                if isinstance(v, str) and ABS_PATH.match(v) else v)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(N76, "supp_tables"))
    args = ap.parse_args()
    out = args.outdir
    os.makedirs(out, exist_ok=True)

    manifest = []
    missing = []
    for book, sheets in merge_parts(TABLES).items():
        rows = []
        frames = []
        for sheet, rel, desc in sheets:
            src = os.path.join(N76, rel)
            if not os.path.exists(src):
                missing.append((book, sheet, rel))
                continue
            df = read_any(src)
            frames.append((sheet, df))
            rows.append({"Sheet": sheet, "Contents": desc,
                         "Source file": rel, "Rows": len(df),
                         "Columns": df.shape[1]})
        if not frames:
            print("  SKIP %s (no inputs found)" % book)
            continue
        dest = os.path.join(out, "Supplementary_Table_%s.xlsx" % book)
        cover = pd.DataFrame(rows)
        with pd.ExcelWriter(dest, engine="openpyxl") as xl:
            cover.to_excel(xl, sheet_name="Cover", index=False)
            for sheet, df in frames:
                df.to_excel(xl, sheet_name=sheet[:31], index=False)
        print("  %-38s %d sheets, %6d rows" % (os.path.basename(dest),
                                               len(frames),
                                               sum(len(d) for _, d in frames)))
        manifest.append({"File": os.path.basename(dest), "Sheets": len(frames),
                         "Total rows": sum(len(d) for _, d in frames)})

    pd.DataFrame(manifest).to_csv(os.path.join(out, "MANIFEST.csv"), index=False)
    if missing:
        print("\n  MISSING INPUTS (%d):" % len(missing))
        for book, sheet, rel in missing:
            print("    %-34s %-26s %s" % (book, sheet, rel))
    print("\nWROTE %s" % out)


if __name__ == "__main__":
    main()
