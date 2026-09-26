#!/usr/bin/env python3
"""Does gene-level capacity predict the measured metabolite?

`ibdmdb_capacity_expression.py` closed the gene-to-transcript step. This closes
the remaining one, against a measured molecule rather than another sequencing
assay: in the same stool, does the DNA abundance of the butyrate enzymes track
the faecal butyrate concentration.

The comparison that decides what our figures can claim is DNA against RNA as
predictors of the same metabolite. We have DNA for our cohort and will never
have RNA for it. If RNA predicts butyrate far better than DNA does, our capacity
measure is a weak stand-in and the manuscript has to say so. If the two are
close, the missing RNA costs us little.

Metabolites come from the HMP2 untargeted LC-MS panel. Only named features are
used, and butyrate, propionate, valerate/isovalerate and caproate all come from
HILIC-neg with pooled-QC coefficients of variation between 5 and 7 percent.
Acetate is not in the named panel: it is too small and polar for these methods,
so no acetate claim is made here.

Feature intensities are not absolute concentrations, so every test is Spearman
on ranks within the cohort, never a slope or an absolute level.

Usage: ibdmdb_metabolite_test.py <ibdmdb_dir> <outdir>
"""
import csv
import gzip
import json
import os
import sys

import numpy as np
from scipy import stats

from ibdmdb_capacity_expression import read_unstratified, load_metadata

# Same enzymes as the expression test, grouped by the acid they make.
PATHWAY = {
    "butyrate": ["2.8.3.8", "2.7.2.7", "2.3.1.19", "1.3.8.1", "1.1.1.157", "4.2.1.17"],
    "propionate": ["5.4.99.2"],
}
# Named LC-MS features. Keys are the metabolite name as the panel reports it.
TARGETS = ["butyrate", "propionate", "valerate/isovalerate", "caproate"]
MIN_PAIR = 30


def load_mbx(path):
    """Named metabolite features -> (sample ids, {name: intensity vector}).

    The biom is sparse triplets. 81,867 features are mostly unnamed mass
    peaks; only the 592 with a `Metabolite` annotation are interpretable, and
    only those are kept.
    """
    d = json.load(gzip.open(path, "rt"))
    rows, cols = d["rows"], d["columns"]
    samples = [c["id"] for c in cols]
    wanted = {}
    for i, r in enumerate(rows):
        md = r.get("metadata") or {}
        name = md.get("Metabolite")
        if name in TARGETS and name not in wanted:
            wanted[i] = (name, md.get("Method"), md.get("Pooled QC sample CV"))
    mat = {i: np.zeros(len(samples)) for i in wanted}
    for entry in d["data"]:
        i, j, v = entry[0], entry[1], entry[2]
        if i in mat:
            mat[i][j] = v
    out = {wanted[i][0]: mat[i] for i in wanted}
    qc = {wanted[i][0]: (wanted[i][1], wanted[i][2]) for i in wanted}
    return samples, out, qc


def collection_index(samples, meta_by_type):
    """sample column index -> stool collection, for one assay."""
    idx = {}
    for j, s in enumerate(samples):
        if s in meta_by_type:
            idx.setdefault(meta_by_type[s][0], j)
    return idx


def pathway_score(tab, ec_list, col):
    """Summed relative abundance of a pathway's enzymes in one sample column.

    A sum, not a mean, because a community carrying several of the enzymes has
    more of the machinery than one carrying a single copy of one, and the
    missing-enzyme case must not be rewarded by averaging over fewer terms.
    """
    present = [ec for ec in ec_list if ec in tab]
    if not present:
        return None, 0
    return np.sum([tab[ec][col] for ec in present], axis=0), len(present)


def spear(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < MIN_PAIR:
        return None
    if np.all(a[ok] == a[ok][0]) or np.all(b[ok] == b[ok][0]):
        return None
    r = stats.spearmanr(a[ok], b[ok])
    return {"rho": float(r.statistic), "p": float(r.pvalue), "n": int(ok.sum())}


def main():
    src, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    mgx_meta, mtx_meta, coll_meta = load_metadata(os.path.join(src, "hmp2_metadata.csv"))
    mbx_meta = {}
    with open(os.path.join(src, "hmp2_metadata.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            if r["data_type"] == "metabolomics":
                mbx_meta[r["External ID"]] = (r["site_sub_coll"], r.get("diagnosis", ""),
                                              r["Participant ID"])

    say("loading tables")
    mgx_samples, dna = read_unstratified(os.path.join(src, "mgx_ecs_relab.tsv.gz"))
    mtx_samples, rna = read_unstratified(os.path.join(src, "mtx_ecs_relab.tsv.gz"))
    mbx_samples, mbx, qc = load_mbx(os.path.join(src, "mbx.biom.gz"))
    say(f"  DNA {len(dna)} ECs x {len(mgx_samples)}")
    say(f"  RNA {len(rna)} ECs x {len(mtx_samples)}")
    say(f"  MBX {len(mbx)} named SCFA features x {len(mbx_samples)}")
    for m, (meth, cv) in sorted(qc.items()):
        say(f"    {m:<24} {meth:<10} pooled-QC CV {cv}")

    d_idx = collection_index(mgx_samples, mgx_meta)
    r_idx = collection_index(mtx_samples, mtx_meta)
    m_idx = collection_index(mbx_samples, mbx_meta)

    dm = sorted(set(d_idx) & set(m_idx))
    rm = sorted(set(r_idx) & set(m_idx))
    say()
    say(f"collections with DNA + metabolome: {len(dm)}")
    say(f"collections with RNA + metabolome: {len(rm)}")

    rows = []
    for label, keep in (("all", None), ("nonIBD", {"nonIBD"}), ("IBD", {"CD", "UC"})):
        dm_k = [c for c in dm if keep is None or coll_meta.get(c, ("", "", ""))[1] in keep]
        rm_k = [c for c in rm if keep is None or coll_meta.get(c, ("", "", ""))[1] in keep]
        if len(dm_k) < MIN_PAIR:
            continue
        say()
        say(f"## {label}: {len(dm_k)} DNA+MBX, {len(rm_k)} RNA+MBX")
        say(f"   {'pathway':<12} {'metabolite':<22} {'assay':<5} {'n':>4} {'rho':>7} {'p':>10}")
        for path, ecs in PATHWAY.items():
            for target in TARGETS:
                if target not in mbx:
                    continue
                for assay, tab, colls, idx in (("DNA", dna, dm_k, d_idx),
                                               ("RNA", rna, rm_k, r_idx)):
                    if len(colls) < MIN_PAIR:
                        continue
                    cols = np.array([idx[c] for c in colls])
                    score, n_ec = pathway_score(tab, ecs, cols)
                    if score is None:
                        continue
                    met = mbx[target][np.array([m_idx[c] for c in colls])]
                    s = spear(score, met)
                    if s is None:
                        continue
                    rows.append({"stratum": label, "pathway": path, "n_ec": n_ec,
                                 "metabolite": target, "assay": assay,
                                 "n": s["n"], "rho": round(s["rho"], 4), "p": s["p"]})
                    say(f"   {path:<12} {target:<22} {assay:<5} {s['n']:>4} "
                        f"{s['rho']:>+7.3f} {s['p']:>10.2e}")

    # Per-enzyme detail on butyrate, so a pathway sum cannot hide one driver.
    say()
    say("## per-enzyme vs measured butyrate (all collections)")
    say(f"   {'EC':<10} {'assay':<5} {'n':>4} {'rho':>7} {'p':>10}")
    met_all = mbx["butyrate"][np.array([m_idx[c] for c in dm])]
    for ec in PATHWAY["butyrate"]:
        for assay, tab, colls, idx, met in (
                ("DNA", dna, dm, d_idx, met_all),
                ("RNA", rna, rm, r_idx, mbx["butyrate"][np.array([m_idx[c] for c in rm])])):
            if ec not in tab:
                continue
            v = tab[ec][np.array([idx[c] for c in colls])]
            s = spear(v, met)
            if s is None:
                continue
            rows.append({"stratum": "all", "pathway": "butyrate_single", "n_ec": 1,
                         "metabolite": f"butyrate|{ec}", "assay": assay,
                         "n": s["n"], "rho": round(s["rho"], 4), "p": s["p"]})
            say(f"   {ec:<10} {assay:<5} {s['n']:>4} {s['rho']:>+7.3f} {s['p']:>10.2e}")

    out = os.path.join(outdir, "ibdmdb_metabolite_tests.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["stratum", "pathway", "n_ec", "metabolite",
                                           "assay", "n", "rho", "p"])
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(outdir, "IBDMDB_METABOLITE_TEST.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    say()
    say(f"WROTE {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
