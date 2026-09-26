#!/usr/bin/env python3
"""Does encoded capacity predict expression in a gut community?

Figures 4e and 4f report abundance-weighted encoded capacity and treat it as
informative about what a community does. The chain is

    gene in genome -> genome abundant in donor -> community has capacity -> ?

The last arrow is the one we cannot test, because no COVID cohort has faecal
RNA: sixty-plus COVID gut studies are deposited and every one is DNA. This
tests that arrow on the cohort that does have both, HMP2 / IBDMDB, where 1,638
metagenomes and 835 metatranscriptomes were produced from the SAME stool
collections, and asks whether the DNA abundance of an enzyme predicts its
transcript abundance in the same sample.

Design point that matters. A positive DNA-RNA correlation is partly guaranteed:
both tables are compositional and both track who is present. So the butyrate
enzymes are never judged on their own rho. They are ranked inside the
distribution of rho over every EC in the table, which absorbs that shared
baseline. The question is not "is rho above zero", it is "is butyrate
machinery any more decoupled from its own DNA than an average gut enzyme".

Three readouts, none of which is a substitute for measuring our own cohort:
  1. rho(DNA, RNA) per EC across paired samples, with the SCFA enzymes placed
     as a percentile of the all-EC distribution.
  2. median RNA/DNA ratio per EC. Regulation, not presence: a ratio far below
     the table median means encoded but comparatively silent.
  3. both of the above split by diagnosis, because our cases are inflamed and
     an association that holds only in nonIBD would not transfer.

Inputs are the published relative-abundance products, unstratified rows only.
Sample columns are joined to stool collections through `site_sub_coll` in the
metadata, which is what makes a DNA and an RNA column the same stool.

Usage: ibdmdb_capacity_expression.py <ibdmdb_dir> <outdir>
"""
import csv
import gzip
import os
import sys

import numpy as np
from scipy import stats

# Enzymes behind the traits in Figure 4e, in the EC namespace HUMAnN reports.
# Both butyrate routes are represented because the catalogue carries both and a
# rule built on only the kinase route was blind to the commoner one.
FOCUS = {
    "2.8.3.8":  ("butyrate", "butyryl-CoA:acetate CoA-transferase (but route)"),
    "2.7.2.7":  ("butyrate", "butyrate kinase (buk, kinase route)"),
    "2.3.1.19": ("butyrate", "phosphate butyryltransferase (ptb, kinase route)"),
    "1.3.8.1":  ("butyrate", "butyryl-CoA dehydrogenase"),
    "1.1.1.157": ("butyrate", "3-hydroxybutyryl-CoA dehydrogenase"),
    "4.2.1.17": ("butyrate", "enoyl-CoA hydratase"),
    "2.8.3.9":  ("butyrate", "butyrate:acetoacetate CoA-transferase"),
    "2.3.1.8":  ("acetate", "phosphate acetyltransferase (pta)"),
    "2.7.2.1":  ("acetate", "acetate kinase (ackA)"),
    "6.2.1.1":  ("acetate", "acetate-CoA ligase"),
    "5.4.99.2": ("propionate", "methylmalonyl-CoA mutase"),
    "1.2.7.1":  ("central", "pyruvate synthase"),
}
MIN_PREV = 0.50   # EC must be non-zero in at least this share of paired samples
MIN_PAIR = 30     # minimum paired samples for a stratum to be reported


def read_unstratified(path):
    """EC id -> {sample: value}, dropping the taxon-stratified rows.

    HUMAnN writes one unstratified row per EC followed by one row per
    contributing species, tagged with '|g__...'. Keeping both would count the
    same signal twice and weight it by how well MetaPhlAn resolved the clade.
    """
    with gzip.open(path, "rt") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        samples = header[1:]
        out = {}
        for row in rdr:
            name = row[0]
            if "|" in name or name.startswith("#"):
                continue
            if name.startswith(("UNMAPPED", "UNGROUPED", "UniRef")):
                continue
            ec = name.split(":", 1)[0].strip()
            out[ec] = np.asarray(row[1:], dtype=float)
    return samples, out


def load_metadata(path):
    """sample id -> (collection, diagnosis, participant), per data type."""
    mgx, mtx, meta = {}, {}, {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            sid = r["External ID"]
            coll = r["site_sub_coll"]
            rec = (coll, r.get("diagnosis", ""), r["Participant ID"])
            dt = r["data_type"]
            if dt == "metagenomics":
                mgx[sid] = rec
            elif dt == "metatranscriptomics":
                mtx[sid] = rec
            meta[coll] = rec
    return mgx, mtx, meta


def pair_columns(mgx_samples, mtx_samples, mgx_meta, mtx_meta):
    """Indices of DNA and RNA columns that came from the same stool collection.

    MGX and MTX carry different sample identifiers for the same physical
    sample, so joining on the identifier directly would silently return
    nothing. `site_sub_coll` is the collection the aliquots came from.
    """
    mgx_by_coll = {}
    for i, s in enumerate(mgx_samples):
        if s in mgx_meta:
            mgx_by_coll.setdefault(mgx_meta[s][0], i)
    pairs = []
    for j, s in enumerate(mtx_samples):
        if s not in mtx_meta:
            continue
        coll = mtx_meta[s][0]
        if coll in mgx_by_coll:
            pairs.append((coll, mgx_by_coll[coll], j))
    return pairs


def per_ec_stats(dna, rna):
    """Spearman rho and median RNA/DNA over samples where the gene is present.

    The ratio is conditioned on DNA > 0 because RNA/0 is undefined, not
    infinite: a sample with no detected gene says nothing about its
    regulation and must not enter the ratio.
    """
    ok = np.isfinite(dna) & np.isfinite(rna)
    d, r = dna[ok], rna[ok]
    present = (d > 0).sum()
    if present < MIN_PAIR or len(d) < MIN_PAIR:
        return None
    if np.all(d == d[0]) or np.all(r == r[0]):
        return None
    rho, p = stats.spearmanr(d, r)
    if not np.isfinite(rho):
        return None
    both = d > 0
    ratio = np.median(r[both] / d[both]) if both.sum() else np.nan
    return {
        "rho": rho,
        "p": p,
        "n": int(len(d)),
        "prev_dna": present / len(d),
        "prev_rna": float((r > 0).sum()) / len(r),
        "rna_dna": ratio,
        "median_dna": float(np.median(d)),
        "median_rna": float(np.median(r)),
    }


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        running = min(running, p[i] * n / (rank + 1))
        q[i] = running
    return q


def analyse(dna_tab, rna_tab, pairs, label, rows):
    """One stratum. Writes into `rows`, returns the all-EC rho distribution."""
    if len(pairs) < MIN_PAIR:
        return None
    di = np.array([p[1] for p in pairs])
    ri = np.array([p[2] for p in pairs])
    shared = sorted(set(dna_tab) & set(rna_tab))
    stats_by_ec, keep = {}, []
    for ec in shared:
        s = per_ec_stats(dna_tab[ec][di], rna_tab[ec][ri])
        if s is None or s["prev_dna"] < MIN_PREV:
            continue
        stats_by_ec[ec] = s
        keep.append(ec)
    if not keep:
        return None
    qs = bh([stats_by_ec[e]["p"] for e in keep])
    rhos = np.array([stats_by_ec[e]["rho"] for e in keep])
    ratios = np.array([stats_by_ec[e]["rna_dna"] for e in keep])
    ratios = ratios[np.isfinite(ratios)]
    for ec, q in zip(keep, qs):
        s = stats_by_ec[ec]
        grp, desc = FOCUS.get(ec, ("", ""))
        rows.append({
            "stratum": label, "ec": ec, "trait_group": grp, "enzyme": desc,
            "n_pairs": s["n"], "rho": round(s["rho"], 4),
            "p": s["p"], "q": q,
            "rho_pctile": round(100.0 * (rhos < s["rho"]).mean(), 1),
            "rna_dna": s["rna_dna"],
            "rna_dna_pctile": round(100.0 * (ratios < s["rna_dna"]).mean(), 1)
                              if np.isfinite(s["rna_dna"]) else "",
            "prev_dna": round(s["prev_dna"], 3),
            "prev_rna": round(s["prev_rna"], 3),
            "median_dna": s["median_dna"], "median_rna": s["median_rna"],
        })
    return {"n_ec": len(keep), "n_pairs": len(pairs), "rhos": rhos,
            "ratios": ratios, "stats": stats_by_ec}


def main():
    src, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    say_lines = []

    def say(s=""):
        print(s, flush=True)
        say_lines.append(s)

    mgx_meta, mtx_meta, coll_meta = load_metadata(os.path.join(src, "hmp2_metadata.csv"))
    say("loading EC tables (unstratified rows only)")
    mgx_samples, dna = read_unstratified(os.path.join(src, "mgx_ecs_relab.tsv.gz"))
    mtx_samples, rna = read_unstratified(os.path.join(src, "mtx_ecs_relab.tsv.gz"))
    say(f"  DNA {len(dna)} ECs x {len(mgx_samples)} samples")
    say(f"  RNA {len(rna)} ECs x {len(mtx_samples)} samples")

    pairs = pair_columns(mgx_samples, mtx_samples, mgx_meta, mtx_meta)
    say(f"paired stool collections with both DNA and RNA: {len(pairs)}")
    if not pairs:
        sys.exit("no paired samples; check site_sub_coll join")

    strata = {"all": pairs}
    for name, keep in (("nonIBD", {"nonIBD"}), ("IBD", {"CD", "UC"})):
        strata[name] = [p for p in pairs if coll_meta.get(p[0], ("", "", ""))[1] in keep]

    rows, summary = [], {}
    for label, pr in strata.items():
        say(f"  {label}: {len(pr)} pairs")
        res = analyse(dna, rna, pr, label, rows)
        if res:
            summary[label] = res

    fields = ["stratum", "ec", "trait_group", "enzyme", "n_pairs", "rho", "p", "q",
              "rho_pctile", "rna_dna", "rna_dna_pctile", "prev_dna", "prev_rna",
              "median_dna", "median_rna"]
    out_csv = os.path.join(outdir, "ibdmdb_ec_dna_rna.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    for label in ("all", "nonIBD", "IBD"):
        if label not in summary:
            continue
        s = summary[label]
        rhos = s["rhos"]
        say()
        say(f"## {label}: {s['n_pairs']} paired samples, {s['n_ec']} ECs "
            f"at >={int(MIN_PREV*100)}% DNA prevalence")
        say(f"   all-EC rho: median {np.median(rhos):+.3f}  "
            f"IQR [{np.percentile(rhos,25):+.3f}, {np.percentile(rhos,75):+.3f}]  "
            f"share above 0: {100*(rhos>0).mean():.1f}%")
        say(f"   all-EC RNA/DNA ratio: median {np.median(s['ratios']):.2f}")
        say(f"   {'EC':<10} {'rho':>7} {'pct':>6} {'q':>9} {'RNA/DNA':>8} {'pct':>6}  enzyme")
        for ec, (grp, desc) in sorted(FOCUS.items(), key=lambda kv: kv[1]):
            hit = [r for r in rows if r["stratum"] == label and r["ec"] == ec]
            if not hit:
                say(f"   {ec:<10} {'below prevalence floor':>39}  {desc}")
                continue
            r = hit[0]
            say(f"   {ec:<10} {r['rho']:>+7.3f} {r['rho_pctile']:>5.0f}% "
                f"{r['q']:>9.2e} {r['rna_dna']:>8.2f} {r['rna_dna_pctile']:>5.0f}%  {desc}")

    with open(os.path.join(outdir, "IBDMDB_CAPACITY_EXPRESSION.txt"), "w") as fh:
        fh.write("\n".join(say_lines) + "\n")
    say()
    say(f"WROTE {out_csv} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
