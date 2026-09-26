#!/usr/bin/env python3
"""Gut functional traits v2: same rules, complete source.

v1 scored traits from DRAM's distillate workbook, which carries only its curated
module sheets. That workbook is a strict subset of what the pipeline annotated:
it misses 2 plant-CAZyme genomes, 6 mucin-GH genomes, and it cannot express the
butyryl-CoA:acetate CoA-transferase route at all, which costs 30 butyrate
genomes. Rationale frozen in gut_traits_v2/TRAIT_METHODS_v3.md.

Every rule is unchanged from v1 except butyrate, which gains the second terminal
route. The acid trait is unchanged but now covers the genomes whose graphs were
extracted after the full catalogue was reconstructed.

Outputs -> result2/n76/metabolic_support/gut_traits_v2/genome_traits_wide.tsv
"""
import csv, gzip, json, os, sys, collections

BASE = os.environ.get("COVID_MAG_BASE", os.getcwd())
MS   = os.path.join(BASE, "result2/n76/metabolic_support")
ANN  = os.path.join(MS, "annotations_full")
OUT  = os.path.join(MS, "gut_traits_v2")
ACIDS = {"lac__D", "lac__L", "succ", "ac"}
POS, NEG, UNK = "positive", "evidence_not_detected", "unknown"

def load(path, key):
    d = collections.defaultdict(set)
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["catalog_id"]].add(r[key])
    return d

def main():
    os.makedirs(OUT, exist_ok=True)
    ko = load(os.path.join(ANN, "genome_ko_counts.tsv.gz"), "ko")
    cz = load(os.path.join(ANN, "genome_cazy_counts.tsv.gz"), "cazy")
    v1 = {r["catalog_id"]: r for r in csv.DictReader(
        open(os.path.join(MS, "gut_traits_v1/genome_traits_wide.tsv")), delimiter="\t")}
    plant = sorted({r["gene_id"] for r in csv.DictReader(
        open(os.path.join(MS, "gut_traits_v1/marker_definitions.tsv")), delimiter="\t")
        if str(r.get("plant_backbone_selected")).strip().upper() == "TRUE"})

    req = {}
    for p in [os.path.join(MS, "graph/model_support_sets.json"),
              os.path.join(MS, "graph_extra/model_support_sets.json")]:
        if os.path.exists(p):
            for cid, v in json.loads(open(p).read()).items():
                req[cid] = set(v.get("requirements") or [])
            print(f"  support sets from {os.path.basename(os.path.dirname(p))}: cumulative {len(req)}")
        else:
            print(f"  MISSING {p}")

    cols = ["catalog_id", "primary", "phylum", "family", "enterobacteriaceae_overlay",
            "ancom_direction", "plant_backbone_cazyme", "butyrate_terminal_pair",
            "acetate_terminal_pair", "strict_acid_uptake", "mucin_gh_repertoire"]
    rows, tally = [], collections.Counter()
    for cid, r in v1.items():
        k, c = ko.get(cid, set()), cz.get(cid, set())
        plant_pos = bool(c & set(plant))
        # Either terminal route; v1 could only express the first.
        but_pos = ("K00634" in k and "K00929" in k) or "K01034" in k or "K01035" in k
        ace_pos = "K00625" in k and "K00925" in k
        muc_pos = "GH33" in c and ({"GH29", "GH95"} & c) and ({"GH20", "GH35"} & c)
        acid = UNK if cid not in req else (POS if (ACIDS & req[cid]) else NEG)
        for name, val in [("plant", plant_pos), ("butyrate", but_pos),
                          ("acetate", ace_pos), ("mucin", bool(muc_pos))]:
            tally[name] += int(val)
        tally["acid_known"] += int(acid != UNK)
        tally["acid_pos"] += int(acid == POS)
        rows.append({
            "catalog_id": cid, "primary": r["primary"], "phylum": r["phylum"],
            "family": r["family"], "enterobacteriaceae_overlay": r["enterobacteriaceae_overlay"],
            "ancom_direction": r["ancom_direction"],
            "plant_backbone_cazyme": POS if plant_pos else NEG,
            "butyrate_terminal_pair": POS if but_pos else NEG,
            "acetate_terminal_pair": POS if ace_pos else NEG,
            "strict_acid_uptake": acid,
            "mucin_gh_repertoire": POS if muc_pos else NEG})

    dest = os.path.join(OUT, "genome_traits_wide.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(sorted(rows, key=lambda x: x["catalog_id"]))
    print(f"\nwrote {len(rows)} genomes -> {dest}")
    print(f"{'trait':26s} {'v1':>6s} {'v2':>6s}")
    for name, col in [("plant CAZyme", "plant_backbone_cazyme"),
                      ("butyrate (either route)", "butyrate_terminal_pair"),
                      ("acetate pair", "acetate_terminal_pair"),
                      ("mucin GH", "mucin_gh_repertoire")]:
        old = sum(1 for r in v1.values() if r[col] == POS)
        key = {"plant_backbone_cazyme": "plant", "butyrate_terminal_pair": "butyrate",
               "acetate_terminal_pair": "acetate", "mucin_gh_repertoire": "mucin"}[col]
        print(f"{name:26s} {old:6d} {tally[key]:6d}")
    old_known = sum(1 for r in v1.values() if r["strict_acid_uptake"] != UNK)
    print(f"{'acid trait known for':26s} {old_known:6d} {tally['acid_known']:6d}")
    print(f"{'  of which positive':26s} "
          f"{sum(1 for r in v1.values() if r['strict_acid_uptake']==POS):6d} {tally['acid_pos']:6d}")

if __name__ == "__main__":
    main()
