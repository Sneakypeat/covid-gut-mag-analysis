#!/usr/bin/env python3
"""Collect the per-donor MICOM results into donor-level tables.

Panel f asks what each genome could form on its own and weights that by abundance.
This asks the next question: with every member competing for one shared substrate
pool under cooperative tradeoff, what does the community actually export, and who
grows. Net exchange is signed, so a metabolite the community consumes on balance
appears negative and is not counted as production.

Usage: micom_collect_n76.py <micom_dir> <outdir>
"""
import csv, glob, json, os, sys
import pandas as pd

PRODUCTS = {"EX_but_m": "butyrate", "EX_ac_m": "acetate", "EX_ppa_m": "propionate",
            "EX_lac__D_m": "lactate_D", "EX_lac__L_m": "lactate_L",
            "EX_succ_m": "succinate", "EX_for_m": "formate", "EX_etoh_m": "ethanol"}

def cliff(a, b):
    a, b = list(a), list(b)
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (len(a) * len(b)) if a and b else float("nan")

def main():
    src, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    status = [json.load(open(f)) for f in sorted(glob.glob(f"{src}/status_*.json"))]
    st = pd.DataFrame(status)
    st.to_csv(f"{outdir}/micom_status_n76.csv", index=False)
    print(f"donors solved: {len(st)} ({(st.group=='Case').sum()} case, {(st.group=='Control').sum()} control)")
    for g in ["Control", "Case"]:
        d = st[st.group == g]
        if len(d):
            print(f"  {g:8s} members median {d.members.median():.0f} | solve median {d.solve_seconds.median():.0f}s "
                  f"| community growth median {d.community_growth.median():.5f}")

    rows = []
    for f in sorted(glob.glob(f"{src}/exchanges_*.csv")):
        d = pd.read_csv(f, index_col=0)
        col = "flux" if "flux" in d.columns else d.columns[0]
        s = d[col]
        for rid, name in PRODUCTS.items():
            rows.append({"sample_id": d["sample_id"].iloc[0], "group": d["group"].iloc[0],
                         "product": name, "flux": float(s.get(rid, 0.0))})
    ex = pd.DataFrame(rows)
    ex.to_csv(f"{outdir}/micom_exchanges_n76.csv", index=False)

    print(f"\n{'product':11s} {'control':>10s} {'case':>10s} {'delta':>7s} {'p':>10s}   (net export, mmol/gDW/h)")
    out = []
    from scipy.stats import mannwhitneyu
    for name in PRODUCTS.values():
        d = ex[ex["product"] == name]
        c = d.loc[d.group == "Control", "flux"]; k = d.loc[d.group == "Case", "flux"]
        if len(c) < 3 or len(k) < 3:
            continue
        p = mannwhitneyu(k, c, alternative="two-sided").pvalue
        out.append({"product": name, "control_median": c.median(), "case_median": k.median(),
                    "delta": cliff(k, c), "p": p})
    res = pd.DataFrame(out)
    if len(res):
        # Benjamini-Hochberg written out rather than pulled from statsmodels, which
        # is not in every environment this has to run in.
        res = res.sort_values("p").reset_index(drop=True)
        m = len(res)
        raw = [min(1.0, p * m / (i + 1)) for i, p in enumerate(res["p"])]
        q, run = [0] * m, 1.0
        for i in range(m - 1, -1, -1):
            run = min(run, raw[i]); q[i] = run
        res["q"] = q
        res = res.sort_values("delta")
        for _, r in res.iterrows():
            print(f"{r['product']:11s} {r.control_median:10.4f} {r.case_median:10.4f} "
                  f"{r.delta:+7.2f} {r.q:10.3g}")
        res.to_csv(f"{outdir}/micom_exchange_tests_n76.csv", index=False)

    gr = []
    for f in sorted(glob.glob(f"{src}/members_*.csv")):
        d = pd.read_csv(f, index_col=0)
        gr.append(d)
    if gr:
        mem = pd.concat(gr)
        mem.to_csv(f"{outdir}/micom_members_n76.csv")
        print(f"\nmember rows: {len(mem)}")

if __name__ == "__main__":
    main()
