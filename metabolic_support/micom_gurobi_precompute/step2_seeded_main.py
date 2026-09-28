"""Run the unchanged hpc_micom_full.py with MICOM's community-growth LP answered by step 1."""
import os, sys, json, runpy
import pandas as pd
import micom.problems as P
MAIN = os.environ.get("MICOM_MAIN_SCRIPT", os.path.join(os.environ["COVID_MAG_PIPELINE"], "metabolic_support/gapfill_v1/scripts/hpc_micom_full.py"))
idx, outdir = int(sys.argv[1]), sys.argv[2]
meta = pd.read_csv(os.path.join(os.environ["COVID_MAG_PIPELINE"], "R_ready/sample_metadata_76.tsv"), sep="\t")
sample = meta["sample"].tolist()[idx]
mu = float(json.load(open(os.path.join(outdir, f"mustar_{sample}.json")))["mu_star"])
original = P.optimize_with_retry
calls = {"seeded": 0}
def seeded(com, message="could not get optimum."):
    if message == "could not get community growth rate.":
        calls["seeded"] += 1
        print(f"[seeded] community growth LP answered from Gurobi: mu*={mu:.8f}", flush=True)
        return mu * com.scale
    return original(com, message)
P.optimize_with_retry = seeded
# Logging only: record whether MICOM fell back to crossover and the tradeoff actually realised.
original_crossover = P.crossover
def logged_crossover(com, sol, fluxes=False):
    print(f"[crossover] QP status not optimal; crossover from QP community growth {sol.growth_rate:.6f}", flush=True)
    out = original_crossover(com, sol, fluxes=fluxes)
    print(f"[crossover] result community growth {out.growth_rate:.6f} = {out.growth_rate / mu:.4f} x mu*", flush=True)
    return out
P.crossover = logged_crossover
sys.argv = [MAIN, str(idx), outdir]
try:
    runpy.run_path(MAIN, run_name="__main__")
finally:
    print(f"[seeded] seeded calls={calls['seeded']} | mu*={mu:.8f} (maximum community growth precomputed with Gurobi)", flush=True)
