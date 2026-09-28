"""Maximum community growth (MICOM stage A) with Gurobi, built exactly as hpc_micom_full.py builds it.

Only the LP optimum is taken from here. The value of an LP optimum does not depend on
the solver that finds it, so handing it to the unchanged pipeline leaves every later
step (the l2 quadratic stage on OSQP, fraction bounds, crossover, output) as it was
for the 71 donors that solved directly.
"""
import os, sys, time, json
import pandas as pd
from micom import Community
import gurobipy
# The WLS academic licence allows two concurrent sessions, and sessions of killed jobs
# linger on the licence server until they time out. Wait for a free slot rather than fail.
for attempt in range(288):
    try:
        gurobipy.Env().dispose()
        break
    except gurobipy.GurobiError as exc:
        if "Too many sessions" not in str(exc):
            raise
        print(f"licence busy ({exc}); retry {attempt + 1} in 300 s", flush=True)
        time.sleep(300)
else:
    raise SystemExit("no Gurobi licence session within 24 h")
from micom.util import _format_min_growth, _apply_min_growth
BASE = os.environ["COVID_MAG_PIPELINE"]   # root of the MAG pipeline tree on the cluster
MODELS = f"{BASE}/metabolic_support/gapfill_v1/mags_all"
FLOOR = 0.01
idx, outdir = int(sys.argv[1]), sys.argv[2]
ab = pd.read_csv(f"{BASE}/R_ready/MAG_relative_abundance_percent.tsv", sep="\t", index_col=0)
meta = pd.read_csv(f"{BASE}/R_ready/sample_metadata_76.tsv", sep="\t")
sample = meta["sample"].tolist()[idx]
have = {d for d in os.listdir(MODELS) if os.path.isfile(os.path.join(MODELS, d, "model.xml"))}
col = ab.loc[[g for g in ab.index if g in have], sample]
det = col[col > FLOOR]
tax = pd.DataFrame({"id": det.index, "abundance": det.values / det.values.sum(),
                    "file": [os.path.join(MODELS, g, "model.xml") for g in det.index]})
t0 = time.time()
com = Community(tax, id=sample, solver="gurobi", progress=False)
cfg = com.solver.configuration
cfg.timeout = 172800
cfg.tolerances.feasibility = 1e-7
cfg.tolerances.optimality = 1e-7
med = pd.read_csv(f"{BASE}/metabolic_support/gapfill_v1/scripts/western_diet_gut_carveme.csv").groupby(
    "reaction", as_index=False).agg({"flux": "sum"})
avail = {r.id for r in com.exchanges}
com.medium = {r: f for r, f in zip(med.reaction, med.flux) if r in avail}
build = time.time() - t0
# identical to micom.problems.cooperative_tradeoff up to its optimize_with_retry call
_apply_min_growth(com, _format_min_growth(0.0, com.taxa))
com.objective = com.scale * com.variables.community_objective
t0 = time.time()
sol = com.optimize()
out = {"sample": sample, "index": idx, "members": len(tax), "build_s": round(build),
       "lp_s": round(time.time() - t0), "status": getattr(sol, "status", None),
       "scale": com.scale, "mu_star": (sol.objective_value / com.scale) if sol is not None else None}
json.dump(out, open(os.path.join(outdir, f"mustar_{sample}.json"), "w"), indent=1)
print(out, flush=True)
if out["mu_star"] is None or out["status"] != "optimal":
    raise SystemExit("LP not optimal")
