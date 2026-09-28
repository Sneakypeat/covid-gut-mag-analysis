#!/bin/bash
# Combined 76-donor MICOM set: 71 from micom_full (unchanged) + the 5 largest controls whose
# maximum community growth was precomputed with Gurobi (micom_seeded_20260927). Links only.
set -euo pipefail
G="${COVID_MAG_PIPELINE:?set COVID_MAG_PIPELINE}/metabolic_support/gapfill_v1"
C=$G/micom_76_20260927
mkdir -p $C
for f in $G/micom_full/{members,exchanges,fluxes,status}_*; do ln -sf "$f" "$C/"; done
for s in PRJEB7949_H2 predict1_MTG_0686 predict1_MTG_0964 predict1_MTG_0959 predict1_MTG_0767; do
  for k in members_$s.csv exchanges_$s.csv fluxes_$s.tsv.gz status_$s.json; do
    test -s $G/micom_seeded_20260927/$k
    ln -sf $G/micom_seeded_20260927/$k $C/$k
  done
done
echo "status files: $(ls $C/status_*.json | wc -l) | fluxes: $(ls $C/fluxes_*.tsv.gz | wc -l)"
PY="${PYTHON:-python}"
mkdir -p $G/substrate_use_76 $G/collect_76
$PY -u $G/scripts/micom_substrate_use.py $C $G/scripts/fam584.tsv $G/substrate_use_76
export C G
$PY - <<'PYEOF'
import glob, json, os
import pandas as pd
C = os.environ["C"]
prod = {"acetate": "EX_ac_m", "butyrate": "EX_but_m", "propionate": "EX_ppa_m", "succinate": "EX_succ_m",
        "formate": "EX_for_m", "ethanol": "EX_etoh_m", "lactate_L": "EX_lac__L_m", "lactate_D": "EX_lac__D_m",
        "co2": "EX_co2_m", "h2": "EX_h2_m"}
rows = []
for st in sorted(glob.glob(C + "/status_*.json")):
    s = json.load(open(st)); ex = pd.read_csv(C + "/exchanges_" + s["sample"] + ".csv", index_col=0)["flux"]
    r = {"sample": s["sample"], "group": s["group"], "growth": s["community_growth"]}
    r.update({k: float(ex.get(v, 0.0)) for k, v in prod.items()})
    rows.append(r)
d = pd.DataFrame(rows); d.to_csv(os.environ["G"] + "/collect_76/micom_net_exchange_76.csv", index=False)
print("net exchange rows", len(d), d.group.value_counts().to_dict())
PYEOF
echo MERGE_DONE
