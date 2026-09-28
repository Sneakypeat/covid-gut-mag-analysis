# MICOM for the five largest communities

Five control communities (250 to 285 modelled members) exceeded the 21,600 s solver
ceiling in the first stage of MICOM's cooperative tradeoff, the linear maximisation of
community growth, on the HiGHS/OSQP hybrid interface. For those five only:

1. `step1_mustar_gurobi.py <donor_index> <outdir>` builds the community exactly as the
   main community script does and solves that linear programme with Gurobi, writing
   `mustar_<donor>.json`.
2. `step2_seeded_main.py <donor_index> <outdir>` runs the unchanged community script in
   the original environment, answering MICOM's single community-growth call from step 1.
   The quadratic stage, crossover and all outputs are MICOM's own. It also logs whether
   crossover ran and the tradeoff actually realised.
3. `merge76.sh` links the 71 directly solved donors and these five into one 76-donor set
   and rebuilds the substrate-use and net-exchange tables.

The value of a linear programme's optimum does not depend on the solver: for a donor
that solved directly, HiGHS and two Gurobi methods gave the same maximum community growth
to ten significant figures. Paths come from `COVID_MAG_PIPELINE`; the Gurobi licence is
read from `GRB_LICENSE_FILE`.
