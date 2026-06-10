# Log

Concise running log of experiment runs and key changes. Newest first.

## 2026-06-10
- Interrupted overnight run: exp 1 (G=5) reached B=22 in ~8h — ~6× slower than paper's Table 1; B=26/30 still ahead (paper ~3.4h/~19h) → not viable as-is.
- Ruled out causes: rounding/clamp change affects 0 of 130 pilot people; MILP formulation gap (`T*ε`, K-based) unchanged since 2023; Gurobi default MIPGap (1e-4) unchanged 2023→13.0.1.
- Action: cap Gurobi `Threads=8`; params at top of experiments.jl.
- MIPGap: tried 0.5% but Table 1 MILP-vs-greedy diffs are 0.004%–0.33% < 0.5%, so a loose gap makes MILP look worse than greedy → meaningless. Set to **1e-4 (0.01%)**, matching paper. So MIPGap can't be the tractability lever — big budgets stay slow; revisit via budget cap / hardware.
- Added `showspeed=true` to progress bar so per-budget solve time (s/it) + elapsed are visible.
- Exp-1 parallelisation: deferred — measure serial Threads=8 + MIPGap 0.5% first; only split Gurobi/MOSEK if still too slow.
- Earlier: MOSEK aborts (TBB) on Julia worker threads on Linux aarch64 → exps 3/4 run serial; exp 5 (Gurobi-only) parallel.
- Synthetic exps n=200 → n=150 (match paper). Pilot budgets capped at 30.
- Infra: on Spot c8g/c7g us-east-2 (capacity flaky); S3 sync for results; setup.sh provisions instance.
