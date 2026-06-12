# Tiny experiment parameters for smoke-testing the full run.jl → store → analyse.jl
# pipeline locally in seconds. Use with:
#
#     POOLED_CONSTANTS=experiments/constants_test.jl julia --project=. -t 2 experiments/run.jl

const UTIL_UPPER_BOUND  = 50               # scale pilot utilities to integers in [1, this]
const PILOT_BUDGETS     = [2,4,6] # {2, 6, ..., 30} — pilot experiments
const SYNTHETIC_BUDGETS = [2,4,6]  # {2, 6, ..., 14} — synthetic experiments
const SYNTHETIC_N       = 50              # synthetic population size
const SYNTHETIC_REPS    = 10               # synthetic populations per experiment
const MILP_K            = 10               # piecewise-linear segments in the MILP
# Total relative error budgets (split between PWL segments and Gurobi MIPGapAbs
# by accuracy_params; override MILP_K when set in an experiment spec).
const PILOT_ACCURACY     = 0.02           # 0.2% for the pilot experiments
const SYNTHETIC_ACCURACY = 0.02           # 0.2% for the synthetic experiments

# Per-experiment specification: (kind, G, seed, accuracy). The pilot and
# synthetic experiments share the standard approx-vs-greedy pair, with the MILP
# driven by the experiment's relative error budget (`accuracy`; see
# accuracy_params). Set accuracy=nothing to use the fixed MILP_K instead.
# Experiment 5 is the disjoint-vs-2-overlap comparison on small populations
# (exact MILPs, so no accuracy parameter applies).
const EXPERIMENT_SPECS = Dict(
    1 => (kind=:pilot,     G=5,  seed=1002, accuracy=PILOT_ACCURACY),
    2 => (kind=:pilot,     G=10, seed=1003, accuracy=PILOT_ACCURACY),
    3 => (kind=:synthetic, G=5,  seed=2001, accuracy=SYNTHETIC_ACCURACY),
    4 => (kind=:synthetic, G=10, seed=2002, accuracy=SYNTHETIC_ACCURACY),
    5 => (kind=:overlap,   G=10, seed=3001, mipgap=1e-3),
)
