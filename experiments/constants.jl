# Experiment parameters. Edit these to change what the experiments run.

const UTIL_UPPER_BOUND  = 50               # scale pilot utilities to integers in [1, this]
const PILOT_BUDGETS     = collect(2:4:30)  # {2, 6, ..., 30} — pilot experiments
const SYNTHETIC_BUDGETS = collect(2:2:12)  # {2, 4, ..., 12} — synthetic experiments
const SYNTHETIC_N       = 150              # synthetic population size
const SYNTHETIC_REPS    = 20               # synthetic populations per experiment
const MILP_K            = 15               # piecewise-linear segments in the MILP

# Per-experiment specification: (kind, G, seed). The pilot and synthetic
# experiments share the standard approx-vs-greedy pair; experiment 5 is the
# disjoint-vs-2-overlap comparison on small populations.
const EXPERIMENT_SPECS = Dict(
    1 => (kind=:pilot,     G=5,  seed=1002),
    2 => (kind=:pilot,     G=10, seed=1003),
    3 => (kind=:synthetic, G=5,  seed=2001),
    4 => (kind=:synthetic, G=10, seed=2002),
    5 => (kind=:overlap,   G=10, seed=3001),
)
