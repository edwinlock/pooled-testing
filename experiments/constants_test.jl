# Tiny experiment parameters for smoke-testing the full run.jl → store → analyse.jl
# pipeline locally in seconds. Use with:
#
#     POOLED_CONSTANTS=experiments/constants_test.jl julia --project=. -t 2 experiments/run.jl

const UTIL_UPPER_BOUND  = 50
const PILOT_BUDGETS     = [2, 4]   # tiny
const SYNTHETIC_BUDGETS = [2, 4]
const SYNTHETIC_N       = 12       # small populations solve fast
const SYNTHETIC_REPS    = 3
const MILP_K            = 5        # few segments ⇒ small, quick MILP

const EXPERIMENT_SPECS = Dict(
    1 => (kind=:pilot,     G=5,  seed=1002),
    2 => (kind=:pilot,     G=10, seed=1003),
    3 => (kind=:synthetic, G=5,  seed=2001),
    4 => (kind=:synthetic, G=10, seed=2002),
    5 => (kind=:overlap,   G=10, seed=3001),
)
