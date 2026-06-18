# Pooled Testing

Welfare-maximizing pooled testing experiments: an exact/approximate MILP
benchmark (`approx`) compared against a fast `greedy` heuristic. Solving and
analysis are separate — `experiments/run.jl` solves and writes results to a
SQLite store, `experiments/analyse.jl` reads that store and produces tables and
plots.

To deploy and run this on a fresh Scaleway/Ubuntu server, see
[`DEPLOY.md`](DEPLOY.md). The rest of this file covers running it anywhere the
tools below are installed.

## Setup

1. Install Julia (from www.julialang.org).
2. Install Gurobi (https://www.gurobi.com/). You may need a Gurobi academic
   account and licence. Store the licence locally, then follow the
   Gurobi/JuMP/Julia install instructions at https://github.com/jump-dev/Gurobi.jl
3. Install MOSEK (https://mosek.com) and get a licence. Academic licences are
   free.

## Activate the project environment

This repository ships its own Julia environment, defined by `Project.toml` and
`Manifest.toml` (the equivalent of a Python virtual environment). Instantiate it
once to download and pin the exact dependency versions:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then always run with `--project=.` so Julia uses this environment rather than
your global one. In VS Code, the Julia extension activates the project
automatically when you open this folder.

## Where data and outputs live: `--rootdir`

Both `run.jl` and `analyse.jl` take a `--rootdir DIR` option that sets the working
directory for results. Under `DIR` they read/write:

- `DIR/data/solves.db` — the SQLite solve store, and `DIR/data/expN-data.csv`
- `DIR/tables/expN-summary.tex` — summary tables
- `DIR/figs/expN-*.pdf` — plots

`--rootdir` defaults to the `experiments/` directory inside the repo, so the
plain commands below operate on `experiments/data/solves.db`. To work against a
store kept elsewhere — a results copy on your laptop, a shared drive, a scratch
dir — pass your own path. Throughout this README, `$ROOTDIR` stands for that
directory; set it once and reuse it:

```sh
ROOTDIR="/path/to/your/experiments"   # the dir containing data/solves.db
julia --project=. experiments/analyse.jl --rootdir "$ROOTDIR"
```

Omit `--rootdir` to use the repo's own `experiments/` directory.

## Run the experiments

Running and analysis are separate. `run.jl` solves every cell and writes each
result to a SQLite store (`$ROOTDIR/data/solves.db`); `analyse.jl` reads that
store and produces summary tables and plots — and can be run while `run.jl` is
still going.

The experiments run several MILP solves in parallel across Julia threads, and
each solve runs Gurobi with `Threads=8`. So launch Julia with about
**physical cores ÷ 8** threads to use the machine fully without oversubscribing —
e.g. `-t 3` on a 24-core machine, `-t 6` on a 48-core one. (Avoid `-t auto`: one
Julia thread per core × 8 Gurobi threads each heavily oversubscribes.)

Run all experiments, or a subset (add `--rootdir "$ROOTDIR"` to target a store
outside the repo):

```sh
julia --project=. -t 3 experiments/run.jl
julia --project=. -t 3 experiments/run.jl --experiments 1,3,5
```

Solves already in the store are skipped, so an interrupted or reclaimed run
simply resumes where it left off — just re-run the same command. Experiment
parameters (budgets, pool sizes, per-experiment `accuracy`) live in
[`experiments/constants.jl`](experiments/constants.jl).

### Re-solving and overwriting

By default a run only fills in missing cells. Two flags let you redo work:

- `--algs a,b` restricts the run to named algorithms (e.g. `greedy`).
- `--overwrite` re-solves cells already in the store and overwrites them in place
  (otherwise they are skipped).

Recorded times are wall-clock, so a fast algorithm timed while heavy MILP solves
saturate the CPU looks slower than it is. To re-time `greedy` alone with no
contention, run it single-threaded after the MILP work is done:

```sh
julia --project=. -t 1 experiments/run.jl --algs greedy --overwrite
```

Populations are generated deterministically from each experiment's seed *before*
any algorithm runs, and `greedy` is deterministic, so a re-run reproduces the
same populations and welfares — only the recorded time changes.

## Analyse the results

`analyse.jl` only reads the store, so it is safe to run any time, including mid
run. It writes, per experiment, `experiments/data/expN-data.csv`,
`experiments/tables/expN-summary.tex` (paper-ready booktabs tables) and
`experiments/figs/expN-*.pdf`.

```sh
julia --project=. experiments/analyse.jl
julia --project=. experiments/analyse.jl --experiments 3              # a subset
julia --project=. experiments/analyse.jl --rootdir "$ROOTDIR"         # a store elsewhere
```

The SQLite store is portable across OS/architecture, so you can copy just the
results elsewhere and analyse them in place with `--rootdir` — the tables and
figs are written back under that same directory.

Notes on correctness:

- **Budgets are filtered by `constants.jl`.** Each experiment is analysed only
  over the budgets in its spec; surplus budgets left in the store from earlier
  runs are ignored (not deleted).
- **Latest parameter regime wins.** If the store holds solves at more than one
  `accuracy`/`K`, analysis uses the most recent regime per algorithm.
- **Incomplete cells show as `--`** rather than breaking the table.
