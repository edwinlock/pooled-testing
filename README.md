# Instructions

## Setup

1. Install Julia (from www.julialang.org).
2. Install Gurobi (https://www.gurobi.com/). You may need to create a Gurobi academic account and get a licence. Follow Gurobi's instructions to store the licence locally on your computer, then follow the installation instructions for Gurobi/JuMP/Julia here: https://github.com/jump-dev/Gurobi.jl
3. Install MOSEK (https://mosek.com) and follow the instructions to get a licence. Academic licences are available for free.

## Activate the project environment

This repository ships its own Julia environment, defined by `Project.toml` and `Manifest.toml` (the equivalent of a Python virtual environment). Before running anything, instantiate it once to download and pin the exact dependency versions:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then always run with `--project=.` so Julia uses this environment rather than your global one:

```sh
julia --project=. run.jl
```

Alternatively, from the Julia REPL press `]` to enter the package manager and run `activate .` followed by `instantiate`. In VS Code, the Julia extension activates the project automatically when you open this folder.

## Run the experiments

Running and analysis are separate. `run.jl` solves every cell and writes each
result to a SQLite store (`data/solves.db`); `analyse.jl` reads that store and
produces summary tables and plots — and can be run while `run.jl` is still going.

The experiments run several MILP solves in parallel across Julia threads, and
each solve runs Gurobi with `Threads=8`. So launch Julia with about **cores ÷ 8**
threads to use the machine fully without oversubscribing — e.g. `-t 8` on a
64-core machine, `-t 12` on a 96-core one. (Avoid `-t auto`: one Julia thread per
core × 8 Gurobi threads each heavily oversubscribes.)

Run all experiments (adjust `-t` to your core count ÷ 8):

```sh
julia --project=. -t 8 run.jl
```

Run only specific experiments (comma-separated):

```sh
julia --project=. -t 8 run.jl --experiments 1,3,5
```

Solves already in the store are skipped, so an interrupted or reclaimed run
simply resumes where it left off — just re-run the same command.

Analyse the results (any time, including mid-run) — writes `data/expN-data.csv`,
`tables/expN-summary.tex` and `figs/expN-*.pdf`:

```sh
julia --project=. analyse.jl
julia --project=. analyse.jl --experiments 3   # a subset
```
