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
julia --project=. experiments.jl
```

Alternatively, from the Julia REPL press `]` to enter the package manager and run `activate .` followed by `instantiate`. In VS Code, the Julia extension activates the project automatically when you open this folder.

## Run the experiments

The experiments run several MILP solves in parallel across Julia threads, and
each solve runs Gurobi with `Threads=8`. So launch Julia with about **cores ÷ 8**
threads to use the machine fully without oversubscribing — e.g. `-t 12` on a
96-core machine, `-t 4` on a 32-core one. (Avoid `-t auto`: one Julia thread per
core × 8 Gurobi threads each heavily oversubscribes.)

Run all experiments (adjust `-t` to your core count ÷ 8):

```sh
julia --project=. -t 12 experiments.jl
```

Run only specific experiments (comma-separated):

```sh
julia --project=. -t 12 experiments.jl --experiments 1,3,5
```

Write outputs (under `data/`, `tables/`, `figs/`) to a different root directory:

```sh
julia --project=. -t 12 experiments.jl --rootdir path/to/output
```

Experiments whose output CSV (`data/expN-data.csv`) already exists are skipped,
so an interrupted run simply resumes where it left off. Pass `--force` to re-run
them anyway:

```sh
julia --project=. -t 12 experiments.jl --force
```
