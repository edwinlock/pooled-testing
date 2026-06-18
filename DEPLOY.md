# Scaleway Ubuntu Deployment

These instructions deploy the project on a Scaleway Ubuntu server with an
AMD EPYC CPU. They intentionally avoid AWS, S3 buckets, IAM roles and persistent
object storage. Results live on the server disk; copy them back to your laptop
periodically if you want an off-server backup.

The commands assume Ubuntu 24.04 or 22.04 on x86_64.

## 1. Create And Connect To The Server

Create a Scaleway instance with:

- Image: Ubuntu 24.04 LTS or Ubuntu 22.04 LTS
- CPU: AMD EPYC / x86_64
- Disk: at least 30 GB
- SSH: your public key installed

The Ubuntu images use the `ubuntu` user. If your image uses a different user,
replace `ubuntu` below.

```sh
ssh ubuntu@<server-ip>
```

Copy the Gurobi and MOSEK licences from your local machine:

```sh
scp ~/gurobi.lic ~/mosek.lic ubuntu@<server-ip>:~/
```

### Record the server specs

Solve timings are hardware-dependent, so record what machine produced them.
Save a spec sheet on the server (it travels back with the results):

```sh
{ echo "=== HOST ==="; hostnamectl; \
  echo "=== CPU ===";  lscpu; \
  echo "=== MEM ===";  free -h; \
  echo "=== DISK ==="; df -h; lsblk; \
  echo "=== OS ===";   cat /etc/os-release; \
  echo "=== KERNEL ==="; uname -a; } > ~/pooled-testing/server-specs.txt 2>&1
```

`lscpu` reports the **physical core count**, which sets the Julia thread count in
section 6 (`cores / 8`). If `inxi` is available (`sudo apt-get install -y inxi`),
`inxi -Fxz` gives a tidier one-shot summary.

## 2. Install System Packages

On the server:

```sh
sudo apt-get update
sudo apt-get install -y \
  bzip2 \
  ca-certificates \
  curl \
  git \
  htop \
  tar \
  tmux \
  wget
```

## 3. Install Julia, Gurobi And MOSEK

These are the x86_64/Linux builds for the AMD EPYC server. Keep Gurobi on 13.x
to match the Julia project environment and the tuning script.

```sh
cd ~

JULIA_VERSION="1.12.6"
JULIA_SERIES="1.12"
GUROBI_VERSION="13.0.1"
GUROBI_SERIES="13.0"
GUROBI_DIR="gurobi1301"
MOSEK_VERSION="11.2.1"
MOSEK_SERIES="11.2"

wget -q "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_SERIES}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
tar zxf "julia-${JULIA_VERSION}-linux-x86_64.tar.gz"

wget -q "https://packages.gurobi.com/${GUROBI_SERIES}/gurobi${GUROBI_VERSION}_linux64.tar.gz"
tar zxf "gurobi${GUROBI_VERSION}_linux64.tar.gz"

wget -q "https://download.mosek.com/stable/${MOSEK_VERSION}/mosektoolslinux64x86.tar.bz2"
tar xjf "mosektoolslinux64x86.tar.bz2"
```

Add the environment variables:

```sh
cat >> ~/.bashrc <<'EOF'

# >>> pooled-testing env >>>
export GUROBI_HOME="$HOME/gurobi1301/linux64"
export GRB_LICENSE_FILE="$HOME/gurobi.lic"
export MOSEKLM_LICENSE_FILE="$HOME/mosek.lic"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:${GUROBI_HOME}/lib"
export PATH="${PATH}:$HOME/mosek/11.2/tools/platform/linux64x86/bin"
export PATH="${PATH}:$HOME/julia-1.12.6/bin"
# <<< pooled-testing env <<<
EOF

source ~/.bashrc
```

Check the tools:

```sh
julia --version
gurobi_cl --version
```

## 4. Clone And Instantiate The Project

```sh
git clone https://github.com/edwinlock/pooled-testing.git
cd pooled-testing

julia --project=. -e '
    using Pkg
    if isempty(Pkg.Registry.reachable_registries())
        Pkg.Registry.add("General")
    else
        Pkg.Registry.update()
    end
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.precompile()
'
```

Confirm that Julia sees Gurobi 13:

```sh
julia --project=. -e 'using Gurobi; println(Gurobi.GRB_VERSION_MAJOR, ".", Gurobi.GRB_VERSION_MINOR, ".", Gurobi.GRB_VERSION_TECHNICAL)'
```

## 5. Run Experiments

Use a separate tmux session so the run survives disconnects:

```sh
tmux new -s experiments
cd ~/pooled-testing
```

Each MILP solve uses Gurobi `Threads=8`, and the experiment runner parallelizes
independent solves across Julia threads. Choose Julia threads as roughly:

```text
Julia threads = physical cores / 8
```

Use the physical core count from `lscpu` (section 1's `server-specs.txt`). For
example, a 24-core EPYC → `-t 3`; 48 physical cores → `-t 6`. Avoid `-t auto`,
which spawns one Julia thread per core and, with 8 Gurobi threads each, badly
oversubscribes.

First run one multithreaded experiment as a sanity check:

```sh
julia --project=. -t 3 experiments/run.jl --experiments 3
```

Then run all experiments, or a subset:

```sh
julia --project=. -t 3 experiments/run.jl

# Or selected experiments
julia --project=. -t 3 experiments/run.jl --experiments 1,3,5
```

Solves already in the store are skipped, so an interrupted or reclaimed run
simply resumes — re-run the same command. Experiment parameters (budgets, pool
sizes, per-experiment `accuracy`) live in `experiments/constants.jl`.

### Re-solving and overwriting (`--algs`, `--overwrite`)

By default a run only fills in missing cells. Two flags let you redo work:

- `--algs a,b` restricts the run to named algorithms (e.g. `greedy`).
- `--overwrite` re-solves cells already in the store and overwrites them in
  place (otherwise they are skipped).

These matter for **clean timings**. Wall-clock time is recorded per solve, so a
fast algorithm timed while heavy MILP solves saturate the CPU looks far slower
than it is. To re-time greedy alone, with no contention, run it single-threaded
after the MILP work is done:

```sh
julia --project=. -t 1 experiments/run.jl --algs greedy --overwrite
```

Populations are generated deterministically from each experiment's seed *before*
any algorithm runs, and greedy is deterministic, so a greedy re-run reproduces
the exact same populations and welfares — only the recorded time changes.

## 6. Analyse Results

Analysis is separate from solving: `analyse.jl` only reads the store, so it is
safe to run while `run.jl` is still solving.

```sh
julia --project=. experiments/analyse.jl

# A subset
julia --project=. experiments/analyse.jl --experiments 3,4
```

It writes, per experiment, `experiments/data/expN-data.csv`,
`experiments/tables/expN-summary.tex` (paper-ready booktabs tables) and
`experiments/figs/expN-*.pdf`. A few things to know so the output is correct:

- **Budgets are filtered by `constants.jl`.** Each experiment is analysed only
  over the budgets in its spec; surplus budgets left in the store from earlier
  runs are ignored (not deleted). Keep `constants.jl` matching the data you want.
- **Latest parameter regime wins.** If the store holds solves at more than one
  `accuracy`/`K`, analysis uses the most recent regime per algorithm and reports
  it; mixing regimes in one table is avoided by design.
- **Incomplete cells show as `--`.** A budget whose MILP solve has not finished
  yet appears blank rather than breaking the table.

### Running analysis against a copied data directory

If you have copied only the results to your laptop (data-only, without running
the server's code), run the repo's `analyse.jl` and point `--rootdir` at the
copied `experiments` directory holding `data/solves.db`:

```sh
julia --project=. experiments/analyse.jl --rootdir /path/to/copied/experiments
```

`run.jl` accepts the same `--rootdir` if you ever re-solve against a copied store.

## 7. Back Up Results To Your Laptop

There is no bucket or persistent object storage in this workflow. The important
outputs are:

- `experiments/data/`  (includes `solves.db` — the SQLite store; portable across
  OS/architecture)
- `experiments/tables/`
- `experiments/figs/`
- `server-specs.txt`   (the spec sheet from section 1)

The SQLite store is safest to copy when nothing is writing to it. If a run may be
active, checkpoint the write-ahead log first so the main file is self-contained:

```sh
sqlite3 ~/pooled-testing/experiments/data/solves.db "PRAGMA wal_checkpoint(TRUNCATE);"
```

From your local machine, copy results down periodically (re-runnable; only
transfers what changed):

```sh
rsync -avz ubuntu@<server-ip>:~/pooled-testing/experiments/data ./experiments/
rsync -avz ubuntu@<server-ip>:~/pooled-testing/experiments/tables ./experiments/
rsync -avz ubuntu@<server-ip>:~/pooled-testing/experiments/figs ./experiments/
rsync -avz ubuntu@<server-ip>:~/pooled-testing/server-specs.txt .
```

Or create one compressed archive on the server:

```sh
cd ~/pooled-testing
tar czf pooled-testing-results.tgz \
  experiments/data \
  experiments/tables \
  experiments/figs \
  server-specs.txt
```

Then download it locally:

```sh
scp ubuntu@<server-ip>:~/pooled-testing/pooled-testing-results.tgz .
```

## 8. Monitor The Server

Open another tmux window and run:

```sh
htop
```

Useful tmux keys:

- `Ctrl-b c`: new window
- `Ctrl-b n`: next window
- `Ctrl-b p`: previous window
- `Ctrl-b d`: detach
- `tmux attach -t experiments`: reattach
