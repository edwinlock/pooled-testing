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

Scaleway images often allow SSH as `root`. If your image uses a different user,
replace `root` below.

```sh
ssh root@<server-ip>
```

Copy the Gurobi and MOSEK licences from your local machine:

```sh
scp ~/gurobi.lic ~/mosek.lic root@<server-ip>:~/
```

## 2. Install System Packages

On the server:

```sh
apt-get update
apt-get install -y \
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

## 5. Tune One Shared Gurobi Profile

Because this server is AMD EPYC, tune on the server rather than relying on a
profile produced on an Apple laptop or ARM instance. The command below creates
one shared `.prm` profile for `B=18`, `B=22` and `B=30`, with an eight-hour total
tuning cap and 45 minutes per trial/model.

```sh
tmux new -s gurobi-tune
cd ~/pooled-testing
julia --project=. experiments/tune_gurobi.jl \
  --budgets 18,22,30 \
  --multi-model \
  --time-limit 2700 \
  --tune-time-limit 28800 \
  --baseline-runs 1 \
  --validate-runs 2
```

Artifacts go to `experiments/gurobi-tuning/`. At the end, the script prints the
baseline mean runtime, tuned validation mean runtime, seconds saved, percentage
saved and speedup for each budget.

Detach from tmux with `Ctrl-b` then `d`. Reattach with:

```sh
tmux attach -t gurobi-tune
```

To validate an existing shared `.prm` without re-tuning:

```sh
julia --project=. experiments/tune_gurobi.jl \
  --budgets 18,22,30 \
  --multi-model \
  --time-limit 2700 \
  --baseline-runs 1 \
  --skip-tune \
  --validate-runs 2
```

The tuning script uses the `.prm` only for validation. The main experiment run
still uses paper-default Gurobi settings unless you explicitly wire a tuned
`param_file` into the experiment code.

## 6. Run Experiments

Use a separate tmux session:

```sh
tmux new -s experiments
cd ~/pooled-testing
```

Each MILP solve uses Gurobi `Threads=8`, and the experiment runner parallelizes
independent solves across Julia threads. Choose Julia threads as roughly:

```text
Julia threads = physical cores / 8
```

For example, if the server has 64 vCPUs and 32 physical cores, start with
`-t 4`; if it has 96 vCPUs and 48 physical cores, start with `-t 6`. Avoid
`-t auto`, which will oversubscribe badly.

First run one multithreaded experiment as a sanity check:

```sh
julia --project=. -t 4 experiments/run.jl --experiments 3
```

Then run all experiments, or a subset:

```sh
julia --project=. -t 4 experiments/run.jl

# Or selected experiments
julia --project=. -t 4 experiments/run.jl --experiments 1,3,5
```

Analyse results:

```sh
julia --project=. experiments/analyse.jl
```

## 7. Back Up Results To Your Laptop

There is no bucket or persistent object storage in this workflow. The important
outputs are:

- `experiments/data/`
- `experiments/tables/`
- `experiments/figs/`
- `experiments/gurobi-tuning/`

From your local machine, copy them down periodically:

```sh
rsync -avz root@<server-ip>:~/pooled-testing/experiments/data ./experiments/
rsync -avz root@<server-ip>:~/pooled-testing/experiments/tables ./experiments/
rsync -avz root@<server-ip>:~/pooled-testing/experiments/figs ./experiments/
rsync -avz root@<server-ip>:~/pooled-testing/experiments/gurobi-tuning ./experiments/
```

Or create one compressed archive on the server:

```sh
cd ~/pooled-testing
tar czf pooled-testing-results.tgz \
  experiments/data \
  experiments/tables \
  experiments/figs \
  experiments/gurobi-tuning
```

Then download it locally:

```sh
scp root@<server-ip>:~/pooled-testing/pooled-testing-results.tgz .
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
