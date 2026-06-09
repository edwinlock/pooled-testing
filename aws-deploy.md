# Instructions to deploy to AWS EC2 instance

**Instance type:** use an ARM (Graviton) compute-optimised instance, e.g.
`c8g.24xlarge` (96 vCPUs) or `c8g.16xlarge` (64 vCPUs). These match the ARM
Julia/Gurobi/MOSEK builds downloaded below. Pick the size to fit the parallel
experiments: experiments 3–6 run their independent solves across Julia threads,
so launch Julia with `-t auto` (see step 7) to use all the cores. `c8g.48xlarge`
(192 vCPUs) is the largest, but is only worth it if you saturate it.

## 0. Connect to the instance

SSH in using your key pair and the instance's **public** DNS name or IP (find it
in the EC2 console under "Public IPv4 DNS" — the internal `ip-…compute.internal`
name is not reachable from your laptop). Amazon Linux logs in as `ec2-user`:

```sh
ssh -i pooled_testing.pem ec2-user@<public-dns-name>
```

Step 1 runs on your local machine; the remaining steps run on the instance.

## 1. Copy over the MOSEK and Gurobi licence files

Run this from your local machine, copying both licences to the instance's home
directory (the env vars in step 4 point here, so no destination setup is needed):

```sh
scp -i pooled-testing.pem ~/mosek.lic ~/gurobi.lic ec2-user@ec2-3-145-200-200.us-east-2.compute.amazonaws.com:~/
```

## 2. Install bzip2

```sh
sudo yum install bzip2 git
```

## 3. Download and unpack Julia, Gurobi and MOSEK

**Julia**

```sh
wget https://julialang-s3.julialang.org/bin/linux/aarch64/1.12/julia-1.12.6-linux-aarch64.tar.gz
tar zxvf julia-1.12.6-linux-aarch64.tar.gz
```

**Gurobi**

```sh
wget https://packages.gurobi.com/13.0/gurobi13.0.1_armlinux64.tar.gz
tar zxvf gurobi13.0.1_armlinux64.tar.gz
```

**MOSEK**

```sh
wget https://download.mosek.com/stable/11.2.1/mosektoolslinuxaarch64.tar.bz2
tar zxjf mosektoolslinuxaarch64.tar.bz2
```

## 4. Add paths to `~/.bashrc`

Append the following to `~/.bashrc` so the variables persist across sessions
(plain `export` in the shell would only last for the current session). Use
`$HOME` rather than `~`, since `~` is not expanded inside double quotes:

```sh
cat >> ~/.bashrc <<'EOF'
export GUROBI_HOME="$HOME/gurobi1301/armlinux64"
export GRB_LICENSE_FILE="$HOME/gurobi.lic"
export MOSEKLM_LICENSE_FILE="$HOME/mosek.lic"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
export PATH="${PATH}:$HOME/mosek/11.2/tools/platform/linuxaarch64/bin"
export PATH="${PATH}:$HOME/julia-1.12.6/bin"
EOF
```

Then reload it (or just reconnect):

```sh
source ~/.bashrc
```

## 5. Clone the pooled-testing repository

```sh
git clone https://github.com/edwinlock/pooled-testing.git
```

## 6. Instantiate the project environment

The repository ships its own Julia environment (`Project.toml` + `Manifest.toml`).
Resolve (in case the manifest was last pinned with a different Julia version)
and download its dependencies once:

```sh
cd pooled-testing
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

## 7. Run experiments inside `tmux`

Running inside `tmux` keeps the experiments alive if your SSH connection drops —
otherwise the run is tied to your SSH session and is killed on disconnect.

Install `tmux` if needed and start a named session:

```sh
sudo yum install -y tmux
tmux new -s experiments
```

Inside the tmux session, run the experiments. Always pass `--project=.` so Julia
uses this environment, and `-t auto` so the parallel experiments (3–6) use all
available cores.

```sh
# Run all experiments
julia --project=. -t auto experiments.jl

# Run only specific experiments (comma-separated)
julia --project=. -t auto experiments.jl --experiments 1,3,6
```

Detach with `Ctrl-b` then `d` (the run keeps going); you can now disconnect SSH.
Reattach later with `tmux attach -t experiments`, or list sessions with `tmux ls`.

## Running on a Spot Instance (cheaper)

Spot instances cost ~60–90% less than on-demand, in exchange for AWS being able
to reclaim them with a 2-minute warning. This job is well suited to Spot because
it resumes cleanly: each experiment writes `data/expN-data.csv` on completion,
and re-running skips any experiment whose CSV already exists (pass `--force` to
re-run regardless). So an interrupted run just needs to be restarted.

1. Request the instance as Spot when launching — in the EC2 console tick
   "Request Spot Instances", or via the CLI:

   ```sh
   aws ec2 run-instances --instance-market-options '{"MarketType":"spot"}' ...
   ```

   Mid-size instances (e.g. `c8g.16xlarge` / `c8g.24xlarge`) usually have lower
   interruption rates than the largest sizes; check the Spot price history and
   interruption frequency for your region first.

2. Put the repository (and therefore the `data/` directory) on a persistent EBS
   volume, or sync results to S3 periodically, so completed CSVs survive an
   interruption and a fresh instance can resume.

3. After an interruption, relaunch and simply re-run the same command — completed
   experiments are skipped automatically:

   ```sh
   julia --project=. -t auto experiments.jl
   ```
