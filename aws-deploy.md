# Instructions to deploy to AWS EC2 instance

**Instance type:** use an ARM (Graviton) compute-optimised instance, e.g.
`c8g.24xlarge` (96 vCPUs) or `c8g.16xlarge` (64 vCPUs). These match the ARM
Julia/Gurobi/MOSEK builds downloaded below. Pick the size to fit the parallel
experiments: experiments 3–5 run their independent solves across Julia threads,
so launch Julia with `-t auto` (see step 3) to use all the cores. `c8g.48xlarge`
(192 vCPUs) is the largest, but is only worth it if you saturate it.

## Launching on a Spot Instance (cheaper)

Spot instances cost ~60–90% less than on-demand, in exchange for AWS being able
to reclaim them with a 2-minute warning. This job is well suited to Spot because
it resumes cleanly: each experiment writes `data/expN-data.csv` on completion,
and re-running skips any experiment whose CSV already exists (pass `--force` to
re-run regardless). So an interrupted run just needs to be restarted.

Decide on Spot before launching, then follow the numbered setup steps below.

1. Request the instance as Spot when launching — in the EC2 console tick
   "Request Spot Instances", or via the CLI:

   ```sh
   aws ec2 run-instances --instance-market-options '{"MarketType":"spot"}' ...
   ```

   Mid-size instances (e.g. `c8g.16xlarge` / `c8g.24xlarge`) usually have lower
   interruption rates than the largest sizes; check the Spot price history and
   interruption frequency for your region first.

2. Persist the experiment outputs to S3 so they survive a termination. The
   repository itself does not need persisting (it is cloned from git); only the
   generated output directories `data/`, `tables/` and `figs/` matter. `data/`
   is what the resume logic reads; `tables/` and `figs/` are the paper outputs.

   On a **fresh instance**, sync any earlier results down *before* running, so
   completed experiments are recognised and skipped:

   ```sh
   cd pooled-testing
   for d in data tables figs; do aws s3 sync s3://YOUR_BUCKET/pooled-testing/$d $d; done
   ```

   While experiments run, sync the outputs **up** periodically so a 2-minute
   spot-termination warning never costs more than one interval. Run this in a
   separate tmux window (`Ctrl-b c`):

   ```sh
   cd pooled-testing
   while true; do
     for d in data tables figs; do aws s3 sync $d s3://YOUR_BUCKET/pooled-testing/$d; done
     sleep 600   # every 10 minutes
   done
   ```

   Note: only *completed* experiments are recoverable — an experiment
   interrupted mid-solve has no CSV yet and re-runs from scratch, so this loses
   at most the in-progress experiment regardless of sync interval.

   For the instance to access S3, attach an **IAM role** to it (see below). This
   gives the AWS CLI temporary, auto-rotating credentials with no keys stored on
   the box — preferable to `aws configure`, especially on a disposable Spot
   instance.

3. After an interruption, relaunch, sync the results down (step 2), and re-run
   the same command — completed experiments are skipped automatically:

   ```sh
   julia --project=. -t auto experiments.jl
   ```

### Setting up the S3 bucket and IAM role (one-time)

These give the instance somewhere to store results and permission to do so. You
only set this up once; afterwards you just attach the role to each instance.

**Create a bucket** (S3 → Create bucket). Pick a globally-unique name and the
same region as your instance (to avoid cross-region transfer costs). This name
replaces `YOUR_BUCKET` in the sync commands above.

**Create an IAM role for EC2** (IAM → Roles → Create role → trusted entity
*AWS service* → *EC2*). Attach a policy scoped to just this bucket rather than
the broad `AmazonS3FullAccess`:

```json
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Action": ["s3:GetObject", "s3:PutObject", "s3:ListBucket"],
    "Resource": [
      "arn:aws:s3:::YOUR_BUCKET",
      "arn:aws:s3:::YOUR_BUCKET/*"
    ]
  }]
}
```

Name the role (e.g. `pooled-testing-s3`) and create it.

**Attach the role to the instance**: either at launch under "IAM instance
profile", or on a running instance via EC2 → Actions → Security → *Modify IAM
role*. The role is per-instance, so re-attach it (or set it at launch) each time
you start a fresh Spot instance. Once attached, the `aws s3 sync` commands above
work with no further credential setup.

## 1. Connect and copy over the licence files

SSH in using your key pair and the instance's **public** DNS name or IP (find it
in the EC2 console under "Public IPv4 DNS" — the internal `ip-…compute.internal`
name is not reachable from your laptop). Amazon Linux logs in as `ec2-user`:

```sh
ssh -i pooled_testing.pem ec2-user@<public-dns-name>
```

Then, **from your local machine**, copy both licences to the instance's home
directory (the setup script expects them there):

```sh
scp -i pooled_testing.pem ~/mosek.lic ~/gurobi.lic ec2-user@<public-dns-name>:~/
```

## 2. Run the setup script

On the instance, install git, clone the repository, and run `setup.sh`. It
installs the OS packages, downloads and unpacks Julia/Gurobi/MOSEK, writes the
environment variables to `~/.bashrc`, and instantiates the Julia project. It
checks that the licences are present and stops with a clear message if not.

```sh
sudo yum install -y git
git clone https://github.com/edwinlock/pooled-testing.git
cd pooled-testing
./setup.sh
source ~/.bashrc   # load the environment into the current shell
```

To upgrade the Julia/Gurobi/MOSEK versions, edit the version variables at the
top of `setup.sh`.

## 3. Run experiments inside `tmux`

Running inside `tmux` keeps the experiments alive if your SSH connection drops —
otherwise the run is tied to your SSH session and is killed on disconnect
(`tmux` is installed by `setup.sh`). Start a named session:

```sh
tmux new -s experiments
```

Inside the tmux session, run the experiments. Always pass `--project=.` so Julia
uses this environment, and `-t auto` so the parallel experiments (3–5) use all
available cores.

First validate that the parallel solves run cleanly on this many-core machine by
running just experiment 3 (the multithreaded MOSEK-heavy one):

```sh
julia --project=. -t auto experiments.jl --experiments 3
```

Once that completes, start the full run:

```sh
# Run all experiments
julia --project=. -t auto experiments.jl

# Or only specific experiments (comma-separated)
julia --project=. -t auto experiments.jl --experiments 1,3,5
```

Detach with `Ctrl-b` then `d` (the run keeps going); you can now disconnect SSH.
Reattach later with `tmux attach -t experiments`, or list sessions with `tmux ls`.

### Monitoring CPU usage with `htop`

To watch CPU/memory load (e.g. to confirm all cores are busy) while the
experiments run, open a second tmux window and start `htop` there, leaving the
experiments running in the first window.

Install `htop` if needed:

```sh
sudo yum install -y htop
```

With the experiments running in window 0:

- Press `Ctrl-b` then `c` to create a new window, then run `htop` in it.
- Switch between windows with `Ctrl-b` then `n` (next) or `p` (previous), or
  jump directly with `Ctrl-b` then the window number (`Ctrl-b 0` for the
  experiments, `Ctrl-b 1` for htop).
- Press `q` to quit `htop`, and `Ctrl-b` then `&` to close its window when done.

`Ctrl-b w` lists all windows so you can pick one from a menu.
