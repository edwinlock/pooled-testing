# Instructions to deploy to AWS EC2 instance

**Instance type:** use an ARM (Graviton) compute-optimised instance, e.g.
`c8g.24xlarge` (96 vCPUs) or `c8g.16xlarge` (64 vCPUs). These match the ARM
Julia/Gurobi/MOSEK builds downloaded below. Pick the size to fit the parallel
experiments: experiments 3–5 run their independent solves across Julia threads,
so launch Julia with about cores÷8 threads (see step 4). `c8g.48xlarge`
(192 vCPUs) is the largest, but is only worth it if you saturate it.

## Launching on a Spot Instance (cheaper)

Spot instances cost ~60–90% less than on-demand, in exchange for AWS being able
to reclaim them with a 2-minute warning. This job is well suited to Spot because
it resumes cleanly: each experiment writes `data/expN-data.csv` on completion,
and re-running skips any experiment whose CSV already exists (pass `--force` to
re-run regardless). So an interrupted run just needs to be restarted.

All commands below use the Ohio region (`us-east-2`, typically the cheapest).

### One-time prerequisites: S3 bucket and IAM role

These give the instance somewhere to store results and permission to do so. Set
them up **once**; they are account-level resources that survive terminations.

**1. Create the S3 bucket.** The commands below use `pooled-testing-bucket`.
Bucket names are globally unique across all AWS accounts, so if that one is
already taken (`BucketAlreadyExists`), pick a more specific name such as
`edwinlock-pooled-testing` and substitute it everywhere below. Outside
`us-east-1` the `LocationConstraint` is required.

```sh
aws s3api create-bucket \
  --bucket pooled-testing-bucket \
  --region us-east-2 \
  --create-bucket-configuration LocationConstraint=us-east-2
```

**2. Create an IAM role and instance profile for EC2**, scoped to just this
bucket (preferred over the broad `AmazonS3FullAccess`). The CLI is the most
reliable way — with the CLI the role, its permission policy, and the *instance
profile* (the container EC2 actually attaches) are separate objects, so there
are several commands. Run them in order:

Step A — write the trust policy that lets EC2 assume the role:

```sh
cat > trust-policy.json <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Principal": {"Service": "ec2.amazonaws.com"},
    "Action": "sts:AssumeRole"
  }]
}
EOF
```

Step B — write the S3 permission policy that grants bucket access:

```sh
cat > s3-policy.json <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Action": ["s3:GetObject", "s3:PutObject", "s3:ListBucket"],
    "Resource": [
      "arn:aws:s3:::pooled-testing-bucket",
      "arn:aws:s3:::pooled-testing-bucket/*"
    ]
  }]
}
EOF
```

Step C — create the role, attach the permission policy, then create the instance
profile and add the role to it. The instance profile name is what the launch
command references:

```sh
aws iam create-role --role-name pooled-testing-s3 --assume-role-policy-document file://trust-policy.json
aws iam put-role-policy --role-name pooled-testing-s3 --policy-name s3-access --policy-document file://s3-policy.json
aws iam create-instance-profile --instance-profile-name pooled-testing-s3
aws iam add-role-to-instance-profile --instance-profile-name pooled-testing-s3 --role-name pooled-testing-s3
```

IAM is global, so these commands need no `--region`. (In the console the
equivalent is fiddlier: you must first create the *policy* under IAM → Policies
→ Create policy → JSON tab — paste the permission policy from (b) — then create
the role under IAM → Roles and attach that policy. The console creates the
matching instance profile for you automatically.)

### Launching the instance

The launch command needs a **key pair** name and a **security group** ID, both
in `us-east-2`. Find them first.

List your key pairs (use the name matching your `.pem`; key pairs are
per-region, so it must exist in `us-east-2`):

```sh
aws ec2 describe-key-pairs --region us-east-2 --query 'KeyPairs[].KeyName' --output text
```

List security groups and their inbound rules, and pick one that allows SSH
(port 22). A console-created `launch-wizard-*` group usually does:

```sh
aws ec2 describe-security-groups --region us-east-2 \
  --query 'SecurityGroups[].{ID:GroupId,Name:GroupName,Ingress:IpPermissions}'
```

If none allows SSH, create one and open port 22 to your IP (find it with
`curl -s ifconfig.me`; use `0.0.0.0/0` to allow any IP, less secure):

```sh
aws ec2 create-security-group --region us-east-2 --group-name pooled-testing-sg --description "SSH for pooled-testing"
aws ec2 authorize-security-group-ingress --region us-east-2 --group-name pooled-testing-sg --protocol tcp --port 22 --cidr YOUR.IP.ADDRESS/32
```

Then launch — an ARM Spot `c8g.16xlarge` using the latest Amazon Linux 2023 ARM
image (resolved automatically, so no region-specific AMI ID is needed), with the
IAM role attached and a 30 GB root disk. Substitute the key pair name and
security group ID from above:

```sh
aws ec2 run-instances \
  --region us-east-2 \
  --image-id resolve:ssm:/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64 \
  --instance-type c8g.16xlarge \
  --instance-market-options '{"MarketType":"spot"}' \
  --key-name YOUR_KEY_NAME \
  --security-group-ids sg-XXXXXXXXXXXXXXXXX \
  --iam-instance-profile Name=pooled-testing-s3 \
  --block-device-mappings '[{"DeviceName":"/dev/xvda","Ebs":{"VolumeSize":30}}]' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=pooled-testing}]'
```

**On-demand instead of Spot.** Spot is cheaper but can be reclaimed (and large
ARM Spot capacity in a region is sometimes unavailable — you may get
`InsufficientInstanceCapacity`). For an uninterrupted run, drop the
`--instance-market-options` line to launch on-demand (guaranteed capacity, full
price, never reclaimed):

```sh
aws ec2 run-instances \
  --region us-east-2 \
  --image-id resolve:ssm:/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64 \
  --instance-type c8g.16xlarge \
  --key-name YOUR_KEY_NAME \
  --security-group-ids sg-XXXXXXXXXXXXXXXXX \
  --iam-instance-profile Name=pooled-testing-s3 \
  --block-device-mappings '[{"DeviceName":"/dev/xvda","Ebs":{"VolumeSize":30}}]' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=pooled-testing}]'
```

On-demand and Spot use separate vCPU quotas, so if Spot fails with
`MaxSpotInstanceCountExceeded` the on-demand launch may still work (and vice
versa). On-demand bills by the hour until you terminate it, so remember to shut
it down when the run finishes:

```sh
aws ec2 terminate-instances --region us-east-2 --instance-ids i-XXXXXXXXXXXXXXXXX
```

Get the public DNS name to SSH into:

```sh
aws ec2 describe-instances --region us-east-2 \
  --filters 'Name=tag:Name,Values=pooled-testing' 'Name=instance-state-name,Values=running' \
  --query 'Reservations[].Instances[].PublicDnsName' --output text
```

Optionally check the current Spot price/availability across zones first:

```sh
aws ec2 describe-spot-price-history --region us-east-2 \
  --instance-types c8g.16xlarge --product-descriptions "Linux/UNIX" \
  --query 'SpotPriceHistory[0:5].[AvailabilityZone,SpotPrice]' --output table
```

Then connect and set up the instance via the numbered steps below.

Results are persisted to S3 (steps 3 and 5 below) so they survive a termination:
the output directories `data/`, `tables/` and `figs/` are synced to the bucket.
The repository itself is not persisted (it is cloned from git); `data/` is what
the resume logic reads, `tables/` and `figs/` are the paper outputs. On-demand
instances are not reclaimed, so syncing is optional there but still guards
against accidental termination or a crash.

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

## 3. Sync any earlier results down from S3

On a fresh instance, pull down results from previous runs *before* running, so
completed experiments are recognised and skipped (skip this on the very first
run, when the bucket is empty):

```sh
cd pooled-testing
for d in data tables figs; do aws s3 sync s3://pooled-testing-bucket/pooled-testing/$d $d; done
```

## 4. Run experiments inside `tmux`

Running inside `tmux` keeps the experiments alive if your SSH connection drops —
otherwise the run is tied to your SSH session and is killed on disconnect
(`tmux` is installed by `setup.sh`). Start a named session:

```sh
tmux new -s experiments
```

Inside the tmux session, run the experiments. Always pass `--project=.` so Julia
uses this environment.

**Choosing the thread count.** Each MILP solve runs Gurobi with `Threads=8`
(its sweet spot — more is slower for these models), and the experiments run
several solves in parallel across Julia threads. To use all cores *without*
oversubscribing, launch Julia with about **cores ÷ 8** threads — e.g. `-t 12`
on a 96-vCPU `c8g.24xlarge`, `-t 8` on a 64-vCPU `c8g.16xlarge`. Do **not** use
`-t auto`: that starts one Julia thread per core, and with 8 Gurobi threads each
it oversubscribes badly.

First validate that the parallel solves run cleanly on this many-core machine by
running just experiment 3 (the multithreaded, MOSEK-heavy one):

```sh
julia --project=. -t 8 experiments.jl --experiments 3
```

Once that completes, start the full run:

```sh
# Run all experiments
julia --project=. -t 8 experiments.jl

# Or only specific experiments (comma-separated)
julia --project=. -t 8 experiments.jl --experiments 1,3,5
```

Detach with `Ctrl-b` then `d` (the run keeps going); you can now disconnect SSH.
Reattach later with `tmux attach -t experiments`, or list sessions with `tmux ls`.

## 5. Sync results back up to S3

While the experiments run, sync the outputs **up** periodically so an
interruption never costs more than one interval. Run this in a separate tmux
window (`Ctrl-b c`):

```sh
cd pooled-testing
while true; do
  for d in data tables figs; do aws s3 sync $d s3://pooled-testing-bucket/pooled-testing/$d; done
  sleep 600   # every 10 minutes
done
```

Only *completed* experiments are recoverable — an experiment interrupted
mid-solve has no CSV yet and re-runs from scratch, so this loses at most the
in-progress experiment regardless of sync interval. After an interruption,
relaunch, repeat from step 1 (the sync-down in step 3 restores finished
experiments), and re-run — completed experiments are skipped automatically.

### Monitoring CPU usage with `htop`

To watch CPU/memory load (e.g. to confirm all cores are busy) while the
experiments run, open a second tmux window and start `htop` there, leaving the
experiments running in the first window (`htop` is installed by `setup.sh`).

With the experiments running in window 0:

- Press `Ctrl-b` then `c` to create a new window, then run `htop` in it.
- Switch between windows with `Ctrl-b` then `n` (next) or `p` (previous), or
  jump directly with `Ctrl-b` then the window number (`Ctrl-b 0` for the
  experiments, `Ctrl-b 1` for htop).
- Press `q` to quit `htop`, and `Ctrl-b` then `&` to close its window when done.

`Ctrl-b w` lists all windows so you can pick one from a menu.
