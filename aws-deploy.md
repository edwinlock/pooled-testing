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

```sh
# (a) Trust policy — lets EC2 assume the role. Write it to a file.
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

# (b) Permission policy — what the role may do (S3 access to the bucket).
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

# (c) Create the role with the trust policy.
aws iam create-role --role-name pooled-testing-s3 \
  --assume-role-policy-document file://trust-policy.json

# (d) Attach the S3 permission policy to the role (this is where the JSON goes).
aws iam put-role-policy --role-name pooled-testing-s3 \
  --policy-name s3-access --policy-document file://s3-policy.json

# (e) Create the instance profile and add the role to it. The instance profile
#     name (pooled-testing-s3) is what the launch command references.
aws iam create-instance-profile --instance-profile-name pooled-testing-s3
aws iam add-role-to-instance-profile \
  --instance-profile-name pooled-testing-s3 --role-name pooled-testing-s3
```

IAM is global, so these commands need no `--region`. (In the console the
equivalent is fiddlier: you must first create the *policy* under IAM → Policies
→ Create policy → JSON tab — paste the permission policy from (b) — then create
the role under IAM → Roles and attach that policy. The console creates the
matching instance profile for you automatically.)

### Launching the instance

The command launches an ARM Spot `c8g.24xlarge` using the latest Amazon Linux
2023 ARM image (resolved automatically, so no region-specific AMI ID is needed),
with the IAM role attached and a 30 GB root disk:

```sh
aws ec2 run-instances \
  --region us-east-2 \
  --image-id resolve:ssm:/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64 \
  --instance-type c8g.24xlarge \
  --instance-market-options '{"MarketType":"spot"}' \
  --key-name YOUR_KEY_NAME \
  --security-group-ids YOUR_SECURITY_GROUP_ID \
  --iam-instance-profile Name=pooled-testing-s3 \
  --block-device-mappings '[{"DeviceName":"/dev/xvda","Ebs":{"VolumeSize":30}}]' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=pooled-testing}]'
```

Replace `YOUR_KEY_NAME` (the EC2 key pair matching your `.pem` file) and
`YOUR_SECURITY_GROUP_ID` (a security group allowing inbound SSH on port 22 from
your IP). The key pair, security group, and IAM role must all exist in
`us-east-2`.

Get the public DNS name to SSH into:

```sh
aws ec2 describe-instances --region us-east-2 \
  --filters 'Name=tag:Name,Values=pooled-testing' 'Name=instance-state-name,Values=running' \
  --query 'Reservations[].Instances[].PublicDnsName' --output text
```

Optionally check the current Spot price/availability across zones first:

```sh
aws ec2 describe-spot-price-history --region us-east-2 \
  --instance-types c8g.24xlarge --product-descriptions "Linux/UNIX" \
  --query 'SpotPriceHistory[0:5].[AvailabilityZone,SpotPrice]' --output table
```

Then connect and set up the instance via the numbered steps below.

### Persisting results to S3

So results survive a termination, sync the output directories `data/`, `tables/`
and `figs/` to the bucket. The repository itself is not persisted (it is cloned
from git); `data/` is what the resume logic reads, `tables/` and `figs/` are the
paper outputs.

On a **fresh instance**, sync any earlier results down *before* running, so
completed experiments are recognised and skipped:

```sh
cd pooled-testing
for d in data tables figs; do aws s3 sync s3://pooled-testing-bucket/pooled-testing/$d $d; done
```

While experiments run, sync the outputs **up** periodically so a 2-minute
spot-termination warning never costs more than one interval. Run this in a
separate tmux window (`Ctrl-b c`):

```sh
cd pooled-testing
while true; do
  for d in data tables figs; do aws s3 sync $d s3://pooled-testing-bucket/pooled-testing/$d; done
  sleep 600   # every 10 minutes
done
```

Only *completed* experiments are recoverable — an experiment interrupted
mid-solve has no CSV yet and re-runs from scratch, so this loses at most the
in-progress experiment regardless of sync interval.

After an interruption, relaunch (the IAM role is set at launch above), sync the
results down, and re-run the same command — completed experiments are skipped
automatically.

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
experiments running in the first window (`htop` is installed by `setup.sh`).

With the experiments running in window 0:

- Press `Ctrl-b` then `c` to create a new window, then run `htop` in it.
- Switch between windows with `Ctrl-b` then `n` (next) or `p` (previous), or
  jump directly with `Ctrl-b` then the window number (`Ctrl-b 0` for the
  experiments, `Ctrl-b 1` for htop).
- Press `q` to quit `htop`, and `Ctrl-b` then `&` to close its window when done.

`Ctrl-b w` lists all windows so you can pick one from a menu.
