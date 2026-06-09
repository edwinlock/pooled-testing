#!/usr/bin/env bash
#
# Set up an AWS EC2 (Amazon Linux, ARM/Graviton) instance to run the
# pooled-testing experiments. Run this once on a fresh instance, after the
# Gurobi and MOSEK licence files have been copied to the home directory.
#
#   ./setup.sh
#
# It installs OS packages, downloads and unpacks Julia/Gurobi/MOSEK, writes the
# required environment variables to ~/.bashrc, and instantiates the Julia
# project. It does NOT run the experiments — start those yourself inside tmux
# (see aws-deploy.md).

set -euo pipefail

# ---- Versions (bump these to upgrade) --------------------------------------
JULIA_VERSION="1.12.6"
JULIA_SERIES="1.12"
GUROBI_VERSION="13.0.1"
GUROBI_SERIES="13.0"
GUROBI_DIR="gurobi1301"        # directory the Gurobi tarball unpacks to
MOSEK_VERSION="11.2.1"
MOSEK_SERIES="11.2"

# Resolve paths relative to this script's location (the repo root).
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "==> Checking for licence files"
missing=0
for lic in "$HOME/gurobi.lic" "$HOME/mosek.lic"; do
    if [[ ! -f "$lic" ]]; then
        echo "ERROR: $lic not found." >&2
        missing=1
    fi
done
if [[ "$missing" -ne 0 ]]; then
    echo "Copy the licences from your local machine first, e.g.:" >&2
    echo "  scp -i KEY.pem ~/mosek.lic ~/gurobi.lic ec2-user@<host>:~/" >&2
    exit 1
fi

echo "==> Installing OS packages"
sudo yum install -y bzip2 git tmux htop wget

echo "==> Downloading and unpacking Julia, Gurobi and MOSEK"
cd "$HOME"

if [[ ! -d "$HOME/julia-$JULIA_VERSION" ]]; then
    wget -q "https://julialang-s3.julialang.org/bin/linux/aarch64/$JULIA_SERIES/julia-$JULIA_VERSION-linux-aarch64.tar.gz"
    tar zxf "julia-$JULIA_VERSION-linux-aarch64.tar.gz"
fi

if [[ ! -d "$HOME/$GUROBI_DIR" ]]; then
    wget -q "https://packages.gurobi.com/$GUROBI_SERIES/gurobi${GUROBI_VERSION}_armlinux64.tar.gz"
    tar zxf "gurobi${GUROBI_VERSION}_armlinux64.tar.gz"
fi

if [[ ! -d "$HOME/mosek" ]]; then
    wget -q "https://download.mosek.com/stable/$MOSEK_VERSION/mosektoolslinuxaarch64.tar.bz2"
    tar xjf "mosektoolslinuxaarch64.tar.bz2"
fi

echo "==> Writing environment variables to ~/.bashrc"
# Idempotent: only append our block once, marked by a sentinel comment.
SENTINEL="# >>> pooled-testing env >>>"
if ! grep -qF "$SENTINEL" "$HOME/.bashrc" 2>/dev/null; then
    cat >> "$HOME/.bashrc" <<EOF

$SENTINEL
export GUROBI_HOME="\$HOME/$GUROBI_DIR/armlinux64"
export GRB_LICENSE_FILE="\$HOME/gurobi.lic"
export MOSEKLM_LICENSE_FILE="\$HOME/mosek.lic"
export PATH="\${PATH}:\${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="\${LD_LIBRARY_PATH:-}:\${GUROBI_HOME}/lib"
export PATH="\${PATH}:\$HOME/mosek/$MOSEK_SERIES/tools/platform/linuxaarch64/bin"
export PATH="\${PATH}:\$HOME/julia-$JULIA_VERSION/bin"
# <<< pooled-testing env <<<
EOF
    echo "    appended environment block"
else
    echo "    environment block already present, skipping"
fi

# Make the new variables available to the rest of this script.
export PATH="$PATH:$HOME/julia-$JULIA_VERSION/bin"

echo "==> Instantiating the Julia project"
cd "$REPO_DIR"
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'

echo
echo "==> Setup complete."
echo "Open a new shell (or run: source ~/.bashrc) so the environment is loaded,"
echo "then start the experiments inside tmux, e.g.:"
echo
echo "    tmux new -s experiments"
echo "    cd $REPO_DIR"
echo "    julia --project=. -t auto experiments.jl --experiments 3   # validate first"
echo "    julia --project=. -t auto experiments.jl                   # full run"
