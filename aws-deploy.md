# Instructions to deploy to AWS EC2 instance

0. Install bzip2

    sudo yum install bzip2

1. Download and unpack Julia, Gurobi and MOSEK

**Julia**

    wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz
    tar zxvf julia-1.8.5-linux-x86_64.tar.gz

**Gurobi**

    wget https://packages.gurobi.com/10.0/gurobi10.0.1_armlinux64.tar.gz
    tar zxvf gurobi10.0.1_armlinux64.tar.gz

**MOSEK**

    wget https://download.mosek.com/stable/10.0.40/mosektoolslinuxaarch64.tar.bz2
    tar zxjf mosektoolslinuxaarch64.tar.bz2

2. Add paths to ~/.bashrc

    export GUROBI_HOME="~/gurobi1001/linux64"
    export GRB_LICENSE_FILE="~/gurobi1001/gurobi.lic"
    export PATH="${PATH}:${GUROBI_HOME}/bin"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
    export PATH="${PATH}:~/mosek/10.0/tools/platform/linux64x86/bin"
    export PATH="${PATH}:~/julia-1.8.5/bin"

3. Copy over MOSEK and Gurobi licence files

    scp -i ~/Downloads/pooled_testing.pem ~/mosek/mosek.lic ec2-user@ec2-52-56-33-36.eu-west-2.compute.amazonaws.com:/home/ec2-user/mosek/
    scp -i ~/Downloads/pooled_testing.pem ~/Downloads/gurobi.lic ec2-user@ec2-52-56-33-36.eu-west-2.compute.amazonaws.com:/home/ec2-user/gurobi1001/

5. Clone pooled-testing repository

    git clone https://github.com/edwinlock/pooled-testing.git

6. Run experiment!

    cd pooled-testing
    julia experiments.jl