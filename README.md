# Instructions
1. Install Julia (from www.julialang.org)
2. Install Gurobi (https://www.gurobi.com/). You may need to create a Gurobi academic account, and get a licence. Follow Gurobi instructions to store the licence locally on your computer. Then follow installation instructions for Gurobi/JuMP/Julia here: https://github.com/jump-dev/Gurobi.jl
3. Install MOSEK (https://mosek.com).
4. Install Julia packages: start Julia in the terminal, press `]` to enter package mode and type `add JuMP, Gurobi, DataStructures, Roots, Combinatorics, CSV, DataFrames, Statistics, Dates, StatsPlots`. Hit the return key to exit package mode.
5. Run the experiments by typing `julia experiments.jl` in the terminal.
