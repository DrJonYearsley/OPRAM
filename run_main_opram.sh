# Script to run the OPRAM models


nNodes=5

# Run the main simulation results
julia --procs "$nNodes" julia/OPRAM_main_program.jl julia/parameters_userdefined.toml


