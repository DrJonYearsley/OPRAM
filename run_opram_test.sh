# Script to run the OPRAM models


nNodes=3

# Run the main simulation results
julia -p "$nNodes" OPRAM_main_program.jl parameters_test.toml


# Create 30yr average
julia  OPRAM_calculate_average.jl parameters_test.toml

# Create final files for the web app
julia  OPRAM_final_results.jl parameters_test.toml
