# Script to run the OPRAM models


nNodes=3
paramFile = "parameters_test.toml"

# Run the main simulation results
julia -p "$nNodes" OPRAM_main_program.jl parameters_test.toml


# Create 30yr average (may not be required)
julia  OPRAM_calculate_average.jl parameters_test.toml

# Create final files for the web app
julia  OPRAM_final_results.jl parameters_test.toml


# Generate geoTIFF file version of final results
R CMD BATCH  "--args parameters_test.toml" ../R/OPRAM_generate_geotiff.R 

