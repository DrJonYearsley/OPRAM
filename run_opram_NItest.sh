# Script to run the OPRAM models for Northern Ireland using Met Eireann data
#
# This is used to compare the results to the UK pest tool data
# The same script can be run using Met Eireann data


nNodes=3

# Run the main simulation results
julia -p "$nNodes" julia/OPRAM_main_program.jl julia/parameters_NI_test.toml



# Create final files
julia  julia/OPRAM_calculate_emergeDOY.jl julia/parameters_NI_test.toml
