# Script to run the OPRAM models and generate geotiffs
# This assumes that the temperature data has been imported
# and saved as a jld2 file using
#     meteo_to_jld2.jl     for historice data
#     translate_to_jld2.jl for translate future scenarios
#
# There are 4 steps
#     1. Run the main OPRAM model to produce results from
#         degree day models (results in jld2 file)
#     2. Create results that average across years (if required)
#         (usually 30yr average 1991-2020)
#     3. Create final results for 12 specific start dates (biofixes)
#         (results stored in CSV files)
#     4. Ceate geotiff versions of the final results (if requrired)
#
# The main program can run on multiple compute nodes
# The other programs run on a single node
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# May 2026
# ***********************************************


# Define parameters ******************

# Set the name of the main parameter file
paramFile="julia/parameters_test.toml"

nNodes=3   # Set the number of compute nodes
average=0  # Set to 1 to calculate 1991-2020 average
geotiff=1  # Set to 1 to save geotiff version of results

# End of parameter definitions *******


# Do not edit below this line ******************

# Check whether julia is installed
if ! command -v julia >/dev/null 2>&1
then
    echo "julia is not installed"
    exit 1
else
    # If julia is installed...
    
    # Run the main simulation results
    julia -p "$nNodes" julia/OPRAM_main_program.jl "$paramFile"


    if [average -eq 1]; then
      # Create 30yr average (may not be required)
      julia  julia/OPRAM_calculate_average.jl "$paramFile"
    fi
   
    # Create final files for the web app
    julia  julia/OPRAM_final_results.jl "$paramFile"

fi




# If the geotiffs are to be created and R is installed
# then run the R script
if [geotiff -eq 1]; then
    if ! command -v R >/dev/null 2>&1
    then
        echo "R is not installed"
        exit 1
     else
        # Generate geoTIFF file version of final results
        R CMD BATCH  "--args $paramFile" R/OPRAM_generate_geotiff.R 
    fi
fi
