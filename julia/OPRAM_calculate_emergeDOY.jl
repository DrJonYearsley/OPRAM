# Use OPRAM model results to calculate dat of emergence with start date of 1st Jan


# The degree day model results are generatde by OPRAM_main_program.jl
#
#
# Author: Jon Yearsley  (jon.yearsley@ucd.ie)
# Date: Sept 2025
#
# =============================================================================
# =============================================================================

# Load the required packages
using CSV
using DataFrames
using Statistics
using Dates
using JLD2
using Distributed
using SharedArrays
using TOML



include("OPRAM_io_functions.jl");
include("OPRAM_processing_functions.jl");


# ============================================================================================
# =============== Import parameter values =======================================================




# Model parameters are stored in a TOML file https://toml.io/en/
if length(ARGS) == 1
    nNodes, run_params, species_setup, paths = import_parameters(ARGS[1], true)

elseif length(ARGS) == 0 & isfile("parameters_test.toml")
    nNodes, run_params, species_setup, paths = import_parameters("parameters_test.toml", true)

else
    @error "No parameter file given"
end


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set species parameters from the parameter file (can be more than one species)
if species_setup.speciesStr[1] == "all"
    @info "Importing all species from the species file" * string(species_setup.speciesStr)
    species_params = import_species(species_setup.speciesFile, species_setup.speciesStr[1])
else
    @info "Importing species from the species file: " * string(species_setup.speciesStr)
    species_params = [import_species(species_setup.speciesFile, species_setup.speciesStr[s]) for s in eachindex(species_setup.speciesStr)]
end


# =============== End of parameter setup =====================================================
# ============================================================================================



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import location data and thin it using thinFactor
grid = read_grid(run_params)




for s in eachindex(species_params)
    # Find directory matching the species name in run_params
    regex = Regex(replace(lowercase(species_params[s].species_name),
        r"\s" => "\\w"))  # Replace spaces with reg expression
    speciesName = filter(x -> occursin(regex, x), readdir(paths.resultsDir))
    if length(speciesName) > 1
        @error "More than one species name found"
    elseif length(speciesName) == 0
        @error "Species not found"
    end

    @info "Importing data for $(speciesName[1])"


    # Create vector of files to import (one file for each year)
    inFiles = Vector{String}(undef, length(run_params.years))
    for y in eachindex(run_params.years)
        inFile = filter(x -> occursin(r"^" * speciesName[1] * "_" * run_params.country * "_" * string(run_params.years[y]) * "_1km.jld2", x),
            readdir(joinpath(paths.resultsDir, speciesName[1])))

        if length(inFile) > 1
            @error "More than one input file found"
        else
            inFiles[y] = joinpath(paths.resultsDir, speciesName[1], inFile[1])
        end
    end

    # Import the data from the jld2 files
    df_1km = read_OPRAM_JLD2(inFiles, run_params.years, grid)

    # Retain only start date of 1st Jan
    filter!(row -> row.startDOY == 1, df_1km)

    # =========================================================
    # Write result to CSV file
    @info "Writing 1km results to CSV files"
    # Add eastings and northings into the 1km data frames (copy the data frames so we can round before saving)
    out_1km = rightjoin(grid[:, [:ID, :east, :north]], df_1km, on=[:ID])


    # Wrie CSV file
    fileout2 = joinpath(paths.outDir, speciesName[1], speciesName[1] * "_" * "01_Jan" * "_" *
                                                      string(run_params.years[y]) *
                                                      "_" * string(run_params.thinFactor) * "km.csv")
    CSV.write(fileout2, out_1km, missingstring="NA")
end
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++