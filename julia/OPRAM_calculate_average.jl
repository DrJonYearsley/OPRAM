#!//opt/homebrew/bin/julia -p 3
#
##!//home/jon/.juliaup/bin/julia -p 7
#
# Julia program to:
#  extract results from model jld2 files 
#  calculate multi-year average of degree day model results
#  export results to CSV files
#
# The degree day model results are generatde by OPRAM_main_program.jl
#
#
# Author: Jon Yearsley  (jon.yearsley@ucd.ie)
# Date: March 2025
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
    nNodes, run_params, species_params, paths = import_parameters(ARGS[1], true)

elseif length(ARGS) == 0 & isfile("parameters.toml")
    nNodes, run_params, species_setup, paths = import_parameters("parameters.toml", true)

else
    @error "No parameter file given"
end

# Files containing the ID system from Granite.ie
# This info is used to package the output into separate files
granite_hectad_ID = "git_repos/OPRAM/data/granite_hectad_defs.csv"
granite_county_ID = "git_repos/OPRAM/data/granite_county_defs.csv"
granite_hectad_county = "git_repos/OPRAM/data/granite_hectad_county_defs.csv"



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
grid = read_grid(run_params.gridFile, run_params.thinFactor, run_params.country)

# Create easting and northings of bottom left of a hectad
grid.east_hectad = convert.(Int32, floor.(grid.east ./ 1e4) .* 1e4);
grid.north_hectad = convert.(Int32, floor.(grid.north ./ 1e4) .* 1e4);







# =========================================================
# =========================================================
# Import data across the years and calculate the average


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



    # =========================================================
    # Calculate average across years (use median to minimise outlier influence)

    @info "Calculating median values across years"


    # # Set missing DOY to be a large number and missing generations to be -9999
    # # This is just for the quantile calculations
    # idx = ismissing.(df_1km.emergeDOY)
    # df_1km.emergeDOY[idx] .= 9999
    # df_1km.nGenerations[idx] .= -9999

    # Group data frame by location and then starting DOY/Date 
    df_group = groupby(df_1km, [:ID, :startMonth])



    # =========================================================
    # ################# Median across the years ##############
    # Calculate the median of nGenerations and emergeDOY
    d_agg = combine(df_group,
        :nGenerations => (x -> if sum(.!isa.(x,Missing))>0 quantile(skipmissing(x), 0.5) else missing end) => :nGenerations_median,
        :emergeDOY => (x -> if sum(.!isa.(x,Missing))>0 quantile(skipmissing(x), 0.5) else missing end) => :emergeDOY_median,
        :emergeDOY => (x -> sum(ismissing.(x))) => :nMissing)



    # # Return nGenerations and emergeDOY to missing
    # df_1km.nGenerations[idx] .= missing
    # df_1km.emergeDOY[idx] .= missing

    # If emergeDOY is greater than maximum possible DOY (366*maxYears) then set
    # nGenerations_median and emergeDOY_median to missing
    idx_agg = d_agg.emergeDOY_median .> 366 * run_params.maxYears
    allowmissing!(d_agg, [:nGenerations_median, :emergeDOY_median])
    d_agg.nGenerations_median[idx_agg] .= missing
    d_agg.emergeDOY_median[idx_agg] .= missing




    # =========================================================
    # =========================================================
    # Export the anomaly and multi year average

    @info "Writing 1km results to CSV files"



    # Add eastings and northings into the 1km data frames (copy the data frames so we can round before saving)
    out_agg1km = rightjoin(grid[:, [:ID, :east, :north]], d_agg, on=[:ID])

    # Round number of generations and emergeDOY
    out_agg1km.nGenerations_median = round.(out_agg1km.nGenerations_median, digits=2)
    out_agg1km.emergeDOY_median = round.(out_agg1km.emergeDOY_median, digits=2)

    # Remove unwanted variables
    select!(out_agg1km, Not(:nMissing))

    # Write out the average across years data frame
    fileout2 = joinpath(paths.outDir, speciesName[1], "average_" * speciesName[1] * "_" *
                                                      string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) *
                                                      "_" * string(run_params.thinFactor) * "km.csv")
    CSV.write(fileout2, out_agg1km, missingstring="NA")

    # Clear variables before moving on to next species
    out_agg1km = nothing
    d_agg = nothing
    df_1km = nothing
    df_group = nothing

    # Print blank line on the screen
    println(" ")
end
