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
#
# Functions checked by JY: 30th Oct 2025
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

elseif length(ARGS) == 0 & isfile("parameters.toml")
    nNodes, run_params, species_setup, paths = import_parameters("parameters_test.toml", true)

else
    @error "No parameter file given"
end


# Check that this is for past climate, not future climate
if run_params.TRANSLATE_future
    @error "Average calculation not available for future climates"
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
        inFile = filter(x -> occursin(r"^" * speciesName[1] * "_" * run_params.country * "_" * string(run_params.years[y]) * "_1km_" * run_params.method * ".jld2", x),
            readdir(joinpath(paths.resultsDir, speciesName[1])))

        if length(inFile) > 1
            @error "More than one input file found"
        else
            inFiles[y] = joinpath(paths.resultsDir, speciesName[1], inFile[1])
        end
    end

    # Import the data from the jld2 files. 
    # This will set locations with no emergence as having emergeDOY=missing and nGenerations=0
    df_1km = read_OPRAM_JLD2(inFiles, run_params.years, grid)



    # =========================================================
    # Calculate average across years (use median to minimise outlier influence)

    @info "Calculating median values across years"

    # Group data frame by location and then starting DOY/Date 
    df_group = groupby(df_1km, [:ID, :startMonth])


    # =========================================================
    # ################# Median across the years ##############
    # Calculate the median of nGenerations and emergeDOY
    d_agg = combine(df_group,
        :nGenerations => (x -> if sum(.!isa.(x,Missing))>0 quantile(skipmissing(x), 0.5) else 0.0 end) => :nGenerations_median,
        :emergeDOY => (x -> if sum(.!isa.(x,Missing))>0 quantile(skipmissing(x), 0.5) else missing end) => :emergeDOY_median,
        :emergeDOY => (x -> sum(ismissing.(x))) => :nMissing)


    # If emergeDOY is greater than maximum possible DOY (366*maxYears) then set
    # nGenerations_median = 0 and emergeDOY_median = missing
    idx_agg = ismissing.(d_agg.emergeDOY_median) .|| d_agg.emergeDOY_median .> 366 * run_params.maxYears

    if any(idx_agg)
        d_agg.nGenerations_median[idx_agg] .= 0.0
        d_agg.emergeDOY_median[idx_agg] .= missing
    end



    # =========================================================
    # =========================================================
    # Export the  multi year average to the resultsDir

    @info "Writing 1km results to CSV files"

    # Add eastings and northings into the 1km data frames (copy the data frames so we can round before saving)
    out_agg1km = rightjoin(grid[:, [:ID, :east, :north]], d_agg, on=[:ID])

    # Round number of generations and emergeDOY
    out_agg1km.nGenerations_median = round.(out_agg1km.nGenerations_median, digits=2)
    out_agg1km.emergeDOY_median = round.(out_agg1km.emergeDOY_median, digits=2)

    # Remove unwanted variables
    select!(out_agg1km, Not(:nMissing))

    # Write out the average across years data frame
    fileout2 = joinpath(paths.resultsDir, speciesName[1], "average_" * speciesName[1] * "_" *
                                                      string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) *
                                                      "_" * string(run_params.thinFactor) * "km_" * run_params.method * ".csv")
    CSV.write(fileout2, out_agg1km, missingstring="NA")

    # Clear variables before moving on to next species
    out_agg1km = nothing
    d_agg = nothing
    df_1km = nothing
    df_group = nothing

    # Print blank line on the screen
    println(" ")
end
