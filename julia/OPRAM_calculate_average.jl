#!//home/jon/.juliaup/bin/julia -p 7
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

# ==============================================================
# Model parameters are stored in a TOML file https://toml.io/en/
if length(ARGS) == 1
    params = TOML.parsefile(ARGS[1])

elseif length(ARGS) == 0 & isfile("parameters.toml")
    params = TOML.parsefile("parameters_userdefined.toml")

else
    @error "No parameter file given"
end



# Perform some checks on file names
# Check grid file exists
if !isfile(params["inputData"]["gridFile"])
    @info "Can't find grid file!"

    if isfile(joinpath(homedir(), "DATA", "OPRAM", params["inputData"]["gridFile"]))  # Guess 1
        params["inputData"]["gridFile"] = joinpath(homedir(), params["inputData"]["gridFile"])
        @info "gridFile set to" * params["inputData"]["gridFile"]

    elseif isfile(joinpath(homedir(), params["inputData"]["gridFile"]))  # Guess 2
        params["inputData"]["gridFile"] = joinpath(homedir(), params["inputData"]["gridFile"])
        @info "gridFile set to" * params["inputData"]["gridFile"]
    else
        @error "No grid file found"
    end
end



# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
# "frugiperda", "duplicatus", "cembrae", "sexdentatus"

run_params = (
    speciesName=params["model"]["speciesList"],         # Name of the species
    years=params["model"]["thirty_years"][1]:params["model"]["thirty_years"][2],                                  # the years collect(year1:year2)
    maxYears=params["model"]["maxYears"],           # Maximum number of years to complete insect development (must correspond to simulation value)
    # lastMeteoYear=params["model"]["lastMeteoYear"], # Maximum number of years to complete insect development
    country=params["model"]["country"],             # Can be "IE", "NI" or "AllIreland"
    thinFactor=params["model"]["thinFactor"],       # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
    gridFile=params["inputData"]["gridFile"]);      # File containing a 1km grid of lats and longs over Ireland 
# (used for daylength calculations as well as importing and thining of meteo data)







# =================================================================================
# Specify directories for import and export of data




if isdir("/media/jon/Seagate_5TB/OPRAM_results")       # Linux workstation
    paths = (outDir="/media/jon/Seagate_5TB/OPRAM_results",
        dataDir="/media/jon/Seagate_5TB/OPRAM_results")
elseif isdir(joinpath(homedir(), "DATA//OPRAM"))       # Linux workstation
    paths = (outDir=joinpath(homedir(), "Desktop//OPRAM//results//"),
        dataDir=joinpath(homedir(), "DATA//OPRAM//"))

elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//R"))   # Mac
    paths = (outDir=joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//results//"),
        dataDir=homedir())

elseif isdir(joinpath(homedir(), "Desktop//OPRAM//"))
    paths = (outDir=joinpath(homedir(), "Desktop//OPRAM//results//"),
        dataDir=joinpath(homedir(), "Desktop//OPRAM//"))
end

include("OPRAM_io_functions.jl");
include("OPRAM_processing_functions.jl");





# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import species info
# Work out unique species to simulate (Remove obvious duplicates)
speciesStr = Vector{String}(undef, 0)
for s in eachindex(params["model"]["speciesList"])
    global speciesStr = unique(vcat(speciesStr, params["model"]["speciesList"][s]))
end
species_params = (speciesFile=joinpath(homedir(), params["inputData"]["speciesFile"]),  # File containing species parameters
    speciesStr=speciesStr)  # A vector of strings to uniquely identify a species name in the speciesFile

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set species parameters from the parameter file (can be more than one species)
if species_params.speciesStr[1] == "all"
    @info "Importing all species from the species file" * string(species_params.speciesStr)
    species_params = import_species(species_params.speciesFile, species_params.speciesStr[1])
else
    @info "Importing species from the species file: " * string(species_params.speciesStr)
    species_params = [import_species(species_params.speciesFile, species_params.speciesStr[s]) for s in eachindex(species_params.speciesStr)]
end



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import location data and thin it using thinFactor
grid = read_grid(run_params.gridFile, run_params.thinFactor, run_params.country)

# Create easting and northings of bottom left of a hectad
grid.east_hectad = convert.(Int32, floor.(grid.east ./ 1e4) .* 1e4)
grid.north_hectad = convert.(Int32, floor.(grid.north ./ 1e4) .* 1e4)







# =========================================================
# =========================================================
# Import data across the years and calculate the average


for s in eachindex(species_params)

    if ismissing(species_params[s].species_name)
        @warn "Skipping an undefined species"
    else
        @info "\n#######  Importing data for species " * species_params[s].species_name * " #######"

        df_1km = doy_results(run_params.years, species_params[s].species_name, run_params, paths)


        # =========================================================
        # Calculate average across years (use median to minimise outlier influence)

        @info "Calculating median values across years"

        # Define DOY for start and emerge dates (and set as integer)
        idx = .!ismissing.(df_1km.emergeDate)
        df_1km.emergeDOY = Vector{Union{Missing,Int64}}(missing, nrow(df_1km))
        df_1km.emergeDOY[idx] = Dates.value.(df_1km.emergeDate[idx] .- df_1km.startDate[idx]) .+
                                Dates.dayofyear.(df_1km.startDate[idx])
        df_1km.startMonth = Dates.month.(df_1km.startDate)  # Use month rather than DOY to avoid leap year problems

        # Set missing DOY to be a large number and missing generations to be -9999
        # This is just for the quantile calculations
        df_1km.emergeDOY[ismissing.(df_1km.emergeDate)] .= 9999
        df_1km.nGenerations[ismissing.(df_1km.emergeDate)] .= -9999

        # Group data frame by location and then starting DOY/Date 
        df_group = groupby(df_1km, [:ID, :startMonth])

        # =========================================================
        # ################# Median across the years ##############
        # Calculate the median of nGenerations and emergeDOY
        d_agg = combine(df_group,
            :nGenerations => (x -> quantile(x, 0.5)) => :nGenerations_median,
            :emergeDOY => (x -> quantile(x, 0.5)) => :emergeDOY_median,
            :emergeDate => (x -> sum(ismissing.(x))) => :nMissing)



        # Return nGenerations and emergeDOY to missing
        allowmissing!(df_1km, [:nGenerations, :emergeDOY])
        df_1km.nGenerations[ismissing.(df_1km.emergeDate)] .= missing
        df_1km.emergeDOY[ismissing.(df_1km.emergeDate)] .= missing

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

        # Round number of generations
        out_agg1km.nGenerations_median = round.(out_agg1km.nGenerations_median, digits=2)

        # Remove unwanted variables
        select!(out_agg1km, Not(:nMissing))

        # Write out the average across years data frame
        fileout2 = joinpath(paths.outDir, species_params[s].species_name, "average_" * species_params[s].species_name * "_" *
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
end
