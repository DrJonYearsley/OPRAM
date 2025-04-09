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




# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
# "frugiperda", "duplicatus", "cembrae", "sexdentatus"

run_params = (speciesName="agrilus_anxius",      # Name of the species
    years=1991:2020,                    # the years collect(year1:year2)
    maxYears = 3,                       # Maximum number of years to complete insect development (must correspond to simulation value)
    country="IE",                       # Country code (IE or NI)
    thinFactor=1,                       # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
    gridFile="IE_grid_locations.csv",   # File containing a 1km grid of lats and longs over Ireland 
    # (used for daylength calculations as well as importing and thining of meteo data)
    save_figs=true)  # If true save figures







# =================================================================================
# Specify directories for import and export of data

if isdir(joinpath(homedir(),"DATA//OPRAM"))       # Linux workstation
    paths = (outDir=joinpath(homedir(),"Desktop//OPRAM//results//"),
        dataDir=joinpath(homedir(),"DATA//OPRAM//"))

elseif isdir(joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//R"))   # Mac
    paths = (outDir=joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//results//"),
        dataDir=joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"))

elseif isdir(joinpath(homedir(),"Desktop//OPRAM//"))
    paths = (outDir=joinpath(homedir(),"Desktop//OPRAM//results//"),
        dataDir=joinpath(homedir(),"Desktop//OPRAM//"))
end

include("OPRAM_io_functions.jl");
include("OPRAM_ddmodel_functions.jl");






  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Import location data and thin it using thinFactor
  grid = prepare_data(joinpath([paths.dataDir, run_params.gridFile]), run_params.thinFactor, run_params.country)

  # Create easting and northings of bottom left of a hectad
  grid.east_hectad = convert.(Int32, floor.(grid.east ./ 1e4) .* 1e4)
  grid.north_hectad = convert.(Int32, floor.(grid.north ./ 1e4) .* 1e4)







# =========================================================
# =========================================================
# Import data across the years and calculate the average

# Find directory matching the species name in run_params
speciesName = filter(x -> occursin(r""*run_params.speciesName, x), readdir(paths.outDir))
if length(speciesName)>1
    @error "More than one species name found"
elseif length(speciesName)==0
    @error "Species not found"
end

dVec = Vector{DataFrame}(undef, length(run_params.years))
for y in eachindex(run_params.years)

    @info "Importing data for year $(run_params.years[y])"

    # Starting dates for output in CSV files
    # The first of every month
    dates = [Date(run_params.years[y], m, 01) for m in 1:12]

    inFile = filter(x -> occursin(r"^" * speciesName[1] * "_" * run_params.country * "_" * string(run_params.years[y]) * "_1km.jld2", x),
        readdir(joinpath(paths.outDir, speciesName[1])))

        if length(inFile)>1
            @error "More than one input file found"
        end
    adult_emerge = load_object(joinpath(paths.outDir, speciesName[1],inFile[1]))

    @info " ---- Generating output for specific starting dates"
    # Create output for specific days of year
    dVec[y] = create_doy_results(dates, adult_emerge)

    # Select columns to work with
    select!(dVec[y], [:ID, :startDate, :emergeDate, :nGenerations])
end

# Put all the model results data together into one Matrix
df_1km = reduce(vcat, dVec)
dVec = nothing


# =========================================================
# Calculate average across years (use median to avoid outliers)

@info "Calculating median values across years"

# Define DOY for start and emerge dates (and set as integer)
idx = .!ismissing.(df_1km.emergeDate)
df_1km.emergeDOY = Vector{Union{Missing,Int64}}(missing, nrow(df_1km))
df_1km.emergeDOY[idx] = Dates.value.(df_1km.emergeDate[idx] .- df_1km.startDate[idx]) .+ 1
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
out_agg1km = rightjoin(grid[:,[:ID, :east, :north]], d_agg,  on = [:ID])

# Round number of generations
out_agg1km.nGenerations_median = round.(out_agg1km.nGenerations_median, digits=2)   

 # Remove unwanted variables
 select!(out_agg1km, Not(:nMissing))
  
 # Write out the average across years data frame
 fileout2 = joinpath(paths.outDir, speciesName[1], "average_" * speciesName[1] * "_" *
                 string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * "_" *string(run_params.thinFactor) * "km.csv")
 CSV.write(fileout2,  out_agg1km, missingstring="NA")
 
