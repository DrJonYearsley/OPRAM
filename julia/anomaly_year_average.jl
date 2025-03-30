# Julia program to:
#  extract results from model jld2 files 
#  calculate multi-year average of degree day model results
#  calculate corresponding anomaly for each year
#  aggregate results across hectads
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
run_params = (speciesName="sex",      # Name of the species
    years=1991:2020,                    # Either a single year or collect(year1:year2)
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

# Put all the temperature data together into one Matrix
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

# Set missing DOY to be a large number and missing generations to be -1
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
# Combine the multi year median with the original data and calculate anomalies

@info "Calculating anomalies"

# Combine origninal data frame with 30 year average
leftjoin!(df_1km, d_agg, on=[:ID, :startMonth])
# Calculate anomalies
transform!(df_1km, [:nGenerations, :nGenerations_median] => ((a, b) -> a .- b ) => :nGenerations_anomaly)
transform!(df_1km, [:emergeDOY, :emergeDOY_median] => ((a, b) -> a .- b ) => :emergeDOY_anomaly)

# Round results for number of generations
df_1km.nGenerations_median = round.(df_1km.nGenerations_median, digits=2)
d_agg.nGenerations_median = round.(d_agg.nGenerations_median, digits=2)



# =========================================================
 # =========================================================
 # Export the anomaly and multi year average

@info "Writing 1km results to CSV files"


# Add eastings and northings into the 1km data frames (copy the data frames so we can round before saving)
out_1km = rightjoin(grid[:,[:ID, :east, :north]], df_1km, on = [:ID])
out_agg1km = rightjoin(grid[:,[:ID, :east, :north]], d_agg,  on = [:ID])

# Round number of generations
out_1km.nGenerations = round.(out_1km.nGenerations, digits=2)
out_1km.nGenerations_median = round.(out_1km.nGenerations_median, digits=2)
out_1km.nGenerations_anomaly = round.(out_1km.nGenerations_anomaly, digits=2)
out_agg1km.nGenerations_median = round.(out_agg1km.nGenerations_median, digits=2)   

 # Remove unwanted variables
 select!(out_1km, Not([:startMonth,:nMissing, :emergeDOY, :emergeDOY_median]))
 select!(out_agg1km, Not(:nMissing))
 
 # Write out the combined data frame
 fileout1 = joinpath(paths.outDir, speciesName[1], "combined_" * speciesName[1] * "_" *
                 string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * "_" *string(run_params.thinFactor) * "km.csv")
 CSV.write(fileout1,  out_1km, missingstring="NA")
 
 
 # Write out the average across years data frame
 fileout2 = joinpath(paths.outDir, speciesName[1], "average_" * speciesName[1] * "_" *
                 string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * "_" *string(run_params.thinFactor) * "km.csv")
 CSV.write(fileout2,  out_agg1km, missingstring="NA")
 
out_1km = nothing
out_agg1km = nothing



# =========================================================
# =========================================================
# Create 10km summary


@info "Aggregating data to 10km resolution"

 # Add bottom left eastings and northings of a hectad into df_1km
leftjoin!(df_1km, grid[:,[:ID, :east, :north]],  on = [:ID])
 df_1km.east_hectad = convert.(Int32, floor.(df_1km.east ./ 1e4) .* 1e4)
 df_1km.north_hectad = convert.(Int32, floor.(df_1km.north ./ 1e4) .* 1e4)


 # Set missing values to be extremes (e.g. generations are small and emergeDOY is large)
 idx1 = ismissing.(df_1km.emergeDate)
 df_1km.emergeDOY[idx1] .= 9999
 df_1km.nGenerations[idx1] .= -9999
 

 idx2 = ismissing.(df_1km.emergeDOY_median)
 df_1km.emergeDOY_median[idx2] .= 9999
 df_1km.nGenerations_median[idx2] .= -9999

 idx3 = ismissing.(df_1km.emergeDOY_anomaly)
 df_1km.emergeDOY_anomaly[idx3] .= 9999
 df_1km.nGenerations_anomaly[idx3] .= -9999


 # Group data frame by location and then starting Date 
 df_group = groupby(df_1km, [:east_hectad, :north_hectad, :startDate])


 
 # Calculate worst case results within each hectad for nGenerations and emergeDOY
 # for each starting date. 
 # These worst case scenarios are not affected by locations where no emergence occurred
 df_10km = combine(df_group,
    :nGenerations  => (x -> quantile(x, 1.0)) => :nGenerations_max,                  # Max generations per hectad
    :nGenerations_median  => (x -> quantile(x, 1.0)) => :nGenerations_median_max,
    :nGenerations_anomaly  => (x -> quantile(x, 1.0)) => :nGenerations_anomaly_max,
    :emergeDOY => (x -> quantile(x, 0.0)) => :emergeDOY_min,                         # Min emergence DOY
    :emergeDOY_median => (x -> quantile(x, 0.0)) => :emergeDOY_median_min,
    :emergeDOY_anomaly => (x -> quantile(x, 1.0)) => :emergeDOY_anomaly_min,
    :nGenerations  => (x -> sum(x.<0)) => :nMissing  )                              # Count number of no emergence events  

# Add in missing values where 10km worst case is out of bounds
# (i.e. nGenerations<0 or emergeDOY>5000)
allowmissing!(df_10km, [:nGenerations_max, :nGenerations_median_max, :nGenerations_anomaly_max, :emergeDOY_min, :emergeDOY_median_min, :emergeDOY_anomaly_min])
df_10km.nGenerations_max[df_10km.nGenerations_max .< 0] .= missing
df_10km.nGenerations_anomaly_max[abs.(df_10km.nGenerations_anomaly_max).>5000] .= missing
df_10km.nGenerations_median_max[df_10km.nGenerations_median_max .< 0] .= missing
df_10km.emergeDOY_min[df_10km.emergeDOY_min .> 5000] .= missing
df_10km.emergeDOY_anomaly_min[abs.(df_10km.emergeDOY_anomaly_min).>5000] .= missing
df_10km.emergeDOY_median_min[df_10km.emergeDOY_median_min .> 5000] .= missing

 # Add in columns with start and emergence dates where we have no missing data
 df_10km.emergeDate_min = Vector{Union{Missing,Date}}(missing, nrow(df_10km))
 df_10km.emergeDate_median_min = Vector{Union{Missing,Date}}(missing, nrow(df_10km))

 idx1 = .!ismissing.(df_10km.emergeDOY_min)
 df_10km.emergeDate_min[idx1] .= Date.(year.(df_10km.startDate[idx1])) .+ Day.(floor.(df_10km.emergeDOY_min[idx1]))

 idx2 = .!ismissing.(df_10km.emergeDOY_median_min)
 df_10km.emergeDate_median_min[idx2] .= Date.(year.(df_10km.startDate[idx2])) .+ Day.(floor.(df_10km.emergeDOY_median_min[idx2]))





# Add in hectad unique identifier from the grid file
out_10km = rightjoin(unique(grid[:,[:east_hectad, :north_hectad, :hectad]]), df_10km, on = [:north_hectad, :east_hectad])



# =========================================================
 # =========================================================
 # Export 10km summary for the anomaly and multi year average

 @info "Writing 10km results to CSV files"

 # Remove unwanted variables
 select!(out_10km, Not(:nMissing))
 
# Round results for number of generations
out_10km.nGenerations_max = round.(out_10km.nGenerations_max, digits=2)
out_10km.nGenerations_median_max = round.(out_10km.nGenerations_median_max, digits=2)
out_10km.nGenerations_anomaly_max = round.(out_10km.nGenerations_anomaly_max, digits=2)


 # Write out the combined data frame
 fileout3 = joinpath(paths.outDir, speciesName[1], "combined_" * speciesName[1] * "_" *
                 string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * "_10km.csv")
 CSV.write(fileout3,  out_10km, missingstring="NA")
 

 

# =========================================================
# =========================================================
# Visualise some results

@info "Producing a couple of visualisations"


using CSV
using Plots
using DataFrames

d_agg = CSV.read(joinpath(paths.outDir,"ips_sexdentatus","average_ips_sexdentatus_1991_2020_1km.csv"), 
        DataFrame, missingstring="NA")

# Pick a starting month
idx2 = d_agg.startMonth .== unique(d_agg.startMonth)[1]

# leftjoin!(d_agg, grid[:,[:ID, :east, :north]],  on = [:ID])

if length(run_params.years) == 1
    year_label = string(run_params.years)
else
    year_label = string(minimum(run_params.years)) * "-" * string(maximum(run_params.years))
end

# Plots a map with a break in colour at 1 generation
colorRange = (floor(minimum([0.999,minimum(d_agg.nGenerations_median[idx2])]),digits=1), 
            ceil(maximum(d_agg.nGenerations_median[idx2]),digits=1))
plot(d_agg.east[idx2],
    d_agg.north[idx2],
    zcolor=d_agg.nGenerations_median[idx2],
    seriestype=:scatter,
    markersize=0.5,
    markerstrokewidth=0,
    colormap=cgrad(:bam,[0,(1-colorRange[1])/(colorRange[2]-colorRange[1]),1]),
    clims=colorRange,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="nGenerations",
    title="Average number of generations per year\n"*year_label * " (" * speciesName[1] * ")",
    dpi=400)

if run_params.save_figs
        savefig(speciesName[1] * "_ngen_" * year_label * ".png")
end

colorRange = (1, 
            maximum([366,maximum(d_agg.emergeDOY_median[idx2])]))

plot(d_agg.east[idx2],
    d_agg.north[idx2],
    zcolor=d_agg.emergeDOY_median[idx2],
    seriestype=:scatter,
    colormap=cgrad(:managua,[0,(366-colorRange[1])/(colorRange[2]-colorRange[1]),1], rev=true),
    clims=colorRange,
    markersize=0.5,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="Emergence DOY",
    title="Average emergence DOY\n"*year_label* " (" * speciesName[1] * ")",
                dpi=400)

if run_params.save_figs
        savefig(speciesName[1] * "_emergeDOY_" * year_label * ".png")
end