# Julia program to calculate multi-year average of degree day model results
# and the corresponding anomaly for each year
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
run_params = (speciesName="ips_sexdentatus",      # Name of the species
    years=1991:2020,                    # Either a single year or collect(year1:year2)
    country="IE",                       # Country code (IE or NI)
    startMonth=01,                      # Month for start of development
    thinFactor=1,                       # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
    gridFile="IE_grid_locations.csv",   # File containing a 1km grid of lats and longs over Ireland 
    # (used for daylength calculations as well as importing and thining of meteo data)
    save_figs=true)  # If true save figures







# =================================================================================
# Specify directories for import and export of data

if isdir(joinpath(homedir(),"Desktop//OPRAM"))       # Linux workstation
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



# =========================================================
# =========================================================
# Import data across the years and calculate the average

# Initialise an empty data frame
d = DataFrame()

for y in eachindex(run_params.years)
    # Get the correct filename
    inFile = filter(x -> occursin(r"result_[A-Z]{2}_" * run_params.speciesName * "_" * string(run_params.years[y]) * "_1km.csv", x),
        readdir(joinpath(paths.outDir, run_params.speciesName)))

    @info "Importing data for year $(run_params.years[y])"

    d_year = CSV.read(joinpath(paths.outDir, run_params.speciesName, inFile[1]), DataFrame, missingstring="NA")

    # Make missing dates 3 years after start date
    idx = ismissing.(d_year.emergeDate)
    d_year.emergeDate[idx] .= Date(run_params.years[y], 1, 1) + Day(365 * 4)

    d = vcat(d, d_year)
end

# Calculate average across years (use median to avoid outliers)

# Define DOY for start and emerge dates (and set as integer)
d.emergeDOY = Dates.value.(d.emergeDate .- d.startDate) .+ 1
d.startMonth = Dates.month.(d.startDate)  # Use month rather than DOY to avoid leap year problems


# Group data frame by location and then starting DOY/Date 
df_group = groupby(d, [:ID, :east, :north, :startMonth])

# Calculate the median of nGenerations and emergeDOY
d_agg = combine(df_group,
    [:nGenerations, :emergeDOY] =>
        ((x, y) -> (nGenerations_median=quantile(x, 0.5),
            emergeDOY_median=quantile(y, 0.5))) =>
            AsTable)

   
d.nGenerations = convert(Vector{Union{Float32,Missing}}, d.nGenerations)
d.emergeDOY = convert(Vector{Union{Float32,Missing}}, d.emergeDOY)

idx = d.emergeDOY .> 365 * 3.2
d.nGenerations[idx] .= missing
d.emergeDOY[idx] .= missing

# =========================================================
# =========================================================
# Combine the multi year median with the original data and calculate anomalies


# Combine origninal data frame with 30 year average
leftjoin!(d, d_agg, on=[:ID, :east, :north, :startMonth])
# Calculate anomalies
transform!(d, [:nGenerations, :nGenerations_median] => ((a, b) -> round.(a .- b, digits=3) ) => :nGenerations_anomaly)
transform!(d, [:emergeDOY, :emergeDOY_median] => ((a, b) -> a .- b ) => :emergeDOY_anomaly)



# =========================================================
# =========================================================
# Create 10km summary


  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Import location data and thin it using thinFactor
  grid = prepare_data(joinpath([paths.dataDir, run_params.gridFile]), run_params.thinFactor, run_params.country)

  grid.east_hectad = convert.(Int32, floor.(grid.east ./ 1e4) .* 1e4)
  grid.north_hectad = convert.(Int32, floor.(grid.north ./ 1e4) .* 1e4)


  # Checj eastings and northing do identified unique hectads
tmp = groupby(grid, [:east_hectad, :north_hectad])
combine(tmp, :hectad => unique)

 # Add bottom left eastings and northings of a hectad into result_1km
 eastList = sort(unique(grid.east))
 northList = sort(unique(grid.north))
 d.east_hectad = convert.(Int32, floor.(d.east ./ 1e4) .* 1e4)
 d.north_hectad = convert.(Int32, floor.(d.north ./ 1e4) .* 1e4)


 # Set missing values to be extremes (e.g. generations are small and emergeDOY is large)
 idx = ismissing.(d.emergeDOY)
 d.emergeDOY[idx] .= 365*10
 d.emergeDOY_median[idx] .= 365*10
 d.emergeDOY_anomaly[idx] .= 365*10

 d.nGenerations[idx] .= 0
 d.nGenerations_median[idx] .= 0
 d.nGenerations_anomaly[idx] .= 0

 # Group data frame by location and then starting Date 
 df_group = groupby(d, [:east_hectad, :north_hectad, :startDate])


 
 # Calculate worst case results within each hectad for nGenerations and emergeDOY
 # for each starting date. 
 # These worst case scenarios are not affected by locations where no emergence occurred
 df_10km = combine(df_group,
   :nGenerations  => (x -> quantile(x, 1.0)) => :nGenerations_max,                  # Max generations per hectad
   :nGenerations_median  => (x -> quantile(x, 1.0)) => :nGenerations_median_max,
   :nGenerations_anomaly  => (x -> quantile(x, 1.0)) => :nGenerations_anomaly_max,
   :emergeDOY => (x -> quantile(x, 0.0)) => :emergeDOY_min,                         # Min emergence DOY
   :emergeDOY_median => (x -> quantile(x, 0.0)) => :emergeDOY_median_min,
   :emergeDOY_anomaly => (x -> quantile(x, 0.0)) => :emergeDOY_anomaly_min,
   :nGenerations  => (x -> sum(x.==0)) => :nMissing  )                              # Count number of no emergence events  


 # Add in columns with start and emergence dates
 df_10km.emergeDate_min = Vector{Union{Missing,Date}}(missing, nrow(df_10km))
 df_10km.emergeDate_median_min = Vector{Union{Missing,Date}}(missing, nrow(df_10km))
 idx = df_10km.nMissing .== 0
 df_10km.emergeDate_min[idx] .= Date.(year.(df_10km.startDate[idx])) .+ Day.(floor.(df_10km.emergeDOY_min[idx]))
 df_10km.emergeDate_median_min[idx] .= Date.(year.(df_10km.startDate[idx])) .+ Day.(floor.(df_10km.emergeDOY_min[idx]))


# Add in hectad unique identifier from the grid file
out = rightjoin(grid[:,[:east_hectad, :north_hectad, :hectad]],df_10km, on = [:north_hectad, :east_hectad])


# scatter(df_nGen2.nMissing, df_nGen2.nGenerations_max,markersize=0.5)
# scatter(df_nGen2.nMissing.==0, df_nGen2.emergeDOY_min, markersize=0.5)
# scatter(df_nGen2.nMissing, df_nGen2.startDate, markersize=0.5)


# histogram(df_nGen2.nMissing)

# idx = df_nGen2.nMissing.>00

# histogram(df_nGen2.emergeDOY_min[idx], xlimits=(0,1e3))
# histogram(df_nGen2.emergeDOY_min[Not(idx)], xlimits=(0,1e3))

# histogram(df_nGen2.nGenerations_max[idx], xlimits=(0,2))
# histogram(df_nGen2.nGenerations_max[Not(idx)], xlimits=(0,2))

# =========================================================
 # =========================================================
 # Export the anomaly and multi year average

 # Remove unwanted variables
 select!(d, Not(:startMonth))
 
 # Write out the combined data frame
 fileout1 = joinpath(paths.outDir, run_params.speciesName, "combined_" * run_params.speciesName *
                 "_years_" * string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * ".csv")
 CSV.write(fileout1,  d)
 
 
 # Write out the average data frame
 fileout2 = joinpath(paths.outDir, run_params.speciesName, "average_" * run_params.speciesName *
                 "_years_" * string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * ".csv")
 CSV.write(fileout2,  d_agg)
 

 

# =========================================================
# =========================================================
# Visualise some results

using Plots
idx = d_agg.startMonth .== unique(d_agg.startMonth)[2]

if length(run_params.years) == 1
    year_label = string(run_params.years)
else
    year_label = string(minimum(run_params.years)) * "-" * string(maximum(run_params.years))
end

plot(d_agg.east[idx],
    d_agg.north[idx],
    zcolor=d_agg.nGenerations_median[idx],
    seriestype=:scatter,
    markersize=0.5,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="nGenerations",
    title="Average number of generations per year\n"*year_label,
    dpi=400)

if run_params.save_figs
        savefig(run_params.speciesName * "_ngen_" * year_label * ".png")
end


plot(d_agg.east[idx],
    d_agg.north[idx],
    zcolor=d_agg.emergeDOY_median[idx],
    seriestype=:scatter,
    markersize=0.5,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="Emergence DOY",
    title="Average emergence DOY\n"*year_label,
                dpi=400)

if run_params.save_figs
        savefig(run_params.speciesName * "_emergeDOY_" * year_label * ".png")
end