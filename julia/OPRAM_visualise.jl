# Visualise OPRAM results (individual year, 30 year average and anomaly)
#
#
#
# Jon yearsley (jon.yearsley@ucd.ie)
# April 2025
#
# ===============================================================================



using DataFrames
using Statistics
using CSV
using Dates
using Plots
using JLD2


# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
# "frugiperda", "duplicatus", "cembrae", "sexdentatus"

run_params = (speciesName="agrilus_anxius",      # Name of the species
    years=2020,                         # A single year
    startMonth = 1,                     # The month to Visualise
    country="IE",                       # Country code (IE or NI)
    thinFactor=1,                       # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
    gridFile="IE_grid_locations.csv",   # File containing a 1km grid of lats and longs over Ireland 
    averagedPeriod="1991_2020",
    # (used for daylength calculations as well as importing and thining of meteo data)
    save_figs=true)  # If true save figures




# =================================================================================
# Specify directories for import and export of data

if isdir(joinpath(homedir(), "DATA//OPRAM"))       # Linux workstation
    paths = (outDir=joinpath(homedir(), "Desktop//OPRAM//results//"),
        dataDir=joinpath(homedir(), "DATA//OPRAM//"))

elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//R"))   # Mac
    paths = (outDir=joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//results//"),
        dataDir=joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"))

elseif isdir(joinpath(homedir(), "Desktop//OPRAM//"))
    paths = (outDir=joinpath(homedir(), "Desktop//OPRAM//results//"),
        dataDir=joinpath(homedir(), "Desktop//OPRAM//"))
end

include("OPRAM_io_functions.jl");
include("OPRAM_ddmodel_functions.jl");






# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import location data and thin it using thinFactor
# (needed to add the hectad codes to the 10km output) 
grid = prepare_data(joinpath([paths.dataDir, run_params.gridFile]), run_params.thinFactor, run_params.country)

# Create easting and northings of bottom left of a hectad
grid.east_hectad = convert.(Int32, floor.(grid.east ./ 1e4) .* 1e4)
grid.north_hectad = convert.(Int32, floor.(grid.north ./ 1e4) .* 1e4)



# =========================================================
# =========================================================
# Import data

# Find directory matching the species name in run_params
speciesName = filter(x -> occursin(r"" * run_params.speciesName, x), readdir(paths.outDir))
if length(speciesName) > 1
    @error "More than one species name found"
elseif length(speciesName) == 0
    @error "Species not found"
end


# Import 30 year average
aggFile = joinpath(paths.outDir, run_params.speciesName, "average_" * speciesName[1] * "_" *
                                                         run_params.averagedPeriod * "_1km.csv")
d_agg = CSV.read(aggFile, DataFrame, missingstring="NA")


# Import single year result
yearFile_1km = joinpath(paths.outDir, speciesName[1], "combined_" * speciesName[1] * "_" *
                                                  string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * "_" * string(run_params.thinFactor) * "km.csv")

d_1km = CSV.read(yearFile_1km, DataFrame, missingstring="NA")


# Import single year 10km result
yearFile_10km = joinpath(paths.outDir, speciesName[1], "combined_" * speciesName[1] * "_" *
                                                  string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * "_10km.csv")

d_10km = CSV.read(yearFile_10km, DataFrame, missingstring="NA")






# =========================================================
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create some maps!

# Pick a starting month
idx2 = d_agg.startMonth .== run_params.startMonth


if length(run_params.years) == 1
    year_label = string(run_params.years)
else
    year_label = string(minimum(run_params.years)) * "-" * string(maximum(run_params.years))
end

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Map of number of generations =================================================
# Plots a map with a break in colour at 1 generation
colorRange1 = (floor(minimum([0.999,minimum(d_agg.nGenerations_median[idx2])]),digits=1), 
            ceil(maximum(d_agg.nGenerations_median[idx2]),digits=1))
plot(d_agg.east[idx2],
    d_agg.north[idx2],
    zcolor=d_agg.nGenerations_median[idx2],
    seriestype=:scatter,
    markersize=0.5,
    markerstrokewidth=0,
    colormap=cgrad(:bam,[0,(1-colorRange1[1])/(colorRange1[2]-colorRange1[1]),1]),
    clims=colorRange1,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="nGenerations",
    title="Average number of generations per year\n"*run_params.averagedPeriod * " (" * speciesName[1] * ")",
    dpi=400)

if run_params.save_figs
        savefig(speciesName[1] * "_ngen_" * run_params.averagedPeriod * ".png")
end



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Map emergence date =================================================
colorRange2 = (1, 
            maximum([366,maximum(d_agg.emergeDOY_median[idx2])]))

plot(d_agg.east[idx2],
    d_agg.north[idx2],
    zcolor=d_agg.emergeDOY_median[idx2],
    seriestype=:scatter,
    colormap=cgrad(:balance,[0,(366-colorRange2[1])/(colorRange2[2]-colorRange2[1])-0.001,1], rev=true),
    clims=colorRange2,
    markersize=0.5,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="Emergence DOY",
    title="Average emergence DOY\n"* run_params.averagedPeriod * " (" * speciesName[1] * ")",
                dpi=400)

if run_params.save_figs
        savefig(speciesName[1] * "_emergeDOY_" * run_params.averagedPeriod * ".png")
end





# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 10km data =================================================

idx = df_10km.startMonth .== run_params.startMonth .&& Not(ismissing.(df_10km.nGenerations_median_max))
plot(df_10km.east_hectad[idx] .+ 5000.0,
    df_10km.north_hectad[idx] .+ 5000.0,
    seriestype=:scatter,
    zcolor=df_10km.nGenerations_median_max[idx],
    colormap=cgrad(:bam, [0, (1 - colorRange1[1]) / (colorRange1[2] - colorRange1[1]), 1]),
    # color=cgrad(:PiYG, categorical=false),
    clims=colorRange1,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    markersize=3,
    markerstrokewidth=0.5,
    colorbar_title="nGenerations",
    title="Average number of generations per year\n"*run_params.averagedPeriod * " (" * speciesName[1] * ")",
    dpi=400)
   