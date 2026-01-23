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
using Distributed;
using TOML;



include("OPRAM_io_functions.jl");
include("OPRAM_processing_functions.jl");


# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
# "frugiperda", "duplicatus", "cembrae", "sexdentatus"

rcp = 26
year = 2055
month = 1
scale = 10     # Either 1 for 1km scale or 10 for hectad scale
county = 5      # Set to 0 to import all counties
paramFile = "parameters_future_halymorpha.toml"

nNodes, run_params, species_setup, paths = import_parameters(paramFile, true)






# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import location data and thin it using thinFactor
# (needed to add the hectad codes to the 10km output) 
grid = read_grid(run_params)

# Create easting and northings of bottom left of a hectad
grid.east_hectad = convert.(Int32, floor.(grid.east ./ 1e4) .* 1e4)
grid.north_hectad = convert.(Int32, floor.(grid.north ./ 1e4) .* 1e4)


hectad_list = unique(grid.hectad)
idx = [findfirst(grid.hectad.==hectad_list[l]) for l in eachindex(hectad_list)]
grid_hectad = grid[idx,:]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set species parameters from the parameter file (can be more than one species)
if species_setup.speciesStr[1] == "all"
    @info "Importing all species from the species file" * string(species_setup.speciesStr)
    species_params = import_species(species_setup.speciesFile, species_setup.speciesStr[1])
else
    @info "Importing species from the species file: " * string(species_setup.speciesStr)
    species_params = [import_species(species_setup.speciesFile, species_setup.speciesStr[s]) for s in eachindex(species_setup.speciesStr)]
end


# =========================================================
# =========================================================
# Import data

# Find directory matching the species name in run_params
regex = Regex(replace(lowercase(species_params[1].species_name),
    r"\s" => "\\w"))  # Replace spaces with reg expression
speciesName = filter(x -> occursin(regex, x), readdir(paths.outDir))
if length(speciesName) > 1
    @error "More than one species name found"
elseif length(speciesName) == 0
    @error "Species not found"
end


# Default values
ms = 0.5
inFile = [speciesName[1] * "_" * string(scale) * "_" * string(county) * "_" * string(year) * ".csv"]

if scale == 10
    inFile = [speciesName[1] * "_100_" * string(year) * ".csv"]
    ms = 3
else
    if county == 0
        # Import mulitple files
        reg = Regex(speciesName[1] * "_" * string(scale) * "_" * "[[:digit:]]{1,2}" * "_" * string(year) * ".csv")
        inFile = filter(x -> occursin(reg, x), readdir(joinpath(paths.outDir, speciesName[1])))
    end
end

d = CSV.read(joinpath(paths.outDir, speciesName[1], inFile[1]), DataFrame, missingstring="NA")
for f in 2:length(inFile)
      append!(d, CSV.read(joinpath(paths.outDir, speciesName[1], inFile[f]), DataFrame, missingstring="NA"))
end


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define starting month for start dates (and set as integer)
idx = .!ismissing.(d.startDate)
d.startMonth = Vector{Union{Missing,Int}}(missing, nrow(d))
d.startMonth[idx] = Dates.month.(d.startDate[idx])  # Use month rather than DOY to avoid leap year problems


# Add in coordinates
if scale==10
   leftjoin!(d, grid_hectad[:,[:hectad, :east_hectad, :north_hectad]], on=[:hectad])
   rename!(d, :east_hectad => "eastings", :north_hectad => "northings")
end


# =========================================================
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create some maps!

# Pick a starting month
idx2 = d.startMonth .== month


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Map of number of generations =================================================
# Plots a map with a break in colour at 1 generation
x=d.eastings[idx2]
y=d.northings[idx2]
z = d.nGen[idx2]

plot_idx = .!ismissing.(z)

plot(x[plot_idx],
    y[plot_idx],
    zcolor=z[plot_idx],
    seriestype=:scatter,
    markersize=ms,
    markerstrokewidth=0,
    # colormap=cgrad(:bam,[0,(1-colorRange1[1])/(colorRange1[2]-colorRange1[1]),1]),
    # clims=colorRange1,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    # colorbar_title="nGenerations",
    # title="Average number of generations per year\n"*run_params.averagedPeriod * " (" * speciesName[1] * ")",
    dpi=400)

    # Add in points with missing data
plot!(x[Not(plot_idx)],
    y[Not(plot_idx)],
    seriestype=:scatter,
    markersize=ms,
    markerstrokewidth=0,
    markercolor="green")


