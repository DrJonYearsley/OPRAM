# Visualise meteo data used as input into OPRAM
#
#
#
# Jon yearsley (jon.yearsley@ucd.ie)
# April 2025
#
# ===============================================================================



using CSV
using DataFrames
using Dates
using JLD2
using Plots
using Statistics
using Distributed;
using TOML;



include("OPRAM_io_functions.jl");
include("OPRAM_processing_functions.jl");


# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
# "frugiperda", "duplicatus", "cembrae", "sexdentatus"

# run_params = (speciesName="halyomorpha_halys",      # Name of the species
#     years=2024,                         # A single year
#     startMonth = 1,                     # The month to Visualise
#     country="IE",                       # Country code (IE or NI)
#     thinFactor=1,                       # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
#     gridFile="IE_grid_locations.csv",   # File containing a 1km grid of lats and longs over Ireland 
#     averagedPeriod="1991_2020",
#     # (used for daylength calculations as well as importing and thining of meteo data)
#     save_figs=true)  # If true save figures

year = 2024
month = 1
scale = 10     # Either 1 for 1km scale or 10 for hectad scale
county = 5      # Set to 0 to import all counties

# Model parameters are stored in a TOML file https://toml.io/en/
if length(ARGS) == 1
    nNodes, run_params, species_setup, paths = import_parameters(ARGS[1], true)

elseif length(ARGS) == 0 & isfile("parameters.toml")
    nNodes, run_params, species_setup, paths = import_parameters("parameters_future.toml", true)

else
    @error "No parameter file given"
end






# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import location data and thin it using thinFactor
# (needed to add the hectad codes to the 10km output) 
grid = read_grid(run_params)

# # Create easting and northings of bottom left of a hectad
# grid.east_hectad = convert.(Int32, floor.(grid.east ./ 1e4) .* 1e4)
# grid.north_hectad = convert.(Int32, floor.(grid.north ./ 1e4) .* 1e4)


# hectad_list = unique(grid.hectad)
# idx = [findfirst(grid.hectad.==hectad_list[l]) for l in eachindex(hectad_list)]
# grid_hectad = grid[idx,:]


# =========================================================
# =========================================================
# Import data

Tmax, Tmaxsd, Tmin, Tminsd, DOY, ID = read_JLD2_translate2(paths.meteoDir_IE, "26", "2021-2050", grid.ID)



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define coordinates
coord_df = DataFrame(ID = ID)

leftjoin!(coord_df, grid[:,[:ID, :east, :north]], on=[:ID])


# =========================================================
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create some maps!

x=coord_df.east
y=coord_df.north
z = Tmin[1,:]
ms=0.5

# x = grid2.east
# y = grid2.north
# z = Tmax_mean[1,:]

plot(x,
    y,
    zcolor=z,
    seriestype=:scatter,
    markersize=ms,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal)
