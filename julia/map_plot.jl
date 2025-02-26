# Plot a map of Ireland using 10km circles and coloured by county or province
# This is a pilot visualisation for the OPRAM web app

using DataFrames
using Statistics
using CSV
using Dates
using Plots
using JLD2
using PlotlyBase

# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
save_figs = true  # If true save figures
hectad = true     # If true aggregate data into 10km squares


# =================================================================================
# Specify directories for import and export of data

if isdir("//home//jon//Desktop//OPRAM")
    outDir = "//home//jon//Desktop//OPRAM//results//"
    dataDir = "//home//jon//DATA//OPRAM//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
    outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
    dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"

elseif isdir("//users//jon//Desktop//OPRAM//")
    outDir = "//users//jon//Desktop//OPRAM//results//"
    dataDir = "//users//jon//Desktop//OPRAM//"
 end





# =================================================================================
# Import coastline data
coast = CSV.read(joinpath([dataDir, "coastline.csv"]), DataFrame,
    types=[Float64, Float64, Int, Int, Int])


# =================================================================================
# Import grid location data
grid = CSV.read(joinpath([dataDir, "IE_grid_locations.csv"]), DataFrame,
    types=[Int, Int, Int, String,  Float64, Float64,String, String,String])

# Calculate east and north of hectad bottom left corner
grid.east_hectad = 10000.0 .* floor.(grid.east/10000)
grid.north_hectad = 10000.0 .* floor.(grid.north/10000)


# =================================================================================
# Plot map

grid_hectads = combine(groupby(grid, :hectad), [:east_hectad, :north_hectad, :county, :province] .=> x -> x[1])


countyList = unique(grid_hectads.county_function)
provinceList = unique(grid_hectads.province_function)

x = sort(unique(coast.idx_east))
y = sort(unique(coast.idx_north))



grid_hectads.county = [findfirst(x .== countyList) for x in grid_hectads.county_function]
grid_hectads.province = [findfirst(x .== provinceList) for x in grid_hectads.province_function]

# Use plotly for figures
# plotly() 
gr()
plot(grid_hectads.east_hectad_function.+5e3,
    grid_hectads.north_hectad_function.+5e3,
        seriestype=:scatter,
        zcolor=grid_hectads.province,
        color=cgrad(:Dark2_4, categorical=true),
        showaxis=false,
        legend=false,
        grid=false,
        cbar=false,
        label=false,
        markersize=4,
        markerstrokewidth=0.5,
        aspect_ratio=:equal,
        dpi=600)

plot!(coast.east,
        coast.north,
        seriestype=:scatter,
        showaxis=false,
        grid=false,
        markersize=1,
        markercolor="black",
        dpi=600)

if save_figs
    savefig("ireland_province_map.png")
end



plot(grid_hectads.east_hectad_function.+5e3,
    grid_hectads.north_hectad_function.+5e3,
        seriestype=:scatter,
        zcolor=grid_hectads.county,
        color=cgrad(:seaborn_bright, categorical=true),
        showaxis=false,
        legend=false,
        grid=false,
        cbar=false,
        label=false,
        markersize=4,
        markerstrokewidth=0.5,
        aspect_ratio=:equal,
        dpi=600)

plot!(coast.east,
        coast.north,
        seriestype=:scatter,
        showaxis=false,
        grid=false,
        markersize=1,
        markercolor="black",
        dpi=600)

if save_figs
    savefig("ireland_county_map.png")
end