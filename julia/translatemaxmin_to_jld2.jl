# Script to import TRANSLATE project meteo files and convert them 
# to jld2 format for more efficient import
# 
# The TRANSLATE files are imported, the median, 10th and 90th quantiles 
# are combined to estimate a mean and standard deviation (that will be used 
# later to produce random temp data). Then the mean and standard deviation
# are interpolated onto a 1km grid of locations in Ireland.
#
# Interpolation uses nearest-neighbour model with 10 neighbours. 
# This avoids extrapolation errors, but may need to be revised
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 23th Feb 2025
#
# ====================================================================



using CSV;
using JLD2;
using DataFrames;
using GeoStats;
using NetCDF;
using SpecialFunctions;
using Distributed;

include("OPRAM_io_functions.jl")  # Includes functions to import meteo data


rcpList = ["26", "45", "85"];              # The RCP scenario (26 = RCP2.6, 45=RCP4.5, 85=RCP8.5)
periodList = ["2021-2050", "2041-2070"];   # The future scenario (2021-2050 or 2041-2070)
gridFile = "IE_grid_locations.csv";        # File containing a 1km grid of lats and longs over Ireland 
quantiles = ["10", "50", "90"];            # The quantiles of the ensemble to import

# rcpList = ["85"];              # The RCP scenario (26 = RCP2.6, 45=RCP4.5, 85=RCP8.5)
# periodList = ["2041-2070"];   # The future scenario (2021-2050 or 2041-2070)


# =========================================================
# =========================================================

# Specify directories for import and export of data
if isdir("//home//jon//DATA//OPRAM")
    outDir = "//home//jon//DATA//OPRAM//Climate_JLD2//"
    dataDir = "//home//jon//DATA//OPRAM//"
    translateDir_max = "//home//jon//DATA//OPRAM//TRANSLATE//tmax"
    translateDir_min = "//home//jon//DATA//OPRAM//TRANSLATE//tmin"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
    outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2//"
    dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
    translateDir_max = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//TRANSLATE//tmax"
    translateDir_min = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//TRANSLATE//tmin"

elseif isdir("//users//jon//Desktop//OPRAM//")
    outDir = "//users//jon//Desktop//OPRAM//Climate_JLD2//"
    dataDir = "//users//jon//Desktop//OPRAM//"
    translateDir_max = "//users//jon//Desktop//OPRAM//TRANSLATE//tmax"
    translateDir_min = "//users//jon//Desktop//OPRAM//TRANSLATE//tmin"
end


# Make directory for the output if one doesn't exist
if !isdir(outDir)
    @info "Making directory for output"
    mkpath(outDir)
end

# =========================================================
# =========================================================
# Import grid of location data 

# Import grid location data
grid = CSV.read(joinpath([dataDir, gridFile]), DataFrame);

# Sort locations in order of IDs
grid = grid[sortperm(grid.ID), :];

# Create a geotable of the grid data
grid_gt = georef((x=grid.east, y=grid.north, ID=grid.ID), ("x", "y"), crs=EPSG{29903})


# =========================================================
# =========================================================
# Check size of data to be imported (using Tmax) and 

# Find correct directory and file
files = filter(x -> occursin(Regex("^tmax_rcp" * string(rcpList[1]) * "_" *string(periodList[1]) * "_" * "\\w+"  *  "_ens50.nc") , x), readdir(translateDir_max))

# Import one file of the TRANSLATE data 
Tmax1 = ncread(joinpath(translateDir_max, files[1]), "tmax");
lat = ncread(joinpath(translateDir_max, files[1]), "lat");
lon = ncread(joinpath(translateDir_max, files[1]), "lon");

# Calculate some dimensions
TSize = size(Tmax1)
nDays = TSize[3];

# Reshape arrays in place
Tmax_long = reshape(Tmax1,:,nDays);
lon_2D = reshape([lon[i] for i in eachindex(lon), j in eachindex(lat)],:);
lat_2D = reshape([lat[j] for i in eachindex(lon), j in eachindex(lat)], :);

# Clear some variables
Tmax1 = nothing;
lat = nothing;
lon = nothing;


# Find data that has unrealsitic max temps
nonzero_idx = dropdims(all(Tmax_long.>-273, dims=2), dims=2);

# Keep only lats and longs with non-zero temps
lon_2D = lon_2D[nonzero_idx];
lat_2D = lat_2D[nonzero_idx];



# =========================================================
# =========================================================
# For each 1km grid square find nearest neighbour TRANSLATE location

# Convert TRANSLATE lat lon into eastings and northings
translate_xy = georef((lat=lat_2D, lon=lon_2D, ID=collect(1:length(lon_2D))),
    ("lat", "lon"),
    crs=EPSG{4326}) |> Proj(EPSG{29903})  # Convert to Irish National Grid TM75

# Define interpolation model
interpolation_model = InterpolateNeighbors(domain(grid_gt),
                   model=NN(Euclidean()),
                   maxneighbors=10) 

# Find nearest neighbour IDs for the TRANSLATE data
translate_NN =  translate_xy |> interpolation_model
grid.translateID = translate_NN.ID


# Array to hold temperature data for the quantiles
Tmax_quantiles = Array{Float32,3}(undef, sum(nonzero_idx), nDays, length(quantiles));
Tmin_quantiles = Array{Float32,3}(undef, sum(nonzero_idx), nDays, length(quantiles));

# Array to hold interpolated temperature data (mean and tandard deviation)
Tmax_interp = Array{Float32,2}(undef, nDays, length(grid.ID));
Tmaxsd_interp = Array{Float32,2}(undef, nDays, length(grid.ID));
Tmin_interp = Array{Float32,2}(undef, nDays, length(grid.ID));
Tminsd_interp = Array{Float32,2}(undef, nDays, length(grid.ID));

# Create array for days of year
DOY = collect(1:nDays);


# Calculate distance between 10th and 90th quantiles in units of standard deviation
# (used to estimate standard deviation of temp))
sd_estimator = (erfinv(2*0.9 - 1) - erfinv(2*0.1 - 1)) * sqrt(2)


# =========================================================
# =========================================================
# Import the data for each RCP and period, interpolate data onto 1km grid 
# and save in a JLD2 file

for r in eachindex(rcpList)
    @time "Converted TRANSLATE data:" for p in eachindex(periodList)
        # +++++++++++++++++++++++++++++++++++++++++++++++++++
        # Import the TRANSLATE data
        @info "Importing data for RCP " * string(rcpList[r]) * " and period " * string(periodList[p])


        # Find the three files for a specific RCP and time period
        # rcpDir = filter(x -> occursin(string(rcpList[r]), x), readdir(joinpath(translateDir)))
        # files_max = filter(x -> occursin(r"_" * string(periodList[p]) * "_ens", x), readdir(joinpath(translateDir, rcpDir[1])))
        files_max = filter(x -> occursin(Regex("^tmax_rcp" * string(rcpList[r]) * "_" * string(periodList[p]) ), x), readdir(translateDir_max))
        files_min = filter(x -> occursin(Regex("^tmin_rcp" * string(rcpList[r]) * "_" * string(periodList[p]) ), x), readdir(translateDir_min))

        for q in eachindex(quantiles)
            @info "Processing quantile " * string(quantiles[q])
            # Read the Tmax and Tmin data and reshape into a 2D matrix
            fileIn = filter(x -> occursin("ens" * quantiles[q], x), files_max)
            Tmax = reshape(ncread(joinpath(translateDir_max, fileIn[1]), "tmax"), : , nDays)

            fileIn = filter(x -> occursin("ens" * quantiles[q], x), files_min)
            Tmin = reshape(ncread(joinpath(translateDir_min, fileIn[1]), "tmin"), : , nDays)

            # Save non-zero temp data
            Tmax_quantiles[:, :, q] .= Tmax[nonzero_idx, :];
            Tmin_quantiles[:, :, q] .= Tmin[nonzero_idx, :];
        end

        # Calculate standard deviation from 10th and 90th quantiles
        Tmaxsd = (Tmax_quantiles[:,:,3] .- Tmax_quantiles[:,:,1]) ./ sd_estimator;
        Tminsd = (Tmin_quantiles[:,:,3] .- Tmin_quantiles[:,:,1]) ./ sd_estimator;

        # Interpolate the median and standard deviations onto the 1km grid
        for doy in eachindex(DOY)
            @info "Interpolating for Day" * string(doy)
            # Interpolate the temperature data onto the 1km grid

            # Max temp interpolation (based on nearest neighbour)
            Tmax_interp[doy, :] = Tmax_quantiles[grid.translateID, doy, 2]
            Tmaxsd_interp[doy, :] = Tmaxsd[grid.translateID, doy]

            # Min temp interpolation
            Tmin_interp[doy, :] = Tmin_quantiles[grid.translateID, doy, 2]
            Tminsd_interp[doy, :] = Tminsd[grid.translateID, doy]
        end

        # An alternative interpolation approach could be to interpolate the data directly
        # This would be better if we were not doing nearest neightbour
        # e.g. 
        # data_xy = georef((lat=lat_2D, lon=lon_2D,
        #         X=Tmax_quantiles[:, doy, 2]), ("lat", "lon"), crs=EPSG{4326}) |> Proj(EPSG{29903})
        # x_interp = data_xy |> interpolation_model

        # Save this to a JLD2 file
        outfile = joinpath([outDir, "TRANSLATE_Tmaxmin_rcp" * string(rcpList[r]) * "_" * string(periodList[p]) * ".jld2"])
        @info "Saving file " * outfile
        jldsave(outfile; Tmax_mean=Tmax_interp, Tmax_sd=Tmaxsd_interp, Tmin_mean=Tmin_interp, Tmin_sd=Tminsd_interp, DOY=DOY, ID=grid.ID)
    end
end

# tmp = Tmax_interp .- 8.0
# tmp[tmp.<0.0] .= 0.0
# z = [sum(Tmaxsd_interp[:,i]) for i in eachindex(grid.ID)]

# Plots.plot(grid.east,
#      grid.north,
#       zcolor=z,
#      seriestype=:scatter,
#      markerstrokewidth=0,
#      markersize=0.5,
#         showaxis=true,
#     grid=false,
#     legend=false,
#     cbar=true,
#     aspect_ratio=:equal)

# z

# =========================================================
# =========================================================
# scatter(grid.longitude, grid.latitude, marker_z=Tavg_interp[1,:,1], markerstrokewidth=0, markersize=1)


# Plots.scatter(lon_2D, lat_2D, marker_z=Tsd[:,1], markerstrokewidth=0, markersize=1)
# Plots.scatter(lon_2D, lat_2D, marker_z=Tavg_quantiles[:,1,2], markerstrokewidth=0, markersize=1)

# ind = ismissing.(Tavg_interp2)
# scatter(grid.longitude[.!ind], grid.latitude[.!ind],
#     marker_z=Tavg_interp2[.!ind],
#     markerstrokewidth=0,
#     markersize=1)



# # Look at symmetry of the 10th and 90th quantiles
# tmp10 = Tavg_interp[:,:,2] -  Tavg_interp[:,:,1]
# tmp90 = Tavg_interp[:,:,3] -  Tavg_interp[:,:,2]
# histogram(tmp90[:,1] .- tmp10[:,1])

# x=[mean(tmp90[:,i] .- tmp10[:,i]) for i in 1:size(tmp10,2)] # Average across year
# y=[mean(tmp90[i,:] .- tmp10[i,:]) for i in 1:size(tmp10,1)] # Average across locations

# rat = (tmp90 .- tmp10) ./ Tavg_interp[:,:,2]


# histogram(x)
# histogram(y)

# DOY_2D = [DOY[i] for i in eachindex(DOY), j in eachindex(grid.ID)]

# Plots.scatter(DOY, y)
# Makie.boxplot(ceil.(DOY_2D[:].*(12/365)), tmp90[:] .- tmp10[:])
# Makie.boxplot(ceil.(DOY_2D[:].*(12/365)), rat[:])

# Plots.scatter(grid.longitude, grid.latitude, marker_z=x, markerstrokewidth=0, markersize=1)


# std(tmp10 .- tmp90)




# # Visualise the TRANSLATE and 1km grid
# # Take a subset of the data
# using Plots
# idx = grid.east .> 20000 .&& grid.east .< 80000 .&& 
#         grid.north .> 25000 .&& grid.north .< 80000

# idx = grid.east .> 3.1e5 .&& grid.east .< 3.5e5 .&& 
#         grid.north .> 2.25e5 .&& grid.north .< 2.4e5

# idx = grid.east .> 3.16e5 .&& grid.east .< 3.27e5 .&& 
#         grid.north .> 2.25e5 .&& grid.north .< 2.36e5
         

# Plots.plot(grid.east[idx],
#      grid.north[idx],
#      seriestype=:scatter,
#      markerstrokewidth=0,
#      markersize=2,
#         showaxis=true,
#     grid=false,
#     legend=false,
#     cbar=true,
#     aspect_ratio=:equal)


# xCoord = [ustrip(translate_xy.geometry.geoms[i].coords.x) for i in eachindex(translate_xy.geometry.geoms)];
# yCoord = [ustrip(translate_xy.geometry.geoms[i].coords.y) for i in eachindex(translate_xy.geometry.geoms)];
    
# Plots.plot!(xCoord[grid.translateID[idx]],
#      yCoord[grid.translateID[idx]],
#      seriestype=:scatter,
#      markerstrokewidth=0,
#      markersize=2,
#      markercolor=:red,
#     showaxis=false,
#     grid=false,
#     legend=false,
#     cbar=false,
#     aspect_ratio=:equal)

#     # using Makie
#     # using GLMakie
# using CairoMakie

# xs = grid.east[idx]
# ys = grid.north[idx]
# xs2 = xCoord[grid.translateID[idx]]
# ys2 = yCoord[grid.translateID[idx]]
# us = xs .- xs2
# vs = ys .- ys2
# strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))

# f = Figure(size = (800, 800))
# Axis(f[1, 1], 
#  xlabel = "Eastings (m)", 
#  ylabel = "Northings (m)",
#  xlabelsize=28,
#     ylabelsize=28,
#  backgroundcolor = "white")


# CairoMakie.scatter!(xs,ys, color = :blue, markersize = 15, marker=:circle)
# CairoMakie.scatter!(xs2,ys2, color = :red, markersize = 15, marker=:xcross)
# CairoMakie.arrows2d!(xs2, ys2, us, vs, lengthscale = 1, color = :red)
# f

 