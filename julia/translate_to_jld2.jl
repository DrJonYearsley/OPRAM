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

include("OPRAM_io_functions.jl")  # Includes functions to import meteo data


rcpList = ["26", "45", "85"];              # The RCP scenario (26 = RCP2.6, 45=RCP4.5, 85=RCP8.5)
periodList = ["2021-2050", "2041-2070"];   # The future scenario (2021-2050 or 2041-2070)
gridFile = "IE_grid_locations.csv";        # File containing a 1km grid of lats and longs over Ireland 
quantiles = ["10", "50", "90"];            # The quantiles of the ensemble to import

# =========================================================
# =========================================================

# Specify directories for import and export of data
if isdir("//home//jon//DATA//OPRAM")
    outDir = "//home//jon//DATA//OPRAM//Climate_JLD2//"
    dataDir = "//home//jon//DATA//OPRAM//"
    translateDir = "//home//jon//DATA//OPRAM//TRANSLATE//tmean"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
    outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2//"
    dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
    translateDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//TRANSLATE//tmean"

elseif isdir("//users//jon//Desktop//OPRAM//")
    outDir = "//users//jon//Desktop//OPRAM//Climate_JLD2//"
    dataDir = "//users//jon//Desktop//OPRAM//"
    translateDir = "//users//jon//Desktop//OPRAM//TRANSLATE//tmean"
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
grid_gt = georef((lon=grid.longitude, lat=grid.latitude), ("lon", "lat"), crs=EPSG{4326})

# Define interpolation model
interpolation_model = InterpolateNeighbors(domain(grid_gt),
                   model=NN(),
                   maxneighbors=10) 

# =========================================================
# =========================================================
# Import the data for each RCP and period and save in a JLD2 file

# Check size of data to be imported and create arrays to hold the results

# Find correct directory and file
rcpDir = filter(x -> occursin(string(rcpList[1]), x), readdir(joinpath(translateDir)))
files = filter(x -> occursin("_" * string(periodList[1]) * "_ens50", x), readdir(joinpath(translateDir, rcpDir[1])))

# Import one file of the TRANSLATE data 
Tavg = ncread(joinpath(translateDir, rcpDir[1], files[1]), "tmean");
lat = ncread(joinpath(translateDir, rcpDir[1], files[1]), "lat");
lon = ncread(joinpath(translateDir, rcpDir[1], files[1]), "lon");

# Calculate some dimensions
TavgSize = size(Tavg)
nDays = TavgSize[3];

# Reshape arrays in place
Tavg_long = reshape(Tavg,:,nDays);

lon_2D = reshape([lon[i] for i in eachindex(lon), j in eachindex(lat)],:);
lat_2D = reshape([lat[j] for i in eachindex(lon), j in eachindex(lat)], :);

# Find data the is non-zero temps
nonzero_idx = dropdims(abs.(sum(Tavg_long, dims=2)) .> 0, dims=2);

# Keep only lats and longs with non-zero temps
lon_2D = lon_2D[nonzero_idx];
lat_2D = lat_2D[nonzero_idx];

# Array to hold temperature data for the quantiles
Tavg_quantiles = Array{Float32,3}(undef, sum(nonzero_idx), nDays, length(quantiles));

# Array to hold interpolated temperature data (mean and tandard deviation)
Tmean_interp = Array{Float32,2}(undef, nDays, length(grid.ID));
Tsd_interp = Array{Float32,2}(undef, nDays, length(grid.ID));

# Create array for days of year
DOY = collect(1:nDays);


# Calculate distance between 10th and 90th quantiles in units of standard deviation
# (used to estimate standard deviation of temp))
sd_estimator = (erfinv(2*0.9 - 1) - erfinv(2*0.1 - 1)) * sqrt(2)


# =========================================================
# For each RCP and time period import data, interpolate it and save in a JLD2 file
for r in eachindex(rcpList)
    @time "Converted TRANSLATE data:" for p in eachindex(periodList)
        # +++++++++++++++++++++++++++++++++++++++++++++++++++
        # Import the TRANSLATE data
        @info "Importing data for RCP " * string(rcpList[r]) * " and period " * string(periodList[p])


        # Find the three files for a specific RCP and time period
        rcpDir = filter(x -> occursin(string(rcpList[r]), x), readdir(joinpath(translateDir)))
        files = filter(x -> occursin(r"_" * string(periodList[p]) * "_ens", x), readdir(joinpath(translateDir, rcpDir[1])))

        for q in eachindex(quantiles)
            @info "Processing quantile " * string(quantiles[q])
            # Read the Tavg data and reshape into a 2D matrix
            fileIn = filter(x -> occursin("ens" * quantiles[q], x), files)
            Tavg = reshape(ncread(joinpath(translateDir, rcpDir[r], fileIn[1]), "tmean"), : , nDays)

            # Save non-zero temp data
            Tavg_quantiles[:, :, q] .= Tavg[nonzero_idx, :];
        end

        # Calculate standard deviation from 10th and 90th quantiles
        Tsd = (Tavg_quantiles[:,:,3] .- Tavg_quantiles[:,:,1]) ./ sd_estimator;

        # Interpolate the median and standard deviations onto the 1km grid
        for doy in eachindex(DOY)
            # Create geotable of mean daily temperature data
            df_median = georef((lon=lon_2D,
                    lat=lat_2D,
                    temp_median=Tavg_quantiles[:, doy, 2]),
                ("lon", "lat"),
                crs=EPSG{4326})

                # Create geotable of mean daily temperature standard deviation
            df_sd = georef((lon=lon_2D,
                    lat=lat_2D,
                    temp_sd=Tsd[:, doy]),
                ("lon", "lat"),
                crs=EPSG{4326})

            # Interpolate the temperature data onto the 1km grid
            tmp = df_median |> interpolation_model
            Tmean_interp[doy, :] = tmp.temp_median

            tmp = df_sd |> interpolation_model
            Tsd_interp[doy, :] = tmp.temp_sd
        end


        # Save this to a JLD2 file
        outfile = joinpath([outDir, "TRANSLATE_Tavg_rcp" * string(rcpList[r]) * "_" * string(periodList[p]) * ".jld2"])
        @info "Saving file " * outfile
        jldsave(outfile; Tavg_mean=Tmean_interp, Tavg_sd=Tsd_interp, DOY=DOY, ID=grid.ID)
    end
end


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
