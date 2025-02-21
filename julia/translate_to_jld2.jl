# Script to import TRANSLATE project meteo files and convert them 
# to jld2 format for more efficient import
# 
# The files have already been interpolated onto the Irish 1km grid using the 
# R script import_translate.R 
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================
 


using CSV;
using JLD2;
using DataFrames;
using Interpolations;
using NetCDF;

include("import_functions.jl")  # Includes functions to import meteo data


rcpList = ["26", "45","85"]              # The RCP scenario (26 = RCP2.6, 45=RCP4.5, 85=RCP8.5)
periodList = ["2021-2050", "2041-2070"]  # The future scenario (2021-2050 or 2041-2070)
gridFile = "IE_grid_locations.csv"  # File containing a 1km grid of lats and longs over Ireland 


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
    translateDir =  "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//TRANSLATE//tmean"

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



# =========================================================
# =========================================================
# Import the data for each RCP and period and save in a JLD2 file

for r in eachindex(rcpList)
    for p in eachindex(periodList)
    # +++++++++++++++++++++++++++++++++++++++++++++++++++
    # Import the TRANSLATE data
    @info "Importing data for RCP " * string(rcpList[r]) * " and period " * string(periodList[p])
    # Tavg, DOY, ID = read_TRANSLATE(translateDir, grid, rcpList[r], periodList[p])


    # Read one meteo file to find locations information and calculate a filter
    rcpDir = filter(x -> occursin(string(rcpList[r]), x), readdir(joinpath(translateDir)))
    files =  filter(x -> occursin(string(periodList[p]), x), readdir(joinpath(translateDir, rcpDir[1])))

    ncinfo(joinpath(translateDir, rcpDir[1], files[1]))

    # Interpolate TRANSLATE data onto 1km grid
    if r==1 && p==1
      lat = ncread(joinpath(translateDir, rcpDir[1], files[1]), "lat")
      lon = ncread(joinpath(translateDir, rcpDir[1], files[1]), "lon")
    end

    # Tavg = Array{Union{Float32, Missing},3}(undef,length(longitude), length(latitude), 365)

    Tavg = ncread(joinpath(translateDir, rcpDir[1], files[1]), "tmean")
    lon_2D = [lon[i] for i in eachindex(lon), j in eachindex(lat)]
    lat_2D = [lat[j] for i in eachindex(lon), j in eachindex(lat)]



    ind = Tavg[:,:,1] .>0
tmp = Tavg[:,:,1]

    df = georef(lon = lon_2D[ind], 
                lat = lat_2D[ind], 
                tmp = tmp[ind], 
                ("lon", "lat"))

    interp_linear2 = linear_interpolation((, ), tmp);

    Tavg_interp = interp_linear.(grid.longitude, grid.latitude)
    Tavg_interp2 = interp_linear2.(grid.longitude, grid.latitude)

    # Save this to a JLD2 file
    outfile = joinpath([outDir, "TRANSLATE_Tavg_" * string(rcpList[r]) * "_" * tring(periodList[p]) * ".jld2"])
    @info "Saving file " * outfile
    jldsave(outfile; Tavg=Tavg_IE, DOY=DOY_IE, ID=ID_IE)   
    end
end



scatter(grid.longitude, grid.latitude, marker_z = Tavg_interp, markerstrokewidth=0, markersize=1)


ind = ismissing.(Tavg_interp2)
scatter(grid.longitude[.!ind], grid.latitude[.!ind], 
marker_z = Tavg_interp2[.!ind], 
markerstrokewidth=0,
markersize=1)