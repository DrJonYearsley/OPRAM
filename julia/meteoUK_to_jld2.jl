# Script to import UK Met Office meteo files and convert them to jld2 format 
# for more efficient import

using Distributed;
using NCDatasets;
using JLD2;
using CSV;
using DataFrames;

include("OPRAM_io_functions.jl")  # Includes functions to import meteo data


years = collect(2023:2024)    # The years to import
gridFile = joinpath(homedir(), "git_repos/OPRAM/data/UK_grid_locations.csv")  # File containing a 1km grid of lats and longs over Ireland 


# using CSV   

# dataset = CSV.read(download("https://mywebsite.edu/ml/machine-learning-databases/my.data"))



# =========================================================
# =========================================================

# Specify directories for import and export of data
if isdir("//home//jon//Desktop//OPRAM")
    outDir = "//home//jon//DATA//OPRAM//ClimateUK_JLD2//"
    dataDir = "//home//jon//DATA//OPRAM//"
    meteoDir = "//home//jon//DATA//OPRAM//UKTemp2024//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM")
    outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//ClimateUK_JLD2//"
    dataDir = "//Users//jon//git_repos//OPRAM//data/"
    meteoDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//UKTemp2024//NCIC//"
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
grid = CSV.read(gridFile, DataFrame);

# Sort locations in order of IDs
grid = grid[sortperm(grid.ID), :];



# =========================================================
# =========================================================
# Import temp data

# =========================================================
# ========================================================
# Import the data for each year and save each year in a JLD2 file

for y in years
    # +++++++++++++++++++++++++++++++++++++++++++++++++++
    # Import the IE data
    @info "Importing data for year " * string(y)
    if !isnothing(meteoDir)
        Tmax, Tmin, DOY, ID = read_netCDF_meteoUK(meteoDir, grid, y)

        # Save this to a JLD2 file

        outfile = joinpath([outDir, "meteoUK_Tminmax_" * string(y) * ".jld2"])
        @info "Saving file " * outfile
        jldsave(outfile; Tmax=Tmax, Tmin=Tmin, DOY=DOY, ID=ID)

        outfile = joinpath([outDir, "meteoUK_Tavg_" * string(y) * ".jld2"])
        @info "Saving file " * outfile
        jldsave(outfile; Tavg=(Tmax.+Tmin)./2, DOY=DOY, ID=ID)

    end
end
