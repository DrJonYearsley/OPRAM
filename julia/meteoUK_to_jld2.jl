# Script to import UK Met Office meteo files and convert them to jld2 format 
# for more efficient import

using Distributed;
using CSV;
using JLD2;
using DataFrames;
using SharedArrays;

include("OPRAM_io_functions.jl")  # Includes functions to import meteo data


years = collect(2023:2024)    # The years to import
gridFile = "IE_grid_locations.csv"  # File containing a 1km grid of lats and longs over Ireland 


# using CSV   

# dataset = CSV.read(download("https://mywebsite.edu/ml/machine-learning-databases/my.data"))



# =========================================================
# =========================================================

# Specify directories for import and export of data
if isdir("//home//jon//Desktop//OPRAM")
    outDir = "//home//jon//DATA//OPRAM//Climate_JLD2//"
    dataDir = "//home//jon//DATA//OPRAM//"
    meteoDir_IE = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"
    meteoDir_NI = nothing

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
    outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2//"
    dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
    meteoDir_IE = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"
    meteoDir_NI = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Northern_Ireland_Climate_Data//"

elseif isdir("//users//jon//Desktop//OPRAM//")
    outDir = "//users//jon//Desktop//OPRAM//Climate_JLD2//"
    dataDir = "//users//jon//Desktop//OPRAM//"
    meteoDir_IE = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
    meteoDir_NI = "//users//jon//Desktop//OPRAM//Northern_Ireland_Climate_Data//"
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
# Import the data for each year and save each year in a JLD2 file

for y in years

    # +++++++++++++++++++++++++++++++++++++++++++++++++++
    # Import the IE data
    @info "Importing data for year " * string(y)
    if !isnothing(meteoDir_IE)
        Tavg_IE, DOY_IE, ID_IE = read_CSV_meteoIE(meteoDir_IE, grid, y)

        # Save this to a JLD2 file

        outfile = joinpath([outDir, "meteoIE_Tavg_" * string(y) * ".jld2"])
        @info "Saving file " * outfile
        jldsave(outfile; Tavg=Tavg_IE, DOY=DOY_IE, ID=ID_IE)
    end
    # +++++++++++++++++++++++++++++++++++++++++++++++++++
    # Import the IE data
    if !isnothing(meteoDir_NI)
        Tavg_NI, DOY_NI, ID_NI = read_CSV_meteoNI(meteoDir_NI, grid, y)

        # Save this to a JLD2 file
        outfile = joinpath([outDir, "meteoNI_Tavg_" * string(y) * ".jld2"])
        @info "Saving file " * outfile
        jldsave(outfile; Tavg=Tavg_NI, DOY=DOY_NI, ID=ID_NI)
    end
end
