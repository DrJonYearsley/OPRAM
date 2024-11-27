# Import two years of data and compare them

using Distributed
using CSV;
using JLD2;
using SharedArrays
using DataFrames

years = 2020:2021
gridFile = "IE_grid_locations.csv"  # File containing a 1km grid of lats and longs over Ireland 




# Specify directories for import and export of data
if isdir("//home//jon//Desktop//OPRAM")
    outDir = "//home//jon//Desktop//OPRAM//results//"
    dataDir = "//home//jon//DATA//OPRAM//"
    meteoDir_IE = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"
    meteoDir_NI = "//home//jon//DATA//OPRAM//Northern_Ireland_Climate_Data//"
  
  elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
    outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
    dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
    meteoDir_IE = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"
    meteoDir_NI = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Northern_Ireland_Climate_Data//"
  
  elseif isdir("//users//jon//Desktop//OPRAM//")
    outDir = "//users//jon//Desktop//OPRAM//results//"
    dataDir = "//users//jon//Desktop//OPRAM//"
    meteoDir_IE = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
    meteoDir_NI = "//users//jon//Desktop//OPRAM//Northern_Ireland_Climate_Data//"
  end


  # Include the functions to import the data 
include("GDD_functions.jl")

# =========================================================
# =========================================================
# Import location data and thin it using thinFactor

# Import grid location data
grid = CSV.read(joinpath([dataDir, gridFile]), DataFrame);

# Sort locations in order of IDs
grid = grid[sortperm(grid.ID), :];

# Thin the locations  using the thinFactor=1
thinInd = findall(mod.(grid.east, (1 * 1e3)) .< 1e-8 .&& mod.(grid.north, (1 * 1e3)) .< 1e-8);

grid_thin = grid[thinInd, :];


meteoIE_1 = read_meteoIE(meteoDir_IE, grid_thin, years[1])
meteoNI_1 = read_meteoNI(meteoDir_NI, grid_thin, years[1])


meteoIE_2 = read_meteoIE(meteoDir_IE, grid_thin, years[2])
meteoNI_2 = read_meteoNI(meteoDir_NI, grid_thin, years[2])


meteoNI_1[1]
meteoNI_2[1]


meteoIE_1[1]
meteoIE_2[1]