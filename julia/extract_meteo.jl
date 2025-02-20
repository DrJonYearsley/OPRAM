# Extract meteo data for BIOL20060 module
# Saves data for 10 years from one random location in each of four counties


using Distributed
using CSV;
using JLD2;
using SharedArrays
using DataFrames
using StatsBase;

# ===========================================================================
# Set parameters ============================================================
years = 2011:2020                                 # Specify the 10 years
county = ["Waterford"]   # Specify the 4 counties
gridFile = "IE_grid_locations.csv"                # File containing a 1km grid
outPrefix = "BIOL20060"                           # Output file prefix
# ===========================================================================
# ===========================================================================



# Specify directories for import and export of data
if isdir("//home//jon//Desktop//OPRAM")
    outDir = "//home//jon//Desktop"
    dataDir = "//home//jon//DATA//OPRAM//"
    meteoDir_IE = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"
    meteoDir_NI = "//home//jon//DATA//OPRAM//Northern_Ireland_Climate_Data//"
  
  elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
    outDir = "//users//jon//Desktop"
    dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
    meteoDir_IE = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"
    meteoDir_NI = nothing
  
  elseif isdir("//users//jon//Desktop//OPRAM//")
    outDir = "//users//jon//Desktop"
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

# Find one location in each county
ID_sample = [sample(grid.ID[grid.county.==x]) for x in county]
df = Vector{DataFrame}(undef, length(county))

for y in eachindex(years)
  @time "Imported meteo data" Tavg, DOY, ID  = read_meteo(years[y], [meteoDir_IE, meteoDir_NI], grid_thin)

  idx = findall([x in ID_sample for x in ID])
  grid_idx = [findfirst(grid_thin.ID.==ID_sample[i]) for i in eachindex(ID_sample)]
  for i in eachindex(county)
    if y==1
       df[i] = DataFrame(ID = ID[idx[i]], 
                         east=grid_thin.east[grid_idx[i]],
                        north=grid_thin.north[grid_idx[i]],                        
                        county=grid_thin.county[grid_idx[i]],
                        year = years[y],
                        DOY=DOY[1:365],
                        Tavg = Tavg[1:365,idx[i]])
    else
        tmp = DataFrame(ID = ID[idx[i]], 
                        east=grid_thin.east[grid_idx[i]],
                        north=grid_thin.north[grid_idx[i]],
                        county=grid_thin.county[grid_idx[i]],
                        year = years[y],
                        DOY=DOY[1:365],
                        Tavg = Tavg[1:365,idx[i]]) 
        df[i] = vcat(df[i],tmp)
    end
  end
end


for i in eachindex(county)
    CSV.write(joinpath([outDir, county[i] * "_" * outPrefix * ".csv"]), df[i])
end



