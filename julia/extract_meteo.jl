# Extract meteo data for BIOL20060 module
# Saves data for 10 years from one random location in each of four counties

using CSV;
using JLD2;
using DataFrames
using StatsBase;

# ===========================================================================
# Set parameters ============================================================
meteoYear = [collect(1971:1975);collect(2016:2020)]                              # Specify the 10 years
county = ["Waterford"]       # Specify the 4 counties
gridFile = "IE_grid_locations.csv"                # File containing a 1km grid
outPrefix = "BIOL20060"                           # Output file prefix
nLocations = 1
# ===========================================================================
# =============================================

# Kate's sites
# Location 1: Wicklow mts, high altitude, low temperature. (53.189607, -6.275902)
# Location 2: Rosslare Harbour, Co. Wexford, high temperature, coastal, port. (52.24747, -6.349449)
# Location 3: Somewhere in the midlands, inland, temperate and altitude not too high or too low, not as extreme as some of the others. (53.465285, -8.094387)
# Location 4: Shannon Foynes Port, Co. Limerick, coastal, port, quite high temperatures, very low altitude. (52.611098, -9.070196)

# The lat and long for each region
regions_kate = [(53.189607, -6.275902),
    (52.24747, -6.349449),
    (53.465285, -8.094387),
    (52.611098, -9.070196)]


regions_saoirse = [(52.242722, -6.493963) ] # Crosstown, Wexford

regions_jack = [(53.292223, -8.879707) ] # Frenchfort, Galway

regions = regions_saoirse

# Specify directories for import and export of data

if isdir("//home//jon//Desktop//OPRAM")
  outDir = "//home//jon//Desktop"
  dataDir = "//home//jon//DATA//OPRAM//"
  meteoDir_IE = "//home//jon//DATA//OPRAM//Climate_JLD2"
  meteoDir_NI = "//home//jon//DATA//OPRAM//Climate_JLD2"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
  outDir = "//users//jon//Google Drive//My Drive//Teaching//BIOL20060//Data/Climate/"
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
include("import_functions.jl")

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


# Pick locations
idx_ID = zeros(Int, length(regions), nLocations)
for r in eachindex(regions)
idx_region = abs.(regions[r][1] .- grid_thin.latitude).<0.02 .&& 
        abs.(regions[r][2] .- grid_thin.longitude).<0.02;

# Pick n locations
idx_ID[r,:] = sample(grid_thin.ID[idx_region], nLocations);       
end

# Initiate a data frame for each region
df = Vector{DataFrame}(undef, length(regions))
for y in eachindex(meteoYear)
@info "Importing data for year " * string(meteoYear[y])

  Tavg, DOY, ID = read_meteo(meteoYear[y], [meteoDir_IE, meteoDir_NI], grid, 1)
  
  for r in eachindex(regions)
    # Indices of locations within 2km of region lat-long 
    # (make sure they are in the meteo data)
    idx = findall([x in idx_ID[r,:] for x in ID])

    for i in eachindex(idx)  # Loop around each idx in the meteo data
      # Find the first row with the correct location ID
      # grid_idx = findfirst(grid_thin.ID .== ID[idx[i]])
      grid_idx = findfirst(grid_thin.ID .== idx_ID[r,i])


      if y == 1 && i == 1
        df[r] = DataFrame(ID=ID[idx[i]],
          east=grid_thin.east[grid_idx],
          north=grid_thin.north[grid_idx],
          longitude=grid_thin.longitude[grid_idx],
          latitude=grid_thin.latitude[grid_idx],
          county=grid_thin.county[grid_idx],
          year=meteoYear[y],
          DOY=DOY[1:365],
          Tavg=Tavg[1:365, idx[i]])
      else
        tmp = DataFrame(ID=ID[idx[i]],
          east=grid_thin.east[grid_idx],
          north=grid_thin.north[grid_idx],
          longitude=grid_thin.longitude[grid_idx],
          latitude=grid_thin.latitude[grid_idx],
          county=grid_thin.county[grid_idx],
          year=meteoYear[y],
          DOY=DOY[1:365],
          Tavg=Tavg[1:365, idx[i]])
          df[r] = vcat(df[r], tmp)
      end
    end
  end
end


for r in eachindex(regions)
    CSV.write(joinpath([outDir, "dummy" * df[r].county[1] * "_" * outPrefix * ".csv"]), df[r])
end






# # ====================================================================================
# # Extract data for different counties

# # Find one location in each county
# ID_sample = [sample(grid.ID[grid.county.==x], 1) for x in county]
# df = Vector{DataFrame}(undef, length(county))

# for y in eachindex(years)
#   @time "Imported meteo data" Tavg, DOY, ID  = read_meteo(years[y], [meteoDir_IE, meteoDir_NI], grid_thin)

#   idx = findall([x in ID_sample for x in ID])
#   grid_idx = [findfirst(grid_thin.ID.==ID_sample[i]) for i in eachindex(ID_sample)]
#   for i in eachindex(county)
#     if y==1
#        df[i] = DataFrame(ID = ID[idx[i]], 
#                          east=grid_thin.east[grid_idx[i]],
#                         north=grid_thin.north[grid_idx[i]],                        
#                         county=grid_thin.county[grid_idx[i]],
#                         year = years[y],
#                         DOY=DOY[1:365],
#                         Tavg = Tavg[1:365,idx[i]])
#     else
#         tmp = DataFrame(ID = ID[idx[i]], 
#                         east=grid_thin.east[grid_idx[i]],
#                         north=grid_thin.north[grid_idx[i]],
#                         county=grid_thin.county[grid_idx[i]],
#                         year = years[y],
#                         DOY=DOY[1:365],
#                         Tavg = Tavg[1:365,idx[i]]) 
#         df[i] = vcat(df[i],tmp)
#     end
#   end
# end


# for i in eachindex(county)
#     CSV.write(joinpath([outDir, county[i] * "_" * outPrefix * ".csv"]), df[i])
# end



