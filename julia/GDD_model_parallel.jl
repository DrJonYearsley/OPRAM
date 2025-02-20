#!/snap/bin/julia -p 8

#!/opt/homebrew/bin/julia -p 3
#
# Julia script to run the degree day development model from Met Eireann's
# gridded data (read gridded data directly from Met Eireann's csv files)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 4th July 2024
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Load packages that will be used
using Distributed;
using CSV;
using JLD2;



nNodes = 3;                         # Number of compute nodes to use (if in interactive)
meteoYear = 1991:2021               # Years to run model
saveToFile = true;                  # If true save the result to a file
gridFile = "IE_grid_locations.csv"  # File containing a 1km grid of lats and longs over Ireland 
# (used for daylength calculations as well as importing and thining of meteo data)

# Factor to think the spatial grid (2 means sample every 2 km, 5 = sample every 5km)
thinFactor = 1;

# Define species parameters
outPrefix = "agrilus_anxius";   # Prefix to use for results files
# Important:
# outPrefix must correspond to part of the variable name 
# for the species parameters. For example, "dummy" if the 
# variable dummy_species is to be used in the model, or "agrius" 
# if the pre-defined parameters for agrilus_anxius are to be used.

# Predefined species are:
#  :agrilus_anxius
#  :halyomorpha_halys
#  :ips_cembrae
#  :ips_duplicatus
#  :ips_sexdentatus
#  :ips_typographus
#  :leptinotarsa_decemlineata
#  :oulema_melanopus
#  :pseudips_mexicanus
#  :spodoptera_frugiperda

# Pre-defined parameters for some species are in species_params.jl
# or you can define your own parameters below in dummy_species
dummy_species = (name="to be defined",
  base_temperature=1.7f0,            # Degrees C
  threshold=1004.0f0,                        # Degrees C
  diapause_photoperiod=missing,              # Hours
  diapause_temperature=missing);             # Degrees C


# =========================================================
# =========================================================

# Specify directories for import and export of data
# Set meteoDir to nothing to stop importing these data
#    e.g. meteoDir_NI = nothing

if isdir("//home//jon//Desktop//OPRAM")
  outDir = "//home//jon//Desktop//OPRAM//results//"
  dataDir = "//home//jon//DATA//OPRAM//"
  meteoDir_IE = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"
  meteoDir_NI = "//home//jon//DATA//OPRAM//Northern_Ireland_Climate_Data//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
  outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
  dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
  meteoDir_IE = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"
  meteoDir_NI = nothing

elseif isdir("//users//jon//Desktop//OPRAM//")
  outDir = "//users//jon//Desktop//OPRAM//results//"
  dataDir = "//users//jon//Desktop//OPRAM//"
  meteoDir_IE = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
  meteoDir_NI = "//users//jon//Desktop//OPRAM//Northern_Ireland_Climate_Data//"
end


# Make directory for the output prefix if one doesn't exist
if !isdir(joinpath(outDir, outPrefix))
  @info "Making directory " * outPrefix
  mkpath(joinpath(outDir, outPrefix))
end
outDir = joinpath(outDir, outPrefix)



# =========================================================
# =========================================================

if isinteractive() & nprocs()==1
  # Enable multiple nodes 
  addprocs(nNodes)
end


@info "Workers : $(workers())"
@info "Interactive Session: $(isinteractive())"


@everywhere using SharedArrays
@everywhere using DataFrames;



# =========================================================
# =========================================================

# Include the functions to import the data and run the degree day model
include("GDD_functions.jl")
include("species_params.jl")
include("import_functions.jl")



# =========================================================
# =========================================================

# Find species data corresponding to outPrefix and set this as the variable params
species_params = filter(x -> occursin(outPrefix, string(x)), names(Main))
if length(species_params) > 0
  @info "Using parameters for species " * string(species_params[1])
  params = eval(species_params[1])    # Define parameters to use
else
  @error "Species parameters not found"
end




# =========================================================
# =========================================================
# Import location data and thin it using thinFactor

# Import grid location data
grid = CSV.read(joinpath([dataDir, gridFile]), DataFrame);

# Remove locations not in the meteo data
if isnothing(meteoDir_IE)
  subset!(grid, :country => c->c.!="IE")
elseif isnothing(meteoDir_NI)
  subset!(grid, :country => c->c.!="NI")
end

# Sort locations in order of IDs
grid = grid[sortperm(grid.ID), :];

# Thin the locations using the thinFactor
subset!(grid, :east=> x-> mod.(x, (thinFactor * 1e3)) .< 1e-8,  :north=> x-> mod.(x, (thinFactor * 1e3)) .< 1e-8 )

# =========================================================
# =========================================================
# Calculate whether daylength allows diapause for each DOY and latitude
# Returns a matrix where rows are unique latitudes and columns are days of year
if !ismissing(params.diapause_photoperiod)
  @info "Calculating DOY threshold for diapause"
  diapause_minDOY = 150
  diapauseDOY_threshold = photoperiod(grid.latitude, diapause_minDOY, params.diapause_photoperiod)
else
  diapauseDOY_threshold = nothing
end


# =========================================================
# =========================================================
# Start the main loop of the program

@time "Total Calculation time: " for year in meteoYear

  # Import meteo data and calculate degree days
  @info "Running model for year " * string(year)
  @time "Total GDD calculation" GDDsh, idx, locInd1, locInd2 = calculate_GDD(year, params, grid, [meteoDir_IE, meteoDir_NI], diapauseDOY_threshold)

  # Create a shared array to hold results for every day when GDD updates
  result = SharedArray{Int16,2}(length(GDDsh), 3)

  # Fill the first column of results
  result[:, 1] = [idx[i][1] for i in eachindex(idx)]  # Calculate day of year
  result[:, 2] .= Int16(-1)
  result[:, 3] .= Int16(-1)


  # Loop over all locations and run the model
  @info "Running the model"
  @time "Model calculations:" begin

    @everywhere thresh = convert(Float32, $params.threshold)
    location_loop!(locInd1, locInd2, result, GDDsh, thresh)

    # Remove rows that have emergeDOY<=0 
    idxKeep = result[:, 2] .> 0

    # Extract location index for every row in result (column coord in idx)
    loc_idx = [convert(Int32, idx[i][2]) for i in eachindex(idx)]  # ID[loc_idx] will give the location's ID 
    
    # Remove rows where the final prediction at the same location doesn't change (they can be recalculated later)
    idxKeep2 = (loc_idx[1:end-1] .!= loc_idx[2:end]) .|| (loc_idx[1:end-1] .== loc_idx[2:end] .&& result[1:end-1, 2] .!= result[2:end, 2])
    push!(idxKeep2, true)    # Add a true value at the end

    # Remove all but the first result with DOY >= 365 (results only for 1 year)
    # Keep data with DOY<=365 or when DOY>365 and the previous result was less than 365
    idxKeep3 = result[2:end, 1] .<= 365 .|| (loc_idx[1:end-1] .== loc_idx[2:end] .&& (result[1:end-1, 1] .< 365 .&& result[2:end, 1] .>= 365))
    pushfirst!(idxKeep3, true)    # Add a true value at the start

    # Combine all 3 indices together
    idxKeep = idxKeep .&& idxKeep2 .&& idxKeep3


    # Create a data frame, using real location ID (second element of meteo)
    adult_emerge = DataFrame(ID=grid.ID[loc_idx[idxKeep]], 
                               DOY=result[idxKeep, 1], 
                               emergeDOY=result[idxKeep, 2])
  end

  @info "Saving the results"
  # Save using jld2 format
  save_object(joinpath([outDir, "result_" * outPrefix * string(year) * "_par_thin" * string(thinFactor) * ".jld2"]), adult_emerge)
  println(" ")
end

if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


