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
using JLD2



nNodes = 2;        # Number of compute nodes to use (if in interactive)
meteoYear = 1961   # Years to run model
saveToFile = true;   # If true save the result to a file
latlonFile = "locations.CSV"  # File containing a grid of lats and longs over Ireland (used for daylength calculations)

# Factor to think the spatial grid (2 means sample every 2 km, 5 = sample every 5km)
thinFactor = 1;

# Define species parameters
outPrefix = "agrilus";   # Prefix to use for results files
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
dummy_species = (base_temperature=1.7f0,            # Degrees C
  threshold=1004.0f0,                        # Degrees C
  diapause_photoperiod=missing,              # Hours
  diapause_temperature=missing);             # Degrees C


# =========================================================
# =========================================================

# Specify directories for import and export of data

if isdir("//home//jon//Desktop//OPRAM")
  outDir = "//home//jon//Desktop//OPRAM//results//"
  meteoDir_IE = "//home//jon//Desktop//OPRAM//Irish_Climate_Data//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
  outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
  meteoDir_IE = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"

elseif isdir("//users//jon//Desktop//OPRAM//")
  outDir = "//users//jon//Desktop//OPRAM//results//"
  meteoDir_IE = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
end


# Make directory for the output prefix if one doesn't exist
if !isdir(joinpath(outDir, outPrefix))
  println("Making directory " * outPrefix)
  mkdir(joinpath(outDir, outPrefix))
end
outDir = joinpath(outDir, outPrefix)


# Make complete path to photoperiod file
latlonFile = joinpath([meteoDir_IE, latlonFile])

# =========================================================
# =========================================================

if isinteractive() & nworkers() == 1
  # Enable multiple nodes 
  addprocs(nNodes)
end


println("Workers : $(workers())")
println("interactive : $(isinteractive())")


@everywhere using SharedArrays
@everywhere using DataFrames;



# =========================================================
# =========================================================

# Include the functions to import the data and run the degree day model
include("GDD_functions.jl")
include("species_params.jl")




# =========================================================
# =========================================================

# Find species data corresponding to outPrefix and set this as the variable params
species_params = filter(x -> occursin(outPrefix, string(x)), names(Main))
if length(species_params) > 0
  println("Using parameters for species " * string(species_params[1]))
  params = eval(species_params[1])    # Define parameters to use
else
  error("Species parameters not found")
end



# =========================================================
# =========================================================
# Start the main loop of the program


for year in meteoYear
  # Import weather data along with location and day of year
  println("Importing meteo data starting from year " * string(year))
  @time meteo = read_meteo(year, meteoDir_IE, thinFactor)



  println("Preparing data for the model")

  # Create ID's for thinned locations
  ID = meteo[2]

  # Calculate average temp (first element of meteo) minus the base temp
  GDD = meteo[1] .- params.base_temperature

  # Calculate days where development updates
  gdd_update = GDD .> 0


  # List ID's for all GDDs above the base temp
  idx = findall(gdd_update)
  IDvec = [convert(Int32, idx[i][2]) for i in eachindex(idx)]       # Location ID index


  # Add in diapause (if necessary)
  if !ismissing(params.diapause_photoperiod)
    # Calulate day of year for all points where GDD > base temp
    DOY = [meteo[5][idx[i][1]] for i in eachindex(idx)]

    if ismissing(params.diapause_temperature)   # diapause determine by photoperiod
      no_overwinter = photoperiod(latlonFile, DOY, IDvec, params.diapause_photoperiod) 
    else # diapause determined by photoperiod (if day length is decreasing) and temp
      no_overwinter = photoperiod(latlonFile, DOY, IDvec, params.diapause_photoperiod).||
                      GDD[gdd_update] .> (params.diapause_temperature - params.base_temperature)
    end

    # Update gdd_update, idx and IDvec to remove overwintering days
    gdd_update[gdd_update] = no_overwinter
    idx = idx[no_overwinter]
    IDvec = IDvec[no_overwinter]

    # Free up some memory
    DOY = nothing
    latitude = nothing
    no_overwinter = nothing
  end

  # Make GDD array a shared array
  GDDsh = SharedArray{Float32,1}(GDD[gdd_update])

  # Find indices separatng different locations
  ind = vec(sum(gdd_update, dims=1))
  locInd2 = accumulate(+, ind)  # End locations
  locInd1 = vcat(1, locInd2[1:end-1] .+ 1)


  # Create a shared array to hold results for every day when GDD updates
  result = SharedArray{Int16,2}(sum(gdd_update), 3)

  # Fill the first 2 columns of results
  result[:, 1] = [idx[i][1] for i in eachindex(idx)]  # Calculate day of year
  result[:, 2] .= Int16(-1)
  result[:, 3] .= Int16(-1)

  # Free up more memory
  meteo = nothing     # Remove the meteo data
  thinInd = nothing
  gdd_update = nothing
  GDD = nothing
  idx = nothing



  # Loop over all locations and run the model
  println("Running the model")

  @everywhere thresh = convert(Float32, $params.threshold)
  @time location_loop!(locInd1, locInd2, result, GDDsh, thresh)


  # Loop over all locations and run the model
  println("Saving the results")


  # Clean up the data

  # Remove rows that have emergeDOY>0 
  idxKeep = result[:, 2] .> 0

  # Remove rows where the final prediction doesn't change (they can be recalculated later)
  idxKeep2 = (IDvec[1:end-1] .!= IDvec[2:end, 1]) .|| (IDvec[1:end-1] .== IDvec[2:end, 1] .&& result[1:end-1, 2] .!= result[2:end, 2])
  push!(idxKeep2, true)    # Add a true value at the end

  # Remove all but the first result with DOY >= 365 (results only for 1 year)
  idxKeep3 = result[2:end, 1] .< 365 .|| (IDvec[1:end-1] .== IDvec[2:end, 1] .&& (result[1:end-1, 1] .< 365 .&& result[2:end, 1] .>= 365))
  pushfirst!(idxKeep3, true)    # Add a true value at the start

  # Combine all 3 indices together
  idxKeep = idxKeep .&& idxKeep2 .&& idxKeep3


  # Create a data frame, using real location ID (second element of meteo)
  tm = DataFrame(ID=ID[IDvec[idxKeep]], DOY=result[idxKeep, 1], emergeDOY=result[idxKeep, 2])


  # # Save the results (replacing ID indices with original ID values)
  # CSV.write(joinpath([outDir,"result_" * outPrefix * string(year) * "_par_thin" * string(thinFactor) * ".csv"]), tm)

  # Save using jld2 format
  save_object(joinpath([outDir, "result_" * outPrefix * string(year) * "_par_thin" * string(thinFactor) * ".jld2"]), tm)


end

if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


