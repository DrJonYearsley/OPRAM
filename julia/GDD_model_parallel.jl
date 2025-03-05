#!/snap/bin/julia -p 8

#!/opt/homebrew/bin/julia -p 3
#
# Julia script to run the degree day development model from Met Eireann's
# gridded data (read gridded data directly from Met Eireann's csv files or
# from JLD2 files)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 4th July 2024
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Load packages that will be used
using Distributed;
using CSV;
using JLD2;
using Dates;
using Statistics;




nNodes = 3;                         # Number of compute nodes to use (if in interactive)
meteoYear = 1991                    # Years to run model
maxYears = 3;                       # Maximum number of years to complete insect development

saveToFile = true;                  # If true save the full result to a JLD2 file
saveSummaryCSV = true;              # If true save the specific results to a CSV file

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

if isinteractive() & nprocs() == 1
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
grid, diapauseDOY_threshold = prepare_data(joinpath([dataDir, gridFile]), thinFactor, params);

# Remove locations not in the meteo data
if isnothing(meteoDir_IE)
  subset!(grid, :country => c -> c .!= "IE")
elseif isnothing(meteoDir_NI)
  subset!(grid, :country => c -> c .!= "NI")
end



# =========================================================
# =========================================================
# Start the main loop of the program

for y in eachindex(meteoYear)
  # Import meteo data and calculate degree days
  @info "Running model for year " * string(meteoYear[y])

  @info "Importing meteo data" 
  Tavg, DOY, ID = read_meteo(meteoYear[y], [meteoDir_IE, meteoDir_NI], grid, maxYears)

  @info "Calculate GDD" 
  GDDsh, idx, locInd1, locInd2 = calculate_GDD_version2(Tavg, DOY, params, diapauseDOY_threshold)

  # Loop over all locations and run the model
  @info "Calculating adult emergence dates"
  @everywhere thresh = convert(Float32, $params.threshold)
  result = location_loop2(locInd1, locInd2, idx, GDDsh, thresh)

  # Simplify the results and put them in a DataFrame
  adult_emerge = cleanup_results(result, idx, grid.ID)

  if saveToFile
    @info "Saving the results to JLD2 file"
    # Save using jld2 format
    outFile = joinpath([outDir, "result_" * outPrefix * string(meteoYear[y]) * "_" * string(thinFactor) * "km.jld2"])
    save_object(outFile, adult_emerge)
  end


  # =========================================================
  # =========================================================
  @info "Creating output for specific days of year"

  # Obtain outputs for starting dates on the first of every month
  dates = [Date(meteoYear[y], m, 01) for m in 1:12]

  # Work out corresponding day of year
  DOY = convert.(Int32, dayofyear.(dates))

  # Create output for specific days of year
  output_1km = create_doy_results(adult_emerge, DOY)

  # Create 10km summary
  output_10km = aggregate_to_hectad(output_1km, grid)

  # Remove unwanted columns
  select!(output_1km, Not([:east_idx, :north_idx, :east_hectad, :north_hectad]))

  # Add eastings and northings to 1km output
  output_1km = rightjoin(grid[:,[:ID, :east, :north]], output_1km, on=:ID)

  if saveSummaryCSV
    @info "Saving the extracted results"
    # Save using csv format
    outFile1km = joinpath([outDir, "result_startdates_" * outPrefix * string(year) * "_" * string(thinFactor) * "km.csv"])
    CSV.write(outFile1km, output_1km)

    outFile10km = joinpath([outDir, "result_startdates_" * outPrefix * string(year) * "_10km.csv"])
    CSV.write(outFile10km, output_10km)
  end

  println(" ")
end

if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


