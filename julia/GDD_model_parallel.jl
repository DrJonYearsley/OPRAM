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


# ============================================================================================
# =============== Set parameter values =======================================================

nNodes = 3;                         # Number of compute nodes to use (if in interactive)
run_params = (years = 1991:1992,               # Years to run model
              maxYears = 3,                        # Maximum number of years to complete insect development
              country = "IE",                      # Can be "IE", "NI" or "AllIreland"
              saveJLDFile = false,                  # If true save the full result to a JLD2 file
              saveSummaryCSV = true,               # If true save the specific results to a CSV file
              thinFactor = 1)                       # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)

              gridFile = "IE_grid_locations.csv";  # File containing a 1km grid of lats and longs over Ireland 
              # (used for daylength calculations as well as importing and thining of meteo data)
              
# Define species 
speciesFile = "/users/jon/git_repos/OPRAM/data/species_parameter.csv";  # File containing species parameters
speciesStr = "typo"  # A string to uniquely identify a species name in the speciesFile
# If speciesStr has no match then user_species is used
user_species = (name="to be defined",
  base_temperature=1.7f0,            # Degrees C
  threshold=1004.0f0,                        # Degrees C
  diapause_daylength=missing,                # Hours
  diapause_temperature=missing);             # Degrees C

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




# =========================================================
# =========================================================

# Specify directories for import and export of data
# Set meteoDir to nothing to stop importing these data
#    e.g. meteoDir_NI = nothing

if isdir("//home//jon//Desktop//OPRAM")       # Linux workstation
  paths = (outDir="//home//jon//Desktop//OPRAM//results//",
    dataDir="//home//jon//DATA//OPRAM//",
    meteoDir_IE="//home//jon//DATA//OPRAM//Climate_JLD2",
    meteoDir_NI="//home//jon//DATA//OPRAM//Climate_JLD2")

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")   # Mac
  paths = (outDir="//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//",
    dataDir="//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//",
    meteoDir_IE="//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2",
    meteoDir_NI=nothing)

elseif isdir("//users//jon//Desktop//OPRAM//")
  paths = (outDir="//users//jon//Desktop//OPRAM//results//",
    dataDir="//users//jon//Desktop//OPRAM//",
    meteoDir_IE="//users//jon//Desktop//OPRAM//Irish_Climate_Data//",
    meteoDir_NI="//users//jon//Desktop//OPRAM//Northern_Ireland_Climate_Data//")
end

# =============== End of parameter setup =====================================================
# ============================================================================================



# Load packages that will be used
using Distributed;
using CSV;
using JLD2;
using Dates;
using Statistics;

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
include("import_functions.jl")
include("GDD_functions.jl")
# include("species_params.jl")


# =========================================================
# =========================================================

# Set species parameters from the parameter file
species_params = import_species(speciesFile, speciesStr)
if ismissing(params.species_name)
  @info "Using '" * user_species.name * "' in user_species as parameters"
  species_params = dummy_species
end

# # Specify the output directory using the species name
# outPrefix = replace(lowercase(params.species_name), " " => "_")


# # Make directory for the output prefix if one doesn't exist
# if !isdir(joinpath(paths.outDir, outPrefix))
#   @info "Making directory " * outPrefix
#   mkpath(joinpath(paths.outDir, outPrefix))
# end
# outDir = joinpath(paths.outDir, outPrefix)


# =========================================================
# =========================================================

# # Find species data corresponding to outPrefix and set this as the variable params
# species_params = filter(x -> occursin(outPrefix, string(x)), names(Main))
# if length(species_params) > 0
#   @info "Using parameters for species " * string(species_params[1])
#   params = eval(species_params[1])    # Define parameters to use
# else
#   @error "Species parameters not found"
# end



# =========================================================
# =========================================================
# Import location data and thin it using thinFactor

# Import grid location data
grid = prepare_data(joinpath([paths.dataDir, gridFile]), run_params.thinFactor, country);


# =========================================================
# =========================================================
# Start the main loop of the program
run_model(run_params, paths, grid, species_params, outPrefix)


# for y in eachindex(run_params.years)
#   # Import meteo data and calculate degree days
#   @info "Running model for year " * string(run_params.years[y])

#   @info "Importing meteo data" 
#   Tavg, DOY, ID = read_meteo(run_params.years[y], [paths.meteoDir_IE, paths.meteoDir_NI], grid, run_params.maxYears)

# # Remove grid points not in the meteo data
#   keep_ID = [in(grid.ID[i],ID) for i in eachindex(grid.ID)]
#   grid_final = grid[keep_ID,:]
#   # diapauseDOY_final = diapauseDOY[keep_ID]

#   @info "Calculate GDD" 
#   GDDsh, idx, locInd1, locInd2 = calculate_GDD_version2(Tavg, DOY, species_params, diapauseDOY)

#   # Loop over all locations and run the model
#   @info "Calculating adult emergence dates"
#   # @everywhere thresh = convert(Float32, $params.threshold)
#   result = location_loop2(locInd1, locInd2, idx, GDDsh, species_params.threshold)

#   # Simplify the results and put them in a DataFrame
#   adult_emerge = cleanup_results(result, idx, grid_final.ID)

#   if run_params.saveToFile
#     @info "Saving the results to JLD2 file"
#     # Save using jld2 format
#     outFile = joinpath([paths.outDir, "result_" * outPrefix * string(run_params.years[y]) * "_" * string(run_params.thinFactor) * "km.jld2"])
#     save_object(outFile, adult_emerge)
#   end


#   # =========================================================
#   # =========================================================
#   # Create summary outputs

#   if run_params.saveSummaryCSV
#     @info "Creating 1km and 10km summary for specific starting dates"
#       # Obtain outputs for starting dates on the first of every month
#     dates = [Date(run_params.years[y], m, 01) for m in 1:12]

#     save_results_as_csv(dates, adult_emerge, grid_final, outPrefix, paths)
#   end
#   # # Work out corresponding day of year
#   # DOY = convert.(Int32, dayofyear.(dates))

#   # Create output for specific days of year
#   output_1km = create_doy_results(adult_emerge, dates)

#   # Create 10km summary
#   output_10km = aggregate_to_hectad(output_1km, grid_final)

#   # Remove unwanted columns
#   select!(output_1km, Not([:east_idx, :north_idx, :east_hectad, :north_hectad, :startDOY, :emergeDOY ]))
#   select!(output_10km, Not([:startDOY, :emergeDOY_min ]))

#   # Add eastings and northings to 1km output
#   output_1km = rightjoin(grid[:,[:ID, :east, :north]], output_1km, on=:ID)

#   # Add bottom left eastings and northings to 10km output

#   if saveSummaryCSV
#     @info "Saving the extracted results"
#     # Save using csv format
#     outFile1km = joinpath([paths.outDir, "result_startdates_" * outPrefix * string(year) * "_" * string(thinFactor) * "km.csv"])
#     CSV.write(outFile1km, output_1km)

#     outFile10km = joinpath([paths.outDir, "result_startdates_" * outPrefix * string(year) * "_10km.csv"])
#     CSV.write(outFile10km, output_10km)
#   end

#   println(" ")
# end

if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


