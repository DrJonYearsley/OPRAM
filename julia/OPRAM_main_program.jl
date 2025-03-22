#!/snap/bin/julia -p 8

#!/opt/homebrew/bin/julia -p 3
#
# This is the main program to run the OPRAM degree day model
#
# The program imports Met Eireann's gridded data 
# (read gridded data directly from Met Eireann's csv files or from JLD2 files)
# and uses it to run a degree day model for a range of insect species
#
# Full model output is saved in a JLD2 file
# Summary results for specific starting dates are saved at raw resolution and 10km resolution in CSV files
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 4th July 2024
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ============================================================================================
# =============== Set parameter values =======================================================

nNodes = 3;                         # Number of compute nodes to use (if in interactive)
run_params = (years = 1991:2021,                   # Years to run model
              maxYears = 3,                        # Maximum number of years to complete insect development
              country = "IE",                      # Can be "IE", "NI" or "AllIreland"
              saveJLDFile = true,                  # If true save the full result to a JLD2 file
              saveSummaryCSV = true,               # If true save the specific results to a CSV file
              thinFactor = 1)                      # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)

              gridFile = "IE_grid_locations.csv";  # File containing a 1km grid of lats and longs over Ireland 
              # (used for daylength calculations as well as importing and thining of meteo data)
              
# Define species 
species_setup = (speciesFile = "/users/jon/git_repos/OPRAM/data/species_parameter.csv",  # File containing species parameters
                  speciesStr = ["sex", "typo"])  # A vector of strings to uniquely identify a species name in the speciesFile

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

# Load additional packages on all compute nodes
@everywhere using SharedArrays
@everywhere using DataFrames;


# =========================================================
# =========================================================

# Include the functions to import the data and run the degree day model
include("OPRAM_io_functions.jl")
include("OPRAM_ddmodel_functions.jl")


# =========================================================
# =========================================================

# Set species parameters from the parameter file (can be more than one species)
species_params = [import_species(species_setup.speciesFile, species_setup.speciesStr[s]) for s in eachindex(species_setup.speciesStr)]

# =========================================================
# =========================================================
# Import location data and thin it using thinFactor
grid = prepare_data(joinpath([paths.dataDir, gridFile]), run_params.thinFactor, run_params.country);


# =========================================================
# =========================================================
# Start the main loop of the program
@time "OPRAM model run complete:" run_model(run_params, species_params, paths, grid)


# =========================================================
# =========================================================
if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


