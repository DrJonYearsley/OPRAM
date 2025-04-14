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


using TOML

# Model parameters are stored in a TOML file https://toml.io/en/
paramFile = "parameters.toml"
params = TOML.parsefile(paramFile)

nNodes = params["runtime"]["nNodes"];           # Number of compute nodes to use (if in interactive)
run_params = (years = 1991:2023,                   # Years to run model
              meteoRCP = "26",                     # Climate scenario to use (2.6, 4.5, 8.5)
              meteoPeriod = "2021-2050",           # Climate period to use (2021-2050, 2041-2070)
              maxYears = 1,                        # Maximum number of years to complete insect development
              country = "IE",                      # Can be "IE", "NI" or "AllIreland"
              saveJLDFile = true,                  # If true save the full result to a JLD2 file
              thinFactor = 1,                      # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
              gridFile = "IE_grid_locations.csv");  # File containing a 1km grid of lats and longs over Ireland 
              # (used for daylength calculations as well as importing and thining of meteo data)
              
# Define species (predefined)
species_setup = (speciesFile = joinpath(homedir(),"git_repos/OPRAM/data/species_parameters.csv"),  # File containing species parameters
                  speciesStr = ["anxius","frugiperda", "duplicatus", "cembrae", "decemlineata", " halys", "typo"])  # A vector of strings to uniquely identify a species name in the speciesFile

species_setup = (speciesFile=joinpath(homedir(), "git_repos/OPRAM/data/species_parameters.csv"),  # File containing species parameters
  speciesStr=["anxius"])  # A vector of strings to uniquely identify a species name in the speciesFile



# # User defined options
# species_setup = (speciesFile = joinpath(homedir(),"git_repos/OPRAM/data/userdefined_parameters.csv"),  # File containing species parameters
#                   speciesStr = ["base0_thresh200", "base5_thresh200", "base10_thresh200", "base15_thresh200",
#                   "base0_thresh400", "base5_thresh400", "base10_thresh400", "base15_thresh400",
#                   "base0_thresh600", "base5_thresh600", "base10_thresh600", "base15_thresh600",
#                   "base0_thresh800", "base5_thresh800", "base10_thresh800", "base15_thresh800",
#                   "base0_thresh1000", "base5_thresh1000", "base10_thresh1000", "base15_thresh1000"])  # A vector of strings to uniquely identify a species name in the speciesFile

# species_setup = (speciesFile = joinpath(homedir(),"git_repos/OPRAM/data/userdefined_parameters.csv"),  # File containing species parameters
#                   speciesStr = ["base0_thresh200"])

# Predefined species are:
#  :agrilus_anxius
#  :halyomorpha_halys
#  :ips_cembrae
#  :ips_duplicatus
#  :ips_sexdentatus
#  :ips_typographus
#  :leptinotarsa_decemlineata
#  :oulema_melanopus
#  :spodoptera_frugiperda




# =========================================================
# =========================================================

# Specify directories for import and export of data
# Set meteoDir to nothing to stop importing these data
#    e.g. meteoDir_NI = nothing

if isdir(joinpath(homedir(),"DATA//OPRAM"))       # Linux workstation
  paths = (outDir=joinpath(homedir(),"Desktop//OPRAM//results//"),
    dataDir=joinpath(homedir(),"DATA//OPRAM//"),
    meteoDir_IE=joinpath(homedir(),"DATA//OPRAM//Climate_JLD2"),
    meteoDir_NI=nothing)

elseif isdir(joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//R"))   # Mac
  paths = (outDir=joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//results//"),
    dataDir=joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"),
    meteoDir_IE=joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"),
    meteoDir_NI=nothing)

elseif isdir(joinpath(homedir(),"Desktop//OPRAM//"))
  paths = (outDir=joinpath(homedir(),"Desktop//OPRAM//results//"),
    dataDir=joinpath(homedir(),"Desktop//OPRAM//"),
    meteoDir_IE=joinpath(homedir(),"Desktop//OPRAM//Irish_Climate_Data//"),
    meteoDir_NI=nothing)
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
# Start the main loop of the program
if isnothing(run_params.meteoRCP) || isnothing(run_params.meteoPeriod)
  # If no climate scenario is specified, run the model for past climates
  @time "OPRAM model run complete:" run_model(run_params, species_setup, paths)
else
  using Distributed;
  using Distributions, Random;
   @time "OPRAM future model run complete:" run_model_futures(run_params, species_setup, paths)
end

# =========================================================
# =========================================================
if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


