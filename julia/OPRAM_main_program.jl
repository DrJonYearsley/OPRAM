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
# =============== Import parameter values =======================================================


using TOML

# Model parameters are stored in a TOML file https://toml.io/en/
if length(ARGS)==1
  params = TOML.parsefile(ARGS[1])

elseif length(ARGS)==0 & isfile("parameters.toml")
  params = TOML.parsefile("parameters.toml")

else
  @error "No parameter file given"
end




nNodes = params["runtime"]["nNodes"];           # Number of compute nodes to use (if in interactive)

# Put parameters into a named tuple
if in("simYears",keys(params["model"]))
  # Run model on past data
  run_TRANSLATE_future = false
  run_params = (years = params["model"]["simYears"], # Years to run model
                maxYears = params["model"]["maxYears"], # Maximum number of years to complete insect development
                lastMeteoYear = params["model"]["lastMeteoYear"], # Maximum number of years to complete insect development
                country = params["model"]["country"],         # Can be "IE", "NI" or "AllIreland"
                saveJLDFile = params["runtime"]["save2file"], # If true save the full result to a JLD2 file
                thinFactor = params["model"]["thinFactor"],   # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
                gridFile = params["inputData"]["gridFile"]);          # File containing a 1km grid of lats and longs over Ireland 
                # (used for daylength calculations as well as importing and thining of meteo data)
elseif in("rcp",keys(params["model"]))
  # Run model on future data
  run_TRANSLATE_future = true
  run_params = (futurePeriod = params["model"]["futurePeriod"], # Years to run model
                rcp = params["model"]["rcp"],                 # Future RCP scenario
                nReps = params["model"]["nReps"],             # Number of replicate climate scenarios
                maxYears = params["model"]["maxYears"],       # Maximum number of years to complete insect development
                country = params["model"]["country"],         # Can be "IE", "NI" or "AllIreland"
                saveJLDFile = params["runtime"]["save2file"], # If true save the full result to a JLD2 file
                thinFactor = params["model"]["thinFactor"],   # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
                gridFile = params["inputData"]["gridFile"]);  # File containing a 1km grid of lats and longs over Ireland 
  # (used for daylength calculations as well as importing and thining of meteo data)
else
  error("parameter file doesn't contain simulation years")  
end


# Define species (predefined)
# species_setup = (speciesFile = joinpath(homedir(),"git_repos/OPRAM/data/species_parameters.csv"),  # File containing species parameters
                  # speciesStr = ["anxius","frugiperda", "duplicatus", "cembrae", "decemlineata", " halys", "typo"])  # A vector of strings to uniquely identify a species name in the speciesFile

species_setup = (speciesFile=joinpath(homedir(), params["inputData"]["speciesFile"]),  # File containing species parameters
  speciesStr=params["model"]["speciesList"])  # A vector of strings to uniquely identify a species name in the speciesFile



# # User defined options
# species_setup = (speciesFile = joinpath(homedir(),"git_repos/OPRAM/data/userdefined_parameters.csv"),  # File containing species parameters
#                   speciesStr = ["base0_thresh200", "base5_thresh200", "base10_thresh200", "base15_thresh200",
#                   "base0_thresh400", "base5_thresh400", "base10_thresh400", "base15_thresh400",
#                   "base0_thresh600", "base5_thresh600", "base10_thresh600", "base15_thresh600",
#                   "base0_thresh800", "base5_thresh800", "base10_thresh800", "base15_thresh800",
#                   "base0_thresh1000", "base5_thresh1000", "base10_thresh1000", "base15_thresh1000"])  # A vector of strings to uniquely identify a species name in the speciesFile


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
# If the given paths for output and data do not exist, try to find alternative paths to the data

# Add on the home directory to the paths
params["paths"]["data"] = joinpath(homedir(), params["paths"]["data"])
params["paths"]["output"] = joinpath(homedir(), params["paths"]["output"])
params["paths"]["meteoIE"] = joinpath(homedir(), params["paths"]["meteoIE"])
params["paths"]["meteoNI"] = joinpath(homedir(), params["paths"]["meteoNI"])

# Check output dir exists
if !isdir(params["paths"]["output"])
  @warn "Can't find output directory!"

  if isdir(joinpath(homedir(), "Desktop//OPRAM//results"))  # Linux guess
    @info "      Setting output directory to Desktop//OPRAM//results"
    params["paths"]["output"] = joinpath(homedir(), "Desktop//OPRAM//results")

  elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//results"))  # Mac guess
    @info "        Setting output directory to Google Drive//My Drive//Projects//DAFM_OPRAM//results"
    params["paths"]["output"] = joinpath(homedir(),"Google Drive//My Drive//Projects//DAFM_OPRAM//results")
  else
    @error "No output directory found"
  end
end

# Check data dir exists
if !isdir(params["paths"]["data"])
  @info "Can't find data directory!"
  if isdir(joinpath(homedir(), "DATA", "OPRAM","Data"))  # Linux guess
    @info "       Setting data directory to DATA//OPRAM//Data"
    params["paths"]["data"] = joinpath(homedir(),"DATA//OPRAM//Data")

  elseif isdir(joinpath(homedir(), "git_repos/OPRAM/data"))  # Mac guess
    @info "       Setting data directory to home directory"
    params["paths"]["data"] = homedir()
  else
    @error "No data directory found"
  end
end

# Check meteo dirs exist (if simulation needs the data)
if in(params["model"]["country"],["IE","AllIreland"]) & !isdir(params["paths"]["meteoIE"])
  @info "Can't find meteoIE directory!"
  if isdir(joinpath(homedir(), "DATA", "OPRAM","Data","Climate_JLD2"))  # Linux guess
    @info "       Setting meteoIE directory to DATA//OPRAM//Data//Climate_JLD2"
    params["paths"]["meteoIE"] = joinpath(homedir(), "DATA//OPRAM//Data//Climate_JLD2")

  elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"))  # Mac guess
    @info "       Setting meteoIE directory to Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"
    params["paths"]["meteoIE"] = joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2")
  else
    @error "No meteoIE directory found"
  end
end

if in(params["model"]["country"],["NI","AllIreland"]) & !isdir(params["paths"]["meteoNI"])
  @info "Can't find meteoNI directory!"
  if isdir(joinpath(homedir(), "DATA", "OPRAM","Data","Climate_JLD2"))  # Linux guess
    @info "       Setting meteoNI directory to DATA//OPRAM//Data//Climate_JLD2"
    params["paths"]["meteoNI"] = joinpath(homedir(), "DATA//OPRAM//Data//Climate_JLD2")

  elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"))  # Mac guess
    @info "       Setting meteoNI directory to Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"
    params["paths"]["meteoNI"] = joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2")
  else
    @error "No meteoNI directory found"
  end
end

# Set folders to nothing if no data requried
if params["model"]["country"] == "IE"
  params["paths"]["meteoNI"] = nothing
elseif params["model"]["country"] == "NI"
  params["paths"]["meteoIE"] = nothing
end  


# Put all the paths together
paths = (outDir=joinpath(homedir(),params["paths"]["output"]),
          dataDir=joinpath(homedir(),params["paths"]["data"]),
          meteoDir_IE = params["paths"]["meteoIE"],
          meteoDir_NI = params["paths"]["meteoNI"])


# =============== End of parameter setup =====================================================
# ============================================================================================





# ===============================================================
# ================ Set up parallel computing======================

using Distributed;
if isinteractive() & nprocs() == 1
  # Enable multiple nodes 
  addprocs(nNodes)
end


@info "Workers : $(workers())"
@info "Interactive Session: $(isinteractive())"



# ===============================================================
# ===============================================================
# Load packages that will be used
using Distributed;
using CSV;
using JLD2;
using Dates;
using Statistics;

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
if run_TRANSLATE_future
  # If a climate scenario is specified, run the model for future climates
  using Distributed;
  using Distributions, Random;

  @info "Running OPRAM model for future climate scenarios"

   @time "OPRAM future model run complete:" run_model_futures(run_params, species_setup, paths)
else
  # If no climate scenario is specified, run the model for past climates

  @info "Running OPRAM model for past climate years" 

  @time "OPRAM model run complete:" run_model(run_params, species_setup, paths)
end

# =========================================================
# =========================================================
if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


