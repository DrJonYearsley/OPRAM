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
# 18th June 2025
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# ===============================================================
# ===============================================================
# Load packages that will be used



using TOML
using Distributed;
using CSV;
using JLD2;
using Dates;
using Statistics;
using DataFrames;
using SharedArrays;



# =========================================================
# =========================================================

# Include the functions to import the data and run the degree day model
include("OPRAM_io_functions.jl")

# ============================================================================================
# =============== Import parameter values =======================================================




# Model parameters are stored in a TOML file https://toml.io/en/
if length(ARGS)==1
  nNodes, run_params, species_params, paths =  process_parameters(ARGS[1])

elseif length(ARGS)==0 & isfile("parameters.toml")
  nNodes, run_params, species_params, paths =  import_parameters("parameters.toml")

else
  @error "No parameter file given"
end

# =============== End of parameter setup =====================================================
# ============================================================================================





# ===============================================================
# ================ Set up parallel computing======================


if isinteractive() & nprocs() == 1
  # Enable multiple nodes 
  addprocs(nNodes)
end

# Load additional packages on all compute nodes
@everywhere using SharedArrays
@everywhere using DataFrames;

# This must be included ater the worker nodes have been assigned
include("OPRAM_ddmodel_functions.jl")


@info "Workers : $(workers())"
@info "Interactive Session: $(isinteractive())"


# =========================================================
# =========================================================
# Start the main loop of the program
if run_params.TRANSLATE_future
  # If a climate scenario is specified, run the model for future climates
  using Distributed;
  using Distributions, Random;

  @info "Running OPRAM model for future climate scenarios"

  # Run future climate degree-day model
  @time "OPRAM future model run complete:" run_model_futures(run_params, species_params, paths)
else
  # If no climate scenario is specified, run the model for past climates

  @info "Running OPRAM model for past climate years" 

  # Run past climate degree-day model
  @time "OPRAM model run complete:" run_model(run_params, species_params, paths)
end

# =========================================================
# =========================================================
if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


