# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This is the main program to run the OPRAM degree day model
#
# The program imports Met Eireann's gridded data 
# (read gridded data directly from Met Eireann's csv files or from JLD2 files)
# and uses it to run a degree day model for a range of insect species
#
# Full model output is saved in a JLD2 file
# Summary results for specific starting dates are saved at raw 
# resolution and 10km resolution in CSV files
#
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#    Copyright (C) 2026  Jon Yearsley  (Jon.Yearsley@ucd.ie)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
  nNodes, run_params, species_params, paths =  import_parameters(ARGS[1])

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


