This directory contains the Julia version of the degree day models

The primary file to run the model is OPRAM_main_program.jl



Filename  | Description
----------| --------------------------------
OPRAM_main_program.jl | The main file to start running the model
OPRAM_ddmodel_functions.jl | The underlying functions for running the degree daymodel
OPRAM_io_functions.jl | Functions to import/export data (meteo, grid and species data)

## Parameter files
Parameters are defined in TOML formatted files (https://toml.io/en/) that are human readable

Filename  | Description
----------| --------------------------------
parameters.toml  |  Model parameter settings for past climates and pre-defined species
parameters_userdefined.toml |  Model parameter settings for past climates and user-defined species
parameters_future.toml   |  Model parameter settings for future climates and pre-defined species
parameters_userdefined_future.toml |  Model parameter settings for future climates and user-defined species

## Helper programs

Filename  | Description
----------| --------------------------------
OPRAM_visualise.jl  | Create some visualisations of the model (Model visualisation is better done in R with OPRAM_visualisations.R)
visualise_GDD_results.jl | Create maps of Ireland showing key model outputs (can also display averages across years)
compare_to_R.jl | Compare output from Julia to the output from the model in R
compare_meteo.jl| Compare two years of meteo data (both ROI and NI gridded data)
matrix_model.jl | Create a metric model that iterates the model and viualises the long-term distribution of emergence dates

