This directory contains the Julia version of the degree day models

The primary file to run the model is GDD_model_parrallel.jl

Other files:


Filename  | Description
----------| --------------------------------
OPRAM_main_program.jl | The main file to start running the model
OPRAM_ddmodel_functions.jl | The underlying functions for running the degree daymodel
OPRAM_io_functions.jl | Functions to import/export data (meteo, grid and species data)
species_params.jl | File containng predefined model parameters for species (deprecated)



## Helper programs

Filename  | Description
----------| --------------------------------
visualise_GDD_results.jl | Create maps of Ireland showing key model outputs (can also display averages across years)
agrilus_ngen_2018_2020.png | Example of a map from visualise_GDD_results.jl
compare_to_R.jl | Compare output from Julia to the output from the model in R
compare_meteo.jl| Compare two years of meteo data (both ROI and NI gridded data)
matrix_model.jl | Create a metric model that iterates the model and viualises the long-term distribution of emergence dates

