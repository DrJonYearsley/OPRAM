This directory contains the Julia version of the degree day models

The primary file to run the model is OPRAM_main_program.jl



Filename  | Description
----------| --------------------------------
OPRAM_main_program.jl | The main file to start running the model
OPRAM_ddmodel_functions.jl | The underlying functions for running the degree daymodel
OPRAM_io_functions.jl | Functions to import/export data (meteo, grid and species data)
OPRAM_processing_functions.jl | Functions used to process the saved model data (e.g. calculate multi-year averages and anomalies)



## Parameter files
Parameters are defined in TOML formatted files (https://toml.io/en/) that are human readable.
These parameter files are used by 

  + OPRAM_main_program.jl  
  + OPRAM_calculate_average.jl  
  + OPRAM_calculate_anomaly.jl  



Filename  | Description
----------| --------------------------------
parameters.toml  |  Model parameter settings for past climates and pre-defined species
parameters_userdefined.toml |  Model parameter settings for past climates and user-defined species
parameters_future.toml   |  Model parameter settings for future climates and pre-defined species
parameters_userdefined_future.toml |  Model parameter settings for future climates and user-defined species


## Data Preparation

Filename  | Description
----------| --------------------------------
meteo_to_jld.jl  | Convert Met Eireann gridded csv files into a jld2 file
translate_to_jld.jl  | Convert Met Eireann TRANSLATE csv file into a jld2 file




## Helper programs

Filename  | Description
----------| --------------------------------
OPRAM_visualise.jl  | Create some visualisations of the model (Model visualisation is better done in R with OPRAM_visualisations.R)
visualise_GDD_results.jl | Create maps of Ireland showing key model outputs (can also display averages across years)
map_plot.jl  | Another program to plot maps of the results
matrix_model.jl | Create a metric model that iterates the model and viualises the long-term distribution of emergence dates
calculate_daylength.jl   |  Calculate day length from latitude and day of year 




