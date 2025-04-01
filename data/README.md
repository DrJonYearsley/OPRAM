Some data files needed to run the degree day model and visualise results

Filename | Description
----------| ----------------
IE_grid_locations.csv | CSV file containing all unique 1km grid locations for ROI and NI.  Variables are ID (unique location ID), east (eastings), north (northings), hectad (10km hectad grid reference), country (one of  IE or NI), province (one of Connacht, Munster, Leinster, Ulster), longitude (degrees), latitude (degrees). Coordinates are for the bottom left corner of each grid square.  
IE_hectad_locations.csv | CSV file containing all unique 10km (hectad) grid locations for ROI and NI.  Variables are hectad (unique hectad grid reference), east (eastings), north (northings), longitude (degrees), latitude (degrees). Coordinates are for the bottom left corner of each hectad.  
coastline.csv | CSV file with points around the coast of Ireland (spacing roughly 500m). These data are used to add a coastline onto maps of results
FR_grid_locations.csv | CSV file containing latitude-longitude grid over all of France (based upon a 0.1 degree spacing used by the EOBS climate data).  Variables are ID (unique location ID), longitude (degrees), latitude (degrees), country code (i.e. FR), region (area within France), department (name of the department within a region). Coordinates are for the bottom left corner of each grid square.    
species_parameters.csv | A file containing degree day model parameters for 10 insect species. Variables are: species (species Latin name), baseline_temperature (the baseline use to calculate degree days), threshold (the accumulated degree day target for insect to develop to adult stage), diapause (the type of overwintering: quiescence, obligate, facultative), diapause_daylength (the day length at which diapause begins, if applicable), diapause_temperature (the mean daily temperature at which diapause begins, if applicable) Reference 1, 2, 3 (DOI to the sources of the data)    
userdefined_parameters.csv | A file containing degree day model parameters for user-defined species. This is used in the OPRAM web app to allow a user to vary model parameters.   

