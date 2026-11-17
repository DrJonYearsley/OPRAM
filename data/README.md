Some data files needed to run the degree day model and visualise results

Filename | Description
----------| ----------------
IE_grid_locations.csv | CSV file containing all unique 1km grid locations for ROI and NI.  Variables are ID (unique location ID), east (eastings), north (northings), hectad (10km hectad grid reference), , longitude (degrees), latitude (degrees), country (one of  IE or NI), province (one of Connacht, Munster, Leinster, Ulster). Coordinates are for the bottom left corner of each grid square (Irish grid CRS EPSG:29903).     
IE_hectad_locations.csv | CSV file containing all unique 10km (hectad) grid locations for ROI and NI.  Variables are hectad (unique hectad grid reference), east (eastings), north (northings), longitude (degrees), latitude (degrees). Coordinates are for the bottom left corner of each hectad (Irish grid CRS EPSG:29903).    
UK_grid_locations.csv | CSV file containing all unique 1km grid locations for UK.  Variables are ID (unique location ID), east (eastings), north (northings), hectad (10km hectad grid reference), longitude (degrees), latitude (degrees), country (always UK), province (always NA), county. Coordinates are for the bottom left corner of each grid square. (eastings and northings are for UK grid reference, CRS EPSG:27700)    
UK_county_defs.csv | Table giving numerical ID's for each county in the UK    
UK_hectad_county_defs.csv | Table giving the relationship between hectads and counties and hectad lat-longs for the UK      
coastline.csv | CSV file with points around the coast of Ireland (spacing roughly 500m). These data are used to add a coastline onto maps of results   
FR_grid_locations.csv | CSV file containing latitude-longitude grid over all of France (based upon a 0.1 degree spacing used by the EOBS climate data).  Variables are ID (unique location ID), longitude (degrees), latitude (degrees), country code (i.e. FR), region (area within France), department (name of the department within a region). Coordinates are for the bottom left corner of each grid square.    
species_parameters.csv | A file containing degree day model parameters for 10 insect species. Variables are: species (species Latin name), baseline_temperature (the baseline use to calculate degree days), threshold (the accumulated degree day target for insect to develop to adult stage), diapause (the type of overwintering: quiescence, obligate, facultative), diapause_daylength (the day length at which diapause begins, if applicable), diapause_temperature (the mean daily temperature at which diapause begins, if applicable) Reference 1, 2, 3 (DOI to the sources of the data)    
userdefined_parameters.csv | A file containing degree day model parameters for user-defined species. This is used in the OPRAM web app to allow a user to vary model parameters.   
granite_hectad_defs.csv | Table giving the relationship between hectads and counties (used by granite.ie for the web page)
granite_county_defs.csv | Table giving numerical ID's for each county (used by granite.ie for the web page)
granite_hectad_county_defs.csv | Table giving the relationship between hectads and counties and hectad lat-longs (used by granite.ie for the web page, lat=long is bottom-left corner of hectad)


