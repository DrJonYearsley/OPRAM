# OPRAM / R

R code for the OPRAM project

## Main R Scripts

Filename   |      Description

----------| ----------------
OPRAM_generate_geotiff.R  | Convert the OPRAM final results into geotiff files. These files are for users to download.
create_locations_file.R | Create the file containing the coordinates (bottom left corner, using Irish Grid, EPSG::29903) for all 1km and 10km squares for the Republic of Ireland and Northern Ireland. These files are call IE_grid_locations.csv, IE_hectad_locations.csv
create_locations_file_UK.R | Similar to create_locations_file.R but for the UK (including Northern Ireland). Coordinates are on UK Grid (EPSG::27700)
create_locations_file_FR.R | Similar to create_locations_file.R but for France with a grid defined on a 0.1 degree resolution (corresponding to EOBS climate data). Coordinates are latitude and longitude (WGS84, EPSG::4326)
eastnorth2os.R     |  Convert eastings and northings to a grid reference (works on Irish grid or UK Grid
OPRAM_visualisations.R  | View output from OPRAM simulations  (see also Julia file OPRAM_visualise.jl)




## Secondary R scripts


Filename   |      Description

----------| ----------------
import_translate.R  | Import Met Eireann's TRANSLATE data (deprecated in favour of Julia program translate_to_jld2.jl)
compareOpram_vs_MetOffice.R    | Import OPRAM results using UK Met Office data and Met Eireann data and quantify the discrepancy /Users/jon//view_records.R   | View GBIF and NBN records (in development)