# OPRAM_generate_geotiff
#
# Import CSV files from the OPRAM online pest tool 
# and convert them into geotiff files using the TM75 CRS (EPSG:29903)
#
# Add option to zip files from same year together
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# June 2025
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(sf)
library(lubridate)
library(terra)
# library(raster)

rm(list=ls())

# Specify folder containing the CSV files
dataFolder = "~/Google Drive/My Drive/Projects/DAFM_OPRAM/results/granite_output"
outFolder = "~/Desktop/gisdata/"
githubFolder = "~/git_repos/OPRAM/data"
speciesfile = "species_parameters.csv"

years = c(1991:2024)
speciesList = c("spodoptera_frugiperda","oulema_melanopus_model_2",
                "oulema_melanopus_model_1","leptinotarsa_decemlineata",
                "ips_sexdentatus","agrilus_anxius","ips_duplicatus",
                "ips_cembrae","halyomorpha_halys")

speciesList = c("oulema_melanopus_model_2",
                "oulema_melanopus_model_1","leptinotarsa_decemlineata",
                "ips_sexdentatus","agrilus_anxius","ips_duplicatus",
                "ips_cembrae","halyomorpha_halys")


zip_output = FALSE
useSpeciesFile = FALSE



# ===================================
# Find and import data --------------

# Import species names
if (useSpeciesFile & file.exists(file.path(githubFolder, speciesfile))) {
  species_info = read.csv(file.path(githubFolder, speciesfile))
  
  speciesList = species_info$species
} 

# Import grid data
granite = read.csv(file.path(githubFolder,"granite_hectad_county_defs.csv"))
grid = read.csv(file.path(githubFolder,"IE_grid_locations.csv"))
grid_hectad = read.csv(file.path(githubFolder,"IE_hectad_locations.csv"))
counties = unique(granite[,c(4,5)])

# Find data Folder for species

for (species in speciesList) {
  for (y in years) {
    
    print(paste0("Processing year ", y, " for species ",species))
    
    # **************************************
    # Import results -------
    
    for (c in 1:nrow(counties)) {
      # Look for CSV files which are at the 1km scale
      files = list.files(path=file.path(dataFolder,species), 
                         pattern=paste0(species,"_1_",counties$Id[c],"_",y,"[[:graph:]]+csv"),
                         recursive = TRUE)
      if (length(files)>0) {
        d = read.csv(file.path(dataFolder,species, files), 
                     header=TRUE,
                     colClasses = c("integer","integer","integer","Date","numeric","integer",
                                    "numeric", "numeric", "numeric", "numeric"))
        idx = match(d$ID, grid$ID)
        if (!exists("d_final")) {
          # Not needed if results contain eastings and northings
          
          # d_final = cbind(eastings=grid$east[idx] + 500, 
          #                 northings = grid$north[idx] + 500,
          #                 county = counties$Id[c],
          #                 d[,-c(1)]) 
          d_final = cbind(d, county = counties$Id[c])
        }    else {
          # tmp = cbind(eastings=grid$east[idx] + 500, 
          #             northings = grid$north[idx] + 500,
          #             county = counties$Id[c],
          #             d[,-c(1)])
          
          d_final = rbind(d_final, 
                          cbind(d, county = counties$Id[c]))
          # rm(list="tmp")
        }
        
        rm(list="d")
      }
    }
    # Put eastings and northings in the middle of the grad square
    d_final$eastings = d_final$eastings + 500
    d_final$northings = d_final$northings + 500
    
    
    # **************************************
    # Create and save geotiffs -------
    starts = unique(d_final$startDate)
    
    # Look for output folder and create dirs if it doesn't exist
    if (!dir.exists(file.path(outFolder,species))) {
      dir.create(file.path(outFolder,species))
    }
    if (zip_output) {
      outPrefix = paste0(species,"_1_",y)
      if (!dir.exists(file.path(outFolder,species,outPrefix))) {
        dir.create(file.path(outFolder,species,outPrefix))
      }
    } 
    
    
    for (m in 1:length(starts)) {
      if (zip_output) {
        fileout = file.path(outFolder,species,outPrefix,
                            paste0(species,"_1_",format(starts[m],"%Y_%m"),"_01.tif"))
      } else {
        fileout = file.path(outFolder,species,
                            paste0(species,"_1_",format(starts[m],"%Y_%m"),"_01.tif"))
        
      }
      
      d_sub = subset(d_final, startDate==starts[m], select=c("eastings","northings","emergeDOY",
                                                             "nGen","emergeDOY_30yr", "nGen_30yr",
                                                             "emergeDOY_anomaly", "nGen_anomaly", "county"))
      d_rast = rast(d_sub, 
                    type="xyz", 
                    crs="EPSG:29903", 
                    extent = ext(min(d_sub$eastings)-500,
                                 max(d_sub$eastings)+500,
                                 min(d_sub$northings)-500,
                                 max(d_sub$northings)+500))
      terra::writeRaster(d_rast, 
                         filename = fileout,
                         filetype="GTiff",
                         datatype = "FLT4S",
                         overwrite=TRUE,
                         gdal=c("COMPRESS=DEFLATE", "TFW=NO"))
    }
    
    
    if (zip_output) {
      # **************************************
      # Zip geotiffs for one complete year -------
      
      # Finally zip files
      setwd(file.path(outFolder,species))
      zip(zipfile=paste0(outPrefix,".zip"),
          files=outPrefix)
      
      # Remove the raw files
      removeFiles = list.files(path=outPrefix, 
                               pattern=outPrefix, 
                               full.names = TRUE)
      file.remove(removeFiles)
      # Remove directory
      file.remove(outPrefix)
    }
    
    if (exists("d_final")) {
      rm(list="d_final")
    } 
    if (exists("d_rast")) {
      rm(list=c("d_rast", "d_sub"))
    } 
  }
}

