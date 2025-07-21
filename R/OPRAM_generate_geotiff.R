# OPRAM_generate_geotiff
#
# Import CSV files from the OPRAM online pest tool 
# and convert them into geotiff files using the TM75 CRS (EPSG:29903)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# June 2025
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(sf)
# library(stars)
library(lubridate)
library(terra)
library(raster)

rm(list=ls())

# Specify folder containing the CSV files
dataFolder = "/media/jon/Seagate_5TB/OPRAM_results"
outFolder = "/media/jon/Seagate_5TB/OPRAM_results/gisdata/"
speciesFile = "userdefined_parameters.csv"
githubFolder = "~/git_repos/OPRAM/data"
years = c(1991:2024)

species2Import = "all"

# speciesList = c("spodoptera_frugiperda","oulema_melanopus_model_2","oulema_melanopus_model_1","leptinotarsa_decemlineata","ips_sexdentatus","agrilus_anxius","ips_duplicatus","ips_cembrae","halyomorpha_halys")
# dataFolder = "~/Google Drive/My Drive/Projects/DAFM_OPRAM/results/granite_output"
# outFolder = "~/Google Drive/My Drive/Projects/DAFM_OPRAM/results/granite_output/gisdata/"

# ===================================
# Find and import data --------------


# Import grid data
granite = read.csv(file.path(githubFolder,"granite_hectad_county_defs.csv"))
grid = read.csv(file.path(githubFolder,"IE_grid_locations.csv"))
grid_hectad = read.csv(file.path(githubFolder,"IE_hectad_locations.csv"))
counties = unique(granite[,c(4,5)])


# Check output folder exists
if (!dir.exists(outFolder)) {
  dir.create(outFolder)
}


# Import parameter list for all species
species = read.csv(file.path(githubFolder,speciesFile))
if (length(species2Import)==1 && species2Import=="all") {
  speciesList = species$species
}  else {
  # Match species2Import against species in the file
  idx = match(species2Import, species$species)
  speciesList = species$species[idx]
}

# Loop through each species
for (species in speciesList) {
  for (y in years) {
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
        if (c==1) {
          d_final = cbind(eastings=grid$east[idx] + 500, 
                          northings = grid$north[idx] + 500,
                          county = counties$Id[c],
                          d[,-c(1)]) 
        }    else {
          tmp = cbind(eastings=grid$east[idx] + 500, 
                      northings = grid$north[idx] + 500,
                      county = counties$Id[c],
                      d[,-c(1)])
          
          d_final = rbind(d_final, tmp)
          rm(list="tmp")
        }
        
        rm(list="d")
      }
    }
    
    
    # **************************************
    # Create and save geotiffs -------
    
    starts = unique(d_final$startDate)
    
    # Look for output folder and create dirs if it doesn't exist
    if (!dir.exists(file.path(outFolder,species))) {
      dir.create(file.path(outFolder,species))
    }
    outPrefix = paste0(species,"_1_",y)
    if (!dir.exists(file.path(outFolder,species,outPrefix))) {
      dir.create(file.path(outFolder,species,outPrefix))
    }
    
    
    for (m in 1:length(starts)) {
      fileout = file.path(outFolder,species,outPrefix,
                          paste0(species,"_1_",format(starts[m],"%Y_%b"),".tif"))
      
      d_sub = subset(d_final, startDate==starts[m], select=-4)
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
}

