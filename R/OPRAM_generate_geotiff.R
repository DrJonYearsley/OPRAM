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


library(RcppTOML)
library(sf)
library(lubridate)
library(terra)
library(data.table)

rm(list=ls())

zip_output = FALSE    # Zip together all files for 1 year

# Import parameter file
args <- commandArgs(trailingOnly = TRUE)

if (length(args > 0)) {
  toml_file = args[1]
 } else {
  toml_file = "../julia/parameters_future.toml"
}

params = parseTOML(toml_file)  

# Specify folder containing the CSV files
outFolder = "~/scratch/gisdata/"    # Folder to hold the geotiff output
dataFolder = params$paths$output
granite_defs = params$inputData$countyDefs
gridFile = params$inputData$gridFile
speciesFile = params$inputData$speciesFile


species = params$model$speciesList

if (!is.null(params$model$rcp)) {
  # Set up years for future climate scenarios
  translate=TRUE
  rcp = rep(params$model$rcp, times=length(params$model$futurePeriod))
  yearStr = paste("rcp", rcp, "_",
                      rep(params$model$futurePeriod, each=length(params$model$rcp)),
                          sep="")
  
  # years = array(NA, dim=length(yearStr))
  # years[grepl("2021-2050",yearStr)] = 2035
  # years[grepl("2041-2070",yearStr)] = 2055
} else {
  # Assume results are for past climates
  translate=FALSE
  yearStr = as.character(params$model$simYears)
}



# Add home directory to start of file paths if necessary
if (!grepl("^[~/\\]{1}",dataFolder)) {
  dataFolder = file.path("~", dataFolder)
}
if (!grepl("^[~/\\]{1}",speciesFile)) {
  speciesFile = file.path("~", speciesFile)
}
if (!grepl("^[~/\\]{1}",granite_defs)) {
  granite_defs = file.path("~", granite_defs)
}
if (!grepl("^[~/\\]{1}",gridFile)) {
  gridFile = file.path("~", gridFile)
}


# ===================================
# Find and import data --------------

# Import species names
species_info = read.csv(file.path(speciesFile))

# Find rows that correspond to species in list
if (any(species=="all")) {
  speciesList = species_info$species
} else {
  speciesList = array(NA, dim=length(species))
  for (s in 1:length(species)) {
    species_name = grep(pattern=gsub("[[:space:][:punct:]]+","[[:space:][:punct:]]*",species[s]),
                        species_info$species,
                        ignore.case = TRUE,
                        value=TRUE)
    # Replace all spaces with underscores and make lower case
    speciesList[s] = gsub(" ","_",tolower(species_name))
  }
}


# Import grid data
granite = read.csv(granite_defs)
grid = read.csv(gridFile)
counties = unique(granite[,c("County","Id")])

# Find data Folder for species

for (species in speciesList) {
  for (y in 1:length(yearStr)) {
    
    print(paste0("Processing year ", yearStr[y], " for species ",species))
    
    # **************************************
    # Import results -------
    
 
   for (c in 1:nrow(counties)) {
      # Look for CSV files which are at the 1km scale
      files = list.files(path=file.path(dataFolder,species), 
                         pattern=paste0(species,"_1_",counties$Id[c],"_",yearStr[y],"[[:graph:]]+csv"),
                         recursive = TRUE)
      if (length(files)>0) {
        # Use fread to speed up import
        d = data.table::fread(file.path(dataFolder,species, files), 
                     header=TRUE,
                     colClasses = c("integer","integer","integer","Date", # Spatial ID, eastings, northings and date
                                    "numeric","numeric", "numeric",       # emergeDOY
                                    "numeric", "numeric", "numeric"))     #nGen

                idx = match(d$ID, grid$ID)
        if (!exists("d_final")) {
          d_final = cbind(d, county = counties$Id[c])
        }    else {
          d_final = rbind(d_final, 
                          cbind(d, county = counties$Id[c]))
        }
        
        rm(list="d")
      }
    }
    # Put eastings and northings in the middle of the 1km grid square
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
      outPrefix = paste0(species,"_1_",yearStr[y])
      if (!dir.exists(file.path(outFolder,species,outPrefix))) {
        dir.create(file.path(outFolder,species,outPrefix))
      }
    } 
    
    
    for (m in 1:length(starts)) {
      if (translate) {
        filename = paste0(species,"_1_",format(starts[m],"%Y_%m"),"_01_rcp",rcp[y],".tif")
       } else {
        filename = paste0(species,"_1_",format(starts[m],"%Y_%m"),"_01.tif")
       }
      
      if (zip_output) {
        fileout = file.path(outFolder,species,outPrefix,filename)
      } else {
        fileout = file.path(outFolder,species,filename)
        
      }
      
      d_sub = subset(d_final, startDate==starts[m], select=c("eastings","northings","emergeDOY",
                                                             "nGen","emergeDOY_30yr", "nGen_30yr",
                                                             "emergeDOY_anomaly", "nGen_anomaly", "county"))
      d_rast = terra::rast(d_sub, 
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

