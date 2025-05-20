# OPRAM_generate_geotiff
#
# Import CSV files form the OPRAM online pest tool 
# and convert them into geotiff files
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# May 2025
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(sf)
# library(stars)
library(lubridate)
library(terra)
library(raster)

rm(list=ls())

# Specify folder containing the CSV files
dataFolder = "~/Google Drive/My Drive/Projects/DAFM_OPRAM/results/granite_output/agrilus_anxius/"
githubFolder = "~/git_repos/OPRAM/data"
years = c(1991:2024)
years=1991
species = "agrilus_anxius"

# ===================================
# Find and import data --------------


# Import grid data
granite = read.csv(file.path(githubFolder,"granite_hectad_county_defs.csv"))
grid = read.csv(file.path(githubFolder,"IE_grid_locations.csv"))
grid_hectad = read.csv(file.path(githubFolder,"IE_hectad_locations.csv"))
counties = unique(granite[,c(4,5)])

for (y in years) {
  for (c in 1:nrow(counties)) {
    # Look for CSV files which are at the 1km scale
    files = list.files(path=dataFolder, 
                       pattern=paste0(species,"_1_",counties$Id[c],"_",y,"[[:graph:]]+csv"),
                       recursive = TRUE)
    if (length(files)>0) {
    d = read.csv(file.path(dataFolder,files), 
                 header=TRUE,
                 colClasses = c("integer","Date","numeric","integer",
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
  
  starts = unique(d_final$startDate)
  for (m in 1:length(starts)) {
    fileout = paste0(species,"_1_",format(starts[m],"%Y_%b"),".tif")
    
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
  
}


