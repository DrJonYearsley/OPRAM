#  Create a file containing all grid locations for and EOBS climate data set
#
#  Jon Yearsley (Jon.Yearsley@ucd.ie)
#  Jan 2025
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd("~/git_repos/OPRAM/R")
require(tidyr, quietly = TRUE)
require(sf, quietly = TRUE)
require(data.table, quietly = TRUE)
require(sf, quietly=TRUE)


# Set the directory containing the data
dataDir = "//users//jon//Research//Data"

# Admin boundary data
regionFile= "GIS//FranceBoundaries//gadm41_FRA_2.shp"

# Use directories containing max daily temps
eobsFile = "Climate//eobs_0.1deg_v30.0_France_1965_1979.csv"

# Output filename
fileout = "FR_grid_locations.csv"



# Import regional data and coastline boundary
p = st_read(file.path(dataDir,regionFile))

# Make sure Ulster name is just Ulster
regionNames = p$NAME_1
deptNames = p$NAME_2

# Import the data
meteo_data = data.table::fread(file.path(dataDir,eobsFile))


locations = unique(meteo_data[,c('longitude','latitude')])
locations$country = "FR"
locations$region=NA
locations$department=NA


###################### Province
# Include department boundary info 
sp_regions <- sf::st_as_sf(locations, 
                      coords=c("longitude","latitude"),
                      crs=4326)

sp_region = st_intersects(p,sp_regions)

locations$region = NA
for (i in 1:length(deptNames)) {
  locations$region[sp_region[[i]]] = regionNames[i]
  locations$department[sp_region[[i]]] = deptNames[i]
}


# Find points not given a province and give them the nearest province.
noProv = which(is.na(locations$region))
if (length(noProv)>0) {
  sp_noProv <- sf::st_as_sf(locations[noProv,], 
                            coords=c("longitude","latitude"),
                            crs=4326)
  closest <- array(NA,dim=length(noProv))
  for(i in seq_len(nrow(sp_noProv))){
    a=p[which.min(st_distance(p, sp_noProv[i,])),]
    closest[i] <- a$OBJECTID
  }
  locations$region[noProv] = regionNames[as.numeric(closest)]
  locations$department[noProv] = deptNames[as.numeric(closest)]
}


# Plot locations to check
plot(locations$longitude, locations$latitude, pch=".")



# Create a unique ID
lon_list = sort(unique(locations$longitude))
lat_list = sort(unique(locations$latitude))

ex = ceiling(log10(max(length(lon_list), length(lat_list))))
ID = match(locations$lon,lon_list)*10^ex + match(locations$latitude,lat_list)



# Save all the locations to a CSV
df_final = data.frame(ID=ID, 
                      longitude=locations$longitude, 
                      latitude=locations$latitude,
                      country=locations$country,
                      region=locations$region,
                      department=locations$department)
write.csv(df_final, 
          file = file.path(dataDir, fileout),
          quote=FALSE,
          row.names = FALSE)



library(ggplot2)
ggplot(data=locations,
       aes(x=longitude,
           y=latitude,
           colour=region)) + 
  geom_point(size=0.05)
