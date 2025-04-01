#  Create a file containing all grid locations for ROI and Northern Ireland
#
#  Jon Yearsley (Jon.Yearsley@ucd.ie)
#  Nov 2024
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
setwd("~/git_repos/OPRAM/R")
require(tidyr, quietly = TRUE)
require(sf, quietly = TRUE)
require(data.table, quietly = TRUE)
require(sf, quietly=TRUE)


# Set the directory containing the data
dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data"

# Admin boundary data
provFile= "GIS//Province.shp"
countyFile_IE = "GIS//counties.shp"

county_NI = c("Antrim","Down","Armagh","Fermanagh","Tyrone","Derry")

# Use directories containing max daily temps
NI_Dir = "Northern_Ireland_Climate_Data/NI_TX_daily_grid/"
IE_Dir = "Irish Climate Data/maxtemp_grids/"

# Output filename
fileout = "IE_grid_locations_test.csv"

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define functions -----------------

convert_to_lonlat <- function(df) {
  # Convert Irish grid eastings and northings to lat long (WGS84)
  sp_ig <- sf::st_as_sf(df, 
                        coords=c("east","north"),
                        crs=29903)
  crds <- sf::st_coordinates(sf::st_transform(sp_ig, crs=4326))
  attr(crds, "dimnames") <- list(NULL, c("longitude", "latitude"))
  return(crds)
}


# ++++++++++++++++++++++End of function definitions++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Import province, county data and coastline boundary
p = st_read(file.path(dataDir,provFile))
county_IE = st_read(file.path(dataDir,countyFile_IE))

# Change Londonderry to Derry
county_IE$NAME_TAG[grepl("London",county_IE$NAME_TAG)] = "Derry"

# Add in country info
county_IE$country="IE"
county_IE$country[county_IE$NAME_TAG %in% county_NI] = "NI"


# Make sure Ulster name is just Ulster
pNames = p$geographic
pNames[grepl("Ulster", pNames)] = "Ulster"

# Remember county names
countyNames = county_IE$NAME_TAG


# Convert to EPSG29903
p = sf::st_transform(p, crs=29903)
county_IE = sf::st_transform(county_IE, crs=29903)


# Add province into county_IE
county_IE$province = NA
for (i in 1:nrow(county_IE)) {
  tmp = st_intersection(county_IE[i,], p)
  if (nrow(tmp)==0) {
    county_IE$province[i] = "Ulster"
  } else {
  provInd = which.max(st_area(tmp))
  county_IE$province[i] = pNames[as.numeric(tmp$OBJECTID[provInd])]
  }
}


# Check country, county and provinces (looks OK)
table(county_IE$province, county_IE$country)
table(county_IE$NAME_TAG, county_IE$province)


# Pick one file to import for NI and ROI
NI_file = list.files(path=file.path(dataDir,NI_Dir), 
                     pattern=".csv",
                     full.names = TRUE)[1]
IE_file = list.files(path=file.path(dataDir,IE_Dir), 
                     pattern=".csv",
                     full.names = TRUE)[1]



# Import the data
data_IE = data.table::fread(IE_file)
data_NI = data.table::fread(NI_file)



# Create a list of all unique locations (with duplicates removed)
locations = unique(rbind(data_IE[,c('east','north')],
                         data_NI[,c('east','north')]))

# Move grid coords to centre of the grid
loc_tmp = locations
loc_tmp$east = loc_tmp$east+500
loc_tmp$north = loc_tmp$north+500
sp_ig <- sf::st_as_sf(loc_tmp, 
                      coords=c("east","north"),
                      crs=29903)


# Assign county names to each grid ============ 
locations$county = NA
sp_county = st_intersects(county_IE, sp_ig)
for (i in 1:length(countyNames)) {
  locations$county[sp_county[[i]]] = countyNames[i]
}
# Find points not given a province and give them the nearest province.
noCounty = which(is.na(locations$county))
for(i in seq_len(length(noCounty))){
  a=which.min(st_distance(county_IE, sp_ig[noCounty[i],]))
  locations$county[noCounty[i]] = county_IE$NAME_TAG[a]
}

# Assign province names =================
locations$province = county_IE$province[match(locations$county, county_IE$NAME_TAG)]

# Assign country names =================
locations$country = county_IE$country[match(locations$county, county_IE$NAME_TAG)]


# Plot points with missing county
ind = is.na(locations$county)


# Include hectad (i.e. 10km square)
source("eastnorth2os.R")

hectad = eastnorth2os(as.matrix(loc_tmp[,c("east","north")]), 
                   system="OSI", 
                   format="hectad")


# Include latitude and longitude of each point
locations_lonlat = convert_to_lonlat(locations)


# Plot locations to check
# plot(locations$east, locations$north, pch=".")
# plot(locations_lonlat[,1], locations_lonlat[,2], pch=".", col=locations$county)

library(ggplot2)
ggplot(data=locations,
       aes(x=east,
           y=north,
           colour=county)) +
  geom_point(size=0.1) + 
  theme(legend.position = "none")



# Create a unique ID
east_list = sort(unique(locations$east))
north_list = sort(unique(locations$north))

ID = match(locations$east,east_list)*1000 + match(locations$north,north_list)



# Save all the locations to a CSV
df_final = data.frame(ID=ID, 
                      east=locations$east, 
                      north=locations$north,
                      hectad=hectad$GR,
                      longitude=locations_lonlat[,1],
                      latitude=locations_lonlat[,2],
                      country=locations$country,
                      province=locations$province,
                      county=locations$county)


library(ggplot2)
ggplot(data=df_final,
       aes(x=east,
           y=north)) +
  geom_point(aes(colour=county),size=0.1) + 
  theme(legend.position = "none")



write.csv(df_final, 
          file = file.path(dataDir, fileout),
          quote=FALSE,
          row.names = FALSE)


table(locations$county)


# =============================================================================
# Write a hectad file --------------------
# This file contains all unique hectads and the coordinates for the bottom left corner
d_hectads = aggregate(cbind(east,north)~hectad, 
                      data=d_final, 
                      FUN=function(x){floor(min(x,na.rm=TRUE)/1e4)*1e4})



# Add in lat long
hectad_lonlat = convert_to_lonlat(d_hectads)

