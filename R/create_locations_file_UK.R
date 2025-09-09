#  Create a file containing all grid locations for UK
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
require(ncdf4)


# Set the directory containing the data
# dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data"
dataDir = "/users/jon/Dropbox/UKTemp2024/NCIC/maxtemp/2024/01/"
outDir = "/users/jon/git_repos/OPRAM/data"
# Admin boundary data
provFile_UK = "/users/jon/Research/Data/GIS/UK/Counties_and_Unitary_Authorities_May_2023_UK_BGC_-5421515367278074247.gpkg"
countryFile_UK = "/users/jon/Research/Data/GIS/UK/Countries_December_2018_FCB_UK_2022_3573476627889399049.gpkg"


# Output filename
fileout_1km = "UK_grid_locations.csv"
fileout_10km = "UK_hectad_locations.csv"
county_defs = "UK_county_defs.csv"

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define functions -----------------

convert_to_lonlat <- function(df) {
  # Convert UK OS grid eastings and northings to lat long (WGS84)
  sp_ig <- sf::st_as_sf(df, 
                        coords=c("east","north"),
                        crs=27700)
  crds <- sf::st_coordinates(sf::st_transform(sp_ig, crs=4326))
  attr(crds, "dimnames") <- list(NULL, c("longitude", "latitude"))
  return(crds)
}


# ++++++++++++++++++++++End of function definitions++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Import on UK met office data file
nc = nc_open(filename=file.path(dataDir, "01.nc"))

temp = ncvar_get(nc, varid="daily_maxtemp")

# Subtract 500 to put coordinate in bottom left of square
north_list=ncvar_get(nc, varid="projection_y_coordinate")-500
east_list=ncvar_get(nc, varid="projection_x_coordinate")-500

nc_close(nc)

uk_coords = data.frame(temp=c(temp),
  east=rep(east_list, times=ncol(temp)),
  north=rep(north_list,each=nrow(temp)))

# Remove points with no data or
# with eastings and northings <100 (so grid refs can be found)
idx = is.finite(uk_coords$temp) & uk_coords$temp<100 & 
  uk_coords$east>100 & uk_coords$north>100

uk_coords = uk_coords[idx,]

# summary(uk_coords)
# id = uk_coords$east<50000
# plot(uk_coords$east[id], uk_coords$north[id],pch=".")



# Import province, county data and coastline boundary
p = st_read(provFile_UK)
country_UK = st_read(countryFile_UK)


# Remove commas in county names
p$CTYUA23NM = gsub(",", "",p$CTYUA23NM)

# Rename countries
idx = grepl("England",country_UK$ctry18nm)
country_UK$ctry18nm[idx] = "EN"
idx = grepl("Scotland",country_UK$ctry18nm)
country_UK$ctry18nm[idx] = "SC"
idx = grepl("Northern Ireland",country_UK$ctry18nm)
country_UK$ctry18nm[idx] = "NI"
idx = grepl("Wales",country_UK$ctry18nm)
country_UK$ctry18nm[idx] = "WL"


# Convert to EPSG27700
p = sf::st_transform(p, crs=27700)
country_UK = sf::st_transform(country_UK, crs=27700)


# Add country_UK into province
p$country = NA
for (i in 1:nrow(p)) {
  tmp = st_intersection(p[i,], country_UK)
  if (nrow(tmp)!=0) {
    p$country[i] = tmp$ctry18nm
  }
}


# Check country, county and provinces (looks OK)
table(p$CTYUA23NM, p$country)






# Create a list of all unique locations (with duplicates removed)
locations = unique(uk_coords[,c('east','north')])

# Move grid coords to centre of the grid before finding intersection with
# polygons
loc_tmp = locations
loc_tmp$east = loc_tmp$east+500
loc_tmp$north = loc_tmp$north+500
sp_ukos <- sf::st_as_sf(loc_tmp, 
                      coords=c("east","north"),
                      crs=27700)


# Assign country names to each grid ============ 
locations$county = NA
locations$country = NA
sp_county = st_intersects(p, sp_ukos)
# Remember county names
countyNames = p$CTYUA23NM

for (i in 1:length(countyNames)) {
  locations$county[sp_county[[i]]] = countyNames[i]
}



# Assign province names =================
locations$country = p$country[match(locations$county, p$CTYUA23NM)]



# Remove points with missing country
ind = is.na(locations$country)
locations = locations[!ind,]


# Include hectad (i.e. 10km square)
source("eastnorth2os.R")

hectad = eastnorth2os(as.matrix(locations[,c("east","north")]+500), 
                   system="OSGB", 
                   format="hectad")


# Include latitude and longitude of each point
locations_lonlat = convert_to_lonlat(locations)


# Plot locations to check
library(ggplot2)
ggplot(data=locations,
       aes(x=east,
           y=north,
           colour=county)) +
  geom_point(size=0.1) + 
  theme(legend.position = "none")



# Create a unique ID
# east_list = sort(unique(locations$east))
# north_list = sort(unique(locations$north))

ID = match(locations$east,east_list)*1000 + match(locations$north,north_list)



# Save all the locations to a CSV
df_final = data.frame(ID=ID, 
                      east=locations$east, 
                      north=locations$north,
                      hectad=hectad$GR,
                      longitude=locations_lonlat[,1],
                      latitude=locations_lonlat[,2],
                      country=locations$country,
                      province= NA,
                      county = locations$county)

df_final$country[df_final$country=="NA"] = NA

library(ggplot2)
ggplot(data=df_final,
       aes(x=east,
           y=north)) +
  geom_point(aes(colour=country),size=0.1) + 
  theme(legend.position = "none")


# Set scipen so that numbers are not in scientific notation
options(scipen=9999)
write.csv(df_final, 
          file = file.path(outDir, fileout_1km),
          quote=FALSE,
          row.names = FALSE)




# =============================================================================
# Write a hectad file --------------------
# This file contains all unique hectads and the coordinates for the bottom left corner
df_hectads = aggregate(cbind(east,north)~hectad+country, 
                      data=df_final, 
                      FUN=function(x){floor(min(x,na.rm=TRUE)/1e4)*1e4})



# Add in lat long
hectad_lonlat = convert_to_lonlat(df_hectads)

df_hectads$longitude = hectad_lonlat[,"longitude"]
df_hectads$latitude = hectad_lonlat[,"latitude"]


write.csv(df_hectads, 
          file = file.path(outDir, fileout_10km),
          quote=FALSE,
          row.names = FALSE)


# =============================================================================
# Write county names to file with an ID --------------------


county = data.frame(County=countyNames,
                    Id=c(1:length(countyNames)))


write.csv(county, 
          file = file.path(outDir, county_defs),
          quote=FALSE,
          row.names = FALSE)
