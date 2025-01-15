#  Create a file containing all grid locations for ROI and Northern Ireland
#
#  Jon Yearsley (Jon.Yearsley@ucd.ie)
#  Nov 2024
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())

require(tidyr, quietly = TRUE)
require(sf, quietly = TRUE)
require(data.table, quietly = TRUE)
require(sf, quietly=TRUE)


# Set the directory containing the data
dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data"

# Admin boundary data
provFile= "GIS//Census2011_Province_generalised20m.shp"

# Use directories containing max daily temps
NI_Dir = "Northern_Ireland_Climate_Data/NI_TX_daily_grid/"
IE_Dir = "Irish Climate Data/maxtemp_grids/"

# Output filename
fileout = "IE_grid_locations.csv"

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



# Import province data and coastline boundary
p = st_read(file.path(dataDir,provFile))

# Make sure Ulster name is just Ulster
pNames = p$PROVNAME
pNames[grepl("Ulster", pNames)] = "Ulster"


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
locations_tmp = rbind(data_IE[,c('east','north')],
                  data_NI[,c('east','north')])

locations = unique(locations_tmp)

# Include hectad (i.e. 10km square)
source("eastnorth2os.R")

hectad = eastnorth2os(as.matrix(locations[,c("east","north")]+500), 
                   system="OSI", 
                   format="hectad")

# Include country and province boundary info 
# (add 500 to move to middle of 1km grid square)
sp_ig <- sf::st_as_sf(locations+500, 
                      coords=c("east","north"),
                      crs=29903)
p = sf::st_transform(p, crs=29903)
sp_province = st_covers(p, sp_ig)

locations$province = "Ulster"
locations$country = "NI"
for (i in 1:length(pNames)) {
  locations$province[sp_province[[i]]] = pNames[i]
  locations$country[sp_province[[i]]] = "IE"
}




# Include latitude and longitude of each point
locations_lonlat = convert_to_lonlat(locations)


# Plot locations to check
plot(locations$east, locations$north, pch=".")
plot(locations_lonlat[,1], locations_lonlat[,2], pch=".")



# Create a unique ID
east_list = sort(unique(locations$east))
north_list = sort(unique(locations$north))

ID = match(locations$east,east_list)*1000 + match(locations$north,north_list)



# Save all the locations to a CSV
df_final = data.frame(ID=ID, 
                      east=locations$east, 
                      north=locations$north,
                      hectad=hectad$GR,
                      country=locations$country,
                      province=locations$province,
                      longitude=locations_lonlat[,1],
                      latitude=locations_lonlat[,2])
write.csv(df_final, 
          file = file.path(dataDir, fileout),
          quote=FALSE,
          row.names = FALSE)
