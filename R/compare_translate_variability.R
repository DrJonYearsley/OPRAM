#  Import TRANSLATE TMAX data and  Met Eireann 2024 data
#  Compare the weekly variability of the two data sets
#
#
#
#  Jon Yearsley (Jon.Yearsley@ucd.ie)
#  April 2026
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
require(tidyr, quietly = TRUE)
require(sf, quietly = TRUE)
require(data.table, quietly = TRUE)
require(sf, quietly=TRUE)
library(ncdf4)

rcpStr = "85"
period = "2021-2050"

pastYear = 2024

fileLoaded = c(FALSE, FALSE)
# Find translate data
translateFile = paste0("translate_tmax_rcp",rcpStr,"_",period,".Rdata")
if (file.exists(file.path("~/MEGAsync/Projects/DAFM_OPRAM/ExampleData",translateFile))) {
  translateDir = "~/MEGAsync/Projects/DAFM_OPRAM/ExampleData"
  load(file.path(translateDir,translateFile))
  fileLoaded[1] = TRUE
} else if (dir.exists("~/DATA/OPRAM/TRANSLATE/")) {
  translateDir = "~/DATA/OPRAM/TRANSLATE/tmax/"
} else {
  translateDir = NA
}


# Find Met Eireann gridded data
meteireannFile = "meteo_AllIreland_1km_2015_2024.RData"
if (file.exists(file.path("~/MEGAsync/Projects/DAFM_OPRAM/ExampleData",meteireannFile))) {
  griddedDir = "~/MEGAsync/Projects/DAFM_OPRAM/ExampleData"
  load(file.path(griddedDir,meteireannFile))
  fileLoaded[2] = TRUE
  
  # Subset to just one year
  d_TN_year = subset(d_TN[[1]], format(Dates, "%yyyy")==pastYear)
  
} else if (dir.exists("~/DATA/OPRAM/Irish_Climate_Data/")) {
  griddedDir = "~/DATA/OPRAM/Irish_Climate_Data/maxtemp_grids/"
} else {
  griddedDir = NA
}



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Open and read ncdf file
if (fileLoaded[1]==FALSE) {
filein = list.files(pattern=paste0(rcpStr,"_",period,"[[:graph:]]+","ens50.nc"), path=translateDir, full.names = TRUE)
nc = nc_open(filename=filein)
# ncatt_get(nc, varid="time")

tmax50 = ncvar_get(nc, varid="tmax")  # Temp data (median)


# Import dimensions
latitude = ncvar_get(nc, varid="lat")
longitude = ncvar_get(nc, varid="lon")
t = ncvar_get(nc, varid="time")        # Time in days since 2005/1/1
nc_close(nc)




# Open and read ncdf file for 10th and 90th quantiles
filein = list.files(pattern=paste0(rcpStr,"_",period,"[[:graph:]]+","ens10.nc"), path=translateDir, full.names = TRUE)
nc = nc_open(filename=filein)
tmax10 = ncvar_get(nc, varid="tmax") # Temp data
nc_close(nc)

filein = list.files(pattern=paste0(rcpStr,"_",period,"[[:graph:]]+","ens90.nc"), path=translateDir, full.names = TRUE)
nc = nc_open(filename=filein)
tmax90 = ncvar_get(nc, varid="tmax") # Temp data
nc_close(nc)

d_long = data.frame(
  longitude = rep(longitude, times=length(latitude)*length(t)),
  latitude = rep(latitude, each=length(longitude),times=length(t)),
  doy = rep(1:length(t), each=length(longitude)*length(latitude)),
  tmax_10 = c(tmax10),
  tmax_50 = c(tmax50),
  tmax_90 = c(tmax90)
)

d_long = d_long[!is.na(d_long$tmax_50),]

save(d_long, file=paste0("~/translate_tmax_rcp",rcpStr,"_",period,".Rdata"))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load Met Eireann gridded data

if (fileLoaded[2]==FALSE) {
  list.files(path=griddedDir, pattern=paste0("TX_",pastYear), full.names = TRUE)
}

tmp = subset(d_long, doy==150 )

library(ggplot2)
ggplot(data=tmp,
       aes(x=longitude,
           y=latitude,
           colour=tmax_50)) + 
  geom_point(size=0.5) + 
  scale_colour_distiller("Daily Max Temp\ndoy=150", palette="OrRd") +
  theme_bw()




