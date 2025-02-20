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
library(ncdf4)

dataDir = "//users//jon//Desktop//OPRAM//TRANSLATE//tmean"
rcpStr = "85"

list.files(pattern=rcpStr, path=dataDir, full.names = TRUE)

nc = nc_open(filename=paste0("tmean_CMIP5_rcp",rcpStr,"_2021-2050_ens50_final_DOY.nc"))
ncatt_get(nc, varid="time")

d_2035=ncvar_get(nc, varid="tmean")

tmp = d_2035[100:102, 100:102,1:3]

# Import dimensions
latitude = ncvar_get(nc, varid="lat")
longitude = ncvar_get(nc, varid="lon")
t_2035 = ncvar_get(nc, varid="time")        # Time in days since 2005/1/1
nc_close(nc)


d_long = data.frame(
  longitude = rep(longitude, times=length(latitude)*length(t_2035)),
  latitude = rep(latitude, each=length(longitude),times=length(t_2035)),
  doy = rep(1:length(t_2035), each=length(longitude)*length(latitude)),
  tmean = c(d_2035)
)

d_long = d_long[d_long$tmean!=0,]

write.csv(d_long[d_long$tmean!=0,], file="~/Desktop/temp_2035.csv",  quote=FALSE, row.names=FALSE)

# Check data is correct
tmp = d_long[sample.int(1000, size=10),]
indlon = match(tmp$longitude, longitude)
indlat = match(tmp$latitude, latitude)

for (i in 1:length(indlon)) {
  print(c(d_2035[indlon[i], indlat[i], tmp$doy[i]], tmp$tmean[i]) )
}

tmp = subset(d_long, doy==150 & tmean!=0)

library(ggplot2)
ggplot(data=tmp,
       aes(x=longitude,
           y=latitude,
           colour=tmean)) + 
  geom_point(size=0.5) + 
  scale_colour_distiller("Mean Temp\ndoy=150", palette="OrRd") +
  theme_bw()




nc = nc_open(filename=paste0("tmean_CMIP5_rcp",rcpStr,"_2041-2070_ens50_final_DOY.nc"))
d_2055=ncvar_get(nc, varid="tmean")
t_2055 = ncvar_get(nc, varid="time")        # Time in days since 2005/1/1
nc_close(nc)



d_long = data.frame(
  longitude = rep(longitude, times=length(latitude)*length(t_2055)),
  latitude = rep(latitude, each=length(longitude),times=length(t_2055)),
  doy = rep(1:length(t_2055), each=length(longitude)*length(latitude)),
  tmean = c(d_2055)
)



write.csv(d_long[d_long$tmean!=0,], file="~/Desktop/temp_2055.csv",  quote=FALSE, row.names=FALSE)


carlow_latlon = c(which.min(abs(latitude-52.836)),which.min(abs(longitude-(-6.9246)))) 

carlow_2035 = data.frame(t=as.Date(t_2035, origin="2005/01/01"),
                         temp50=d_2035[carlow_latlon[2], carlow_latlon[1], ])
carlow_2055 = data.frame(t=as.Date(t_2055, origin="2005/01/01"),
                         temp50=d_2055[carlow_latlon[2], carlow_latlon[1], ])



nc = nc_open(filename=paste0("tmean_CMIP5_rcp",rcpStr,"_2021-2050_ens10_final_DOY.nc"))
d_2035=ncvar_get(nc, varid="tmean")
nc_close(nc)
nc = nc_open(filename=paste0("tmean_CMIP5_rcp",rcpStr,"_2041-2070_ens10_final_DOY.nc"))
d_2055=ncvar_get(nc, varid="tmean")
nc_close(nc)

carlow_2035$temp10 = d_2035[carlow_latlon[2], carlow_latlon[1], ]
carlow_2055$temp10 = d_2055[carlow_latlon[2], carlow_latlon[1], ]

nc = nc_open(filename=paste0("tmean_CMIP5_rcp",rcpStr,"_2021-2050_ens90_final_DOY.nc"))
d_2035=ncvar_get(nc, varid="tmean")
nc_close(nc)
nc = nc_open(filename=paste0("tmean_CMIP5_rcp",rcpStr,"_2041-2070_ens90_final_DOY.nc"))
d_2055=ncvar_get(nc, varid="tmean")
nc_close(nc)

carlow_2035$temp90 = d_2035[carlow_latlon[2], carlow_latlon[1], ]
carlow_2055$temp90 = d_2055[carlow_latlon[2], carlow_latlon[1], ]

carlow = rbind(carlow_2035, carlow_2055)
carlow$year = format(carlow$t,"%Y")

library(ggplot2)
library(tidyr)


carlow_df = pivot_longer(carlow, 
                         cols=c(2:4),
                         names_to="quantile",
                         names_prefix="temp",
                         values_to="tmean")
carlow_df$doy = as.numeric(format(carlow_df$t,"%j"))

ggplot(subset(carlow_df,year=="2035"),
       aes(x=doy,
           y=tmean,
           linewidth=quantile,
           linetype=quantile)) + 
  geom_path() +
  labs(x="Date",
       y="Mean daily temperature (deg C)",
       title=paste0("2035 Daily Temperature for Carlow\n(Met Éireann TRASLATE data, RCP ",rcpStr,")")) +
  scale_linetype_manual("Ensemble\nQuantile",values=c(2,1,2)) +
  scale_linewidth_manual("Ensemble\nQuantile",values=c(0.5,1.5,0.5)) +
  scale_colour_brewer("Year", palette="Set2") +
  theme_bw() + 
  theme(text=element_text(size=18))


ggplot(carlow_df,
       aes(x=doy,
           y=tmean,
           linewidth=quantile,
           linetype=quantile,
           colour=year)) + 
  geom_path() +
  labs(x="Date",
       y="Mean daily temperature (deg C)",
       title=paste0("2035 Daily Temperature for Carlow\n(Met Éireann TRASLATE data, RCP ",rcpStr,")")) +
  scale_linetype_manual("Ensemble\nQuantiles",values=c(2,1,2),labels=c("10th","50th","90th")) +
  scale_linewidth_manual("Ensemble\nQuantiles",values=c(0.5,1.5,0.5),labels=c("10th","50th","90th")) +
  scale_colour_brewer("Year", palette="Dark2") +
  theme_bw() + 
  theme(text=element_text(size=18))
