#  Import individual ensemble data for TRANSLATE TMAX data 
#  and an average across 30 years for the TRANSLATE data
#
#  Compare the weekly variability of the two data sets
#
#
#
#  Jon Yearsley (Jon.Yearsley@ucd.ie)
#  April 2026
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
gc()
require(tidyr, quietly = TRUE)
require(sf, quietly = TRUE)
require(data.table, quietly = TRUE)
require(sf, quietly=TRUE)
library(ncdf4)

rcpStr = "85"
period = "2071-2100"
# translateDir = "/Volumes/CLIMATEDATA/Translate/Ensemble_member_sample_Representative_year_DOY_values/"
translateDir = "/Volumes/CLIMATEDATA/Translate/Ensemble_member_sample_30_years_daily_data/"
varStr = "tmax"


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Open and read ncdf file
filein = list.files(pattern=paste0(varStr,"[[:graph:]]+",  period,"[[:graph:]]+",".nc"), path=translateDir, full.names = TRUE)


for (i in 1:length(filein)) {
  nc = nc_open(filename=filein[[i]])
  
  # Import temp data
  tmax = ncvar_get(nc, varid="tmax")  # Temp data (median)
  
  # Import dimensions
  latitude = ncvar_get(nc, varid="lat")
  longitude = ncvar_get(nc, varid="lon")
  t = ncvar_get(nc, varid="time")        # Time in days since 2071-1-1 00:00:00
  nc_close(nc)
  
  
  
  
  if (i==1) {
    d_long = data.frame(
      ensemble = i,
      longitude = rep(longitude, times=length(latitude)*length(t)),
      latitude = rep(latitude, each=length(longitude),times=length(t)),
      doy = rep(1:length(t), each=length(longitude)*length(latitude)),
      tmax = c(tmax)) 
  } else {
    d_long = rbind(d_long, data.frame(
      ensemble = i,
      longitude = rep(longitude, times=length(latitude)*length(t)),
      latitude = rep(latitude, each=length(longitude),times=length(t)),
      doy = rep(1:length(t), each=length(longitude)*length(latitude)),
      tmax = c(tmax)) ) 
  }
  
  # d_long = d_long[!is.na(d_long$tmax),]
}




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compare TRANSLATE and historic variability




library(dplyr)
library(ggplot2)
library(tidyr)


# Calculate mean and sd f data at each DOY
tmp = d_long |> 
  group_by(doy, ensemble) |>
  summarise(tmax_mean = mean(tmax, na.rm=TRUE),
            tmax_sd = sd(tmax,na.rm=TRUE),
            tmax_50 = quantile(tmax,prob=0.5,na.rm=TRUE),
            tmax_10 = quantile(tmax,prob=0.1,na.rm=TRUE),
            tmax_90 = quantile(tmax,prob=0.9,na.rm=TRUE))
head(tmp)

ggplot(data=tmp,
       aes(x=doy,
           y=tmax_50,
           ymin=tmax_10,
           ymax=tmax_90,
           colour=factor(ensemble),
           fill=factor(ensemble))) + 
  geom_ribbon(aes(colour=NULL),alpha=0.2) +
  geom_path(size=0.5) + 
  labs(x="Day of Year",
       y="Tmax") +
  theme_bw()




