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

# Give file for single ensemble run
translateEnsembleFile = "~/git_repos/OPRAM/translate_ensemble_tmax_aggregated.csv"

# Find translate data
translateFile = paste0("translate_tmax_rcp",rcpStr,"_",period,".Rdata")
if (file.exists(file.path("~/MEGA/Projects/DAFM_OPRAM/ExampleData",translateFile))) {
  translateDir = "~/MEGA/Projects/DAFM_OPRAM/ExampleData"
  load(file.path(translateDir,translateFile))
  fileLoaded[1] = TRUE
} else if (dir.exists("~/DATA/OPRAM/TRANSLATE/")) {
  translateDir = "~/DATA/OPRAM/TRANSLATE/tmax/"
} else {
  translateDir = NA
}


# Find Met Eireann gridded data
meteireannFile = "meteo_AllIreland_1km_2015_2024.RData"
if (file.exists(file.path("~/MEGA/Projects/DAFM_OPRAM/ExampleData",meteireannFile))) {
  griddedDir = "~/MEGA/Projects/DAFM_OPRAM/ExampleData"
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




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compare TRANSLATE and historic variability


# Calculate average 10th and 90th qunatile difference

d_long$diff = d_long$tmax_90 - d_long$tmax_10
a=aggregate(diff~doy, data=d_long, FUN=quantile, probs=c(0.1, 0.5, 0.9))


# Calculate difference between 10th and 90th quantiles for each day of year
# across all locations in Ireland
d_past = data.frame(doy=c(1:365), TX_diff_10=NA,TX_diff_50=NA,TX_diff_90=NA, TX_diff=NA, 
                    TX_future_10=a$diff[,1], 
                    TX_future_50=a$diff[,2], 
                    TX_future_90=a$diff[,3])


rm("d_long")
rm("d_TN")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load one ensemble run from TRANSLATE (30 yr run)

if (file.exists(translateEnsembleFile)) {
  tmax_ensemble = read.csv(translateEnsembleFile)
  tmax_ensemble$tmax_diff = tmax_ensemble$tmax_90-tmax_ensemble$tmax_10
  d_past = cbind(d_past, tmax_ensemble)
}



for (j in 1:length(d_TX)) {
  print(paste0("File ",j))
  d_TX[[j]]$doy = as.numeric(format(d_TX[[j]]$Dates,"%j"))
}

for (i in 1:nrow(d_past)) {
  print(paste0("DOY ",i))
  for (j in 1:length(d_TX)) {
    if (i %in% d_TX[[j]]$doy){
      if (exists("tmp")) {
        tmp = rbind(tmp, subset(d_TX[[j]], doy==i))
      } else {
        tmp = subset(d_TX[[j]], doy==i)
      }
    }
  }
  a = aggregate(TX~east+north, data=tmp, FUN=function(x) {diff(quantile(x, probs=c(0.1, 0.9)))})
  d_past$TX_diff_10[i] = quantile(a$TX, 0.1)
  d_past$TX_diff_50[i] = quantile(a$TX, 0.5)
  d_past$TX_diff_90[i] = quantile(a$TX, 0.9)
  d_past$TX_diff[i] = diff(quantile(tmp$TX, probs=c(0.1, 0.9)))
  rm("tmp")
}

save(d_past, file="TX_quantile_difference.RData")



library(dplyr)
library(ggplot2)
library(tidyr)

# tmp = rename(d_past[,-9], "2015-2024"=TX_diff)

tmp = d_past[,-9] |>
    select(doy, TX_diff_10,  TX_diff_50,  TX_diff_90, 
           TX_future_10, TX_future_50, TX_future_90,
           tmax_diff) |>
    rename("2015-2024 Q10"= TX_diff_10, 
           "2015-2024 Q50"= TX_diff_50,
           "2015-2024 Q90"= TX_diff_90,
           "2021-2050 Q10"= TX_future_10, 
           "2021-2050 Q50"= TX_future_50, 
           "2021-2050 Q90"= TX_future_90,
           "2071-2100 Q50"= tmax_diff) |>
    select(doy, "2015-2024 Q50", "2021-2050 Q50", "2071-2100 Q50") |>
    pivot_longer(cols=-1, names_to="Data", values_to="T_max")

ggplot(data=tmp,
       aes(x=doy,
           y=T_max,
           colour=Data)) + 
  geom_point(size=1) + 
  scale_colour_brewer(name="",palette="Dark2") +
  labs(x="Day of Year",
       y="Tmax (Q90-Q10)",
    title="Tmax 1km Q90-Q10,\nTRANSLATE RCP8.5 versus historic") +
  theme_bw()

ggsave("Tmax_variation.png")



