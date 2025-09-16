# Compare the UK Met Office emergence date with the OPRAM result

rm(list=ls())
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)

metoffice_resultFile = "/Users/jon/OPRAM/resultsUK_MetOffice/emDay_bt10_dd277/emDay_bt10_dd277_2023.tif"

opramUKtriangle_resultFile = "/Users/jon/OPRAM/resultsUK_triangle/agrilus_anxius/agrilus_anxius_01_Jan_2023_1km.csv"
opramUKsine_resultFile = "/Users/jon/OPRAM/resultsUK_sine/agrilus_anxius/agrilus_anxius_01_Jan_2023_1km.csv"
opramUKaverage_resultFile = "/Users/jon/OPRAM/resultsUK_average/agrilus_anxius/agrilus_anxius_01_Jan_2023_1km.csv"
opramIE_resultFile = "/Users/jon/OPRAM/resultsNI_triangle/agrilus_anxius/agrilus_anxius_01_Jan_2023_1km.csv"


IE_CRS = 29903  # Irish grid is 29903, UK is 27700
githubFolder = "~/git_repos/OPRAM/data"


# Import grid data
grid_UK = read.csv(file.path(githubFolder,"UK_grid_locations.csv"))
grid_IE = read.csv(file.path(githubFolder,"IE_grid_locations.csv"))

# ==================================================================
# Import UK Met Office result
metOffice = rast(metoffice_resultFile)




# ==================================================================
# Import OPRAM result



opram_UK_triangle = read.csv(opramUKtriangle_resultFile, 
                    header=TRUE,
                    colClasses = c("integer","integer","integer","integer", "Date","numeric","integer",
                                   "integer"))
opram_UK_sine = read.csv(opramUKsine_resultFile, 
                             header=TRUE,
                             colClasses = c("integer","integer","integer","integer", "Date","numeric","integer",
                                            "integer"))
opram_UK_average = read.csv(opramUKaverage_resultFile, 
                             header=TRUE,
                             colClasses = c("integer","integer","integer","integer", "Date","numeric","integer",
                                            "integer"))
opram_IE = read.csv(opramIE_resultFile, 
                            header=TRUE,
                            colClasses = c("integer","integer","integer","integer", "Date","numeric","integer",
                                           "integer"))

idxUK = match(opram_UK_triangle$ID, grid_UK$ID)
idxIE = match(opram_IE$ID, grid_IE$ID)


# ==================================================================
# Convert OPRAM IE result to a raster

# Put eastings and northings in the middle of the grad square
opram_IE$east = opram_IE$east + 500
opram_IE$north = opram_IE$north + 500

opram_sub = subset(opram_IE,select=c("east","north","emergeDOY"))
opram_rast = rast(opram_sub, 
              type="xyz", 
              crs=paste0("EPSG:29903"), 
              extent = ext(min(opram_sub$east)-500,
                           max(opram_sub$east)+500,
                           min(opram_sub$north)-500,
                           max(opram_sub$north)+500))

opram_UK = project(opram_rast, metOffice)


# Crop both opram_UK and metOffice to NI extent
NI_ext = ext(min(opram_UK_average$east),max(opram_UK_average$east)+1000,min(opram_UK_average$north),max(opram_UK_average$north)+1000)

opram_UK_crop = crop(opram_UK, NI_ext)
met_UK_crop = crop(metOffice, NI_ext)



# fileout = paste0(strsplit(opram_resultFile,".csv")[[1]],".tif")
# 
# terra::writeRaster(opram_rast, 
#                    filename = fileout,
#                    filetype="GTiff",
#                    datatype = "FLT4S",
#                    overwrite=TRUE,
#                    gdal=c("COMPRESS=DEFLATE", "TFW=NO"))


# Turn ratsers back into data frames
opram_UK_df = as.data.frame(opram_UK_crop, xy=TRUE)
met_UK_df = as.data.frame(met_UK_crop, xy=TRUE)


df_final = merge(opram_UK_df, 
            met_UK_df, 
            by=c("x","y"), 
            all=FALSE)
names(df_final) = c("east","north","opram_IE","metoffice")

df_final$east = df_final$east-500
df_final$north = df_final$north-500

df_final = merge(df_final, 
           opram_UK_average[,c("east","north","emergeDOY")], 
           by=c("east","north"), 
           all=FALSE)
names(df_final) = c("east","north","opram_IE","metoffice", "opram_UK_average")


df_final = merge(df_final, 
                 opram_UK_triangle[,c("east","north","emergeDOY")], 
                 by=c("east","north"), 
                 all=FALSE)
names(df_final) = c("east","north","opram_IE","metoffice", "opram_UK_average", "opram_UK_triangle")



df_final = merge(df_final, 
                 opram_UK_sine[,c("east","north","emergeDOY")], 
                 by=c("east","north"), 
                 all=FALSE)
names(df_final) = c("east","north","opram_IE","metoffice", "opram_UK_average", "opram_UK_triangle", "opram_UK_sine")




# Try merging rasters
# tmp2 = merge(opram_rast, metOffice2, algo=1)
# names(tmp2) = c("east","north","opram","metoffice")

ggplot(data=df_final,
       aes(x=east,
           y=north,
           colour=opram_UK_triangle-metoffice)) +
  geom_point() +
  scale_colour_distiller(palette = "PuOr", limits=c(-30,30)) +
  theme_bw()

# ggplot(data=df_final,
#        aes(x=east,
#            y=north,
#            colour=opram_IE-metoffice)) + 
#   geom_point() + 
#   scale_colour_distiller(palette = "PuOr", limits=c(-50,50)) +
#   theme_bw()
# 
# 
#   
# ggplot(data=df_final,
#        aes(x=east,
#            y=north,
#            colour=opram_UK_sine-metoffice)) + 
#   geom_point() + 
#   scale_colour_distiller(palette = "PuOr", limits=c(-30,30)) +
#   theme_bw()
# 
# ggplot(data=df_final,
#        aes(x=opram_UK_average-metoffice)) + 
#   geom_histogram() +
#   lims(x=c(-10,50)) + 
#   theme_bw()
# 
# ggplot(data=df_final,
#          aes(x=opram_UK_triangle-metoffice)) + 
#     geom_histogram() +
#     lims(x=c(-10,20)) + 
#   theme_bw()

# ggplot(data=df_final,
#        aes(x=opram_UK_sine - metoffice)) + 
#   geom_histogram() +
#   lims(x=c(-10,20)) + 
#   theme_bw()



library(tidyr)

df_long = pivot_longer(df_final, cols=-c(1:2), 
                       names_to="simulation", 
                       values_to="emerge_DOY")
df_long2 = pivot_longer(df_final, cols=-c(1:2,4), 
                       names_to="simulation", 
                       values_to="emerge_DOY")

ggplot(data=df_long,
       aes(x=simulation, 
           y=emerge_DOY)) +
  geom_boxplot() + 
  lims(y=c(150,300)) +
  theme_bw()

ggplot(data=df_long2, 
       aes(x=simulation,
           y=emerge_DOY-metoffice)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0, colour="darkred")+
  lims(y=c(-20,100)) +
  labs(x="OPRAM Simulation",
       y="Discrepancy (days)") +
  scale_x_discrete(labels=c("Met Eireann data\n (single triangle method)",
                            "UK Met Office data, \n (simple average method)",
                            "UK  Met Office data, \n(single sine method)",
                            "UK  Met Office data, \n(single triangle method)")) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=17),
        axis.text = element_text(size=17,
                                 angle=45,
                                 hjust=1))
ggsave("OPRAM_NI_discrepancy_boxplot.png")


# Convert OPRAM result (run on UK data) into raster
# ==================================================================


# Put eastings and northings in the middle of the grad square
opram_UK_triangle$east = opram_UK_triangle$east + 500
opram_UK_triangle$north = opram_UK_triangle$north + 500

opram_sub = subset(opram_UK_triangle,select=c("east","north","emergeDOY"))
opram_rast = rast(opram_sub, 
                  type="xyz", 
                  crs=crs(metOffice), 
                  extent = ext(min(opram_sub$east)-500,
                               max(opram_sub$east)+500,
                               min(opram_sub$north)-500,
                               max(opram_sub$north)+500))
tmp = crop(metOffice, opram_rast) 
opram_diff = opram_rast - tmp

UK = st_read(dsn = "~/Research/Data/GIS/UK/gadm41_GBR_1.shp")
IE = st_read(dsn = "~/Research/Data/GIS/country.shp")
UK_crop = st_crop(st_transform(UK, crs(opram_diff)), opram_diff)
IE_crop = st_crop(st_transform(IE, crs(opram_diff)), opram_diff)

ggplot() +
  geom_spatraster(data=opram_diff) + 
  geom_sf(data=UK_crop, fill=NA, linewidth=1, colour="black") +
  # geom_sf(data=IE_crop, fill=NA, linewidth=1, colour="black") +
  scale_fill_fermenter("Discrepancy \n(days)",
                       palette="RdYlBu", 
                       limits=c(-50,50), 
                       breaks=c(-50,-25,-5,5,25,50),
                       na.value=rgb(0,0,0,0)) + 
  theme_bw() + 
  theme(axis.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))

ggsave("OPRAM_NI_discrepancy_map.png")


# Produce similar map but with OPRAm results run on Met Eireann data
ggplot() +
  geom_spatraster(data=opram_UK_crop - met_UK_crop) + 
  geom_sf(data=UK_crop, fill=NA, linewidth=1, colour="black") +
  # geom_sf(data=IE_crop, fill=NA, linewidth=1, colour="black") +
  scale_fill_fermenter("Discrepancy \n(days)",
                       palette="RdYlBu", 
                       limits=c(-50,50), 
                       breaks=c(-50,-25,-5,5,25,50),
                       na.value=rgb(0,0,0,0)) + 
  theme_bw() + 
  theme(axis.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))

ggsave("OPRAM_NI_discrepancy_map2.png")



# Generate quantile
quantile(df_final$opram_UK_average-df_final$metoffice, c(0.025,0.25, 0.5, 0.75, 0.975))
quantile(df_final$opram_UK_sine-df_final$metoffice, c(0.025,0.25, 0.5, 0.75, 0.975))
quantile(df_final$opram_UK_triangle-df_final$metoffice, c(0.025,0.25, 0.5, 0.75, 0.975))
quantile(df_final$opram_IE-df_final$metoffice, c(0.025,0.25, 0.5, 0.75, 0.975))
