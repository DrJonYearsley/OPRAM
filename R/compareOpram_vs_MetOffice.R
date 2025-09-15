# Compare the UK Met Office emergence date with the OPRAM result

rm(list=ls())
library(sf)
library(terra)
library(ggplot2)

metoffice_resultFile = "/Users/jon/OPRAM/resultsUK_MetOffice/emDay_bt10_dd277/emDay_bt10_dd277_2023.tif"

opramUK_resultFile = "/Users/jon/OPRAM/resultsUK_triangle/agrilus_anxius/agrilus_anxius_01_Jan_2023_1km.csv"
opramUK_resultFile = "/Users/jon/OPRAM/resultsUK_sine/agrilus_anxius/agrilus_anxius_01_Jan_2023_1km.csv"
opramIE_resultFile = "/Users/jon/OPRAM/resultsNI_triangle/agrilus_anxius/agrilus_anxius_01_Jan_2023_1km.csv"


IE_CRS = 27700  # Irish grid is 29903, UK is 27700
githubFolder = "~/git_repos/OPRAM/data"


# Import grid data
grid_UK = read.csv(file.path(githubFolder,"UK_grid_locations.csv"))

# ==================================================================
# Import UK Met Office result
metOffice = rast(metoffice_resultFile)
metOffice2 = project(metOffice, "epsg:27700")


met_crds = crds(metOffice)



# ==================================================================
# Import OPRAM result
opram_UK = read.csv(opramUK_resultFile, 
             header=TRUE,
             colClasses = c("integer","integer","integer","integer", "Date","numeric","integer",
                             "integer"))
idx = match(opram$ID, grid_UK$ID)



opram_UK = read.csv(opramUK_resultFile, 
                    header=TRUE,
                    colClasses = c("integer","integer","integer","integer", "Date","numeric","integer",
                                   "integer"))


# ==================================================================
# Convert OPRAM result to a raster

# Put eastings and northings in the middle of the grad square
opram$east = opram$east + 500
opram$north = opram$north + 500

opram_sub = subset(opram,select=c("east","north","emergeDOY"))
opram_rast = rast(opram_sub, 
              type="xyz", 
              crs=paste0("EPSG:",IE_CRS), 
              extent = ext(min(opram_sub$east)-500,
                           max(opram_sub$east)+500,
                           min(opram_sub$north)-500,
                           max(opram_sub$north)+500))

fileout = paste0(strsplit(opram_resultFile,".csv")[[1]],".tif")

terra::writeRaster(opram_rast, 
                   filename = fileout,
                   filetype="GTiff",
                   datatype = "FLT4S",
                   overwrite=TRUE,
                   gdal=c("COMPRESS=DEFLATE", "TFW=NO"))

opram_crds = crds(opram_rast)


image(is.na(opram_rast))


head(opram_sub)
head(metOffice)


tmp = as.data.frame(metOffice)
met_df = cbind(met_crds, tmp)


tmp = merge(opram_sub, 
            met_df, 
            by.x=c("east","north"), 
            by.y=c("x","y"), 
            all=FALSE)
names(tmp) = c("east","north","opram","metoffice")


# Try merging rasters
# tmp2 = merge(opram_rast, metOffice2, algo=1)
# names(tmp2) = c("east","north","opram","metoffice")

ggplot(data=tmp,
       aes(x=east,
           y=north,
           colour=opram-metoffice)) + 
  geom_point() + 
  scale_colour_distiller(palette = "PuOr", limits=c(-10,100)) +
  theme_bw()

  
  
ggplot(data=tmp,
         aes(x=opram-metoffice)) + 
    geom_histogram() +
    lims(x=c(-10,100)) + 
  theme_bw()
  