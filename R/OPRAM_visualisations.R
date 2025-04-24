# An example of visualisations for the OPRAM web app

library(ggplot2)
library(sf)
library(lubridate)
rm(list=ls())
species_name = "ips_sexdentatus"
year = 2010
county = 6  # 6=cork
dataFolder = "~/Google Drive/My Drive/Projects/DAFM_OPRAM/results/granite_output"
githubFolder = "~/git_repos/OPRAM/data"
start = paste0(year,"-01-01")



# ===================================

# Import data
filein = paste0(species_name,"_100_",year,".csv")
d = read.csv(file.path(dataFolder,species_name,filein))

filein_1km = paste0(species_name,"_1_",county,"_",year,".csv")
d1km = read.csv(file.path(dataFolder,species_name,filein_1km))



granite = read.csv(file.path(githubFolder,"granite_hectad_county_defs.csv"))
grid = read.csv(file.path(githubFolder,"IE_grid_locations.csv"))
grid_hectad = read.csv(file.path(githubFolder,"IE_hectad_locations.csv"))


d_start01 = subset(d, startDate==as.Date(start))
d1km_start01 = subset(d1km, startDate==as.Date(start))

# Add in coords
d_start01$east = grid_hectad$east[match(d_start01$hectad, grid_hectad$hectad)]
d_start01$north = grid_hectad$north[match(d_start01$hectad, grid_hectad$hectad)]


d1km_start01$east = grid$east[match(d1km_start01$ID, grid$ID)]
d1km_start01$north = grid$north[match(d1km_start01$ID, grid$ID)]


################################################################################
# All Ireland


# Make GIS object for data
sf_start01 = st_as_sf(d_start01, coords = c("east","north"), crs=29903)
sf_start01

br = c(1, 60, 90,
       105, 120, 130, 140, 151,
       161, 171, 181, 191, 201, 212,
       222, 232, 243, 258, 273, 365,
       731, 1096, 2000)

cols = c("#a1d99b", "#238b45",
  "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c",
  "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506",
  "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58",
  "#ef00ef", "#ffffbf","#B2B2B2")

labs_emerge = c("Jan-Feb",
                format(as.Date(start) + br[c(2)],"%b"),
                format(as.Date(start) + br[-c(1,2,21:23)],"%d %b"),
                format(as.Date(start) + br[c(21:22)],"%Y"),
                "No emergence")

# cols = c(  "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506")
# br = c(1,30,60,90, 120,250,350)

ggplot() +
  geom_sf(data=sf_start01, aes(colour=emergeDOY), size=3) + 
  scale_color_stepsn("Date of\nEmergence",
    colours = cols,
    breaks = br,
    values = (br-1)/(2000-1),
    limits=c(1,2000),
    # labels = ~ format(as.Date(start, format="%Y-%m-%d")+pmax(., 1,na.rm=TRUE), format = "%d %b %Y")
    labels = labs_emerge,
    right=FALSE
  ) +
  labs(title=species_name) +
  theme_void() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches")) +
  guides(colour = guide_coloursteps(show.limits = TRUE))



# Bin width is 7 days
ggplot(data=d_start01,
       aes(x=emergeDOY,
           fill=emergeDOY,
           group=emergeDOY)) + 
  geom_histogram(boundary=1,binwidth=7) + 
  labs(x="Date of emergence",
       y = "Number of 10km squares",
       title = species_name) + 
  scale_x_continuous(
    labels = ~ format(as.Date(start, format="%Y-%m-%d")+pmax(., 1,na.rm=TRUE), format = "%d %b %Y")
  ) +
  scale_fill_stepsn("Date of\nEmergence",
                    aesthetics="fill",
                          colours = cols,
                          breaks = br,
                          values = (br-1)/(2000-1),
                          limits=c(1,2000),
                          labels = ~ format(as.Date(start, format="%Y-%m-%d")+pmax(., 1,na.rm=TRUE), format = "%d %b %Y")
  ) +
  theme_bw() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        axis.text = element_text(size=14)) +
  guides(colour = guide_coloursteps(show.limits = FALSE))
  




# Number of generations -------


cols_nGen = c("#7f3b08","#b35806","#e08214",  "#fdb863","#fee0b6","#f7f7f7","#bcbddc","#9e9ac8",
"#807dba","#6a51a3","#54278f","#54278f","#3f007d")

br_nGen = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5,3, 3.5, 4,10)

lab_nGen = as.character(br_nGen)
lab_nGen[length(lab_nGen)] = ""
lab_nGen[length(lab_nGen)-1] = "4+"

ggplot() +
  geom_sf(data=sf_start01, aes(colour=nGen), size=3) + 
  scale_color_stepsn("Number of\nGenerations",
                     colours = cols_nGen,
                     breaks = br_nGen,
                     values = (br_nGen-0)/(10-0),
                     limits=c(0,10),
                     labels=lab_nGen
  ) +
  labs(title=paste(species_name,"starting date:",start)) +
  theme_void() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches")) +
  guides(colour = guide_coloursteps(show.limits = FALSE))




ggplot(data=d_start01,
       aes(x=nGen,
           fill=nGen,
           group=nGen)) + 
  geom_histogram(boundary=0, binwidth=0.1) + 
  labs(x="Number of Generations",
       y = "Number of 10km squares across Ireland",
       title=paste(species_name,"starting date:",start)) + 
  scale_fill_stepsn("Number of\nGenerations",
                    aesthetics="fill",
                    colours = cols_nGen,
                    breaks = br_nGen,
                    values = (br_nGen-0)/(10-0),
                    limits=c(0,10),
                    labels=lab_nGen
  ) +
  theme_bw() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        axis.text = element_text(size=14)) +
  guides(colour = guide_coloursteps(show.limits = FALSE))



#################################################################
# Plots for Cork



# Make GIS object for data
sf_1km = st_as_sf(d1km_start01, coords = c("east","north"), crs=29903)
sf_1km



ggplot() +
  geom_sf(data=sf_1km, aes(colour=emergeDOY), size=0.8) + 
  scale_color_stepsn("Date of\nEmergence",
                     colours = cols,
                     breaks = br,
                     values = (br-1)/(2000-1),
                     limits=c(1,2000),
                     labels = ~ format(as.Date(start, format="%Y-%m-%d")+pmax(., 1,na.rm=TRUE), format = "%d %b %Y")
  ) +
  labs(title=species_name) +
  theme_void() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches")) +
  guides(colour = guide_coloursteps(show.limits = FALSE))






# Bin width is 7 days
ggplot(data=d1km_start01,
       aes(x=emergeDOY,
           fill=emergeDOY,
           group=emergeDOY)) + 
  geom_histogram(boundary=1,binwidth=7) + 
  labs(x="Date of emergence",
       y = "Number of 1km squares",
       title = species_name) + 
  scale_x_continuous(
    labels = ~ format(as.Date(start, format="%Y-%m-%d")+pmax(., 1,na.rm=TRUE), format = "%d %b %Y")
  ) +
  scale_fill_stepsn("Date of\nEmergence",
                    aesthetics="fill",
                    colours = cols,
                    breaks = br,
                    values = (br-1)/(2000-1),
                    limits=c(1,2000),
                    labels = ~ format(as.Date(start, format="%Y-%m-%d")+pmax(., 1,na.rm=TRUE), format = "%d %b %Y")
  ) +
  theme_bw() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        axis.text = element_text(size=14)) +
  guides(colour = guide_coloursteps(show.limits = FALSE))





# Number of generations -------


cols_nGen = c("#7f3b08","#b35806","#e08214",  "#fdb863","#fee0b6","#f7f7f7","#bcbddc","#9e9ac8",
              "#807dba","#6a51a3","#54278f","#54278f","#3f007d")

br_nGen = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5,3, 3.5, 4,10)

lab_nGen = as.character(br_nGen)
lab_nGen[length(lab_nGen)] = ""
lab_nGen[length(lab_nGen)-1] = "4+"

ggplot() +
  geom_sf(data=sf_1km, aes(colour=nGen), size=0.8) + 
  scale_color_stepsn("Number of\nGenerations",
                     colours = cols_nGen,
                     breaks = br_nGen,
                     values = (br_nGen-0)/(10-0),
                     limits=c(0,10),
                     labels=lab_nGen
  ) +
  labs(title=paste(species_name,"starting date:",start)) +
  theme_void() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches")) +
  guides(colour = guide_coloursteps(show.limits = FALSE))




ggplot(data=d1km_start01,
       aes(x=nGen,
           fill=nGen,
           group=nGen)) + 
  geom_histogram(boundary=0, binwidth=0.1) + 
  labs(x="Number of Generations",
       y = "Number of 1km squares across the county",
       title=paste(species_name,"starting date:",start)) + 
  scale_fill_stepsn("Number of\nGenerations",
                    aesthetics="fill",
                    colours = cols_nGen,
                    breaks = br_nGen,
                    values = (br_nGen-0)/(10-0),
                    limits=c(0,10),
                    labels=lab_nGen
  ) +
  theme_bw() + 
  theme(legend.key.height = unit(dev.size()[2]/10 , "inches"),
        panel.grid = element_blank(),
        axis.title = element_text(size=20),
        axis.text = element_text(size=14)) +
  guides(colour = guide_coloursteps(show.limits = FALSE))


