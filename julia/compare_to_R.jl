# Compare julia model results withe R calculations by paul
#
#
# Jon Yearsley#
# 17th Oct 2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


using Distributed;
using DataFrames
using Dates
using CSV
using Plots
using JLD2
using SharedArrays
using RData
import CodecBzip2

# =================================================================================
# Filename to import

paulFile = "Code that is under development//Changing Start Dates//Updated_Code__Removal_Of_Years//first_occurrence_threshold_with_date_range_1_1961_Agrilus_anxius_GDD_multiple_dates.rds"
paulFile2 = "WebDeveloper//Web Developer Quote - Shiny App//Model Outputs//Agrilus anxius//1st Occurrence - Select Individual Years//first_occurrence_1961_Agrilus_anxius.rds"

juliaFile = "result_agrilus1961_par_thin1.jld2"
year = 1961
thinFactor = 1



# =================================================================================
# Specify directories for import and export of data
if isdir("//home//jon//DATA//OPRAM")
        outDir = "//home//jon//DATA//OPRAM//results//"
        paulDir = "//home//jon//DATA//OPRAM//results//"
        meteoDir = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
        outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
        paulDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R"
        meteoDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"

elseif isdir("//users//jon//Desktop//OPRAM//")
        outDir = "//users//jon//Desktop//OPRAM//results//"
        paulDir = "//users//jon//Desktop//OPRAM//results//"
        meteoDir = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
end

# Include the functions to import the data and run the degree day model
include("GDD_functions.jl")
include("species_params.jl")

# ==============================================================================
# Import location data
latlongs = CSV.read(joinpath(meteoDir,"locations.CSV"), DataFrame,
            types=[Int32, Int64, Int64, Float64, Float64])

# Import meteo data
meteo = read_meteo(year, meteoDir, thinFactor)


# Import Paul data
# paul = CSV.read(normpath(joinpath(paulDir, paulFile)),DataFrame)
paul = load(normpath(joinpath(paulDir, paulFile)))

# This has starting day on 1st Jan
paul2 = load(normpath(joinpath(paulDir,"..", paulFile2)))

# Import julia data
julia = load(joinpath(outDir, "agrilus",juliaFile), "single_stored_object")



# Sort the latlongs dataframe
sort!(latlongs, [:ID])

# Remove data not aligned with a lat longs
idx = [searchsortedfirst(latlongs.ID, julia.ID[id]) for id in eachindex(julia.ID)]
keep = findall(x->x!=(nrow(latlongs)+1), idx)

julia.east .= 0
julia.north .= 0

julia.east[keep] .= latlongs.east[idx[keep]]
julia.north[keep] = latlongs.north[idx[keep]]


julia.ID2 = string.(julia.east) .* "_" .*string.(julia.north)
paul.ID2 = string.(paul.east) .* "_" .*string.(paul.north)
paul2.ID2 = string.(paul2.East) .* "_" .*string.(paul2.North)

# Add in starting day of year
paul.DOY = dayofyear.(paul.Start_Date)

# Subset julia data to contain locations in paul's data
# and then sort by ID2 and then start date
d_jul = julia[in.(julia.ID2,Ref(unique(paul.ID2))),:]
sort!(d_jul, [:ID2, :DOY])

# Subset paul2 data to contain locations in paul's data
# and then sort by ID2 and then start date
d_paul2 = paul2[in.(paul2.ID2,Ref(unique(paul.ID2))),:]
sort!(d_paul2, [:ID2])


# Subset meteo data and make it into a data frame
idx_meteo = in.(meteo[2], Ref(unique(d_jul.ID)))
unique(meteo[2][idx_meteo])

Tavg = meteo[1][:,idx_meteo]
ID = meteo[2][idx_meteo]
DOY = meteo[5]
d_meteo = DataFrame(ID=repeat(ID, inner=length(DOY)), 
                     DOY = repeat(DOY, sum(idx_meteo)),
                    Tavg = reshape(Tavg, prod(size(Tavg))))


# Pick specific columns from paul
d_paul = select(paul, :DOY, :Day, :ID2, :east, :north, :Date)
sort!(d_paul, [:ID2, :DOY])


# Check d and paul contain the same locations
unique(d_jul.ID2)
unique(d_paul.ID2)

idx = findall(diff(d_paul.Day).==1)
d_paul[idx,:]
d_jul


plot(d_paul.DOY, d_paul.Day, seriestype=:scatter, ms=1, legend=false);
xlims!(0,366);
ylims!(0,366)
plot(d_jul.DOY, d_jul.emergeDOY, seriestype=:scatter, ms=1, legend=false);
xlims!(0,366);
ylims!(0,366)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Check dates of first emergence
# Pull out first emergence date

# Julia data
idx = [searchsortedfirst(d_jul.ID2, id2)  for id2 in unique(d_paul.ID2)];
d_jul[idx,:]

# Paul's changing dstart date data
idx = [searchsortedfirst(d_paul.ID2, id2)  for id2 in unique(d_paul.ID2)];
d_paul[idx,:]

d_paul2


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Check DD calculation for one location 

ID = 68001    # The location ID
d_sub = d_meteo[in.(d_meteo.ID, ID),:]

d_sub.dT = d_sub.Tavg .- agrilus_anxius.base_temperature
d_sub.dT[d_sub.dT.<0] .= 0

d_sub.cumT = cumsum(d_sub.dT)



d_sub


d_sub[d_sub.cumT.>agrilus_anxius.threshold, :]