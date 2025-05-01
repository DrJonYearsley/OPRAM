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
# using SharedArrays
using RData
import CodecBzip2

# =================================================================================
# Filename to import

run_params = (speciesName = "decem",
  year = 1992,
  thinFactor = 1)

# Files containing the ID system from Granite.ie
# This info is used to package the output into separate files
granite_hectad_ID = "granite_hectad_defs.csv"
granite_county_ID = "granite_county_defs.csv"
granite_hectad_county = "granite_hectad_county_defs.csv"


# Put all the paths together
paths = (outDir=joinpath(homedir(),"Google Drive/My Drive/Projects/DAFM_OPRAM/results"),
          dataDir=joinpath(homedir(),"git_repos/OPRAM/data"))


# Find directory matching the species name in run_params
speciesName_jl = filter(x -> occursin(r"" * run_params.speciesName, x), readdir(joinpath(paths.outDir,"granite_output")))
if length(speciesName_jl) > 1
    @error "More than one species name found"
elseif length(speciesName_jl) == 0
    @error "Species not found"
end

speciesName_R = filter(x -> occursin(r"" * run_params.speciesName, x), readdir(joinpath(paths.outDir,"Changing_Start_Dates")))
if length(speciesName_R) > 1
    @error "More than one species name found"
elseif length(speciesName_R) == 0
    @error "Species not found"
end





# =================================================================================
# Specify directories for import and export of data

# if isdir("//home//jon//Desktop//OPRAM")
#         outDir = "//home//jon//Desktop//OPRAM//results//"
#         dataDir = "//home//jon//DATA//OPRAM//"
#                 paulDir = "//home//jon//DATA//OPRAM//results//"
#         meteoDir_IE = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"
#         meteoDir_NI = "//home//jon//DATA//OPRAM//Northern_Ireland_Climate_Data//"
      
#       elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
#         outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
#         dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
#         paulDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R"        
#         meteoDir_IE = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"
#         meteoDir_NI = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Northern_Ireland_Climate_Data//"
      
#       elseif isdir("//users//jon//Desktop//OPRAM//")
#         outDir = "//users//jon//Desktop//OPRAM//results//"
#         dataDir = "//users//jon//Desktop//OPRAM//"
#         paulDir = "//users//jon//Desktop//OPRAM//results//"        
#         meteoDir_IE = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
#         meteoDir_NI = "//users//jon//Desktop//OPRAM//Northern_Ireland_Climate_Data//"
#       end

# Include the functions to import the data and run the degree day model
# include("GDD_functions.jl")
include("OPRAM_io_functions.jl")

# ==============================================================================
# Import location data
grid = read_grid(joinpath(paths.dataDir,"IE_grid_locations.csv"), run_params.thinFactor)

# Import granite's county ID's
county_defs = CSV.read(joinpath([paths.dataDir, granite_county_ID]), DataFrame);


# # Import meteo data
# meteo = read_meteo(year, meteoDir_IE, meteoDir_NI, grid_thin)


# Import R model data
RFile1 = joinpath(paths.outDir,"Changing_Start_Dates", speciesName_R[1],
          "first_occurrence_threshold_with_date_range_1_" * string(run_params.year) * 
          "_" * speciesName_R[1] * ".rds")

RFile2 = joinpath(paths.outDir,"Changing_Start_Dates", speciesName_R[1],
          "number_of_thresholds_" * string(run_params.year) * 
          "_" * speciesName_R[1] * ".rds")


# paul = CSV.read(normpath(joinpath(paulDir, paulFile)),DataFrame)
R1_out = load(normpath(RFile1))
R2_out = load(normpath(RFile2))

R_out = leftjoin(select(R1_out,[:Start_Date, :End_Date, :north, :east,:Day, :Threshold_Exceeded_Count]), 
                select(R2_out,[:Start_Date, :End_Date, :north, :east, :Day, :thresholds_exceeded ]),
                 on = [:Start_Date, :End_Date, :north, :east], makeunique=true, renamecols = "_R1" => "_R2")


                #  R_out=R2_out

# Import julia data
# Find county corresponding to R_out data
locations = unique(R_out[:,[:east, :north]])

idx = grid.east.==locations.east  .&& grid.north .== locations.north;

county = unique(grid.county[idx])[1]
countyID = county_defs.Id[county_defs.County.==county][1]

juliaFile = joinpath(paths.outDir, "granite_output",speciesName_jl[1], 
                    speciesName_jl[1] * "_1_" * string(countyID) * 
                    "_" * string(run_params.year)  * ".csv")
jl_all = CSV.read(juliaFile, DataFrame)


jl_out = jl_all[jl_all.ID .== grid.ID[idx],:]


jl_out[:,[2,3,4]]
R_out[:,[1,5,6,8]]

# # Add in starting day of year
# paul.DOY = dayofyear.(paul.Start_Date)

# # # Subset julia data to contain locations in paul's data
# # # and then sort by ID2 and then start date
# # d_jul = julia[in.(julia.ID2,Ref(unique(paul.ID2))),:]
# # sort!(d_jul, [:ID2, :DOY])

# # # Subset paul2 data to contain locations in paul's data
# # # and then sort by ID2 and then start date
# # d_paul2 = paul2[in.(paul2.ID2,Ref(unique(paul.ID2))),:]
# # sort!(d_paul2, [:ID2])


# # Subset meteo data and make it into a data frame
# idx_meteo = in.(meteo[3], Tuple(unique(d_jul.ID)))
# unique(meteo[3][idx_meteo])

# Tavg = meteo[1][:,idx_meteo]
# ID = meteo[3][idx_meteo]
# DOY = meteo[2]
# d_meteo = DataFrame(ID=repeat(ID, inner=length(DOY)), 
#                      DOY = repeat(DOY, sum(idx_meteo)),
#                     Tavg = reshape(Tavg, prod(size(Tavg))))


# # Pick specific columns from paul
# d_paul = select(paul, :DOY, :Day, :ID2, :east, :north, :Date)
# sort!(d_paul, [:ID2, :DOY])


# # Check d_jul and paul contain the same locations
# unique(d_jul.ID2)
# unique(d_paul.ID2)

# # Compare all the d_jul results to d_paul results
# idx = findall(diff(d_paul.Day).>0)
# d_paul[idx,:]
# d_jul


# # Plot Paul's results alongside the julia output
# plot(d_paul.DOY, d_paul.Day, seriestype=:scatter, ms=1, legend=false);
# xlims!(0,366);
# # ylims!(0,366)
# plot!(d_jul.DOY, d_jul.emergeDOY, seriestype=:scatter, markercolor=:red, ms=2, legend=false)
# # xlims!(0,366);
# # ylims!(0,366)

# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # Check dates of first emergence
# # Pull out first emergence date

# # Julia data
# idx = [searchsortedfirst(d_jul.ID2, id2)  for id2 in unique(d_paul.ID2)];
# d_jul[idx,:]

# # Paul's changing dstart date data
# idx = [searchsortedfirst(d_paul.ID2, id2)  for id2 in unique(d_paul.ID2)];
# d_paul[idx,:]

# d_paul2
