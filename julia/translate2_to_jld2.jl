# Script to import raw ensemble level TRANSLATE project meteo files and 
# convert them  to jld2 format for more efficient import
# 
# The TRANSLATE files for maximum and minimum daily temperatures are imported.
# These files give 30yr of simulations  for an ensemble member 
# The data are interpolated onto the 1km grid of Ireland and 
# saved in a jld2 file for each RCP and time period.
#
# Interpolation uses nearest-neighbour model with 10 neighbours. 
# This avoids extrapolation errors, but may need to be revised.
#
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 15th July 2026
#
# ====================================================================



using CSV;
using JLD2;
using DataFrames;
using GeoStats;
using NCDatasets;
using Distributed;

include("OPRAM_io_functions.jl")  # Includes functions to import meteo data


rcpList = ["26", "45", "85"];              # The RCP scenario (26 = RCP2.6, 45=RCP4.5, 85=RCP8.5)
periodList = ["2021-2050", "2041-2070"];   # The future scenario (2021-2050 or 2041-2070)
gridFile = "IE_grid_locations.csv";        # File containing a 1km grid of lats and longs over Ireland 

rcpList = ["85"];              # The RCP scenario (26 = RCP2.6, 45=RCP4.5, 85=RCP8.5)
periodList = ["2071-2100"];   # The future scenario (2021-2050 or 2041-2070)

nDays = 365;  # Number of days in the data (assumes no leap years)

# =========================================================
# =========================================================

# Specify directories for import and export of data
if isdir("//home//jon//DATA//OPRAM")
    outDir = "//home//jon//DATA//OPRAM//Climate_JLD2//"
    dataDir = "//home//jon//DATA//OPRAM//"
    translateDir = "//home//jon//DATA//OPRAM//TRANSLATE"

elseif isdir("/Users/jon/Research/Data/Climate/Translate/ensemble/")
    outDir = "/Users/jon/Research/OPRAM/Climate_JLD2"
    dataDir = "/Users/jon/Library/CloudStorage/GoogleDrive-jon.yearsley@ucd.ie/My Drive//Projects//DAFM_OPRAM//Data//"
    translateDir = "/Users/jon/Research/Data/Climate/Translate/ensemble"

elseif isdir("//users//jon//Desktop//OPRAM//")
    outDir = "//users//jon//Desktop//OPRAM//Climate_JLD2//"
    dataDir = "//users//jon//Desktop//OPRAM//"
    translateDir = "//users//jon//Desktop//OPRAM//TRANSLATE"
end


# Make directory for the output if one doesn't exist
if !isdir(outDir)
    @info "Making directory for output"
    mkpath(outDir)
end

# =========================================================
# =========================================================
# Import grid of location data 

# Import grid location data
grid = CSV.read(joinpath([dataDir, gridFile]), DataFrame);

# Sort locations in order of IDs
grid = grid[sortperm(grid.ID), :];

# Create a geotable of the grid data
grid_gt = georef((x=grid.east, y=grid.north, ID=grid.ID), ("x", "y"), crs=EPSG{29903})


# =========================================================
# =========================================================
# Check size of data to be imported (using Tmax) and 

# Find correct directory and file
reg = Regex("^tmax_" * "[-\\w]+" * "_rcp" * string(rcpList[1]) * "_" * string(periodList[1]) * "[-\\w]+" * ".nc")
files = filter(x -> occursin(reg, x), readdir(translateDir))

# Extract the name of the ensemble member from the file name
ensemble_member = [match(Regex("tmax_([-\\w]+).nc"), file).captures[1] for file in files]

# Import one file of the TRANSLATE data 
ds = NCDataset(joinpath(translateDir, files[1]), "r");
Tmax1 = ds["tmax"];
lat = ds["lat"];
lon = ds["lon"];


# Indices to remove data that has unrealsitic max temps
nonzero_idx = findall(dropdims(.!all(ismissing.(Tmax1[:, :, 1:3]), dims=3), dims=3));

# Calculate some dimensions
TSize = size(Tmax1)
nYears = floor(Int, TSize[3] / nDays);

# Reshape lat and lon arrays into 2D arrays 
lon_2D = [lon[i] for i in eachindex(lon), j in eachindex(lat)];
lat_2D = [lat[j] for i in eachindex(lon), j in eachindex(lat)];


# Remove locations with missing values
lon = lon_2D[nonzero_idx]
lat = lat_2D[nonzero_idx]


# =========================================================
# =========================================================
# For each 1km grid square find nearest neighbour TRANSLATE location

# Convert TRANSLATE lat lon into eastings and northings for the full grid of TRANSLATE data
translate_xy = georef((lat=lat, lon=lon, ID=collect(1:length(nonzero_idx))),
    ("lat", "lon"),
    crs=EPSG{4326}) |> Proj(EPSG{29903})  # Convert to Irish National Grid TM75

# Define interpolation model
interpolation_model = InterpolateNeighbors(domain(grid_gt),
    model=NN(Euclidean()),
    maxneighbors=10)

# Find nearest neighbour IDs for the TRANSLATE data
translate_NN = translate_xy |> interpolation_model
grid.translateID = translate_NN.ID

# Create array for days of year
DOY = collect(1:nDays);

# Clear memory
translate_NN =  nothing;
translate_xy = nothing;

# =========================================================
# =========================================================
# Import the data for each RCP and period, interpolate data onto 1km grid 
# and save each individual year in a JLD2 file

for r in eachindex(rcpList)
    @time "Converted TRANSLATE data:" for p in eachindex(periodList)
        for ens in eachindex(ensemble_member)
            # +++++++++++++++++++++++++++++++++++++++++++++++++++
            # Import the TRANSLATE data
            @info "Importing data for ensemble member " * string(ens) * ": " * ensemble_member[ens]

            # Extract data from the netCDF and put nearest TRANSLATE data point onto the 1km grid

            # Tmax data first
            @info "Importing Tmax data for ensemble member " * string(ens)
            # Find tmax files for a specific RCP, time period and ensemble member
            files_max = filter(x -> occursin(Regex("^tmax_" * ensemble_member[ens]), x), readdir(translateDir))
            ds = NCDataset(joinpath(translateDir, files_max[1]), "r");
            Tmax = convert.(Float32, ds["tmax"][nonzero_idx[grid.translateID], :]);
            close(ds)

            # Write data to files (one file per year)
            for year in 1:nYears
                @info "Saving data for year " * string(year)

                # Save this to a JLD2 file
                outfile = joinpath([outDir, "TRANSLATE_" * ensemble_member[ens] * "_year" * string(year) * ".jld2"])
                fh = jldopen(outfile, "w");
                # Save Tmax... transpose so that spatial locations are in columns and days are in rows
                fh["Tmax"]=transpose(Tmax[:, (year-1)*nDays+1:year*nDays]);
                fh["DOY"]=DOY;
                fh["ID"]=grid.ID;
                close(fh)
            end
            Tmax = nothing  # Clear memory

            # Now Tmin data
            @info "Importing Tmin data for ensemble member " * string(ens)
            # Find tmin files for a specific RCP, time period and ensemble member
            files_min = filter(x -> occursin(Regex("^tmin_" * ensemble_member[ens]), x), readdir(translateDir))

            ds = NCDataset(joinpath(translateDir, files_min[1]), "r");
            Tmin = convert.(Float32, ds["tmin"][nonzero_idx[grid.translateID], :]);
            close(ds)

            # Write data to files (one file per year)
            for year in 1:nYears
                @info "Appending Tmin data for year " * string(year)

                # Append Tmin to a JLD2 file
                outfile = joinpath([outDir, "TRANSLATE_" * ensemble_member[ens] * "_year" * string(year) * ".jld2"])
                fh = jldopen(outfile, "a");
                # Save Tmin... transpose so that spatial locations are in columns and days are in rows
                fh["Tmin"]=transpose(Tmin[:, (year-1)*nDays+1:year*nDays]);
                close(fh)
            end
            Tmin = nothing  # Clear memory
        end
   end
end





# # PLot data to test

# using Plots
# plot(grid.east,
#     grid.north,
#     zcolor=tmin[:, 200], 
#     st=:scatter, 
#     ms=1,
#     markerstrokewidth=0,
#     c=:thermal,
#     clims=(0, 30),
#         grid=false,
#     legend=false,
#     cbar=true,
#     aspect_ratio=:equal,
#     xlabel="Eastings",
#     ylabel="Northings",
#     title="Max Temp Day 1")



# # Try loading a JLD2 file to check it has been saved correctly
# tmp = load(joinpath([outDir, "TRANSLATE_" * ensemble_member[1] * "_year" * string(1) * ".jld2"]));

# tmin = load(joinpath([outDir, "TRANSLATE_" * ensemble_member[1] * "_year" * string(1) * ".jld2"]), "Tmin");