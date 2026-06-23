# Import individual ensemble TRANSLATE data and 
# visualise variability in tmax
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# June 2026
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# using CSV;
using JLD2;

using DataFrames;
using NetCDF;
using Plots
using Statistics;
using CSV;


translateDir_max = "/Volumes/CLIMATEDATA/Translate/Ensemble_member_sample_30_years_daily_data";
varStr = "tmax"

# =========================================================
# =========================================================



# =========================================================
# =========================================================
# Check size of data to be imported (using Tmax) and 

# Find correct directory and file
files = filter(x -> occursin(Regex("^" * varStr) , x), readdir(translateDir_max))

# Import one file of the TRANSLATE data 
Tmax = ncread(joinpath(translateDir_max, files[1]), "tmax");
lat = ncread(joinpath(translateDir_max, files[1]), "lat");
lon = ncread(joinpath(translateDir_max, files[1]), "lon");
doy = ncread(joinpath(translateDir_max, files[1]), "time");

# Calculate some dimensions
TSize = size(Tmax)
nDays = TSize[3];


# Reshape arrays in place
Tmax_long = reshape(Tmax,:,nDays);
lon_2D = reshape([lon[i] for i in eachindex(lon), j in eachindex(lat)],:);
lat_2D = reshape([lat[j] for i in eachindex(lon), j in eachindex(lat)], :);

Tmax = nothing
lat = nothing
lon = nothing

# Find data that has  max temps more than -273 (i.e. not missing data) 
# and keep only those lats and longs
nonzero_idx = dropdims(all(Tmax_long.>-273, dims=2), dims=2);


# Keep only lats and longs with non-zero temps
lon_2D = lon_2D[nonzero_idx];
lat_2D = lat_2D[nonzero_idx];


Tmax_long2 = reshape(Tmax_long[nonzero_idx, :], :, 1);

Tmax_long = nothing
nLocation = length(lon_2D)

size(Tmax_long2)
size(lon_2D)
size(lat_2D)

tmax_df = DataFrame(
    # lat=repeat(lat_2D,outer=nDays), 
    #                 lon=repeat(lon_2D, outer=nDays), 
                    doy = mod1.(repeat(doy, inner=nLocation),365),
                    tmax=Tmax_long2[:])

jldsave("translate_ensemble_tmax.jld2"; tmax_df=tmax_df)

lat_2D=nothing
lon_2D=nothing
Tmax_long2=nothing

 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Group the data by doy
tmax_group = groupby(tmax_df, :doy)
tmax_df = nothing
tmax_agg = combine(tmax_group,
            :tmax => (x ->  quantile(skipmissing(x), 0.5)) => :tmax_50,
            :tmax => (x ->  quantile(skipmissing(x), 0.9)) => :tmax_90,
            :tmax => (x ->  quantile(skipmissing(x), 0.1)) => :tmax_10,
            )


# Write results to CSV file
CSV.write("translate_ensemble_tmax_aggregated.csv", tmax_agg)

ind = 1

x=mod1.(doy, 365)
y = Tmax_2D[ind,:]
ms=0.5


plot(x[1:365],
    y[1:365],
    seriestype=:line,
    markersize=ms,
    markerstrokewidth=0,
    showaxis=true,
    grid=false,
    legend=false,
    cbar=false)


x=lon_2D[:,ind]
y=lat_2D[:, ind]
z = Tmax_2D[:, ind]

    plot(x,
    y,
    zcolor=z,
    seriestype=:scatter,
    markersize=ms,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal)