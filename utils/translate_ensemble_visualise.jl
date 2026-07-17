# Import individual ensemble TRANSLATE data and 
# visualise variability in tmax
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# June 2026
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# using CSV;
using JLD2;
using DataFrames;
using NCDatasets;
using Plots
using Statistics;
using CSV;

if (isdir("/Volumes/CLIMATEDATA"))
    translateDir_max = "/Volumes/CLIMATEDATA/Translate/Ensemble_member_sample_30_years_daily_data";
else
    translateDir_max = "/Volumes/CLIMATEDATA/Translate/Ensemble_member_sample_30_years_daily_data";
end
varStr = "tmax"

# =========================================================
# =========================================================



# =========================================================
# =========================================================
# Check size of data to be imported (using Tmax) and 
if (isfile("translate_ensemble_tmax.jld2"))
    println("File translate_ensemble_tmax.jld2 already exists, skipping import of Tmax data")

    tmax_df, lat, lon, doy = load("translate_ensemble_tmax.jld2", "tmax_df", "lat", "lon", "doy")
else
    println("Importing Tmax data from TRANSLATE ensemble")
    # Find correct directory and file
    files = filter(x -> occursin(Regex("^" * varStr), x), readdir(translateDir_max))

    # Import one file of the TRANSLATE data 
    Tmax = ncread(joinpath(translateDir_max, files[1]), "tmax");
    lat = ncread(joinpath(translateDir_max, files[1]), "lat");
    lon = ncread(joinpath(translateDir_max, files[1]), "lon");
    doy = ncread(joinpath(translateDir_max, files[1]), "time");

    # Calculate some dimensions
    TSize = size(Tmax)
    nDays = TSize[3];


    # Find data that has  max temps more than -273 (i.e. not missing data) for 
    # every time point 
    nonzero_idx = findall(view(Tmax,:,:,1) .> -273)

    # Create data frame for every location with non-zero data
    tmax_df = reduce(vcat, [DataFrame(ID=Int32(x), tmax=Tmax[nonzero_idx[x], :]) for x in eachindex(nonzero_idx)]);
    # Keep only lats and longs with non-zero temps

    # Clear memory
    Tmax = nothing


    jldsave("translate_ensemble_tmax.jld2"; tmax_df=tmax_df, lat=lat, lon=lon, doy=doy)
end


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Filter one day of year and calculate Q90 -Q10 for tmax from every spatial location 
# and then calculate quantiles of that for every day of year

# Create array to hold the final results d
delta_tmax = Array{Float64}(undef, (365, 4))

for d in 1:365
    @info "Processing day $d of 365"

    # Filter the data frame for the d day of year
    tmp = filter(df->mod1(rownumber(df), 365)==d, tmax_df)
    # Group by spatial location
    tmax_group = groupby(tmp, :ID)
    tmp = nothing

    # Calculate the median and 90-10 range for each spatial location
    tmax_agg = combine(tmax_group,
        :tmax => (x -> quantile(skipmissing(x), 0.5)) => :tmax_50,
        :tmax => (x -> quantile(skipmissing(x), 0.9) - quantile(skipmissing(x), 0.1)) => :tmax_90_10,
    )

    # Calculate the median and quantiles of the 90-10 range for all spatial locations
    delta_tmax[d, :] = Matrix(combine(tmax_agg,
        :tmax_50 => (x -> quantile(skipmissing(x), 0.5)) => :tmax_50,
        :tmax_90_10 => (x -> quantile(skipmissing(x), 0.1)) => :diff_10,
        :tmax_90_10 => (x -> quantile(skipmissing(x), 0.5)) => :diff_50,
        :tmax_90_10 => (x -> quantile(skipmissing(x), 0.9)) => :diff_90,
    ))
end

# Write results to CSV file
df = DataFrame([:doy => collect(1:365), 
    :tmax_50 => delta_tmax[1:365, 1], 
    :diff_10 => delta_tmax[1:365, 2], 
    :diff_50 => delta_tmax[1:365, 3], 
    :diff_90 => delta_tmax[1:365, 4]])
CSV.write("translate_ensemble_tmax_aggregated.csv", df, writeheader=true)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Visualise the results
x=df[:, :doy]
y = df[:, :diff_90]
ms=0.5


plot(x,
    y,
    seriestype=:line,
    markersize=ms,
    markerstrokewidth=0,
    showaxis=true,
    grid=false,
    legend=false,
    cbar=false)

