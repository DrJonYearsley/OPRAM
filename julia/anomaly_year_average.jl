# Julia program to calculate multi-year average of degree day model results
# and the corresponding anomaly for each year
#
# The degree day model results are generatde by OPRAM_main_program.jl
#
#
# Author: Jon Yearsley  (jon.yearsley@ucd.ie)
# Date: March 2025
#
# =============================================================================
# =============================================================================

# Load the required packages
using CSV
using DataFrames
using Statistics
using Dates
using JLD2
using Distributed
using SharedArrays




# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
run_params = (speciesName="ips_sexdentatus",      # Name of the species
    years=1991:2020,                    # Either a single year or collect(year1:year2)
    country="IE",                       # Country code (IE or NI)
    startMonth=01,                      # Month for start of development
    thinFactor=1,                       # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
    gridFile="IE_grid_locations.csv",   # File containing a 1km grid of lats and longs over Ireland 
    # (used for daylength calculations as well as importing and thining of meteo data)
    save_figs=true)  # If true save figures







# =================================================================================
# Specify directories for import and export of data

if isdir("//home//jon//Desktop//OPRAM")       # Linux workstation
    paths = (outDir="//home//jon//Desktop//OPRAM//results//",
        dataDir="//home//jon//DATA//OPRAM//")

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")   # Mac
    paths = (outDir="//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//",
        dataDir="//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//")

elseif isdir("//users//jon//Desktop//OPRAM//")
    paths = (outDir="//users//jon//Desktop//OPRAM//results//",
        dataDir="//users//jon//Desktop//OPRAM//")
end

include("OPRAM_io_functions.jl");
include("OPRAM_ddmodel_functions.jl");



# =========================================================
# =========================================================
# Import data across the years and calculate the average

# Initialise an empty data frame
d = DataFrame()

for y in eachindex(run_params.years)
    # Get the correct filename
    inFile = filter(x -> occursin(r"result_[A-Z]{2}_" * run_params.speciesName * "_" * string(run_params.years[y]) * "_1km.csv", x),
        readdir(joinpath(paths.outDir, run_params.speciesName)))

    @info "Importing data for year $(run_params.years[y])"

    d_year = CSV.read(joinpath(paths.outDir, run_params.speciesName, inFile[1]), DataFrame, missingstring="NA")

    # Make missing dates 3 years after start date
    idx = ismissing.(d_year.emergeDate)
    d_year.emergeDate[idx] .= Date(run_params.years[y], 1, 1) + Day(365 * 4)

    d = vcat(d, d_year)
end

# Calculate average across years (use median to avoid outliers)

# Define DOY for start and emerge dates (and set as integer)
d.emergeDOY = Dates.value.(d.emergeDate .- d.startDate) .+ 1
d.startMonth = Dates.month.(d.startDate)  # Use month rather than DOY to avoid leap year problems


# Group data frame by location and then starting DOY/Date 
df_group = groupby(d, [:ID, :east, :north, :startMonth])

# Calculate the median of nGenerations and emergeDOY
d_agg = combine(df_group,
    [:nGenerations, :emergeDOY] =>
        ((x, y) -> (nGenerations_median=quantile(x, 0.5),
            emergeDOY_median=quantile(y, 0.5))) =>
            AsTable)

   
d.nGenerations = convert(Vector{Union{Float32,Missing}}, d.nGenerations)
d.emergeDOY = convert(Vector{Union{Float32,Missing}}, d.emergeDOY)

idx = d.emergeDOY .> 365 * 3.2
d.nGenerations[idx] .= missing
d.emergeDOY[idx] .= missing

# COmbine origninal data frame with 30 year average
leftjoin!(d, d_agg, on=[:ID, :east, :north, :startMonth])
# Calculate anomalies
transform!(d, [:nGenerations, :nGenerations_median] => ((a, b) -> round.(a .- b, digits=3) ) => :nGenerations_anomaly)
transform!(d, [:emergeDOY, :emergeDOY_median] => ((a, b) -> a .- b ) => :emergeDOY_anomaly)

select!(d, Not(:startMonth))



# Write out the combined data frame
fileout1 = joinpath(paths.outDir, run_params.speciesName, "combined_" * run_params.speciesName *
                "_years_" * string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * ".csv")
CSV.write(fileout1,  d)


# Write out the average data frame
fileout2 = joinpath(paths.outDir, run_params.speciesName, "average_" * run_params.speciesName *
                "_years_" * string(minimum(run_params.years)) * "_" * string(maximum(run_params.years)) * ".csv")
CSV.write(fileout2,  d_agg)



using Plots
idx = d_agg.startMonth .== unique(d_agg.startMonth)[2]

if length(run_params.years) == 1
    year_label = string(run_params.years)
else
    year_label = string(minimum(run_params.years)) * "-" * string(maximum(run_params.years))
end

plot(d_agg.east[idx],
    d_agg.north[idx],
    zcolor=d_agg.nGenerations_median[idx],
    seriestype=:scatter,
    markersize=0.5,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="nGenerations",
    title="Average number of generations per year\n"*year_label,
    dpi=400)

if run_params.save_figs
        savefig(run_params.speciesName * "_ngen_" * year_label * ".png")
end


plot(d_agg.east[idx],
    d_agg.north[idx],
    zcolor=d_agg.emergeDOY_median[idx],
    seriestype=:scatter,
    markersize=0.5,
    markerstrokewidth=0,
    showaxis=false,
    grid=false,
    legend=false,
    cbar=true,
    aspect_ratio=:equal,
    colorbar_title="Emergence DOY",
    title="Average emergence DOY\n"*year_label,
                dpi=400)

if run_params.save_figs
        savefig(run_params.speciesName * "_emergeDOY_" * year_label * ".png")
end