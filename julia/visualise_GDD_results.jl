#!/usr/bin/julia

# A script to plot the output of the degree day models
#
# Script plots spatial distribution of first emergence DOY and 
# average number of generations per year
# If more than one year is give an average over the years is plotted

using DataFrames
using Statistics
using CSV
using Dates
using Plots
using JLD2
using Distributed
using SharedArrays


# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
speciesName = "ips_sexdentatus"
years = collect(2018:2021)  # Either a single year or collect(year1:year2)
thin = 1         # Spatial thining (thin=1 is 1km scale, thin=10 is 10km scale)
doy::Int32 = 1   # Day of year development started
save_figs = true  # If true save figures
hectad = true     # If true aggregate data into 10km squares



# =================================================================================
# Specify directories for import and export of data

if isdir("//home//jon//Desktop//OPRAM")
        outDir = "//home//jon//Desktop//OPRAM//results//"
        dataDir = "//home//jon//DATA//OPRAM//"
        meteoDir_IE = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"
        meteoDir_NI = "//home//jon//DATA//OPRAM//Northern_Ireland_Climate_Data//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
        outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
        dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
        meteoDir_IE = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"
        meteoDir_NI = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Northern_Ireland_Climate_Data//"

elseif isdir("//users//jon//Desktop//OPRAM//")
        outDir = "//users//jon//Desktop//OPRAM//results//"
        dataDir = "//users//jon//Desktop//OPRAM//"
        meteoDir_IE = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
        meteoDir_NI = "//users//jon//Desktop//OPRAM//Northern_Ireland_Climate_Data//"
end

include("import_functions.jl")
include("GDD_functions.jl")


# =================================================================================
# Import coastline data
coast = CSV.read(joinpath([dataDir, "coastline.csv"]), DataFrame,
        types=[Float64, Float64, Int, Int, Int])


# =================================================================================
# Import grid location data
grid_df = CSV.read(joinpath([dataDir, "IE_grid_locations.csv"]), DataFrame,
        types=[Int, Int, Int, String,  Float64, Float64,String, String,String])

# Calculate east and north of hectad bottom left corner
grid_df.east_hectad = 10000.0 .* floor.(grid_df.east/10000)
grid_df.north_hectad = 10000.0 .* floor.(grid_df.north/10000)

# =================================================================================
# Start of function definitions



function year_average(years::Union{Vector{Int64},Int64},
        outDir::String,
        thin::Int64,
        speciesName::String, 
        doy::Int32,
        hectad::Bool)
        # Import model data for several years, extract results for a specific starting day and
        # average the results across all years
        #
        #    years       The years to import model results and average over
        #    outDir      Directory for saving results
        #    thin        The spatial thining parameter
        #    speciesName The species for which results will be calculated
        #    doy         Day of year for development to start
        #    hectad      Average data over 10km squares
        #
        #   Output:
        #     a data frame with columns:
        #       ID            unique ID for each soatial location
        #       east          eastings of the spatial locations in result
        #       north         northings of the spatial location in result
        #       emergeDOY     adult emergence day of year (average over years)
        #       nGen          number of generations within the year (average over years)
        #                     Can be fractional for each year representing partial development 
        #       nGenInt       number of generations within the year (average over years)
        #                     This is nGen rounded down to nearest whole number 
        #       emergeDOYSD   adult emergence day of year (std. dev. average over years)
        #       nGenSD        number of generations within the year (std. dev. over years)
        #       N             number of years being averaged over
        #############################################################

        out = []

        for y in eachindex(years)
                @info "Processing data for " * string(years[y])

                fileName = "result_" * speciesName * string(years[y]) * "_" * string(thin) * "km.jld2"
                resultFile = joinpath([outDir, speciesName, fileName])
                result = load(resultFile, "single_stored_object")

                # Average across years
                d = extract_results(doy, result)
                if y == 1
                        out = d
                        out.nGenInt = floor.(d.nGenerations)
                        out.N = ones(size(out.nGenerations))
                        out.nGenSD = d.nGenerations .^ 2
                        out.emergeSD = d.emergeDOY .^ 2
                else
                        idx = indexin(d.ID, out.ID)    # Equate spatial IDs in d to spatial IDs in out
                        use = idx .!= nothing .&& d.emergeDOY .!= 0  # Ignore locations with emerge DOY of zero
                        out.nGenerations[idx[use]] .+= d.nGenerations[use]
                        out.nGenInt[idx[use]] += floor.(d.nGenerations[use])
                        out.nGenSD[idx[use]] .+= d.nGenerations[use] .^ 2
                        out.emergeDOY[idx[use]] .+= d.emergeDOY[use]
                        out.emergeSD[idx[use]] .+= d.emergeDOY[use] .^ 2
                        out.N[idx[use]] .+= 1
                end
        end

        # Calculate average across the years
        out.nGenerations = out.nGenerations ./ out.N
        out.nGenInt = out.nGenInt ./ out.N
        tmp = Float64.(out.emergeDOY)
        out.emergeDOY = tmp ./ out.N

        # Calculate standard deviation across years
        out.nGenSD = (out.nGenSD ./ out.N - out.nGenerations .^ 2) .^ 0.5
        tmp = Float64.(out.emergeSD)
        out.emergeSD = (tmp ./ out.N - out.emergeDOY .^ 2) .^ 0.5


        return out
end



# ============================================================


function plot_map(varStr::String,
        data::DataFrame,
        coast::Union{Nothing,DataFrame},
        titleStr::String="NoTitle",
        cscale=:matter, 
        discrete=false,
        size=0.5)

if isnothing(coast)
        x=sort(unique(data.east))
        y=sort(unique(data.north))
else
        x = sort(unique(coast.idx_east))
        y = sort(unique(coast.idx_north))

end

        C = fill(NaN, maximum(y), maximum(x)) # make a 2D C with NaNs
        for i in 1:nrow(data)
                C[data.north[i], data.east[i]] = data[i, varStr]
        end


        heatmap(C,
                showaxis=false,
                grid=false,
                axis_ratio=:equal,
                color=cgrad(cscale, categorical=discrete),
                legend=false,
                cbar=true,
                title=titleStr)



       if coast!=nothing

        plot!(coast.idx_east,
                coast.idx_north,
                seriestype=:scatter,
                showaxis=false,
                grid=false,
                markersize=size,
                markercolor="black",
                dpi=600)
        end

end

# =================================================================================
# End of function definitions
# =================================================================================



# Create data to visualised
# number of generations and first day of emergence
@time d = year_average(years, outDir, thin, speciesName, doy, hectad)

# Combine grid and results
d2 = leftjoin(d, grid_df[:,[1,2,3,4,5,6,8]], on=:ID)

# Aggregate by hectad
d3a = combine(groupby(d2, :hectad), [:nGenerations] .=> x -> quantile(x,0.95))
d3b = combine(groupby(d2, :hectad), [:emergeDOY] .=> x -> quantile(x,0.05))
d3c = leftjoin(d3a, d3b, on=:hectad)

d3 = leftjoin(d3c, grid_df, on=:hectad)



# Create subset for Munster
d_munster = subset(d2, :province => x-> x.=="Munster")


# ===================================================================
# Visualise output

if length(years) == 1
        year_label = string(years)
else
        year_label = string(minimum(years)) * " - " * string(maximum(years))
end

titleStr = "Emergence DOY (Start DOY = " * string(doy) * ",  " * year_label * ")"
d_subset = d2[d2.emergeDOY.!=0,:]
plot_map("emergeDOY", d_subset, nothing, titleStr, :matter)
# Save output from the last plot
if save_figs
        savefig(speciesName * "_emergence_" * year_label * ".png")
end



titleStr = "Number of Generations (" * year_label * ")";
plot_map("nGen", d, coast, titleStr, :PiYG)
# Save output from the last plot
if save_figs
        savefig(speciesName * "_ngen_" * year_label * ".png")
end

titleStr = "Number of Whole Generations (" * year_label * ")";
plot_map("nGenInt", d, coast, titleStr, :Dark2_3, true)


titleStr = "Number of Generations SD (" * year_label * ")";
plot_map("nGenSD", d, coast, titleStr, :BuPu)


titleStr = "Emergence DOY  SD (" * year_label * ")";
plot_map("emergeSD", d, coast, titleStr, :BuPu)


titleStr = "Number of Years of Data (" * year_label * ")";
plot_map("N", d, coast, titleStr, :YlGnBu_5, true)
# Save output from the last plot
if save_figs
        savefig(speciesName * "_ndatayears_" * year_label * ".png")
end


# =================================================
# Try plotting just hectads
plot(d3.east_hectad.+5000.0,
     d3.north_hectad.+5000.0,
        seriestype=:scatter,
        zcolor=d3.nGen_function,
        color=cgrad(:PiYG, categorical=false),
        showaxis=false,
        legend=false,
        grid=false,
        cbar=true,
        clims=(minimum(d.nGen),maximum(d.nGen)),
        markersize=3,
        markerstrokewidth=0.5,
        aspect_ratio=:equal,
        dpi=600,
        title="90th Percentile in\n Number of Generations (" * year_label * ")");

plot!(coast.east,
        coast.north,
                seriestype=:scatter,
                showaxis=false,
                grid=false,
                markersize=1,
                markercolor="black",
                dpi=600)
if save_figs
        savefig(speciesName * "_hectads_ngen_" * year_label * ".png")
end


plot(d3.east_hectad,
     d3.north_hectad,
        seriestype=:scatter,
        zcolor=d3.emergeDOY_function,
        color=cgrad(:matter, categorical=false),
        showaxis=false,
        legend=false,
        grid=false,
        cbar=true,
        clims=(minimum(d.emergeDOY[d.emergeDOY.>0]),maximum(d.emergeDOY)),
        markersize=3,
        markerstrokewidth=0.5,
        aspect_ratio=:equal,
        dpi=600,
        title="10th Percentile in \n Emergence DOY (Start DOY = " * string(doy) * ",  " * year_label * ")");

plot!(coast.east,
        coast.north,
                seriestype=:scatter,
                showaxis=false,
                grid=false,
                markersize=1,
                markercolor="black",
                dpi=600)
if save_figs
        savefig(speciesName * "_hectads_emergence_" * year_label * ".png")
end


# Plot Munster region

coast_munster = subset(coast, 
:idx_east => x->x.>=minimum(d_munster.east),
:idx_east => x->x.<=maximum(d_munster.east), 
:idx_north => x->x.>=minimum(d_munster.north),
:idx_north => x->x.<=maximum(d_munster.north))


titleStr = "Number of Generations \n(" * year_label * ", Munster)";
plot_map("nGen", d_munster, nothing, titleStr, :PiYG, false, 1)


# Save output from the last plot
if save_figs
        savefig(speciesName * "_ngen_munster_" * year_label * ".png")
end


titleStr = "Emergence DOY \n(Start DOY = " * string(doy) * ",  " * year_label * ", Munster)"
plot_map("emergeDOY", 
        subset(d_munster,:emergeDOY=> x->x.!=0),
        coast_munster, 
        titleStr, 
        :matter,
        false,
        1)

# Save output from the last plot
if save_figs
        savefig(speciesName * "_emergence_munster" * year_label * ".png")
end