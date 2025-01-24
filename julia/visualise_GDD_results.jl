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

# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
speciesName = "ips_cembrae"
years = collect(2000:2002)  # Either a single year or collect(year1:year2)
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





# =================================================================================
# Import coastline data
coast = CSV.read(joinpath([dataDir, "coastline.csv"]), DataFrame,
        types=[Float64, Float64, Int, Int, Int])


# =================================================================================
# Import grid location data
grid = CSV.read(joinpath([dataDir, "IE_grid_locations.csv"]), DataFrame,
        types=[Int, Int, Int, String, String, String, Float64, Float64])

# Calculate east and north of hectad bottom left corner
grid.east_hectad = 10000.0 .* floor.(grid.east/10000)
grid.north_hectad = 10000.0 .* floor.(grid.north/10000)

# =================================================================================
# Start of function definitions
function extract_results(doy::Int32, result::DataFrame)
        # Function to calculate number of generations per year, and first
        # day of adult emergence based on a starting development on doy
        #
        #   doy    day of year when larval development begins
        #   result data frame containing the results from the model
        #  
        #  Output:
        #   a data frame with the following columns:
        #       ID            unique ID for each soatial location        
        #       emergeDOY     adult emergence day of year 
        #       nGen          the number of generations within the year
        #       east          index for eastings of the unique spatial locations in result
        #       north         index for northings of the unique spatial location in result
        ########################################################################

        # Find start and end indicies for each location
        # @time idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)]
        # @time idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)]
        @time idx1 = indexin(unique(result.ID), result.ID)
        @time idx2 = vcat(idx1[2:end] .- 1, length(result.ID))


        # Create data frame to hold number of generations
        out_res = DataFrame(ID=result.ID[idx1], nGen=0.0, emergeDOY=0)
        out_res.north = mod.(out_res.ID, 1000)
        out_res.east = div.((out_res.ID .- out_res.north), 1000)

        # Count number of generations per year
        startDOY = zeros(Int32, length(idx1)) .+ doy
        first_pass = true

        # Set startDOY to zero for a location if no more generations can be completed in the year
        while any(startDOY .> 0)
                # Find index of result that corresponds to desired day of year (ie doy>=result.DOY)
                idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], startDOY[i]) - 1 + idx1[i] for i = eachindex(idx1)]

                println(maximum(idx3))

                # Find out if end of current developmental generation is in the results (i.e. occurs before data on next location)
                development_complete = startDOY .> 0 .&& idx2 .+ 1 .> idx3

                # Find out if development happens within the first year
                if (idx3[end] == nrow(result) + 1)
                        # If final entry resulted in no developement (i.e. index is nrow+1) then set within year to false
                        within_a_year = result.emergeDOY[idx3[1:end-1]] .<= 365   # Is the generation complete within the year?
                        push!(within_a_year, false)
                else
                        within_a_year = result.emergeDOY[idx3] .<= 365   # Is the generation complete within the year?
                end

                # Increment generations
                full_generation = development_complete .& within_a_year
                partial_generation = development_complete .& .!within_a_year

                # Add on one full generation
                out_res.nGen[full_generation] .+= 1

                # Calculate end of year fraction of time towards next generation
                out_res.nGen[partial_generation] .+= (365 .- startDOY[partial_generation]) ./
                                                     (result.emergeDOY[idx3[partial_generation]] - startDOY[partial_generation])

                # If this is the first time through the loop, record the emergance day
                if first_pass
                        out_res.emergeDOY[full_generation] .= result.emergeDOY[idx3[full_generation]]
                        first_pass = false
                end

                # Reset starting DOY to zero
                startDOY = zeros(Int32, length(idx1))

                # If a full gneration completed within the year, update starting doy
                startDOY[full_generation] .= result.emergeDOY[idx3[full_generation]] .+ 1

                # Remove startDOY >= 365
                startDOY[startDOY.>=365] .= 0

                # println([sum(startDOY .> 0), maximum(result.emergeDOY[idx3])])
        end


        return (out_res)
end



# ============================================================



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

                fileName = "result_" * speciesName * string(years[y]) * "_par_thin" * string(thin) * ".jld2"
                resultFile = joinpath([outDir, speciesName, fileName])
                result = load(resultFile, "single_stored_object")

                # Average across years
                d = extract_results(doy, result)
                if y == 1
                        out = d
                        out.nGenInt = floor.(d.nGen)
                        out.N = ones(size(out.nGen))
                        out.nGenSD = d.nGen .^ 2
                        out.emergeSD = d.emergeDOY .^ 2
                else
                        idx = indexin(d.ID, out.ID)    # Equate spatial IDs in d to spatial IDs in out
                        use = idx .!= nothing .&& d.emergeDOY .!= 0  # Ignore locations with emerge DOY of zero
                        out.nGen[idx[use]] .+= d.nGen[use]
                        out.nGenInt[idx[use]] += floor.(d.nGen[use])
                        out.nGenSD[idx[use]] .+= d.nGen[use] .^ 2
                        out.emergeDOY[idx[use]] .+= d.emergeDOY[use]
                        out.emergeSD[idx[use]] .+= d.emergeDOY[use] .^ 2
                        out.N[idx[use]] .+= 1
                end
        end

        # Calculate average across the years
        out.nGen = out.nGen ./ out.N
        out.nGenInt = out.nGenInt ./ out.N
        tmp = Float64.(out.emergeDOY)
        out.emergeDOY = tmp ./ out.N

        # Calculate standard deviation across years
        out.nGenSD = (out.nGenSD ./ out.N - out.nGen .^ 2) .^ 0.5
        tmp = Float64.(out.emergeSD)
        out.emergeSD = (tmp ./ out.N - out.emergeDOY .^ 2) .^ 0.5


        return out
end



# ============================================================


function plot_map(varStr::String,
        data::DataFrame,
        coast::DataFrame,
        titleStr::String="NoTitle",
        cscale=:matter, 
        discrete=false)

        x = sort(unique(coast.idx_east))
        y = sort(unique(coast.idx_north))





        C = fill(NaN, maximum(y), maximum(x)) # make a 2D C with NaNs
        for i in 1:nrow(data)
                C[data.north[i], data.east[i]] = data[i, varStr]
        end


        heatmap(C,
                showaxis=false,
                grid=true,
                axis_ratio=:equal,
                color=cgrad(cscale, categorical=discrete),
                legend=false,
                cbar=true,
                title=titleStr)




        plot!(coast.idx_east,
                coast.idx_north,
                seriestype=:scatter,
                showaxis=false,
                grid=false,
                markersize=0.5,
                markercolor="black",
                dpi=600)

end

# =================================================================================
# End of function definitions
# =================================================================================



# Create data to visualised
# number of generations and first day of emergence
@time d = year_average(years, outDir, thin, speciesName, doy, hectad)

# Combine grid and results
d2 = leftjoin(d, grid[:,[1,4,5,6]], on=:ID)

# Aggregate by hectad
d3a = combine(groupby(d2, :hectad), [:nGen] .=> maximum)
d3b = combine(groupby(d2, :hectad), [:emergeDOY] .=> minimum)
d3c = leftjoin(d3a, d3b, on=:hectad)

d3 = leftjoin(d3c, grid, on=:hectad)

# ===================================================================
# Visualise output

if length(years) == 1
        year_label = string(years)
else
        year_label = string(minimum(years)) * " - " * string(maximum(years))
end

titleStr = "Emergence DOY (Start DOY = " * string(doy) * ",  " * year_label * ")"
d_subset = d[d.emergeDOY.!=0,:]
plot_map("emergeDOY", d_subset, coast, titleStr, :matter)
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



# Try plotting just hectads
plot(d3.east_hectad,
     d3.north_hectad,
        seriestype=:scatter,
        zcolor=d3.nGen_maximum,
        color=cgrad(:PiYG, categorical=false),
        showaxis=false,
        legend=false,
        grid=false,
        cbar=true,
        clims=(minimum(d.nGen),maximum(d.nGen)),
        markersize=3,
        aspect_ratio=:equal,
        dpi=600,
        title="Maximum Number of Generations (" * year_label * ")");

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
        zcolor=d3.emergeDOY_minimum,
        color=cgrad(:matter, categorical=false),
        showaxis=false,
        legend=false,
        grid=false,
        cbar=true,
        clims=(minimum(d.emergeDOY[d.emergeDOY.>0]),maximum(d.emergeDOY)),
        markersize=3,
        aspect_ratio=:equal,
        dpi=600,
        title="Minimum Emergence DOY (Start DOY = " * string(doy) * ",  " * year_label * ")");

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