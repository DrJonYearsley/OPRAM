#!/usr/bin/julia

# A script to plot the output of the degree day models
#
# Script plots spatial distribution of first emergence DOY and 
# average number of generations per year
# If more than one year is give an average over the years is plotted

using DataFrames
using CSV
using Dates
using Plots
using JLD2

# =================================================================================
# Set parameters for the visualisation
# If more than one year then take average across years
speciesName = "pseudips"
years = 2018:2021
thin = 1         # Spatial thining (thin=1 is 1km scale, thin=10 is 10km scale)
doy::Int32 = 1   # Day of year development started
save_figs = true  # If true save figures




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
        @time idx1 = indexin(unique(result.ID),result.ID) 
        @time idx2 = nrow(result) + 1 .- indexin(unique(result.ID),reverse(result.ID)) 


        # Create data frame to hold number of generations
        out = DataFrame(ID=result.ID[idx1], nGen=0.0, emergeDOY=0)
        out.north = mod.(out.ID, 1000)
        out.east = div.((out.ID .- out.north), 1000)

        # Count number of generations per year
        startDOY = zeros(Int32, length(idx1)) .+ doy

        while any(startDOY .> 0)
                # Find index of result that corresponds to desired day of year (ie doy>=result.DOY)
                idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], startDOY[i]) - 1 + idx1[i] for i = eachindex(idx1)]

                println(maximum(idx3))

                # Set any missing matches to zero
                match = startDOY .> 0 .&& idx2 + ones(Int64, size(idx1)) .> idx3
                complete = result.emergeDOY[idx3] .<= 365   # Is the generation complete within the year?

                # Increment generations
                out.nGen[match.&complete] .+= 1
                out.nGen[match.&.!complete] .+= (365 .- startDOY[match.&.!complete]) ./
                                                (result.emergeDOY[idx3[match.&.!complete]] - startDOY[match.&.!complete])
                if startDOY[1] == doy
                        out.emergeDOY[match.&complete] .+= result.emergeDOY[idx3[match.&complete]]
                end

                # Reset starting DOY
                startDOY = zeros(Int32, length(idx1))
                startDOY[match .& complete] .+= result.emergeDOY[idx3[match .& complete]] .+ 1

                println([sum(startDOY .> 0), maximum(result.emergeDOY[idx3])])
        end


        return (out)
end



# ============================================================



function year_average(years::Union{Vector{Int64}, Int64}, 
                      outDir::String, 
                      thin::Int64,
                      speciesName::String, doy::Int32)
        # Import model data for several years, extract results for a specific starting day and
        # average the results across all years
        #
        #    years       The years to import model results and average over
        #    outDir      Directory for saving results
        #    thin        The spatial thining parameter
        #    speciesName The species for which results will be calculated
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
                        # idx = [findfirst(x -> x == ID, d.ID) for ID in out.ID]
                        idx = indexin(d.ID, out.ID)
                        use = idx.!=nothing
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
                  cscale=:matter, discrete=false)

        x = sort(unique(coast.idx_east))
        y = sort(unique(coast.idx_north))
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




# Make sure years is a vector
years = collect(years)

# Create data to visualised
# number of generations and first day of emergence
@time d = year_average(years, outDir, thin, speciesName, doy)





# ===================================================================
# Visualise output

if length(years)==1
        year_label = string(years)
else
  year_label = string(minimum(years)) * " - " * string(maximum(years))
end

titleStr = "Emergence DOY (Start DOY = " * string(doy) * ",  " * year_label * ")"
plot_map("emergeDOY", d, coast, titleStr, :PiYG)
# Save output from the last plot
if save_figs
        savefig(speciesName * "_emergence_" * year_label * ".png")
end

titleStr="Number of Generations (" * year_label * ")";
plot_map("nGen", d, coast, titleStr, :PiYG)
# Save output from the last plot
if save_figs
        savefig(speciesName * "_ngen_" * year_label * ".png")
end

plot_map("nGenInt", d, coast, titleStr, :Dark2_3, true)
# Save output from the last plot

titleStr="Number of Generations SD (" * year_label * ")";
plot_map("nGenSD", d, coast, titleStr, :Accent)





# # +++++++++++++++++++++++++++++++++++++++++++
# Plot emergence for all locations across Ireland
# plot(result.DOY,
#         result.emergeDOY,
#         group=result.ID,
#         aspect_ratio=:equal,
#         lc=:black,
#         alpha=RGBA(0, 0, 0, 0.05),
#         xaxis="Start DOY",
#         yaxis="Emergence DOY",
#         legend=false)
# plot!([1:365], [1:365])


# # Same result visualised as a 2D histogram
# histogram2d(result.DOY,
#         result.emergeDOY,
#         nbinsx=200,
#         nbinsy=100,
#         xaxis="Start DOY",
#         yaxis="Emergence DOY")




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find a specific location

# Specify eastings and northings
location = [95000, 21000] 
latlonFile = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//locations.CSV"
latlonFile = joinpath(dataDir,"IE_grid_locations.csv")
paulFile = "first_occurrence_threshold_with_date_range_2_1961_Agrilus_anxius_GDD_multiple_dates.csv"

# Import location data
latlongs = CSV.read(latlonFile, DataFrame)

# paul = CSV.read(joinpath(outDir, "Code that is under development//Changing_Start_Dates",paulFile),DataFrame)

idx = latlongs.east.==location[1] .&& latlongs.north.==location[2]

d[d.ID.==latlongs.ID[idx],:]



# julia[julia.ID.==16090,:]