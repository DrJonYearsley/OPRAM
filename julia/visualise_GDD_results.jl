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
speciesName = "agrilus"
years = 1961
thin = 1         # Spatial thining (thin=1 is 1km scale, thin=10 is 10km scale)
doy::Int32 = 10   # Day of year development started




# =================================================================================
# Specify directories for import and export of data
if isdir("//home//jon//DATA//OPRAM")
        outDir = "//home//jon//DATA//OPRAM//results//"
        meteoDir = "//home//jon//DATA//OPRAM//Irish_Climate_Data//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
        outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
        meteoDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"

elseif isdir("//users//jon//Desktop//OPRAM//")
        outDir = "//users//jon//Desktop//OPRAM//results//"
        meteoDir = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
end






# =================================================================================
# Import coastline data
coast = CSV.read(joinpath([meteoDir, "..", "coastline.csv"]), DataFrame,
        types=[Float64, Float64, Int, Int, Int])


# =================================================================================
# Start of function definitions
function extract_results(doy::Int32, result::DataFrame)
        # Function to calculate number of generations per year, and first
        # day of adult emergence based on a starting development on doy

        # Find start and end indicies for each location
        idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)]
        idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)]

        # Create data frame to hold number of generations
        res = DataFrame(ID=result.ID[idx1], nGen=0.0, emergeDOY=0)
        res.north = mod.(res.ID, 1000)
        res.east = div.((res.ID .- res.north), 1000)

        # Count number of generations per year
        startDOY = zeros(Int32, length(idx1)) .+ doy

        while any(startDOY .> 0)
                # Find index of result that corresponds to desired day of year (ie doy>=result.DOY)
                idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], startDOY[i]) - 1 + idx1[i] for i = eachindex(idx1)]

                # println(maximum(idx3))

                # Set any missing matches to zero
                match = startDOY .> 0 .&& idx2 + ones(Int64, size(idx1)) .> idx3
                complete = result.emergeDOY[idx3] .< 365

                # Increment generations
                res.nGen[match.&complete] .+= 1
                res.nGen[match.&.!complete] .+= (365 .- startDOY[match.&.!complete]) ./
                                                (result.emergeDOY[idx3[match.&.!complete]] - startDOY[match.&.!complete])
                if startDOY[1] == doy
                        res.emergeDOY[match.&complete] .+= result.emergeDOY[idx3[match.&complete]]
                end

                # Reset starting DOY
                startDOY = zeros(Int32, length(idx1))
                startDOY[match .& complete] .+= result.emergeDOY[idx3[match .& complete]] .+ 1

                # println([sum(startDOY .> 0), maximum(result.emergeDOY[idx3])])

        end


        return (res)
end



# ============================================================



function year_average(years::Union{Vector{Int64}, Int64}, outDir::String, thin::Int64,
        speciesName::String, doy::Int32)

        out = []
        for y in eachindex(years)
                println("Processing data for " * string(years[y]))

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
                        idx = [findfirst(x -> x == ID, d.ID) for ID in out.ID]
                        out.nGen[idx] .+= d.nGen
                        out.nGenInt[idx] += floor.(d.nGen)   
                        out.nGenSD[idx] .+= d.nGen .^ 2
                        out.emergeDOY[idx] .+= d.emergeDOY
                        out.emergeSD[idx] .+= d.emergeDOY .^ 2
                        out.N[idx] .+= 1
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


function plot_map(varStr::String, data::DataFrame, coast::DataFrame, 
        titleStr::String="NoTitle", cscale=:matter, discrete=false)

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







# Create data to visualised
# number of generations and first day of emergence
@time d = year_average(years, outDir, thin, speciesName, doy)





# ===================================================================
# Visualise output

year_label = string(minimum(years)) * " - " * string(maximum(years))

titleStr = "Emergence DOY (Start DOY = " * string(doy) * ",  " * year_label * ")"
plot_map("emergeDOY", d, coast, titleStr, :RdYlGn)


titleStr="Number of Generations (" * year_label * ")";
plot_map("nGen", d, coast, titleStr, :PiYG)
plot_map("nGenInt", d, coast, titleStr, :Dark2_3, true)

titleStr="Number of Generations SD (" * year_label * ")";
plot_map("nGenSD", d, coast, titleStr, :Accent)



# Save output from the last plot
savefig("agrilus_ngen_2018_2020.png")


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

# Import location data
latlongs = CSV.read(latlonFile, DataFrame)

idx = latlongs.east.==location[1] .&& latlongs.north.==location[2]

d[d.ID.==latlongs.ID[idx],:]