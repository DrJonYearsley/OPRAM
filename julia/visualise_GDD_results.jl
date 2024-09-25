#!/usr/bin/julia

# A script to plot the output of the degree day models

using DataFrames
using CSV
using Dates
using Plots
using JLD2

# Location of results file
speciesName = "agrilus"
years = [2018,2019,2020]
thin = 1

doy::Int32 = 1   # Day of year to be visualised




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






# Import results data
coast = CSV.read(joinpath([meteoDir,"..","coastline.csv"]), DataFrame,
        types=[Float64, Float64, Int64, Int64, Int64])
fileName = "result_" * speciesName * string(years[1]) * "_par_thin" * string(thin) * ".jld2"
resultFile = joinpath([outDir, speciesName, fileName])
result = load(resultFile, "single_stored_object")


# ============================================================
# ============================================================

# function doy_idx(doy::Int32, result::DataFrame)
#         # Produce a dataframe corresponding to a day of year (doy)

#         # Find start and end indicies for each location
#         idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)]
#         idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)]


#         # Find index of result that corresponds to desired day of year (ie doy>=result.DOY)
#         idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], doy) - 1 + idx1[i] for i = eachindex(idx1)]

#         # Set any missing matches to zero
#         match = idx2 + ones(Int64, size(idx1)) .> idx3

#         res = DataFrame(ID=result.ID[idx1[match]],
#                 emergeDOY=result.emergeDOY[idx3[match]])
#         res.north = mod.(res.ID, 1000)
#         res.east = div.((res.ID .- res.north), 1000)

#         return (res)
# end


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

                # Set any missing matches to zero
                match = startDOY .> 0 .&& idx2 + ones(Int64, size(idx1)) .> idx3
                complete = result.emergeDOY[idx3] .<= 365

                # Increment generations
                res.nGen[match.&complete] .+= 1
                res.nGen[match.&.!complete] .+= (365 .- startDOY[match.&.!complete]) ./
                                                (result.emergeDOY[idx3[match.&.!complete]] - startDOY[match.&.!complete])
                if startDOY[1]==doy
                        res.emergeDOY[match .& complete] .+= result.emergeDOY[idx3[match.&complete]]
                end

                # Reset starting DOY
                startDOY = zeros(Int32, length(idx1))
                startDOY[match.&complete] .+= result.emergeDOY[idx3[match.&complete]] .+ 1

                sum(startDOY .> 0)

        end


        return (res)
end



# ============================================================



function year_average(years::Vector{Int64}, outDir::String, thin::Int64, 
        speciesName::String, doy::Int32, result::DataFrame)

        out = [];
        for y in eachindex(years)
                fileName = "result_" * speciesName * string(years[y]) * "_par_thin" * string(thin) * ".jld2"
                resultFile = joinpath([outDir, speciesName, fileName])
                result = load(resultFile, "single_stored_object")

                # Average across years
                d = extract_results(doy, result)
                if y == 1
                        out = d
                        out.N = ones(size(out.nGen))
                else
                        idx = [findfirst(x -> x == ID, d.ID) for ID in out.ID]
                        out.nGen[idx] .+= d.nGen
                        out.emergeDOY[idx] .+= d.emergeDOY
                        out.N[idx] .+= 1
                end
        end

        # Calculate average across the years
        out.nGen = out.nGen ./ out.N
        tmp = Float64.(out.emergeDOY)
        out.emergeDOY = tmp ./ out.N

        return out
end


# ============================================================

# ============================================================




# Create number of generations and first day of emergence
@time d = year_average(years, outDir, thin, speciesName, doy, result)

# ===================================================================
# Visualise output

# plot map of emergence day of year
plot(d.east,
        d.north,
        seriestype=:scatter,
        zcolor=d.emergeDOY,
        axis_ratio=:equal,
        color=:matter,
        cbar=true,
        legend=false,
        showaxis=false,
        markersize=0.5,
        markerstrokewidth=0,
        markershape=:square,
        title="doy = " * string(doy),
        dpi=600);

plot!(coast.idx_east, 
        coast.idx_north,
        seriestype=:scatter,
        showaxis=false,
        grid = false,
        markersize = 0.5,
        markercolor="black",
        dpi=600)
# savefig("test.png")

plot(d.east,
        d.north,
        zcolor=d.nGen,
        seriestype=:scatter,
        color=:Accent_6,
        axis_ratio=:equal,
        clims=:auto,
        cbar=true,
        legend=false,
        showaxis=false,
        markersize=0.5,
        markerstrokewidth=0,
        markershape=:square,
        dpi=600,
        title="Number of Generations");

        plot!(coast.idx_east, 
        coast.idx_north,
        seriestype=:scatter,
        showaxis=false,
        grid = false,
        markersize = 0.5,
        markercolor="black",
        dpi=600)

# +++++++++++++++++++++++++++++++++++++++++++
# Plot emergence for all locations across Ireland
plot(result.DOY,
        result.emergeDOY,
        group=result.ID,
        aspect_ratio=:equal,
        lc=:black,
        alpha=RGBA(0, 0, 0, 0.05),
        xaxis="Start DOY",
        yaxis="Emergence DOY",
        legend=false)
plot!([1:365], [1:365])


# Same result visualised as a 2D histogram
histogram2d(result.DOY,
        result.emergeDOY,
        nbinsx=200,
        nbinsy=100,
        xaxis="Start DOY",
        yaxis="Emergence DOY")