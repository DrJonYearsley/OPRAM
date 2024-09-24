#!/usr/bin/julia

# A script to plot the output of the degree day models

using DataFrames
using CSV
using Dates
using Plots
using JLD2

# Location of results file
resultsFile = "ips_typographus//result_ips_typographus1961_par_thin1.jld2"
doy::Int32 = 1   # Day of year to be visualised




# Specify directories for import and export of data
if isdir("//home//jon//Desktop//OPRAM")
        outDir = "//home//jon//Desktop//OPRAM//results//"
        meteoDir = "//home//jon//Desktop//OPRAM//Irish_Climate_Data//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
        outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
        meteoDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Irish Climate Data//"

elseif isdir("//users//jon//Desktop//OPRAM//")
        outDir = "//users//jon//Desktop//OPRAM//results//"
        meteoDir = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
end






# Import results data
# result = CSV.read(joinpath([outDir,resultsFile]), DataFrame,
#                    types = [Int32,Int32,Union{Missing, Int32}])

result = load(joinpath([outDir, resultsFile]), "single_stored_object")


# ============================================================
# ============================================================

function doy_idx(doy::Int32, result::DataFrame)
        # Produce a dataframe corresponding to a day of year (doy)

        # Find start and end indicies for each location
        idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)]
        idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)]


        # Find index of result that corresponds to desired day of year (ie doy>=result.DOY)
        idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], doy) - 1 + idx1[i] for i = eachindex(idx1)]

        # Set any missing matches to zero
        match = idx2 + ones(Int64, size(idx1)) .> idx3

        res = DataFrame(ID=result.ID[idx1[match]],
                emergeDOY=result.emergeDOY[idx3[match]])
        res.north = mod.(res.ID, 1000)
        res.east = div.((res.ID .- res.north), 1000)

        return (res)
end


function num_generations(doy::Int32, result::DataFrame)
        # Function to calculate number of generations per year, starting
        # development on doy

        # Find start and end indicies for each location
        idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)]
        idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)]

        # Create data frame to hold number of generations
        res = DataFrame(ID=result.ID[idx1], nGen=0.0)
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

                # Reset starting DOY
                startDOY = zeros(Int32, length(idx1))
                startDOY[match.&complete] .+= result.emergeDOY[idx3[match.&complete]] .+ 1

                sum(startDOY .> 0)

        end


        return (res)
end

# ============================================================

# ============================================================



# Create result for a start day of year
@time d = doy_idx(doy, result)

# Create number of generations
@time gen = num_generations(doy, result)


# ===================================================================
# Visualise output

# plot map of emergence day of year
scatter(d.east,
        d.north,
        zcolor=d.emergeDOY,
        color=:matter,
        cbar=true,
        legend=false,
        showaxis=false,
        markersize=1,
        markerstrokewidth=0,
        title="doy = " * string(doy))



scatter(gen.east,
        gen.north,
        zcolor=gen.nGen,
        color=:Accent_6,
        clims=:auto,
        cbar=true,
        legend=false,
        showaxis=false,
        markersize=1,
        markerstrokewidth=0,
        title="Number of Generations")


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