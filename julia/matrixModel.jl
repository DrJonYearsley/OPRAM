#!/usr/bin/julia

# Create a matrix to calculate pest dynamics from degree day models
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 20th July 2024
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

using DataFrames
using CSV
using Plots
using LinearAlgebra
using Profile
using JLD2
using StatsBase

# Location of results file
outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R"

outDir = "//home//jon//Desktop//OPRAM//results/pseudips"

importFilePrefix = "result_pseudips"
importYears = 1971:1990
# importYears = 2001:2020
thin = 10
idx_ID = 200


# Import results data

# resultsFile = "result_pseudips2019_par_thin10.jld2"
# result = load(joinpath([outDir,resultsFile]),"single_stored_object")

# result = CSV.read(joinpath([outDir,resultsFile]), DataFrame,
#                    types = [Int32,Int16,Union{Missing, Int16}])



function development_matrix!(m, res)
    # Add one location's development
    # Rows are emerge days
    # Columns are starting days
    # Sum down one column should equal
    for i in eachindex(res.emergeDOY)
        iStart = res.DOY[i]

        if i==length(res.emergeDOY)
            iEnd = 365
        else
            iEnd = min(365,res.DOY[i+1]-1)
        end

        m[mod1(res.emergeDOY[i],365), iStart:iEnd] .+= 1
    end
end


# Create matrix columns are doy that development starts
# Rows are doy development ends

# m[i,j] is the probability that an insect starting development on day j emerges on day i
m = zeros(Int64,365,365)
ID = unique(result.ID)

# Make maximum DOY 365
result.DOY = min.(result.DOY, 365)

for y in importYears
    resultsFile = importFilePrefix * string(y) * "_par_thin" * string(thin) * ".jld2"

    @show resultsFile
    result = load(joinpath([outDir,resultsFile]),"single_stored_object")

    res = filter(x -> x.ID==ID[idx_ID], result)
    development_matrix!(m, res)
end

mnorm = convert.(Float64, m);
for i in 1:size(mnorm,1)
    # Rows are emerge days
    # Columns are starting days
    # Sum down one column should equal one when normalised

    s = sum(mnorm[:,i])
    if s>0
        mnorm[:,i] = mnorm[:,i] ./ s
    end
end




heatmap(log10.(mnorm))


ev = eigen(mnorm);

emerge_dist = abs.(real(ev.vectors[:,end]));

t2 = mnorm * emerge_dist;


# Display long-term emergence dates
findall(abs.(real(ev.vectors[:,end])).>1e-3)
ev.values[end]

# Plot long term distribution of emergence dates
plot(emerge_dist)


# Produce plot of development
d = transpose(reinterpret(reshape,Int64, findall(!iszero, mnorm)));
plot(d[:,2], [d[:,1], d[:,2]], 
    aspect_ratio=:equal,
    ylimits=(0,365),
    xlimits=(0, 365),
    color=RGBA{Float64}(0,0,0,0.2))


    describe(emerge_dist)