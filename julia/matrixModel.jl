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

# Location of results file
outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R"

outDir = "//home//jon//Desktop//OPRAM//results/pseudips"


resultsFile = "result_pseudips2019_par_thin10.jld2"


# Import results data
result = load(joinpath([outDir,resultsFile]),"single_stored_object")

# result = CSV.read(joinpath([outDir,resultsFile]), DataFrame,
#                    types = [Int32,Int16,Union{Missing, Int16}])



function development_matrix!(m, res)
    # Add one location's development
    for i in eachindex(res.emergeDOY)
        iStart = res.DOY[i]

        if i==length(res.emergeDOY)
            iEnd = 365
        else
            iEnd = res.DOY[i+1]-1
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

for i in 1:100
    @show i
    res = filter(x -> x.ID==ID[i], result)
    development_matrix!(m, res)
end

mnorm = convert.(Float64, m)
for i in 1:365
    s = sum(mnorm[i,:])
    if s>0
        mnorm[i,:] = mnorm[i,:] ./ s
    end
end




heatmap(m)


ev = eigen(mnorm)

emerge_dist = abs.(real(ev.vectors[:,end]))

t2 = mnorm * emerge_dist

m2 = mnorm^10000
heatmap(m2)


findall(real(ev.vectors[:,end]).>1e-3)

plot(emerge_dist)