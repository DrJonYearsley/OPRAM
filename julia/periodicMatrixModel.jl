#!/usr/bin/julia

# Develop a periodic matrix model to simulate the population dynamics
# of an insect based upon degree-day development times and 
# estimates of reproductive output and survival.
#
# The broad approach was suggested by Jenouvrier and Visser (2011)(DOI 10.1007/s00484-011-0458-x)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 26th May 2025
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

outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//agrilus_anxius/"

importFilePrefix = "agrilus_anxius_IE"
importYears = 1971:1990
importYears = 2015:2020
thin = 1
adultLife = 14
idx_ID = 200    # ID of the location to use for development matrix


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create a development matrix for a single location
function development_matrix!(m, res, adultLifespan)
    # Add one location's development
    # Rows are emergence days
    # Columns are starting days
    # Sum down one column should equal
    for i in eachindex(res.emergeDOY)
        iStart = res.DOY[i]

        if i==length(res.emergeDOY)
            iEnd = 365
        else
            iEnd = min(365,res.DOY[i+1]-1)
        end

        # Find emergence day of year (i.e. modulo 365)
        adult = mod1.(res.emergeDOY[i] .+ collect(Int16, 0:adultLifespan-1),365)

        m[adult, iStart:iEnd] .+= 1
    end


end



function leslie_matrix(doy, nStages, K, r, s)
# Function to create a simple Leslie matrix for a single day
# Arguments:
#   doy: the day of year that the matrix corresponds to
#   nStages: the number of stages being tracked by the matrix
#   K: the development kernel. K[i,j] is the probability that an egg 
#      laid on day j develops to become reproductively competent on day i
#   r: the intrinsic growth rate of the population
#   s: the survival rate of the population


    # Create a Leslie matrix for the given day of year
    L = zeros(Float64, nStages, nStages)

    # Fill in the first row with the reproduction probabilities
    for i in 1:nStages
        if i <= size(K, 1)
            L[1, i] = r * K[doy, mod1(doy-i,size(K,2))]
        end
    end

    # Fill in the sub-diagonal with survival probabilities
    for i in 2:nStages
        L[i, i-1] = s
    end

    return L
end



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Create matrix 
# Columns are doy that development starts
# Rows are doy development ends
# mnorm[i,j] is the probability that an adult insect lays eggs on day j (start of development) 
# and the offspring becomes a competent adult on day i 
# The sum of mnorm over all i will be 1 
m = zeros(Int64,365,365);

for y in importYears
    resultsFile = importFilePrefix * "_" * string(y) * "_" * string(thin) * "km.jld2"

    @show resultsFile
    result = load(joinpath([outDir,resultsFile]),"single_stored_object")

    # Make maximum DOY 365
    result.DOY = min.(result.DOY, 365)

    # Identify unique IDs
    ID = unique(result.ID)

    
    res = filter(x -> x.ID==ID[idx_ID], result)
    development_matrix!(m, res, adultLife)
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


L = Vector{Matrix{Float64}}(undef, 365);
for i in 1:365
    L[i] = leslie_matrix(i, 365, mnorm, 4.5, 0.999)
end

a =  reduce(*,reverse(L))

ev = eigen(a);

emerge_dist = abs.(real(ev.vectors[:,end]));

# Display long-term growth rate
ev.values[end]

# heatmap(transpose(log10.(a)), yflip=true, aspect_ratio=:equal)
heatmap(a, yflip=true, aspect_ratio=:equal)


# Plot long term distribution of emergence dates
plot(emerge_dist)
