#!/usr/bin/julia

# A script to plot the output of the degree day models

using DataFrames
using CSV
using Dates
using Plots

# Location of results file
#outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R"

outDir = "//home//jon//Desktop//OPRAM//"


resultsFile = "result_halyomorpha1974_par_thin1.csv"
doy::Int32 = 1   # Day of year to be visualised







# Import results data
result = CSV.read(joinpath([outDir,resultsFile]), DataFrame,
                   types = [Int32,Int32,Union{Missing, Int32}])



function doy_idx(doy, result)
# Produce a dataframe corresponding to a day of year (doy)

        # Find start and end indicies for each location
        idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)];
        idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)];


        # Find index of result that corresponds to desired day of year (ie doy>=result.DOY)
        idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], doy)-1+idx1[i] for i = eachindex(idx1)];

        # Set any missing matches to zero
        match = idx2+ ones(Int64,size(idx1)).> idx3

        res = DataFrame(ID = result.ID[idx1[match]], 
                        emergeDOY = result.emergeDOY[idx3[match]])
        res.north = mod.(res.ID,1000)
        res.east = div.((res.ID .- res.north), 1000)

        return(res)
end





# Create result for a start day of year
@time d = doy_idx(doy, result)



# ===================================================================
# Visualise output

# plot map of emergence day of year
scatter(d.east, 
        d.north,
        zcolor=d.emergeDOY,
        color=:matter,
        cbar=true,
        legend=false,
        showaxis = false,
        markersize=1,
        markerstrokewidth=0,
        title= "doy = " * string(doy))



# +++++++++++++++++++++++++++++++++++++++++++
# Plot emergence for all locations across Ireland
plot(result.DOY, 
     result.emergeDOY, 
     group = result.ID,
     aspect_ratio=:equal,
     lc=:black,
     alpha=RGBA(0,0,0,0.05),
     xaxis = "Start DOY",
     yaxis = "Emergence DOY",
     legend=false)



# Same result visualised as a 2D histogram
histogram2d(result.DOY, 
                result.emergeDOY, 
                nbinsx = 200,
                nbinsy = 100,
                xaxis = "Start DOY",
                yaxis = "Emergence DOY")    