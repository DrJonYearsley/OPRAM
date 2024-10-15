# Compare julia model results withe R calculations by paul
#
#
# Jon Yearsley#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


using DataFrames
using CSV
using Plots
using JLD2

# =================================================================================
# Filename to import

paulFile = "first_occurrence_threshold_with_date_range_2_1961_Agrilus_anxius_GDD_multiple_dates.csv"
juliaFile = "result_agrilus1961_par_thin1.jld2"



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


# ==============================================================================
# Import location data
latlongs = CSV.read(joinpath(meteoDir,"locations.CSV"), DataFrame,
            types=[Int32, Int64, Int64, Float64, Float64])

# Import Paul data
paul = CSV.read(joinpath(outDir, "Changing_Start_Dates",paulFile),DataFrame)

# Import julia data
julia = load(joinpath(outDir, "agrilus",juliaFile), "single_stored_object")


[ latlongs.east[searchsortedfirst(latlongs.ID, id)]for id in eachindex(julia.ID)]

sort!(latlongs, [:ID])
idx = [searchsortedfirst(latlongs.ID, id) for id in eachindex(julia.ID)]

# Remove data not aligned with a lat longs
keep = findall(x->x!=(nrow(latlongs)+1), idx)

julia.east .= 0
julia.north .= 0

julia.east[keep] .= latlongs.east[idx[keep]]
julia.north[keep] = latlongs.north[idx[keep]]

[string(east) * "_" *string(north) for (east,north) in zip(eachindex(julia.east), eachindex(julia.north))]

julia.ID2 = string.(julia.east) .* "_" .*string.(julia.north)
paul.ID2 = string.(paul.east) .* "_" .*string.(paul.north)


# Subset julia data
d = julia[in.(julia.ID2,Ref(unique(paul.ID2))),:]


paul