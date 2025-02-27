#!/snap/bin/julia -p 8

#!/opt/homebrew/bin/julia -p 3
#
# Julia script to run the degree day development model for future climate scenarios
# using  Met Eireann's TRANSLATE data (https://www.met.ie/science/translate)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 22nd Jan 2025
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Load packages that will be used
using Distributed;
using Distributions, Random;
using CSV;
using JLD2;
using Dates;



nNodes = 3;                 # Number of compute nodes to use (if in interactive)
meteoPeriod = "2041-2070"   # Years to run model
meteoRCP = "85"             # RCP scenario to use
nReps = 30;                 # Number of times to simulate future climate
maxYears = 3;               # Maximum number of years to complete insect development
saveRepsToFile = false;     # If true save every replicate result to a file
saveRepSummaryToFile = true;# If true save the summary over replicates to a file
gridFile = "IE_grid_locations.csv"  # File containing a 1km grid of lats and longs over Ireland 
# (used for daylength calculations as well as importing and thining of meteo data)

# Factor to think the spatial grid (2 means sample every 2 km, 5 = sample every 5km)
thinFactor = 1;

# Define species parameters
outPrefix = "oulema_melanopus";   # Prefix to use for results files
# Important:
# outPrefix must correspond to part of the variable name 
# for the species parameters. For example, "dummy" if the 
# variable dummy_species is to be used in the model, or "agrius" 
# if the pre-defined parameters for agrilus_anxius are to be used.

# Predefined species are:
#  :agrilus_anxius
#  :halyomorpha_halys
#  :ips_cembrae
#  :ips_duplicatus
#  :ips_sexdentatus
#  :ips_typographus
#  :leptinotarsa_decemlineata
#  :oulema_melanopus
#  :pseudips_mexicanus
#  :spodoptera_frugiperda

# Pre-defined parameters for some species are in species_params.jl
# or you can define your own parameters below in dummy_species
dummy_species = (base_temperature=1.7f0,            # Degrees C
  threshold=1004.0f0,                        # Degrees C
  diapause_photoperiod=missing,              # Hours
  diapause_temperature=missing);             # Degrees C


# =========================================================
# =========================================================

# Specify directories for import and export of data
if isdir("//home//jon//Desktop//OPRAM")
  outDir = "//home//jon//Desktop//OPRAM//results//"
  dataDir = "//home//jon//DATA//OPRAM//"
  meteoDir = "//home//jon//DATA//OPRAM//Climate_JLD2//"

elseif isdir("//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//R")
  outDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//results//"
  dataDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//"
  meteoDir = "//users//jon//Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"

elseif isdir("//users//jon//Desktop//OPRAM//")
  outDir = "//users//jon//Desktop//OPRAM//results//"
  dataDir = "//users//jon//Desktop//OPRAM//"
  meteoDir = "//users//jon//Desktop//OPRAM//Irish_Climate_Data//"
end


# Make directory for the output prefix if one doesn't exist
if !isdir(joinpath(outDir, outPrefix))
  @info "Making directory " * outPrefix
  mkpath(joinpath(outDir, outPrefix))
end
outDir = joinpath(outDir, outPrefix)



# =========================================================
# =========================================================

if isinteractive() & nworkers() == 1
  # Enable multiple nodes 
  addprocs(nNodes)
end


@info "Workers : $(workers())"
@info "Interactive Session: $(isinteractive())"


@everywhere using SharedArrays
@everywhere using DataFrames;



# =========================================================
# =========================================================

# Include the functions to import the data and run the degree day model
include("GDD_functions.jl")
include("import_functions.jl")
include("species_params.jl")




# =========================================================
# =========================================================
# Set up model variables (grid, species parameters and temperature data


# grid, Tavg, params = setup_model


# Find species data corresponding to outPrefix and set this as the variable params
species_params = filter(x -> occursin(outPrefix, string(x)), names(Main))
if length(species_params) > 0
  @info "Using parameters for species " * string(species_params[1])
  params = eval(species_params[1])    # Define parameters to use
else
  @error "Species parameters not found"
end




# =========================================================
# =========================================================
# Import location data and thin it using thinFactor

# Import grid location data
grid = read_grid(joinpath([dataDir, gridFile]), thinFactor);


# =========================================================
# =========================================================
# Import TRANSLATE climate data
Tavg_mean, Tavg_sd, DOY, ID = read_JLD2_translate(meteoDir, meteoRCP, meteoPeriod, grid.ID)

# =========================================================
# =========================================================
# Calculate whether daylength criteria is sufficient for diapause (Temp criteria may also be needed)
# Returns a matrix where rows are unique latitudes and columns are days of year
if !ismissing(params.diapause_photoperiod)
  @info "Calculating DOY threshold for diapause"
  diapause_minDOY = 150
  diapauseDOY_threshold = photoperiod(grid.latitude, diapause_minDOY, params.diapause_photoperiod)
else
  diapauseDOY_threshold = nothing
end


# =========================================================
# =========================================================
# Start the main loop of the program

adult_emerge = Vector{DataFrame}(undef, nReps);

for r in 1:nReps
  @info "====== Replicate " * string(r) * " ==========="

  # Generate daily temperature for maxYears
  @info "Generating climate data"
  TavgVec = Vector{Array{Float32,2}}(undef, maxYears)
  for y in 1:maxYears
    TavgVec[y] = Tavg_mean .+ Tavg_sd .* rand(Normal(0, 1), size(Tavg_sd))
  end
  Tavg = reduce(vcat, TavgVec)
  DOY = collect(1:size(Tavg, 1))

  # Calculate GDD
  @info "Calculating GDD"
  @time GDDsh, idx, locInd1, locInd2 = calculate_GDD_version2(Tavg, DOY, params, diapauseDOY_threshold)

  # Create a shared array to hold results for every day when GDD updates
  result = SharedArray{Int16,2}(length(GDDsh), 3)

  # Fill the first column of results
  result[:, 1] = [idx[i][1] for i in eachindex(idx)]  # Calculate day of year
  result[:, 2] .= Int16(-1)
  result[:, 3] .= Int16(-1)


  # Loop over all locations and run the model
  @info "Running the model"
  @time "Model calculations:" begin

    @everywhere thresh = convert(Float32, $params.threshold)
    location_loop!(locInd1, locInd2, result, GDDsh, thresh)

    # Remove rows that have emergeDOY<=0 
    idxKeep = result[:, 2] .> 0

    # Extract location index for every row in result (column coord in idx)
    loc_idx = [convert(Int32, idx[i][2]) for i in eachindex(idx)]  # ID[loc_idx] will give the location's ID 

    # Remove rows where the final prediction at the same location doesn't change (they can be recalculated later)
    idxKeep2 = (loc_idx[1:end-1] .!= loc_idx[2:end]) .|| (loc_idx[1:end-1] .== loc_idx[2:end] .&& result[1:end-1, 2] .!= result[2:end, 2])
    push!(idxKeep2, true)    # Add a true value at the end

    # Remove all but the first result with DOY >= 365 (results only for 1 year)
    # Keep data with DOY<=365 or when DOY>365 and the previous result was less than 365
    idxKeep3 = result[2:end, 1] .<= 365 .|| (loc_idx[1:end-1] .== loc_idx[2:end] .&& (result[1:end-1, 1] .< 365 .&& result[2:end, 1] .>= 365))
    pushfirst!(idxKeep3, true)    # Add a true value at the start

    # Combine all 3 indices together
    idxKeep = idxKeep .&& idxKeep2 .&& idxKeep3


    # Create a data frame, using real location ID (second element of meteo)
    adult_emerge[r] = DataFrame(rep=r,
                          ID=grid.ID[loc_idx[idxKeep]],
                          DOY=result[idxKeep, 1],
                          emergeDOY=result[idxKeep, 2])

  end
  println(" ================================== \n")
end



# =========================================================
# =========================================================
# Combine the reps into one Data Frame and save to a file

adult_emerge_all = reduce(DataFrames.vcat, adult_emerge)

if saveRepsToFile
  @info "Saving the results"
  # Save using jld2 format
  outFile = joinpath([outDir, "result_" * outPrefix * "_rcp" * meteoRCP * "_" * meteoPeriod * "_par_thin" * string(thinFactor) * ".jld2"])
  save_object(outFile, adult_emerge_all)
end





# =========================================================
# =========================================================
# Convert results into output for key dates (e.g. first day of each month)

# Obtain outputs for starting dates on the first of every month
dates = [Date(0000,m,01) for m in 1:12];

# Work out corresponding day of year
DOY = convert.(Int32,dayofyear.(dates));

# Extract results for each starting date and each location
df_tmp = DataFrame();
for doy in eachindex(DOY)
  @info "Extracting results for DOY " * string(DOY[doy])
    for r in 1:nReps
    # Extract result for the day of year
    tmp = extract_results(DOY[doy], adult_emerge_all[adult_emerge_all.rep.==r, :]);

    # Add in replicate number
    insertcols!(tmp, 1, :rep=> r, after=false)

    # Add these results to the other replicates
    append!(df_tmp, tmp)
  end
end

# Group the df_tmp by ID and startDOY and calculate the percentiles
df_group1 = groupby(df_tmp, [:ID, :startDOY])


# =========================================================
# Calculate 10, 50 and 90 percentiles of results for nGen and emergeDOY
# for each location and starting DOY


# Calculate quantiles across the reps
out_nGen = combine(df_group1,
  :nGen => (x -> [quantile(x, [0.1, 0.5, 0.9])]) =>
    [:nGen_10, :nGen_50, :nGen_90])

out_emergeDOY = combine(df_group1,
  :emergeDOY => (x -> [quantile(x, [0.1, 0.5, 0.9])]) =>
    [:emergeDOY_10, :emergeDOY_50, :emergeDOY_90])

# Combine these results into one data frame and include grid info
a = innerjoin(out_nGen, out_emergeDOY, on=[:ID, :startDOY])
output_1km = rightjoin(grid[:,[:ID,:east,:north]], a, on=:ID)


if saveRepSummaryToFile
  @info "Saving the extracted results"
  # Save using csv format
  outFile = joinpath([outDir, "result_startdates_" * outPrefix * "_rcp" * meteoRCP * "_" * meteoPeriod * "_par_thin" * string(thinFactor) * "_nRep" * string(nReps) * ".csv"])
  CSV.write(outFile, output_1km)
end




# =========================================================
# =========================================================
# Summarise output at the 10km (hectad) scale


# Add hectad info into the df_tmp data frame
eastList = sort(unique(grid.east))
northList = sort(unique(grid.north))
df_tmp.east_hectad = convert.(Int32,floor.(eastList[df_tmp.east_idx]./1e4))
df_tmp.north_hectad = convert.(Int32, floor.(northList[df_tmp.north_idx]./1e4))


# Calculate worst case results within each hectad for nGen and emergeDOY
# for each starting DOY and replicate
df_group2 = groupby(df_tmp, [:east_hectad, :north_hectad, :startDOY, :rep])
df_nGen = combine(df_group2,
  :nGen => (x -> quantile(x, 1.0)) => :nGen_max)    # Maximum num generations per hectad

df_emergeDOY = combine(df_group2,
  :emergeDOY => (x -> quantile(x, 0.0)) => :emergeDOY_min)

# Average over the replicates
out_nGen2 = combine(groupby(df_nGen, [:east_hectad, :north_hectad, :startDOY]),
  :nGen_max => (x -> [quantile(x, [0.1, 0.5, 0.9])]) =>
    [:nGen_max10, :nGen_max50, :nGen_max90])

out_emergeDOY2 = combine(groupby(df_emergeDOY, [:east_hectad, :north_hectad, :startDOY]),
  :emergeDOY_min => (x -> [quantile(x, [0.1, 0.5, 0.9])]) =>
    [:emergeDOY_min10, :emergeDOY_min50, :emergeDOY_min90])


# Create a code that corresponds to hectad in the grid data frame
h_idx = [findfirst(grid.hectad .==h) for h in unique(grid.hectad)]
h_nGen_idx = [findfirst(out_nGen2.east_hectad[i] .== floor.(grid.east[h_idx]./1e4) .&& out_nGen2.north_hectad[i] .== floor.(grid.north[h_idx]./1e4)) for i in 1:nrow(out_nGen2)]
h_emergeDOY_idx = [findfirst(out_emergeDOY2.east_hectad[i] .== floor.(grid.east[h_idx]./1e4) .&& out_emergeDOY2.north_hectad[i] .== floor.(grid.north[h_idx]./1e4)) for i in 1:nrow(out_emergeDOY2)]

out_nGen2.hectad = grid.hectad[h_idx[h_nGen_idx]]
out_emergeDOY2.hectad = grid.hectad[h_idx[h_emergeDOY_idx]]

# Combine these results into one data frame and include grid info
output_10km = innerjoin(grid[h_idx,[:hectad,:country,:province,:county]],
          innerjoin(out_nGen2, out_emergeDOY2, on=[:east_hectad, :north_hectad, :hectad, :startDOY]),
           on=:hectad)

if saveRepSummaryToFile
  @info "Saving the extracted results"
  # Save using csv format
  outFile = joinpath([outDir, "result_startdates_" * outPrefix * "_rcp" * meteoRCP * "_" * meteoPeriod * "_hectads" * "_nRep" * string(nReps) * ".csv"])
  CSV.write(outFile, output_10km)
end




# =========================================================
# =========================================================
# Stop parallel nodes
if isinteractive()
  # Remove parallel nodes
  a = workers()
  rmprocs(a)
end


# ===========  And Breath =================================
# =========================================================