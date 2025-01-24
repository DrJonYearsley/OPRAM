# GDD_finctions.jl
# A file containing the functions to run the degree day development models
# on Met Eireann daily temperature data
#
# This file contains the following functions
# function photoperiod(latlonFile::String, DOY::Vector{Int16}, ID::Vector{Int32})
# function read_meteo(meteoYear::Int64, meteoDirs::Vector{String}, grid_thin::DataFrame)
# function read_JLD2_meteo(meteoDir::String, years::Vector{Int64}, IDgrid::Vector{Int64}, country::String)
# function read_CSV_meteoIE(meteoDir_IE::String, grid_thin::DataFrame, years::Vector{Int64})
# function read_CSV_meteoNI(meteoDir_NI::String, grid_thin::DataFrame, years::Vector{Int64})
# function emergence(cumGDD::Vector{Float32}, cumGDD_doy::Vector{Int16}, threshold::Float32)
# function location_loop(locInd1::Vector{Int64}, locInd2::Vector{Int64}, result::SharedMatrix{Int16}, 
#                        GDD::SharedVector{Float32}, threshold::Float32)
# function calculate_GDD(year::Int64, params::NamedTuple, grid_thin::DataFrame, meteoDirs::Vector)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================

# =========================================================
# ============= Defne functions ===========================

function photoperiod2(latitude::Vector{Float64}, DOY::Vector{Int16}, threshold)
  # Calculate daylength in hours using the algorithm from the R package geosphere
  # This package uses the algorithm in 
  # Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and 
  # Robert M. Schoolfield, 1995. A model comparison for daylength as a function of 
  # latitude and day of the year. Ecological Modeling 80:87-95.

  # Check latitude and doy have the same length
  if length(latitude) != length(DOY)
    @warn "function photoperiod: latitude and DOY are not the same length"
  end

  # daylength only affects overwintering at end of year (i.e entering diapause)
  doy150_idx = findall(DOY .> 150)

  # Return true if daylength allows diapause to happen 
  # (diapause not relevant before day 150 of year)
  overwinterBool = Vector{Bool}(undef, length(DOY))  # True if overwintering 
  for i in eachindex(doy150_idx)
    P = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (convert(Float64, DOY[doy150_idx[i]]) - 186.0)))))
    a = (sind(0.8333) + sind(latitude[doy150_idx[i]]) * sin(P)) / (cosd(latitude[doy150_idx[i]]) * cos(P))
    a = min(max(a, -1), 1)

    # Overwintering if daylength less than threshold     
    overwinterBool[doy150_idx[i]] = 24.0 - (24.0 / pi) * acos(a) < threshold
  end

  # Return true if daylength allows development to happen 
  # (daylength not relevant before day 150 of year)
  return .!overwinterBool
end


function photoperiod(latitude::Vector{Float64}, minDOY::Int64, threshold::Float64)
  # Calculate daylength in hours using the algorithm from the R package geosphere
  # This package uses the algorithm in 
  # Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and 
  # Robert M. Schoolfield, 1995. A model comparison for daylength as a function of 
  # latitude and day of the year. Ecological Modeling 80:87-95.
  #
  # Output:
  #   A vector giving the DOY on which daylength diapause threshold is reached for each latitude
  #   in the input
  #
  # ========================================================================================



  # overwinterBool is true if daylength allows diapause to happen 
  # (other conditions may be needed to start diapause)
  DOY = convert.(Float64,collect(minDOY:366))
  overwinterBool = Array{Bool}(undef, (length(latitude), length(DOY)))  # True if overwintering 

  for i in eachindex(DOY)
    P = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (DOY[i] - 186.0)))))

    a = (sind.(0.8333) .+ sind.(latitude) .* sin(P)) ./ (cosd.(latitude) .* cos(P))
    a = min.(max.(a, -1), 1)
    
    # Overwintering if daylength less than threshold     
    overwinterBool[:, i] = 24.0 .- (24.0 / pi) .* acos.(a) .< threshold
  end

  # Return DOY on which diapause can start for each latitude 
  diapause_DOY = [findfirst(overwinterBool[i, :]) + minDOY - 1 for i in eachindex(latitude)]
  # If diapause threshold cannot be met this will give diapause_DOY = 367 (meaning no diapause)

  return diapause_DOY
end

# ------------------------------------------------------------------------------------------



function read_meteo(meteoYear::Int64, meteoDirs::Vector, grid_thin::DataFrame)
  # Import multiple years of daily min and max temperature for Republic of Ireland and 
  # Northern Ireland
  # Then creates the daily average temp for each eastings and northings of spatial locations
  # Locations are also given a unique ID
  #
  # Arguments:
  #   meteoYear   the starting year to be imported
  #   meteoDirs[1]the directory containing the data for the Republic of Ireland 
  #               (the daily maximum/minimum temps are assumed to be 
  #                in folders maxtemp_grids and mintemp_grids)
  #   meteoDirs[2] the directory containing the data for Northern Ireland 
  #               (the daily maximum/minimum temps are assumed to be 
  #                in folders NI_TX_daily_grid and NI_TN_daily_grid)
  #   grid_thin   the grid of locations to use
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of average daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:   A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Third entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************

  # Number of years of meteo data to import   
  nYears = 3

  # An array containing the years to be imported
  years = collect(meteoYear:meteoYear+nYears-1)



  # Import the weather data

  # ROI data
  if !isnothing(meteoDirs[1])   
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[1]))) > 0
      # Check for jld2 files and import them
      Tavg_IE, DOY_IE, ID_IE = read_JLD2_meteo(meteoDirs[1], years, grid_thin.ID, "IE");
    else
      # Otherwise import csv files
      Tavg_IE, DOY_IE, ID_IE = read_CSV_meteoIE(meteoDirs[1], grid_thin, years);
    end
  end


  # Northern Ireland data
  if !isnothing(meteoDirs[2])
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[2]))) > 0
      # Check for jld2 files and import them
      Tavg_NI, DOY_NI, ID_NI = read_JLD2_meteo(meteoDirs[2], years, grid_thin.ID, "NI");
    else
      # Otherwise import csv files
    Tavg_NI, DOY_NI, ID_NI = read_CSV_meteoNI(meteoDirs[2], grid_thin, years);
    end
  end

  if isnothing(meteoDirs[2])
    # Order meteo data by location ID in grid_thin
    sort_idx = sortperm(ID_IE)

    # Return sorted data
    return Tavg_IE[:, sort_idx], DOY_IE, ID_IE[sort_idx]

  elseif isnothing(meteoDirs[1])
    # Order meteo data by location ID in grid_thin
    sort_idx = sortperm(ID_NI)

    # Return sorted data
    return Tavg_NI[:, sort_idx], DOY_NI, ID_NI[sort_idx]

  elseif isnothing(meteoDirs[1]) & isnothing(meteoDirs[2])
    @error "No meteo directories have been given"

  else
    # Find duplicate locations and remove NI data
    NI_keep = [isdisjoint(ID_NI[i], ID_IE) for i in eachindex(ID_NI)]

    # Combine meteo data from NI and IE
    Tavg = hcat(Tavg_IE, Tavg_NI[:, NI_keep])

    # Combine location data 
    ID = vcat(ID_IE, ID_NI[NI_keep])

    # Check days of year are the same for both data sets
    if (DOY_IE != DOY_NI)
      @error "IE and NI meteo files have different number of days!"
    end

    # Order meteo data by location ID in grid_thin
    sort_idx = sortperm(ID)

    # Return sorted data
    return Tavg[:, sort_idx], DOY_IE, ID[sort_idx]
  end
end


# ------------------------------------------------------------------------------------------

function read_JLD2_meteo(meteoDir::String, years::Vector{Int64}, IDgrid::Vector{Int64}, country::String)
  # Import multiple years of daily min and max temperature from a jld2 file 
  # for Republic of Ireland and create the daily average temp, and 
  # the ID, eastings and northings of spatial locations
  #
  # Arguments:
  #   meteoDir    the directory containing the meteo data in JLD2 format 
  #   years       years to be read
  #   IDgrid      ID of locations in grid_thin
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of average daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:   A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Third entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # Create empty array of arrays to hold average temps 
  TavgVec = Vector{Array{Float32,2}}(undef, length(years))
  DOYVec = Vector{Array{Int16,1}}(undef, length(years))
  IDVec = Vector{Array{Int64,1}}(undef, length(years))

  for y in eachindex(years)
      # Get the correct filename
      meteoFile = filter(x -> occursin("meteo" * country * "_Tavg_" * string(years[y]), x), 
                        readdir(meteoDir))

      # Import the data for Tavg
      f = jldopen(joinpath(meteoDir,meteoFile[1]), "r");
      TavgVec[y] = read(f, "Tavg")
      DOYVec[y] = read(f, "DOY") 
      IDVec[y] = read(f, "ID")
      close(f)
  end

  # Check all ID's are equivalent
  if any(IDVec[1].!=IDVec[2]) || any(IDVec[2].!=IDVec[3])
    @error "Locations in different meteo files are not the same"
  end

  # Find ID's that are in grid_thin
  ID = intersect(IDgrid, IDVec[1])

  # Put all the temperature data together into one Matrix
  Tavg = reduce(vcat, TavgVec)
  DOY = reduce(vcat, DOYVec)

  # Put location ID's (columns of Tavg) in order of ID
  idx = [findfirst(IDVec[1] .== id) for id in ID]
  Tavg = Tavg[:,idx]
  return Tavg, DOY, ID
end

# ------------------------------------------------------------------------------------------



function read_CSV_meteoIE(meteoDir_IE::String, grid_thin::DataFrame, years)
  # Import multiple years of daily min and max temperature from a CSV file
  # for Republic of Ireland and create the daily average temp, and 
  # the ID, eastings and northings of spatial locations
  #
  # Arguments:
  #   meteoDir_IE the directory containing the data for the Republic of Ireland 
  #               (the daily maximum/minimum temps are assumed to be 
  #                in folders maxtemp_grids and mintemp_grids)
  #   grid_thin   the grid of locations to use
  #   years       years to be read (could be a vector or a scalar)
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of average daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:   A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Third entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # Create empty array of arrays to hold average temps 
  TavgVec = Vector{Array{Float32,2}}(undef, length(years) * 12)

  # Read one meteo file to find locations information and calculate a filter
  coordFile = filter(x -> occursin("TX_" * string(years[1]) * "01", x), readdir(joinpath(meteoDir_IE, "maxtemp_grids")))
  meteoCoords = CSV.read(joinpath([meteoDir_IE, "maxtemp_grids", coordFile[1]]), DataFrame, select=[1, 2], types=Int32)

  # Find indices of meteo data to use 
  # (needed in case meteo file has missing grid points or different order compared to grid_thin)
  IDidx = [findfirst(abs.(grid_thin.east .- meteoCoords.east[i]) .== 0 .&& abs.(grid_thin.north .- meteoCoords.north[i]) .== 0) for i in 1:nrow(meteoCoords)]
  idx = findall(IDidx .!= nothing)

  for y in eachindex(years)
    for month in 1:12
      # Import data for month and year (y)

      # Get the correct filename
      meteoFileTX = filter(x -> occursin("TX_" * string(years[y]) * string(month, pad=2), x), readdir(joinpath(meteoDir_IE, "maxtemp_grids")))
      meteoFileTN = filter(x -> occursin("TN_" * string(years[y]) * string(month, pad=2), x), readdir(joinpath(meteoDir_IE, "mintemp_grids")))

      # Import the data for min and max daily temperature (removing first two coordinate columns)
      meteoTX = CSV.read(joinpath([meteoDir_IE, "maxtemp_grids", meteoFileTX[1]]), DataFrame, drop=[1, 2], types=Float32)
      meteoTN = CSV.read(joinpath([meteoDir_IE, "mintemp_grids", meteoFileTN[1]]), DataFrame, drop=[1, 2], types=Float32)

      # Calculate daily average temp with locations as columns, days as rows
      TavgVec[(y-1)*12+month] = permutedims(Array{Float32,2}((meteoTX[idx, :] .+ meteoTN[idx, :]) ./ 2.0))
    end
  end

  # Put all the temperature data together into one Matrix
  Tavg = reduce(vcat, TavgVec)

  # Calculate day of year (used for photoperiod calculations)
  DOM = size.(TavgVec, 1)      # Day of month
  local DOY = reduce(vcat, [collect(1:sum(DOM[(1+12*(y-1)):12*y])) for y in 1:length(years)])  # Day of year
  DOY = convert.(Int16, DOY)

  # Return mean daily temp, day of year and location ID
  return Tavg, DOY, grid_thin.ID[IDidx[idx]]
end

# ------------------------------------------------------------------------------------------

function read_CSV_meteoNI(meteoDir_NI::String, grid_thin::DataFrame, years)
  # Import multiple years of daily mean temperature for Northern Ireland
  # create the daily average temp, and the ID, eastings and northings of spatial locations
  #
  # The data for the entire year are stored in one file
  #
  # Arguments:
  #   meteoDir_NI the directory containing the data for Northern Ireland 
  #               (the daily maximum/minimum temps are assumed to be 
  #                in folders NI_TX_daily_grid and NI_TN_daily_grid)
  #   grid_thin   the grid of locations to use
  #   years       years to be read
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of average daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:   A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Third entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  # *************************************************************

  # Create empty array of arrays to hold average temps 
  TavgVec = Vector{Array{Float32,2}}(undef, length(years))

  # Read one meteo file to find locations information and calculate a filter
  coordFileNI = filter(x -> occursin("TX_daily_grid_" * string(years[1]), x),
    readdir(joinpath(meteoDir_NI, "NI_TX_daily_grid")))
  meteoCoords = CSV.read(joinpath([meteoDir_NI, "NI_TX_daily_grid", coordFileNI[1]]), DataFrame, select=[1, 2], types=Int32)

  # Find indices of meteo data to use
  # (needed incase meteo file has missing grid points or different order compared to grid_thin)
  IDidx = [findfirst(abs.(grid_thin.east .- meteoCoords.east[i]) .== 0 .&& abs.(grid_thin.north .- meteoCoords.north[i]) .== 0) for i in 1:nrow(meteoCoords)]
  idx = findall(IDidx .!= nothing)


  for y in eachindex(years)
    # Import data for year (y)

    # Get the correct filename
    meteoFileTX = filter(x -> occursin("TX_daily_grid_" * string(years[y]), x), readdir(joinpath(meteoDir_NI, "NI_TX_daily_grid")))
    meteoFileTN = filter(x -> occursin("TN_daily_grid_" * string(years[y]), x), readdir(joinpath(meteoDir_NI, "NI_TN_daily_grid")))

    # Import the data for min and max daily temperature (removing first two coordinate columns)
    meteoTX = CSV.read(joinpath([meteoDir_NI, "NI_TX_daily_grid", meteoFileTX[1]]), DataFrame, drop=[1, 2], types=Float32)
    meteoTN = CSV.read(joinpath([meteoDir_NI, "NI_TN_daily_grid", meteoFileTN[1]]), DataFrame, drop=[1, 2], types=Float32)

    # Calculate daily average temp with locations as columns, days as rows
    TavgVec[y] = permutedims(Array{Float32,2}((meteoTX[idx, :] .+ meteoTN[idx, :]) ./ 2.0))
  end

  # Put all the temperature data together into one Matrix
  Tavg = reduce(vcat, TavgVec)


  # Calculate day of year (used for photoperiod calculations) 
  local DOY = reduce(vcat, [collect(1:size(TavgVec[y], 1)) for y in 1:length(years)])
  DOY = convert.(Int16, DOY)  # Day of year

  # Return mean daily temp, day of year and location ID
  return Tavg, DOY, grid_thin.ID[IDidx[idx]]
end

# ------------------------------------------------------------------------------------------

@everywhere function emergence(cumGDD::Vector{Float32}, cumGDD_doy::Vector{Int16}, threshold::Float32)
  # Function to calculate emergence for all possible starting days at one location
  #
  # Arguments:
  #   cumGDD      cumulative degree growing days (only values following an increment are given)
  #   cumGDD_doy  days of year for each cumGDD
  #   threshold   the threshold cumulative GDD for development
  #
  # Output:
  #  emergeDOY    the day of year for adult emergence 
  #               (zero indicates no emergence within the period)
  # *************************************************************

  # Initialise emergeDOY
  emergeDOY = zeros(Int16, length(cumGDD))


  # Calculate emergence DOY for first starting day
  # (note starting GDD is zero for first day)
  i = 1    # First day of year
  emerge = findfirst(cumGDD .>= threshold)


  # Loop around all other starting days until emergence cannot occur
  while !isnothing(emerge)
    # Calculate day of year at emergence
    emergeDOY[i] = cumGDD_doy[emerge]

    # Calculate emergence for next start doy
    # (note remove the starting GDD)
    i += 1
    emerge = findfirst(cumGDD .>= cumGDD[i-1] + threshold)
  end

  return (emergeDOY)
end



# ------------------------------------------------------------------------------------------



function location_loop!(locInd1::Vector{Int64}, locInd2::Vector{Int64}, result::SharedMatrix{Int16},
  GDD::SharedVector{Float32}, threshold::Float32)
  # Function to loop around locations and run development model emergence()
  #
  # Arguments:
  #   locInd1     vector of starting indices for spatial locations
  #   locInd2     vector of ending indices for spatial locations
  #   GDD         average daily temp minus base temp for each spatial location in locs
  #                (each location must be ordered chronologically)
  #   result      a DataFrame to store the results (this is pre-defined for all locations)
  #               this dataframe will be updated (in place) with the results
  #               Col 1: starting day of year, Col 2: emergence day of year, Col 3: cpu number
  #   threshold   an integer threshold for the emergence model
  #
  # *************************************************************

  # Loop around all spatial locations (distributed across nodes)
  @sync @distributed for l in eachindex(locInd1)
    # Calculate cumulative degree growing days for one location
    cumGDD = cumsum(GDD[locInd1[l]:locInd2[l]])

    # Calculate emergence time
    result[locInd1[l]:locInd2[l], 2] .= emergence(cumGDD, result[locInd1[l]:locInd2[l], 1], threshold)
    result[locInd1[l]:locInd2[l], 3] .= myid()
  end

end




# ------------------------------------------------------------------------------------------


function prepare_data(grid_path::String, thin_factor::Float64, params::NamedTuple)
# Function to prepare the data for the degree day model
#
# Arguments:
#   grid_path   path to the file containing the grid of locations
#   thin_factor factor for thining spatial grid 
#      (thin_factor = 10, every 10th grid point, i.e. 10km spacing)
#   params      parameters of the degree day model
#
# Output:
#   GDDsh       shared array of growing degree days (only giving days when 
#               GDD are accumulated)
#   idx         row, column indices of TRUE elements in gdd_update (where GDD accumulates)
#   locInd1     Start index in GDDsh for each spatial location
#   locInd2     End index in GDDsh for each spatial location
# *************************************************************

# =========================================================
# =========================================================
# Import location data and thin it using thinFactor

@info "Importing spatial grid"
# Import grid location data
grid = CSV.read(grid_path, DataFrame);

# Remove locations not in the meteo data
if isnothing(meteoDirs[1])
  subset!(grid, :country => c->c.!="IE")
elseif isnothing(meteoDirs[2])
  subset!(grid, :country => c->c.!="NI")
end

# Sort locations in order of IDs
grid = grid[sortperm(grid.ID), :];

# Thin the locations  using the thinFactor
subset!(grid, :east=> x-> mod.(x, (thinFactor * 1e3)) .< 1e-8,  :north=> x-> mod.(x, (thinFactor * 1e3)) .< 1e-8 )


# =========================================================
# =========================================================
# Calculate whether daylength allows diapause for each DOY and latitude
# Returns a matrix where rows are unique latitudes and columns are days of year
if !ismissing(params.diapause_photoperiod)
  @info "Calculating DOY threshold for diapause"
  diapause_minDOY = 150
  diapauseDOY_threshold = photoperiod(grid.latitude, diapause_minDOY, params.diapause_photoperiod)
else
  diapauseDOY_threshold = nothing
end




  return grid, diapauseDOY_threshold
end

# ------------------------------------------------------------------------------------------


function calculate_GDD(year::Int64, params::NamedTuple, grid_thin::DataFrame, 
                        meteoDirs::Vector, diapauseDOY::Union{Vector{Int64}, Nothing}=nothing)
  # Function to calculate growing degree days from meteo data for a range of years
  #
  # Arguments:
  #   year        starting year for calculations
  #   params      parameters of the degree day model
  #   grid_thin   data frame giving spatial grid of locations
  #   meteoDirs   a 2-element array giving the paths to directories containing
  #               ROI and NI meteo data  (in this order)
  #
  # Output:
  #   GDDsh       shared array of growing degree days (only giving days when 
  #               GDD are accumulated)
  #   idx         row, column indices of TRUE elements in gdd_update (where GDD accumulates)
  #   locInd1     Start index in GDDsh for each spatial location
  #   locInd2     End index in GDDsh for each spatial location
  # *************************************************************


  # Import weather data along with location and day of year
  @info "Calculating for starting year " * string(year)

  @time "Imported meteo data" Tavg, DOY, ID  = read_meteo(year, meteoDirs, grid_thin)

  # Make grid_thin consisent with ID's from meteo
  idx=[id in ID for id in grid_thin.ID]
  grid_thin = grid_thin[idx,:]

  # Calculate days where development updates (Tavg>baseline)
  gdd_update = Tavg .> params.base_temperature

  # Add in diapause (if necessary)
  if !ismissing(params.diapause_photoperiod)
    @info "Removing days when insect is in diapause"
    @time "Diapause" begin
      # Make grid_thin consisent with ID's from meteo
      diapauseDOY = diapauseDOY[idx]

      # Find unique DOY in diapauseDOY
      uniqueDiapause = unique(diapauseDOY) # Unique DOY when diapause starts

      # Create overwinteringBool
      overwinterBool = falses(size(Tavg))  # True if overwintering 

      # Identify when diapuse will occur
      if ismissing(params.diapause_temperature)   # diapause determined only by photoperiod
        # Is DOY past the diapause threshold DOY?
        for i in eachindex(uniqueDiapause)
          locidx = findall(x-> x==uniqueDiapause[i], diapauseDOY)
          DOYidx = findall(x -> x >= uniqueDiapause[i], DOY)

          overwinterBool[DOYidx, locidx] .= true
        end
      else # diapause determined by photoperiod AND temp
        for i in eachindex(uniqueDiapause)
          locidx = findall(x-> x==uniqueDiapause[i], diapauseDOY)
          DOYidx = findall(x -> x >= uniqueDiapause[i], DOY)

          overwinterBool[DOYidx, locidx] = Tavg[DOYidx, locidx] .<= params.diapause_temperature
        end
      end

      # Remove overwintering from gdd_update
      gdd_update[gdd_update .&& overwinteringBool] .= false
    end
  end

  # List ID's for all GDDs above the base temp
  idx = findall(gdd_update)    # Gives row, column coords of non-zero elements in gdd_update
  # IDvec = [convert(Int32, idx[i][2]) for i in eachindex(idx)]  # Location ID index (column coord in idx)

  # Calculate GDD as a shared array
  GDDsh = SharedArray{Float32,1}(Tavg[gdd_update] .- params.base_temperature)


  # Find indices separatng different locations
  ind = vec(sum(gdd_update, dims=1))
  locInd2 = accumulate(+, ind)  # End locations
  locInd1 = vcat(1, locInd2[1:end-1] .+ 1)

  return GDDsh, idx, locInd1, locInd2

end

# ================ End of Function Definitions ======================
# ===================================================================
