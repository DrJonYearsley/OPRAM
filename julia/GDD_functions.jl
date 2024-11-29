# GDD_finctions.jl
# A file containing the functions to run the degree day development models
# on Met Eireann daily temperature data
#
# This file contains the following functions
# function photoperiod(latlonFile::String, DOY::Vector{Int16}, ID::Vector{Int32})
# function read_meteo(meteoYear, meteoDir, thinFactor)
# function emergence(cumGDD::Vector{Float32}, cumGDD_doy::Vector{Int16}, threshold::Float32)
# function location_loop(locInd1::Vector{Int64}, locInd2::Vector{Int64}, result::SharedMatrix{Int16}, 
#                        GDD::SharedVector{Float32}, threshold::Float32)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================

# =========================================================
# ============= Defne functions ===========================

function photoperiod(latitude::Vector{Float64}, DOY::Vector{Int16}, threshold)
  # Calculate daylength in hours using the algorithm from the R package geosphere
  # This package uses the algorithm in 
  # Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and 
  # Robert M. Schoolfield, 1995. A model comparison for daylength as a function of 
  # latitude and day of the year. Ecological Modeling 80:87-95.

  # Check latitude and doy have the same length
  if length(latitude) != length(DOY) 
    @warn "function photoperiod: latitude and DOY are not the same length"
  end

  daylength = zeros(Float64, length(DOY))        # Calculate day length
  DOY = convert.(Float64, DOY)     # Convert to Float64 for the calculations

  for i in eachindex(DOY)
    P = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (DOY[i] - 186.0)))))
    a = (sind(0.8333) + sind(latitude[i]) * sin(P)) / (cosd(latitude[i]) * cos(P))
    a = min(max(a, -1), 1)
    daylength[i] = 24.0 - (24.0 / pi) * acos(a)
  end

  # Return true if daylength allows development to happen 
  # (photoperiod not relevant before day 150 of year)
  return daylength>threshold || DOY.<150
end

# ------------------------------------------------------------------------------------------



function read_meteo(meteoYear, meteoDir_IE, meteoDir_NI, grid_thin)
  # Import multiple years of daily min and max temperature for Republic of Ireland and 
  # Northern Ireland
  # Then creates the daily average temp for each eastings and northings of spatial locations
  # Locations are also given a unique ID
  #
  # Arguments:
  #   meteoYear   the starting year to be imported
  #   meteoDirIE  the directory containing the data for the Republic of Ireland 
  #               (the daily maximum/minimum temps are assumed to be 
  #                in folders maxtemp_grids and mintemp_grids)
  #   meteoDirNI  the directory containing the data for Northern Ireland 
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
  meteoIE = read_meteoIE(meteoDir_IE, grid_thin, years)
  meteoNI = read_meteoNI(meteoDir_NI, grid_thin, years)

  
  # Find duplicate locations and remove NI data
  NI_keep = [isdisjoint(meteoNI[2][i], meteoIE[2]) for i in eachindex(meteoNI[2])]

  # Combine meteo data from NI and IE
  Tavg = hcat(meteoIE[1], meteoNI[1][:,NI_keep])

  # Combine location data 
  ID = vcat(meteoIE[2],meteoNI[2][NI_keep])

  # Check days of year are the same for both data sets
  if (meteoIE[3]==meteoNI[3])
    local DOY = meteoIE[3]
  else
    @error "IE and NI meteo files have different number of days!"
  end
  
  # Order meteo data by location ID in grid_thin
  sort_idx = sortperm(ID)

  # Return sorted data
  return Tavg[:,sort_idx], DOY, ID[sort_idx]

end
# ------------------------------------------------------------------------------------------



function read_meteoIE(meteoDir_IE, grid_thin, years)
  # Import multiple years of daily min and max temperature for Republic of Ireland and 
  # create the daily average temp, and the ID, eastings and northings of spatial locations
  #
  # Arguments:
  #   meteoDir_IE the directory containing the data for the Republic of Ireland 
  #               (the daily maximum/minimum temps are assumed to be 
  #                in folders maxtemp_grids and mintemp_grids)
  #   grid_thin   the grid of locations to use
  #   years       years to be read
  #
  # *************************************************************

  
  # Create empty array of arrays to hold average temps 
  TavgVec = Vector{Array{Float32,2}}(undef, length(years) * 12)

  # Read one meteo file to find locations information and calculate a filter
  coordFile = filter(x -> occursin("TX_" * string(years[1]) * "01", x), readdir(joinpath(meteoDir_IE, "maxtemp_grids")))
  meteoCoords = CSV.read(joinpath([meteoDir_IE, "maxtemp_grids", coordFile[1]]), DataFrame, select=[1, 2], types=Int32)

  # Find indices of meteo data to use 
  # (needed incase meteo file has missing grid points or different order compared to grid_thin)
  IDidx = [findfirst(abs.(grid_thin.east.-meteoCoords.east[i]).==0 .&& abs.(grid_thin.north.-meteoCoords.north[i]).==0) for i in 1:nrow(meteoCoords)]
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

  return Tavg, grid_thin.ID[IDidx[idx]], DOY
end
# ------------------------------------------------------------------------------------------





function read_meteoNI(meteoDir_NI, grid_thin, years)
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
  # *************************************************************

  # Create empty array of arrays to hold average temps 
  TavgVec = Vector{Array{Float32,2}}(undef, length(years))

  # Read one meteo file to find locations information and calculate a filter
  coordFileNI = filter(x -> occursin("TX_daily_grid_" * string(years[1]), x), 
                     readdir(joinpath(meteoDir_NI, "NI_TX_daily_grid")))
  meteoCoords = CSV.read(joinpath([meteoDir_NI, "NI_TX_daily_grid", coordFileNI[1]]), DataFrame, select=[1, 2], types=Int32)

  # Find indices of meteo data to use
  # (needed incase meteo file has missing grid points or different order compared to grid_thin)
  IDidx = [findfirst(abs.(grid_thin.east.-meteoCoords.east[i]).==0 .&& abs.(grid_thin.north.-meteoCoords.north[i]).==0) for i in 1:nrow(meteoCoords)]
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
  local DOY = reduce(vcat, [collect(1:size(TavgVec[y],1)) for y in 1:length(years)]) 
  DOY = convert.(Int16, DOY)  # Day of year

  return Tavg, grid_thin.ID[IDidx[idx]], DOY
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

# ================ End of Function Definitions ======================
# ===================================================================
