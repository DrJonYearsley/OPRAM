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

function photoperiod(latlonFile::String, DOY::Vector{Int16}, ID::Vector{Int32})
    # Calculate daylength in hours using the algorithm from the R package geosphere
    # This package uses the algorithm in 
    # Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and 
    # Robert M. Schoolfield, 1995. A model comparison for daylength as a function of 
    # latitude and day of the year. Ecological Modeling 80:87-95.
  
    # Import location data
    latlongs = CSV.read(latlonFile, DataFrame)
    # Sort latlongs
    latlongs = latlongs[sortperm(latlongs.ID),:]
    
  
    lat = [latlongs.latitude[searchsortedfirst(latlongs.ID, ID[i])] for i in eachindex(ID)]
  
    
    daylength = zeros(Float64,length(DOY))
    DOY = convert.(Float64, DOY)     # Convert to Float64 for the calculations
  
    for i in eachindex(DOY) 
      P = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (DOY[i] - 186.0)))))
      a = (sind(0.8333) + sind(lat[i]) * sin(P))/(cosd(lat[i]) * cos(P))
      a = min(max(a, -1), 1)
      daylength[i] = 24.0 - (24.0/pi) * acos(a)
    end
  
    return  daylength
  end
  
  
  # ------------------------------------------------------------------------------------------
  
  
  
  function read_meteo(meteoYear, meteoDir, thinFactor)
    # Import multiple years of daily min and max temperature and 
    # create the daily average temp, and the ID, eastings and northings of spatial locations
    #
    # Arguments:
    #   meteoYear   the starting year to be imported
    #   meteoDIR    the directory containing the data (the daily maximum/minimum temps are
    #               assumed to be in folders maxtemp_grids and mintemp_grids)
    #
    # *************************************************************
  
    # Number of years of meteo data to import   
    nYears = 3
  
    # An array containing the years to be imported
    years = collect(meteoYear:meteoYear+nYears-1)
  
    # Create empty array of arrays to hold average temps 
    TavgVec = Vector{Array{Float32,2}}(undef,length(years)*12);
  
  
    # Read one meteo file to find locations information and calculate a filter
    coordFile = filter(x->occursin("TX_" * string(years[1]) * "01", x), readdir(joinpath(meteoDir,"maxtemp_grids")))
    meteoCoords = CSV.read(joinpath([meteoDir, "maxtemp_grids", coordFile[1]]), DataFrame, select = [1,2], types=Int32);
  
    # Create unique eastings and northings values
    east = sort(unique(meteoCoords.east))
    north = sort(unique(meteoCoords.north))
  
    # Thin the meteo data  using the thinFactor
    thinInd = findall(mod.(meteoCoords.east,(thinFactor*1e3)).< 1e-8 .&&  mod.(meteoCoords.north,(thinFactor*1e3)).< 1e-8) ;
  
  
    # Create ID for each location in meteo file
    idxEast = [searchsortedfirst(east, e) for e in meteoCoords.east[thinInd]]
    idxNorth = [searchsortedfirst(north, n) for n in meteoCoords.north[thinInd]]
  
    ID = [Int32(idxEast[e]*1000 + idxNorth[e]) for e in eachindex(idxEast)]
  
    # Make sure output is in increasing ID order
    sortIdx = sortperm(ID)
  
    # Put filtered locations in ID order
    thinInd = thinInd[sortIdx]
  
    for y in eachindex(years)
        for month in 1:12
            # Import data for month and year (y)
  
            # Get the correct filename
            meteoFileTX = filter(x->occursin("TX_" * string(years[y]) * string(month, pad=2),x), readdir(joinpath(meteoDir,"maxtemp_grids")))
            meteoFileTN = filter(x->occursin("TN_" * string(years[y]) * string(month, pad=2),x), readdir(joinpath(meteoDir,"mintemp_grids")))
    
            # Import the data for min and max daily temperature (removing first two coordinate columns)
            meteoTX = CSV.read(joinpath([meteoDir, "maxtemp_grids", meteoFileTX[1]]), DataFrame, drop=[1,2], types=Float32);
            meteoTN = CSV.read(joinpath([meteoDir, "mintemp_grids", meteoFileTN[1]]), DataFrame, drop=[1,2], types=Float32);
  
            # Calculate daily average temp with locations as columns, days as rows
            TavgVec[(y-1)*12+month] = permutedims(Array{Float32,2}((meteoTX[thinInd,:] .+ meteoTN[thinInd,:])./2.0))
        end
    end
  
    # Put all the temperature data together into one Matrix
    Tavg = reduce(vcat, TavgVec)
  
    # Calculate day of year (used for photoperiod calculations)
    DOM = size.(TavgVec,1)      # Day of month
    DOY = reduce(vcat,[collect(1:sum(DOM[(1+12*(y-1)):12*y])) for y in 1:nYears])  # Day of year
    DOY = convert.(Int16, DOY)
  
    return Tavg, ID[sortIdx], meteoCoords.east[thinInd], meteoCoords.north[thinInd], DOY
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
    # *************************************************************
    
    # Initialise emergeDOY
    emergeDOY = zeros(Int16,length(cumGDD))
  
   
    # Calculate emergence DOY for first starting day
    # (note starting GDD is zero for first day)
    i=1;    # First day of year
    emerge = findfirst(cumGDD  .>= threshold);
  
    
    # Loop around all other starting days until emergence cannot occur
    while !isnothing(emerge)
      # Calculate day of year at emergence
      emergeDOY[i] = cumGDD_doy[emerge];
  
      # Calculate emergence for next start doy
      # (note remove the starting GDD)
      i += 1
      emerge = findfirst(cumGDD  .>= cumGDD[i-1] + threshold);
    end
  
    return(emergeDOY)
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
    #   threshold   an integer threshold for the emergence model
    #
    # *************************************************************
  
  # Loop around all spatial locations (distributed across nodes)
  @sync @distributed for l in eachindex(locInd1)
      # Calculate cumulative degree growing days for one location
      cumGDD = cumsum(GDD[locInd1[l]:locInd2[l]])
      
      # Calculate emergence time
      result[locInd1[l]:locInd2[l],2] .= emergence(cumGDD, result[locInd1[l]:locInd2[l],1], threshold)
      result[locInd1[l]:locInd2[l],3] .= myid()
    end
  
  end
  
  # ================ End of Function Definitions ======================
  # ===================================================================
  