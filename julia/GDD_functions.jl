# GDD_finctions.jl
# A file containing the functions to run the degree day development models
# on Met Eireann daily temperature data
#
# This file contains the following functions
# function photoperiod(latlonFile::String, DOY::Vector{Int16}, ID::Vector{Int32})
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
  DOY = convert.(Float64, collect(minDOY:366))
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
  grid = CSV.read(grid_path, DataFrame)

  # Remove locations not in the meteo data
  if isnothing(meteoDirs[1])
    subset!(grid, :country => c -> c .!= "IE")
  elseif isnothing(meteoDirs[2])
    subset!(grid, :country => c -> c .!= "NI")
  end

  # Sort locations in order of IDs
  grid = grid[sortperm(grid.ID), :]

  # Thin the locations  using the thinFactor
  subset!(grid, :east => x -> mod.(x, (thinFactor * 1e3)) .< 1e-8, :north => x -> mod.(x, (thinFactor * 1e3)) .< 1e-8)


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
  meteoDirs::Vector, diapauseDOY::Union{Vector{Int64},Nothing}=nothing)
  # Function to calculate growing degree days from meteo data for a range of years
  #
  # Arguments:
  #   year        starting year for calculations
  #   params      parameters of the degree day model
  #   grid_thin   data frame giving spatial grid of locations
  #   meteoDirs   a 2-element array giving the paths to directories containing
  #               ROI and NI meteo data  (in this order)
  #   diapauseDOY the day of year when diapause will start (nothing if no diapause)
  #
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

  @time "Imported meteo data" Tavg, DOY, ID = read_meteo(year, meteoDirs, grid_thin)

  # Make grid_thin consisent with ID's from meteo
  idx = [id in ID for id in grid_thin.ID]
  grid_thin = grid_thin[idx, :]

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
          locidx = findall(x -> x == uniqueDiapause[i], diapauseDOY)
          DOYidx = findall(x -> x >= uniqueDiapause[i], DOY)

          overwinterBool[DOYidx, locidx] .= true
        end
      else # diapause determined by photoperiod AND temp
        for i in eachindex(uniqueDiapause)
          locidx = findall(x -> x == uniqueDiapause[i], diapauseDOY)
          DOYidx = findall(x -> x >= uniqueDiapause[i], DOY)

          overwinterBool[DOYidx, locidx] = Tavg[DOYidx, locidx] .<= params.diapause_temperature
        end
      end

      # Remove overwintering from gdd_update
      gdd_update[gdd_update.&&overwinterBool] .= false
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






# ------------------------------------------------------------------------------------------



function calculate_GDD_version2(Tavg::Matrix{Float32}, DOY::Vector{Int64},
  params::NamedTuple, diapauseDOY::Union{Vector{Int64},Nothing}=nothing)
  # Function to calculate growing degree days from meteo data for a range of years
  #
  # Arguments:
  #   Tavg        Temperature data
  #   DOY         Day of year corresponding to rows in Tavg
  #   params      parameters of the degree day model
  #   diapauseDOY the day of year when diapause will start (nothing if no diapause)
  #
  # Output:
  #   GDDsh       shared array of growing degree days (only giving days when 
  #               GDD are accumulated)
  #   idx         row, column indices of TRUE elements in gdd_update (where GDD accumulates)
  #   locInd1     Start index in GDDsh for each spatial location
  #   locInd2     End index in GDDsh for each spatial location
  # *************************************************************


  # Calculate days where development updates (Tavg>baseline)
  gdd_update = Tavg .> params.base_temperature

  # Add in diapause (if necessary)
  if !ismissing(params.diapause_photoperiod)
    @info "Removing days when insect is in diapause"
    @time "Diapause" begin

      # Find unique DOY in diapauseDOY
      uniqueDiapause = unique(diapauseDOY) # Unique DOY when diapause starts

      # Create overwinteringBool
      overwinterBool = falses(size(Tavg))  # True if overwintering 

      # Identify when diapuse will occur
      if ismissing(params.diapause_temperature)   # diapause determined only by photoperiod
        # Is DOY past the diapause threshold DOY?
        for i in eachindex(uniqueDiapause)
          locidx = findall(x -> x == uniqueDiapause[i], diapauseDOY)
          DOYidx = findall(x -> x >= uniqueDiapause[i], DOY)

          overwinterBool[DOYidx, locidx] .= true
        end
      else # diapause determined by photoperiod AND temp
        for i in eachindex(uniqueDiapause)
          # Is DOY past the diapause threshold DOY?
          locidx = findall(x -> x == uniqueDiapause[i], diapauseDOY)
          DOYidx = findall(x -> x >= uniqueDiapause[i], DOY)

          # Is Tavg past the diapause threshold temperature?
          overwinterBool[DOYidx, locidx] = Tavg[DOYidx, locidx] .<= params.diapause_temperature
        end
      end

      # Remove overwintering from gdd_update
      gdd_update[gdd_update.&&overwinterBool] .= false
    end
  end

  # List ID's for all GDDs above the base temp
  idx = findall(gdd_update)    # Gives row, column coords of non-zero elements in gdd_update
  
  # Calculate GDD as a shared array
  GDDsh = SharedArray{Float32,1}(Tavg[gdd_update] .- params.base_temperature)


  # Find indices separatng different locations
  ind = vec(sum(gdd_update, dims=1))
  locInd2 = accumulate(+, ind)  # End locations
  locInd1 = vcat(1, locInd2[1:end-1] .+ 1)

  return GDDsh, idx, locInd1, locInd2

end







# ------------------------------------------------------------------------------------------




function extract_results(doy::Int32, result::DataFrame)
        # Function to calculate number of generations per year, and first
        # day of adult emergence based on a starting development on doy
        #
        #   doy    day of year when larval development begins
        #   result data frame containing the results from the model
        #  
        #  Output:
        #   a data frame with the following columns:
        #       ID            unique ID for each soatial location        
        #       emergeDOY     adult emergence day of year 
        #       nGen          the number of generations within the year
        #       east          index for eastings of the unique spatial locations in result
        #       north         index for northings of the unique spatial location in result
        ########################################################################

        # Find start and end indicies for each location
        # @time idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)]
        # @time idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)]
        @time idx1 = indexin(unique(result.ID), result.ID)
        @time idx2 = vcat(idx1[2:end] .- 1, length(result.ID))


        # Create data frame to hold number of generations
        out_res = DataFrame(ID=result.ID[idx1], nGen=0.0, emergeDOY=0)
        out_res.north = mod.(out_res.ID, 1000)
        out_res.east = div.((out_res.ID .- out_res.north), 1000)

        # Count number of generations per year
        startDOY = zeros(Int32, length(idx1)) .+ doy
        first_pass = true

        # Set startDOY to zero for a location if no more generations can be completed in the year
        while any(startDOY .> 0)
                # Find index of result that corresponds to desired day of year (ie doy>=result.DOY)
                idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], startDOY[i]) - 1 + idx1[i] for i = eachindex(idx1)]

                println(maximum(idx3))

                # Find out if end of current developmental generation is in the results (i.e. occurs before data on next location)
                development_complete = startDOY .> 0 .&& idx2 .+ 1 .> idx3

                # Find out if development happens within the first year
                if (idx3[end] == nrow(result) + 1)
                        # If final entry resulted in no developement (i.e. index is nrow+1) then set within year to false
                        within_a_year = result.emergeDOY[idx3[1:end-1]] .<= 365   # Is the generation complete within the year?
                        push!(within_a_year, false)
                else
                        within_a_year = result.emergeDOY[idx3] .<= 365   # Is the generation complete within the year?
                end

                # Increment generations
                full_generation = development_complete .& within_a_year
                partial_generation = development_complete .& .!within_a_year

                # Add on one full generation
                out_res.nGen[full_generation] .+= 1

                # Calculate end of year fraction of time towards next generation
                out_res.nGen[partial_generation] .+= (365 .- startDOY[partial_generation]) ./
                                                     (result.emergeDOY[idx3[partial_generation]] - startDOY[partial_generation])

                # If this is the first time through the loop, record the emergance day
                if first_pass
                        out_res.emergeDOY[full_generation] .= result.emergeDOY[idx3[full_generation]]
                        first_pass = false
                end

                # Reset starting DOY to zero
                startDOY = zeros(Int32, length(idx1))

                # If a full gneration completed within the year, update starting doy
                startDOY[full_generation] .= result.emergeDOY[idx3[full_generation]] .+ 1

                # Remove startDOY >= 365
                startDOY[startDOY.>=365] .= 0

                # println([sum(startDOY .> 0), maximum(result.emergeDOY[idx3])])
        end


        return (out_res)
end


# ================ End of Function Definitions ======================
# ===================================================================





