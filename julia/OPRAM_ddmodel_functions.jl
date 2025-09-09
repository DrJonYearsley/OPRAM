# OPRAM_ddmodel_functions.jl
# A file containing the functions to run the degree day development models
# on Met Eireann daily temperature data
#
# This file contains the following functions
# function run_model(run_params::NamedTuple, species_params::NamedTuple, paths::NamedTuple)
# function run_model_futures(run_params::NamedTuple, species_params::NamedTuple, paths::NamedTuple)
# function location_loop(locInd1::Vector{Int64}, locInd2::Vector{Int64},idx::Vector{CartesianIndex{2}}, 
#                        GDD::SharedVector{Float32}, threshold::Float32)
# function emergence(cumGDD::Vector{Float32}, cumGDD_doy::Vector{Int16}, threshold::Float32)
# function calculate_GDD(Tavg::Matrix{Float32}, grid::DataFrame, DOY::Vector{Int16}, params::parameters)
# function photoperiod(latlonFile::String, DOY::Vector{Int16}, ID::Vector{Int32})
# function extract_results(doy::Int32, result::DataFrame)
# function create_doy_results(adult_emerge::DataFrame, DOY::Vector{Int32})
# function aggregate_to_hectad(result_1km::DataFrame, grid::DataFrame)

#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================

# =========================================================
# ============= Defne functions ===========================


function run_model(run_params::NamedTuple, species_params::NamedTuple, paths::NamedTuple)
  # Function to run the degree day model for a range of years and for a range of species
  #
  # Arguments:
  #   run_params      a named tuple containing the parameters for the model
  #   species_params   a named tuple containing the parameters for the species
  #   paths           a named tuple containing the paths to the data
  #
  # Output:
  #   a data frame containing the adult emergence dates that is saved in JLD2 format
  #                   Columns are: 
  #                     ID         unique ID of spatial location, 
  #                     DOY        starting day of year for development
  #                     emergeDOY  day of year for adult emergence
  #
  #   two CSV files containing the adult emergence dates at specific starting dates
  #
  # ====================================================================
  # ====================================================================


  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Set species parameters from the parameter file (can be more than one species)
  if species_params.speciesStr[1] == "all"
    @info "Importing all species from the species file" * string(species_params.speciesStr)
    species_params = import_species(species_params.speciesFile, species_params.speciesStr[1])
  else
    @info "Importing species from the species file: " * string(species_params.speciesStr)
    species_params = [import_species(species_params.speciesFile, species_params.speciesStr[s]) for s in eachindex(species_params.speciesStr)]
  end

  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Import location data and thin it using thinFactor
  grid = read_grid(run_params)



  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Start the main loop of the program
  for y in eachindex(run_params.years)

    # Import meteo data and calculate degree days
    @info "======Running model for year " * string(run_params.years[y]) * "============="

    if run_params.years[y] > run_params.lastMeteoYear
      @error "Year " * string(run_params.years[y]) * " is greater than the last year of meteo data. Skipping this year."
    end

    @info "Importing " * string(run_params.maxYears) * " years of meteo data, starting at year " * string(run_params.years[y])
    Tavg, DOY, ID = read_meteo(run_params.years[y], [paths.meteoDir_IE, paths.meteoDir_NI], grid, run_params)

    # Remove grid points not in the meteo data
    keep_ID = [in(grid.ID[i], ID) for i in eachindex(grid.ID)]
    grid_final = grid[keep_ID, :]

    # Loop through all the species
    for s in eachindex(species_params)
      if ismissing(species_params[s].species_name)
        @warn "Skipping an undefined species"
      else
        @info "\n#######  Running model for species " * species_params[s].species_name * " #######"


        @info "        Calculate GDD"
        GDDsh, idx, locInd1, locInd2 = calculate_GDD(Tavg, grid_final, DOY, species_params[s])

        # Calculate model at each location
        @info "        Calculating adult emergence dates"
        result = location_loop(locInd1, locInd2, idx, GDDsh, species_params[s].threshold)


        # Specify the output directory using the species name
        outPrefix = replace(lowercase(species_params[s].species_name), " " => "_")

        # Make directory for the output prefix if one doesn't exist
        if !isdir(joinpath(paths.outDir, outPrefix))
          @info "        Making directory " * outPrefix
          mkpath(joinpath(paths.outDir, outPrefix))
        end

        if run_params.saveJLDFile
          @info "        Saving the results to JLD2 file"

          # Simplify the results and put them in a DataFrame
          adult_emerge = cleanup_results(result, idx, grid_final.ID, run_params.save1year)

          # Save using jld2 format
          outFile = joinpath([paths.outDir, outPrefix,
            outPrefix * "_" * run_params.country * "_" * string(run_params.years[y]) * "_" * string(run_params.thinFactor) * "km.jld2"])
          save_object(outFile, adult_emerge)
        end


        # # =========================================================
        # # =========================================================
        # # Create summary outputs
        # if run_params.saveSummaryCSV
        #   @info "        Creating 1km CSV summary for specific starting dates"
        #   # Starting dates for output in CSV files
        #   # The first of every month
        #   dates = [Date(run_params.years[y], m, 01) for m in 1:12]


        #   save_to_csv(dates, adult_emerge, grid_final, outPrefix, paths, run_params)
        # end

        # Clear memory of model results before starting next species
        # NOTE: Don't clear the CLimate data because it will be reused
        adult_emerge = nothing
        result = nothing
        @everywhere GDDsh = nothing   # Make sure the shared array is cleared everywhere
        @everywhere GC.gc()           # Clean up memory (this call may not be needed)
      end
    end
    println(" ")    # Print a blank line
  end
end


# ------------------------------------------------------------------------------------------


function run_model_futures(run_params::NamedTuple, species_params::NamedTuple, paths::NamedTuple)
  # Function to run the degree day model for a range of future climate predictions and for a range of species
  #
  # The differences to run_model are:
  #      + imports data from the TRANSLATE project which gives mean and standard deviation of daily temps
  #      + generates many yearly weather scenario replicates (e.g. 50) based upon TRANSLATE data
  #      
  #
  # Arguments:
  #   run_params      a named tuple containing the parameters for the model
  #   species_params   a named tuple containing the parameters for the species
  #   paths           a named tuple containing the paths to the data
  #
  # Output:
  #   a data frame containing the adult emergence dates that is saved in JLD2 format
  #                   Columns are: 
  #                     ID         unique ID of spatial location, 
  #                     DOY        starting day of year for development
  #                     emergeDOY  day of year for adult emergence
  #
  #   two CSV files containing the adult emergence dates at specific starting dates
  #
  # ====================================================================
  # ====================================================================


  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Set species parameters from the parameter file (can be more than one species)
  if species_params.speciesStr[1] == "all"
    @info "Importing all species from the species file" * string(species_params.speciesStr)
    species_params = import_species(species_params.speciesFile, species_params.speciesStr[1])
  else
    @info "Importing species from the species file" * string(species_params.speciesStr)
    species_params = [import_species(species_params.speciesFile, species_params.speciesStr[s]) for s in eachindex(species_params.speciesStr)]
  end

  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Import location data and thin it using thinFactor
  grid = read_grid(run_params)


  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Import TRANSLATE climate data

  for r in eachindex(run_params.rcp)
    @info "======Running model for RCP " * run_params.rcp[r] * "============="

    for p in eachindex(run_params.futurePeriod)
      @info "======Running model for future period " * run_params.futurePeriod[p] * "============="

      Tavg_mean, Tavg_sd, DOY, ID = read_JLD2_translate(paths.meteoDir_IE, run_params.rcp[r], run_params.futurePeriod[p], grid.ID)

      if run_params.nReps==1
        # If only one replicate then use the mean temperature
        Tavg_sd .= 0
      end

      # Remove grid points not in the meteo data
      keep_ID = [in(grid.ID[i], ID) for i in eachindex(grid.ID)]
      grid_final = grid[keep_ID, :]


      # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Start the main loop of the program

      # Loop through all the species
      for s in eachindex(species_params)
        if ismissing(species_params[s].species_name)
          @warn "Skipping an undefined species"
        else
          @info "\n#######  Running model for species " * species_params[s].species_name * " #######"

          # Initialise data frame
          adult_emerge = Vector{DataFrame}(undef, run_params.nReps)

          for r in 1:run_params.nReps
            @info "====== Replicate " * string(r) * " ==========="

            # Generate daily temperature for maxYears
            @info "        Generating climate data"
            TavgVec = Vector{Array{Float32,2}}(undef, run_params.maxYears)
            for y in 1:run_params.maxYears
              TavgVec[y] = Tavg_mean .+ Tavg_sd .* rand(Normal(0, 1), size(Tavg_sd))
            end
            Tavg = reduce(vcat, TavgVec)
            DOY = convert.(Int16, collect(1:size(Tavg, 1)))

            @info "        Calculate GDD"
            GDDsh, idx, locInd1, locInd2 = calculate_GDD(Tavg, grid_final, DOY, species_params[s])

            # Calculate model at each location
            @info "        Calculating adult emergence dates"
            result = location_loop(locInd1, locInd2, idx, GDDsh, species_params[s].threshold)

            # Simplify the results and put them in a DataFrame
            adult_emerge[r] = cleanup_results(result, idx, grid_final.ID, run_params.save1year)


            # Add in a column for the replicate number
            insertcols!(adult_emerge[r], 1, :rep => Int8(r))

            result = nothing
            @everywhere GDDsh = nothing   # Make sure the shared array is cleared everywhere
            @everywhere GC.gc()           # Clean up memory (this call may not be needed)
          end


          # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          # ========== Save output to a JLD2 file ====================================
          # Specify the output directory using the species name
          outPrefix = replace(lowercase(species_params[s].species_name), " " => "_")

          # Make directory for the output prefix if one doesn't exist
          if !isdir(joinpath(paths.outDir, outPrefix))
            @info "        Making directory " * outPrefix
            mkpath(joinpath(paths.outDir, outPrefix))
          end

          if run_params.saveJLDFile
            @info "        Saving the results to JLD2 file"
            # Save using jld2 format
            outFile = joinpath([paths.outDir, outPrefix,
              outPrefix * "_" * run_params.country * "_rcp" * run_params.rcp[r] * "_" *
              run_params.futurePeriod[p] * "_" * string(run_params.thinFactor) * "km.jld2"])
            save_object(outFile, adult_emerge)
          end


          # # Clear memory of model results before starting next species
          # # NOTE: Don't clear the CLimate data because it will be reused
          # adult_emerge = nothing
        end
        println(" ")    # Print a blank line
      end
    end
  end
end





# ------------------------------------------------------------------------------------------





function location_loop(locInd1::Vector{Int64}, locInd2::Vector{Int64},
  idx::Vector{CartesianIndex{2}},
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

  # Create a shared array to hold results for every day when GDD updates
  result = SharedArray{Int16,2}(length(GDD), 3)

  # Fill the first column of results
  result[:, 1] = [idx[i][1] for i in eachindex(idx)]  # Calculate day of year
  result[:, 2] .= Int16(-1)    # Stores emergence DOY
  result[:, 3] .= Int16(-1)    # Stores worker number 

  # Loop around all spatial locations (distributed across nodes)
  # Don't move on until all calculations are complete
  @sync @distributed for l in eachindex(locInd1)
    # Calculate cumulative degree growing days for one location
    cumGDD = cumsum(GDD[locInd1[l]:locInd2[l]])

    # Calculate emergence time
    result[locInd1[l]:locInd2[l], 2] .= emergence(cumGDD, result[locInd1[l]:locInd2[l], 1], threshold)
    result[locInd1[l]:locInd2[l], 3] .= myid()    # Record ID of the worker node
  end

  return result
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



function calculate_GDD(Tavg::Matrix{Float32}, grid::DataFrame, DOY::Vector{Int16}, params::parameters)
  # Function to calculate growing degree days from meteo data for a range of years
  #
  # Arguments:
  #   Tavg        Temperature data
  #   grid        a data frame containing the spatial grid
  #   DOY         Day of year corresponding to rows in Tavg
  #   params      parameters of the degree day model
  #
  # Output:
  #   GDDsh       shared array of growing degree days (only giving days when 
  #               GDD are accumulated) [Float32]
  #   idx         row, column indices of TRUE elements in gdd_update (where GDD accumulates)
  #   locInd1     Start index in GDDsh for each spatial location
  #   locInd2     End index in GDDsh for each spatial location
  # *************************************************************


  # Calculate days where development updates (Tavg>baseline)
  gdd_update = Tavg .> params.base_temperature


  # =========================================================
  # =========================================================
  # Add in diapause (if necessary)
  if !ismissing(params.diapause_daylength)
    @info "        Diapause: Removing GDD for days when insect is in diapause"
    # Create overwinteringBool
    overwinterBool = falses(size(Tavg))  # True if overwintering 

    # Calculate whether daylength allows diapause for each DOY and latitude
    # Returns a matrix where rows are unique latitudes and columns are days of year
    diapause_minDOY = 150
    diapauseDOY = photoperiod(grid.latitude, diapause_minDOY, params.diapause_daylength)

    # Find unique DOY in diapauseDOY
    uniqueDiapause = unique(diapauseDOY) # Unique DOY when diapause starts

    # Identify when diapause will occur
    if ismissing(params.diapause_temperature)   # diapause determined only by photoperiod
      # Is DOY past the diapause threshold DOY?
      for i in eachindex(uniqueDiapause)
        locidx = findall(x -> x == uniqueDiapause[i], diapauseDOY)  # Find uniqueDiapause index for each location
        DOYidx = findall(x -> x >= uniqueDiapause[i], DOY) # Find DOY indices >= uniqueDiapause

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

  # ================= End of diapause code ==================
  # =========================================================


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



function calculate_GDD2(Tmin::Matrix{Float32}, Tmax::Matrix{Float32}, grid::DataFrame, 
                        DOY::Vector{Int16}, params::parameters)
  # Function to calculate growing degree days from meteo data for a range of years
  # Several algorithms for calculating degree days can be used
  #
  # Arguments:
  #   Tmin        Minimum daily temperature data
  #   Tmax        Maximum daily temperature data
  #   grid        a data frame containing the spatial grid
  #   DOY         Day of year corresponding to rows in Tavg
  #   params      parameters of the degree day model
  #
  # Output:
  #   GDDsh       shared array of growing degree days (only giving days when 
  #               GDD are accumulated) [Float32]
  #   idx         row, column indices of TRUE elements in gdd_update (where GDD accumulates)
  #   locInd1     Start index in GDDsh for each spatial location
  #   locInd2     End index in GDDsh for each spatial location
  # *************************************************************


  # Calculate days where development updates (Tavg>baseline)
  gdd_update = Tavg .> params.base_temperature


  # =========================================================
  # =========================================================
  # Add in diapause (if necessary)
  if !ismissing(params.diapause_daylength)
    @info "        Diapause: Removing GDD for days when insect is in diapause"
    # Create overwinteringBool
    overwinterBool = falses(size(Tavg))  # True if overwintering 

    # Calculate whether daylength allows diapause for each DOY and latitude
    # Returns a matrix where rows are unique latitudes and columns are days of year
    diapause_minDOY = 150
    diapauseDOY = photoperiod(grid.latitude, diapause_minDOY, params.diapause_daylength)

    # Find unique DOY in diapauseDOY
    uniqueDiapause = unique(diapauseDOY) # Unique DOY when diapause starts

    # Identify when diapause will occur
    if ismissing(params.diapause_temperature)   # diapause determined only by photoperiod
      # Is DOY past the diapause threshold DOY?
      for i in eachindex(uniqueDiapause)
        locidx = findall(x -> x == uniqueDiapause[i], diapauseDOY)  # Find uniqueDiapause index for each location
        DOYidx = findall(x -> x >= uniqueDiapause[i], DOY) # Find DOY indices >= uniqueDiapause

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

  # ================= End of diapause code ==================
  # =========================================================


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

function average_degreeday(Tmin, Tmax, T_baseline)
# Calculate degree day from daily max and min temperatures using the daily average method
# This method tends to underestimate degree-days
# 
# Arguments:
#   Tmin         minimum daily temperature
#   Tmax         maximum daily temperature
#   T_baseline   baseline temperature above which degree days start accumulating
#
# Output:
#   DD           degree days
#   DOY          the day of year corresponding to DD
#
#
# ========================================================================================

end


# ------------------------------------------------------------------------------------------

function triangle_degreeday(Tmin, Tmax, T_baseline)
# Calculate degree day from daily max and min temperatures using the single traingle method
# 
# Arguments:
#   Tmin         minimum daily temperature
#   Tmax         maximum daily temperature
#   T_baseline   baseline temperature above which degree days start accumulating
#
# Output:
#   DD           degree days
#   DOY          the day of year corresponding to DD
#
#
# ========================================================================================

end



# ------------------------------------------------------------------------------------------

function sine_degreeday(Tmin, Tmax, T_baseline)
# Calculate degree day from daily max and min temperatures using the single sine method
# 
# Arguments:
#   Tmin         minimum daily temperature
#   Tmax         maximum daily temperature
#   T_baseline   baseline temperature above which degree days start accumulating
#
# Output:
#   DD           degree days
#   DOY          the day of year corresponding to DD
#
#
# ========================================================================================



end


# ------------------------------------------------------------------------------------------


function photoperiod(latitude::Vector{Float64}, minDOY::Int64, threshold::Float32)
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
function save_to_csv(dates::Vector{Date}, adult_emerge::DataFrame, grid::DataFrame,
  outPrefix::String, paths::NamedTuple, run_params::NamedTuple)
  # Function to take raw results and save them to a csv file for specific starting dates
  #
  # Arguments:
  #   dates         a vector of starting dates to extract results for
  #   adult_emerge  a data frame containing the adult emergence dates
  #   grid          a data frame containing the spatial grid
  #   outPrefix     a string to add to the start of the output file
  #   paths         a named tuple containing the paths to the data
  #
  # Output:
  #   Two csv files containing the adult emergence dates at 
  #      the native spatial resolution
  #
  # ====================================================================
  # ====================================================================


  # Create output for specific days of year
  output_1km = create_doy_results(dates, adult_emerge)

  # Select columns to save in CSV file
  select!(output_1km, [:ID, :startDate, :emergeDate, :nGenerations])

  # Set the digits for nGenerations to 2 decimal places 
  # Note: must be done after creating hectad summary
  output_1km.nGenerations = round.(output_1km.nGenerations, digits=2)

  # Add eastings and northings to 1km output
  output_1km = rightjoin(grid[:, [:ID, :east, :north]], output_1km, on=:ID)

  # Extract the year from the dates
  yearStr = string(Dates.year(dates[1]))

  # Save using csv format
  outFile1km = joinpath([paths.outDir, outPrefix,
    "result_" * run_params.country * "_" * outPrefix * "_" * yearStr * "_" * string(run_params.thinFactor) * "km.csv"])
  CSV.write(outFile1km, output_1km, missingstring="NA")
end


# ------------------------------------------------------------------------------------------

function cleanup_results(result::SharedMatrix{Int16}, idx::Vector{CartesianIndex{2}},
  ID::Vector{Int64}, save1year::Bool)
  # A function to removing missing or unneeded results from the degree day model
  #
  # Arguments:
  #   result       A dataframe of results generated by location_loop()
  #   idx          row, column indices of TRUE elements in gdd_update (where GDD accumulates)
  #   ID           IDs of spatial locations
  #   save1year    If true only keep results for the first year
  #
  # Outputs:
  #   A data frame with a reduced number of rows that contains only necessary results
  #
  ###################################################################################

  # Index to remove rows that have emergeDOY<=0 
  idxKeep1 = result[:, 2] .> 0

  # Extract location index for every row in result (column coord in idx)
  loc_idx = [convert(Int32, idx[i][2]) for i in eachindex(idx)]
  # ID[loc_idx] will give the location's ID 

  # Index to remove rows where the final prediction at the same location 
  # doesn't change (they can be recalculated later)
  idxKeep2 = (loc_idx[1:end-1] .!= loc_idx[2:end]) .||
             (loc_idx[1:end-1] .== loc_idx[2:end] .&&
              result[1:end-1, 2] .!= result[2:end, 2])
  push!(idxKeep2, true)    # Add a true value at the end


  # Combine idxKeep1 and idxKeep2
  idxKeep = idxKeep1 .&& idxKeep2


  if save1year
    # Remove all but the first result with DOY >= 365 (results only for 1 year)
    # Keep data with DOY<=365 or when DOY>365 and the previous result was less than 365
    idxKeep3 = result[2:end, 1] .<= 365 .||
               (loc_idx[1:end-1] .== loc_idx[2:end] .&&
                (result[1:end-1, 1] .< 365 .&& result[2:end, 1] .>= 365))
    pushfirst!(idxKeep3, true)    # Add a true value at the start

    # Update idxKeep
    idxKeep = idxKeep .&& idxKeep3
  end



  # Create a data frame, using real location ID (second element of meteo)
  return DataFrame(ID=ID[loc_idx[idxKeep]],
    DOY=result[idxKeep, 1],
    emergeDOY=result[idxKeep, 2])

end





# ------------------------------------------------------------------------------------------




# function extract_results(doy::Int32, result::DataFrame)
#   # Function to calculate number of generations per year, and first
#   # day of adult emergence based on a starting development on doy
#   #
#   #   If development doesn't complete by the end of the simulation (usually 3 years)
#   #  then the emergeDOY and nGenerations is set to missing
#   #
#   #   doy    day of year when larval development begins
#   #   result data frame containing the results from the model
#   #  
#   #  Output:
#   #   a data frame with the following columns:
#   #       ID            unique ID for each spatial location   
#   #       startDOY      day of year to start development     
#   #       nGenerations  the number of generations within the year
#   #       emergeDOY     adult emergence day of year 
#   #       east_idx      index for eastings of the unique spatial locations in result
#   #       north_idx     index for northings of the unique spatial location in result
#   ########################################################################

#   # Find start and end indicies for each location
#   #@time idx1 = [searchsortedfirst(result.ID, loc) for loc in unique(result.ID)]
#   # @time idx2 = [searchsortedlast(result.ID, loc) for loc in unique(result.ID)]
#   idx1 = indexin(unique(result.ID), result.ID)
#   idx2 = vcat(idx1[2:end] .- 1, length(result.ID))


#   # Create data frame to hold results
#   out_res = DataFrame(ID=result.ID[idx1],
#     startDOY=doy,
#     nGenerations=0.0,      # Must be Float because it can be fractional 
#     emergeDOY=Vector{Union{Missing,Int32}}(missing, length(idx1)))

#   out_res.north_idx = mod.(out_res.ID, 1000)
#   out_res.east_idx = div.((out_res.ID .- out_res.north_idx), 1000)

#   allowmissing!(out_res, :nGenerations)

#   # Set up starting DOY for each location
#   startDOY = zeros(Union{Missing,Int32}, length(idx1)) .+ doy  # Start DOY for each location
#   first_pass = true                                            # True for the first generation

#   # Set startDOY to zero for a location if no more generations can be completed in the year
#   # when all startDOY are zero then stop
#   while any(startDOY .> 0)
#     # Find index of result output that corresponds to desired day of year (ie doy>=result.DOY)
#     # this assumes results are ordered by DOY for each location
#     idx3 = [searchsortedfirst(result.DOY[idx1[i]:idx2[i]], startDOY[i]) - 1 + idx1[i] for i in eachindex(idx1)]

#     # Find out if end of current developmental generation is in the results (i.e. occurs before data on next location)
#     # idx2[i] is the last index for the present location
#     # idx2[i]+1 is the first index of the next location   = idx1[i+1]
#     development_complete = startDOY .> 0 .&& idx2 .+ 1 .> idx3

#     # Find out if development happens within the first year
#     if (idx3[end] == nrow(result) + 1)
#       # If final idx3 gives no developement (i.e. idx3[end] is nrow+1) then set its within_a_year to false
#       within_a_year = result.emergeDOY[idx3[1:end-1]] .<= 365   # Is the generation complete within the year?
#       push!(within_a_year, false)
#     else
#       # Work out all values of within_a_year
#       within_a_year = result.emergeDOY[idx3] .<= 365   # Is the generation complete within the year?
#     end

#     # Create booleans that will increment generation count
#     full_generation = development_complete .& within_a_year
#     partial_generation = development_complete .& .!within_a_year

#     # Add on one full generation
#     out_res.nGenerations[full_generation] .+= 1

#     # Calculate end of year fraction of time towards next generation
#     out_res.nGenerations[partial_generation] .+= (365 .- startDOY[partial_generation]) ./
#                                                  (result.emergeDOY[idx3[partial_generation]] - startDOY[partial_generation])

#     # If this is the first time through the loop, record the emergance day
#     if first_pass
#       # Only update adult emergence DOY if emergence is complete (but could take multiple years)   
#       out_res.emergeDOY[development_complete] .= result.emergeDOY[idx3[development_complete]]
#       first_pass = false
#     end

#     # Reset starting DOY to zero
#     startDOY = zeros(Int32, length(idx1))

#     # If a full generation completed within the year, update starting doy to be the
#     # day after generation completed
#     startDOY[full_generation] .= result.emergeDOY[idx3[full_generation]] .+ 1

#     # Remove startDOY >= 365
#     startDOY[startDOY.>=365] .= 0
#   end

#   # Set nGenerations to missing if no emergence occurs within the simulation
#   out_res.nGenerations[ismissing.(out_res.emergeDOY)] .= missing

#   return (out_res)
# end





# # ------------------------------------------------------------------------------------------

# function create_doy_results(dates::Vector{Date}, adult_emerge::DataFrame)
#   # Function to calculate number of generations per year, and first
#   # day of adult emergence for specific starting days of year.
#   # If emergence date exceeds simulation timeframe then set to missing
#   #
#   #   dates           a vector of dates when larval development begins
#   #   run_params      a tuple of runtime parameters
#   #   adult_emerge    data frame containing the results from the model
#   #  
#   #  Output:
#   #   a data frame with the following columns:
#   #       ID            unique ID for each spatial location  
#   #       startDOY      starting DOY for larval development  
#   #       startDate     starting date for larval development  
#   #       nGenerations  number of generations in a year    
#   #       emergeDOY     adult emergence day of year 
#   #       emergeDate    adult emergence date 
#   #       east_idx      index for eastings of the unique spatial locations in result
#   #       north_idx     index for northings of the unique spatial location in result
#   ########################################################################

#   out = DataFrame()    # Initialise dataframe for outputs


#   for d in eachindex(dates)
#     # Produce results for a specific DOY
#     result = extract_results(convert(Int32, dayofyear(dates[d])), adult_emerge)

#     # Add in a column with starting dates and emergence dates rather than DOY
#     result.startDate .= dates[d]

#     # Initialise column with missing values
#     result.emergeDate = Vector{Union{Missing,Date}}(missing, length(result.ID))

#     # Find results with a non-zero emergence DOY
#     idx = .!ismissing.(result.emergeDOY)
#     # Calculate date of emergence from day of year
#     result.emergeDate[idx] = Date(year(dates[d])) .+ Day.(result.emergeDOY[idx] .- 1)

#     # Add these results to the other outputs
#     append!(out, result)
#   end

#   return out
# end




# # ------------------------------------------------------------------------------------------



# function aggregate_to_hectad(result_1km::DataFrame, grid::DataFrame)
#   # Function to take results on a 1km grid an aggregate them to a 10km grid
#   #
#   # Arguments:
#   #   grid              information about the spatial grid
#   #   result_1km        results from the model on 1km grid
#   #  
#   #  Output:
#   #   a data frame with the following columns:
#   #       hectad            unique hectad grid reference
#   #       east              easting of bottom left of hectad
#   #       north             northing of bottom left of hectad
#   #       startDOY          starting DOY for larval development  
#   #       startDate         starting date for larval development  
#   #       nGenerations_max  maximum number of generations in a year in hectad
#   #       emergeDOY_min     minimum adult emergence day of year in a hectad
#   #       emergeDate_min    minimum adult emergence date in a hectad
#   ########################################################################

#   # Add bottom left eastings and northings of a hectad into result_1km
#   eastList = sort(unique(grid.east))
#   northList = sort(unique(grid.north))
#   result_1km.east_hectad = convert.(Int32, floor.(eastList[result_1km.east_idx] ./ 1e4) .* 1e4)
#   result_1km.north_hectad = convert.(Int32, floor.(northList[result_1km.north_idx] ./ 1e4) .* 1e4)

#   # Group data frame by location and then starting DOY/Date 
#   df_group = groupby(result_1km, [:east_hectad, :north_hectad, :startDOY, :startDate])

#   # Calculate worst case results within each hectad for nGenerations and emergeDOY
#   # for each starting DOY. 
#   # These worst case scenarios are not affected by locations where no emergence occurred
#   df_agg = combine(df_group,
#     :nGenerations => (x -> quantile(x, 1.0)) => :nGenerations_max,
#     :emergeDOY => (x -> quantile(x, 0.0)) => :emergeDOY_min)   # Max generations per hectad

#   # df_emergeDOY = combine(df_group,
#   #   :emergeDOY => (x -> quantile(x, 0.0)) => :emergeDOY_min)         # Minimum emergence date

#   # Add in columns with start and emergence dates
#   df_agg.emergeDate_min = Vector{Union{Missing,Date}}(missing, nrow(df_agg))
#   idx = df_agg.emergeDOY_min .> 0
#   df_agg.emergeDate_min[idx] .= Date.(year.(df_agg.startDate[idx])) .+ Day.(floor.(df_agg.emergeDOY_min[idx]))


#   # Create a code that corresponds to hectad in the grid data frame
#   h_idx = [findfirst(grid.hectad .== h) for h in unique(grid.hectad)]
#   h_nGen_idx = [findfirst(df_agg.east_hectad[i] .== floor.(grid.east[h_idx] ./ 1e4) .* 1e4 .&& df_agg.north_hectad[i] .== floor.(grid.north[h_idx] ./ 1e4) .* 1e4) for i in 1:nrow(df_agg)]
#   # h_emergeDOY_idx = [findfirst(df_emergeDOY.east_hectad[i] .== floor.(grid.east[h_idx] ./ 1e4) .* 1e4 .&& df_emergeDOY.north_hectad[i] .== floor.(grid.north[h_idx] ./ 1e4) .* 1e4) for i in 1:nrow(df_emergeDOY)]

#   # Add hectad code into both data frames
#   insertcols!(df_agg, 1, :hectad => grid.hectad[h_idx[h_nGen_idx]], after=false)
#   # insertcols!(df_emergeDOY, 1, :hectad => grid.hectad[h_idx[h_emergeDOY_idx]], after=false)


#   # Combine these results into one data frame and include grid info
#   # return innerjoin(df_nGen, df_emergeDOY,
#   # on=[:hectad, :east_hectad, :north_hectad, :startDOY, :startDate])

#   # Return final dataframe
#   return df_agg

# end


# ================ End of Function Definitions ======================
# ===================================================================


