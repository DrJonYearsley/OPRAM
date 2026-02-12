# Functions to import and export data for the OPRAM model
#
# function import_parameters(tomlFile::String)
# function import_species(speciesFile::String, speciesName::String)
# function read_meteo(meteoYear::Int64, meteoDirs::Vector{String}, grid_thin::DataFrame, run_params::NamedTuple)
# function read_meteo2(meteoYear::Int64, meteoDirs::Vector, grid_thin::DataFrame, run_params::NamedTuple)
# function read_JLD2_meteo(meteoDir::String, years::Vector{Int64}, IDgrid::Vector{Int64}, country::String)
# function read_JLD2_meteo2(meteoDir::String, years::Vector{Int64}, IDgrid::Vector{Int64},
# function read_JLD2_translate(meteoDir::String, rcp::String, period::String, IDgrid::Vector{Int64})
# function read_JLD2_translate2(meteoDir::String, rcp::String, period::String, IDgrid::Vector{Int64})
# function read_netCDF_meteoUK(meteoDir::String, grid::DataFrame, years)
# function read_CSV_meteoIE(meteoDir_IE::String, grid_thin::DataFrame, years::Vector{Int64})
# function read_CSV_meteoIE2(meteoDir_IE::String, grid_thin::DataFrame, years)
# function read_CSV_meteoNI(meteoDir_NI::String, grid_thin::DataFrame, years::Vector{Int64})
# function read_grid(gridFilePath::String, thinFactor::Int, countryStr::String="IE")
# function read_grid(run_params::NamedTuple)
# function read_OPRAM_JLD2(inFiles::Union{String,Vector{String}}, years::Union{Int,Vector{Int}}, grid::DataFrame)
# function save_OPRAM_1km_CSV(df_1km::DataFrame, species_name::String, yearStr::String, paths::NamedTuple, saveDates::Bool=false)
# function save_OPRAM_10km_CSV(out_10km::DataFrame, species_name::String, yearStr::String, paths::NamedTuple, saveDates::Bool=false)
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================

# Define a structure (with default values) to hold model parameters
@everywhere @kwdef struct parameters
  species_name::Union{String,Missing} = missing
  base_temperature::Float32 = 0.0
  threshold::Float32 = 0.0
  diapause_daylength::Union{Float32,Missing} = missing
  diapause_temperature::Union{Float32,Missing} = missing
end


# =========================================================
# ============= Defne functions ===========================

function import_parameters(tomlFile::String, calculate_average::Bool=false)
  # Import parameters from a TOML file
  # The file should contain a section called "parameters" with the following fields:
  #   species_name, base_temperature, threshold, diapause_daylength, diapause_temperature
  #
  # Returns a parameters object with the values from the TOML file

  # Read TOML parameter file
  params = TOML.parsefile(tomlFile)

  nNodes = params["runtime"]["nNodes"]           # Number of compute nodes to use (if in interactive)


  # =========================================================
  # =========================================================
  # Check country name ======================================

  # Perform check on country name
  if params["model"]["country"] ∉ ["IE", "NI", "AllIreland", "UK", "EN", "SC", "WL"]
    @error "Country must be one of IE, NI, AllIreland, UK, EN, SC or WL"
  end





  # =========================================================
  # =========================================================
  # Check paths =============================================
  # If the given paths for output and data do not exist, try to find alternative paths to the data

  # Add on the home directory to the paths if required
  if !occursin(r"^/", params["paths"]["data"])
    params["paths"]["data"] = joinpath(homedir(), params["paths"]["data"])
  end

  if !occursin(r"^/", params["paths"]["results"])
    params["paths"]["results"] = joinpath(homedir(), params["paths"]["results"])
  end

  if !occursin(r"^/", params["paths"]["output"])
    params["paths"]["output"] = joinpath(homedir(), params["paths"]["output"])
  end

  if !occursin(r"^/", params["paths"]["meteoIE"])
    params["paths"]["meteoIE"] = joinpath(homedir(), params["paths"]["meteoIE"])
  end

  if !occursin(r"^/", params["paths"]["meteoNI"])
    params["paths"]["meteoNI"] = joinpath(homedir(), params["paths"]["meteoNI"])
  end

  # Check output dir exists
  if !isdir(params["paths"]["output"])
    @warn "Can't find output directory!"

    if isdir(joinpath(homedir(), params["paths"]["output"]))  # guess 1
      params["paths"]["output"] = joinpath(homedir(), params["paths"]["output"])
      @info "      Setting output directory to " * params["paths"]["output"]

    elseif isdir(joinpath(homedir(), "Desktop//OPRAM//results"))  # Linux guess
      params["paths"]["output"] = joinpath(homedir(), "Desktop//OPRAM//results")
      @info "      Setting output directory to " * params["paths"]["output"]

    elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//results"))  # Mac guess
      params["paths"]["output"] = joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//results")
      @info "      Setting output directory to " * params["paths"]["output"]

    else
      @error "No output directory found"
    end
  end

  # Check output dir exists
  if !isdir(params["paths"]["results"])
    @warn "Can't find model results directory!"

    if isdir(joinpath(homedir(), params["paths"]["results"]))  # guess
      params["paths"]["results"] = joinpath(homedir(), params["paths"]["results"])
      @info "      Setting output directory to " * params["paths"]["results"]

    elseif isdir(joinpath(homedir(), "Desktop//OPRAM//results"))  # Linux guess
      params["paths"]["results"] = joinpath(homedir(), "Desktop//OPRAM//results")
      @info "      Setting output directory to " * params["paths"]["results"]

    elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//results"))  # Mac guess
      params["paths"]["results"] = joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//results")
      @info "      Setting output directory to " * params["paths"]["results"]

    else
      @error "No results directory found"
    end
  end

  # Check data dir exists
  if !isdir(params["paths"]["data"])
    @info "Can't find data directory!"
    if isdir(joinpath(homedir(), params["paths"]["data"]))  # guess
      params["paths"]["data"] = joinpath(homedir(), params["paths"]["data"])
      @info "       Setting data directory to " * params["paths"]["data"]

    elseif isdir(joinpath(homedir(), "DATA", "OPRAM", "Data"))  # Linux guess
      params["paths"]["data"] = joinpath(homedir(), "DATA//OPRAM//Data")
      @info "       Setting data directory to " * params["paths"]["data"]

    elseif isdir(joinpath(homedir(), "git_repos/OPRAM/data"))  # Mac guess
      @info "       Setting data directory to home directory"
      params["paths"]["data"] = joinpath(homedir(), "git_repos/OPRAM/data")
    else
      @error "No data directory found"
    end
  end



  # Check meteo dirs exist (if simulation needs the data)
  if in(params["model"]["country"], ["IE", "AllIreland"]) &
     !isdir(params["paths"]["meteoIE"])


    @info "Can't find meteoIE directory!"
    if isdir(joinpath(homedir(), "DATA", "OPRAM", "Data", "Climate_JLD2"))  # Linux guess
      params["paths"]["meteoIE"] = joinpath(homedir(), "DATA//OPRAM//Data//Climate_JLD2")
      @info "       Setting meteoIE directory to " * params["paths"]["meteoIE"]

    elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"))  # Mac guess
      params["paths"]["meteoIE"] = joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2")
      @info "       Setting meteoIE directory to " * params["paths"]["meteoIE"]

    else
      @error "No meteoIE directory found"
    end
  end

  if in(params["model"]["country"], ["NI", "AllIreland"]) & !isdir(params["paths"]["meteoNI"])
    @info "Can't find meteoNI directory!"
    if isdir(joinpath(homedir(), "DATA", "OPRAM", "Data", "Climate_JLD2"))  # Linux guess
      params["paths"]["meteoNI"] = joinpath(homedir(), "DATA//OPRAM//Data//Climate_JLD2")
      @info "       Setting meteoNI directory to " * params["paths"]["meteoNI"]

    elseif isdir(joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2"))  # Mac guess
      params["paths"]["meteoNI"] = joinpath(homedir(), "Google Drive//My Drive//Projects//DAFM_OPRAM//Data//Climate_JLD2")
      @info "       Setting meteoNI directory to " * params["paths"]["meteoNI"]

    else
      @error "No meteoNI directory found"
    end
  end




  # Set folders to nothing if no data requried
  if params["model"]["country"] == "IE"
    params["paths"]["meteoNI"] = nothing
  elseif in(params["model"]["country"], ["NI", "EN", "SC", "WL", "UK"])
    params["paths"]["meteoIE"] = nothing
  end


  # Put all the paths together
  paths = (outDir=joinpath(homedir(), params["paths"]["output"]),
    resultsDir=joinpath(homedir(), params["paths"]["results"]),
    dataDir=joinpath(homedir(), params["paths"]["data"]),
    meteoDir_IE=params["paths"]["meteoIE"],
    meteoDir_NI=params["paths"]["meteoNI"])





  # ==================================================================  
  # ==================================================================  
  # Check filenames ==================================================


  # Check grid file exists
  if isfile(params["inputData"]["gridFile"])
    # Do nothing
  elseif isfile(joinpath(params["paths"]["data"], params["inputData"]["gridFile"]))
    params["inputData"]["gridFile"] = joinpath(params["paths"]["data"], params["inputData"]["gridFile"])
  else
    @info "Can't find grid file!"

    if isfile(joinpath(homedir(), params["inputData"]["gridFile"]))  # Guess 1
      params["inputData"]["gridFile"] = joinpath(homedir(), params["inputData"]["gridFile"])
      @info "gridFile set to" * params["inputData"]["gridFile"]

    elseif isfile(joinpath(homedir(), "DATA", "OPRAM", params["inputData"]["gridFile"]))  # Guess 1
      params["inputData"]["gridFile"] = joinpath(homedir(), params["inputData"]["gridFile"])
      @info "gridFile set to" * params["inputData"]["gridFile"]

    else
      @error "No grid file found"
    end
  end



  # Check grid file exists
  if isfile(params["inputData"]["countyDefs"])
    # Do nothing
  elseif isfile(joinpath(homedir(), params["inputData"]["countyDefs"]))
    params["inputData"]["countyDefs"] = joinpath(homedir(), params["inputData"]["countyDefs"])

  elseif isfile(joinpath(params["paths"]["data"], params["inputData"]["countyDefs"]))
    params["inputData"]["countyDefs"] = joinpath(params["paths"]["data"], params["inputData"]["countyDefs"])

  else
    @info "Can't find countyDefs file!"

    if isfile(joinpath(homedir(), params["inputData"]["countyDefs"]))  # Guess 1
      params["inputData"]["countyDefs"] = joinpath(homedir(), params["inputData"]["countyDefs"])
      @info "countyDefs file set to" * params["inputData"]["countyDefs"]

    elseif isfile(joinpath(homedir(), "DATA", "OPRAM", params["inputData"]["countyDefs"]))  # Guess 1
      params["inputData"]["countyDefs"] = joinpath(homedir(), params["inputData"]["countyDefs"])
      @info "countyDefs file set to" * params["inputData"]["countyDefs"]

    else
      @error "No country defs file found"
    end
  end



  # =========================================================
  # =========================================================
  # Put all parameters in named tuple =======================

  thirty_years_str = nothing
  if params["model"]["thirty_years"]!="nothing"
    thirty_years_str = string(params["model"]["thirty_years"][1]) *
                   "_" * string(params["model"]["thirty_years"][2]) 
  end

  # Put parameters into a named tuple
  if in("simYears", keys(params["model"]))
    if calculate_average & !isnothing(thirty_years_str)
      # An average across years is being calculated
      yearVec = collect(params["model"]["thirty_years"][1]:params["model"]["thirty_years"][2])
    else
      # A single year is being calculated (maybe across several years)
      yearVec = params["model"]["simYears"]
    end

    # Run model on past data
    run_params = (
      TRANSLATE_future=false,
      method=params["model"]["method"],          # The degree-day algotrithm to use
      saveJLDFile=params["runtime"]["save2file"],  # If true save the full result to a JLD2 file
      save1year=params["runtime"]["save1year"],    # If true save only the first year's results
      missing2zero=params["runtime"]["missing2zero"], # If true missing meteo data are set to zero
      years=yearVec,                               # Years to run model
      maxYears=params["model"]["maxYears"],        # Maximum number of years to complete insect development
      lastMeteoYear=params["model"]["lastMeteoYear"], # Maximum number of years to complete insect development
      country=params["model"]["country"],          # Can be "IE", "NI" or "AllIreland"
      thinFactor=params["model"]["thinFactor"],    # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
      gridFile=params["inputData"]["gridFile"],    # File containing a 1km grid of lats and longs over Ireland 
      countyFile=params["inputData"]["countyDefs"],# File containing names and codes of counties
      thirty_years= thirty_years_str)               # 30 years to create averaged results


  # (gridfile used for daylength calculations as well as importing and thining of meteo data)
  elseif in("rcp", keys(params["model"]))
    # Run model on future data
    run_params = (
      TRANSLATE_future=true,
      method=params["model"]["method"],         # The degree-day algotrithm to use
      saveJLDFile=params["runtime"]["save2file"], # If true save the full result to a JLD2 file
      save1year=params["runtime"]["save1year"],   # If true save only the first year's results
      missing2zero=params["runtime"]["missing2zero"],# If true missing meteo data are set to zero
      futurePeriod=params["model"]["futurePeriod"], # Years to run model
      rcp=params["model"]["rcp"],                 # Future RCP scenario
      nReps=params["model"]["nReps"],             # Number of replicate climate scenarios
      maxYears=params["model"]["maxYears"],       # Maximum number of years to complete insect development
      country=params["model"]["country"],         # Can be "IE", "NI" or "AllIreland"
      thinFactor=params["model"]["thinFactor"],   # Factor to thin grid (2 = sample every 2 km, 5 = sample every 5km)
      gridFile=params["inputData"]["gridFile"],   # File containing a 1km grid of lats and longs over Ireland 
      countyFile=params["inputData"]["countyDefs"],  # File containing names and codes of counties
      thirty_years= thirty_years_str)               # 30 years to create averaged results

  # (gridfile used for daylength calculations as well as importing and thining of meteo data)
  else
    error("parameter file doesn't contain simulation years")
  end

  # Work out unique species to simulate (Remove obvious duplicates)
  speciesStr = Vector{String}(undef, 0)
  for s in eachindex(params["model"]["speciesList"])
    speciesStr = unique(vcat(speciesStr, params["model"]["speciesList"][s]))
  end
  species_params = (speciesFile=joinpath(homedir(), params["inputData"]["speciesFile"]),  # File containing species parameters
    speciesStr=speciesStr)  # A vector of strings to uniquely identify a species name in the speciesFile


  return nNodes, run_params, species_params, paths
end



# ------------------------------------------------------------------------------------------




function import_species(speciesFile::String, speciesStr::String)
  # Import species parameters from a CSV file and find species 
  # that correspond to the speciesStr

  # Import set species parameters
  params = CSV.read(speciesFile, DataFrame, missingstring="NA",
    types=[String, Float32, Float32, String, Union{Float32,Missing}, Union{Float32,Missing}, String, String, String, String, String])

  # Create a list of species names in lowercase
  spList = lowercase.(params.species)

  if lowercase(speciesStr) == "all"
    # If speciesStr is "all", return all parameters for all species in the file
    return parameters.(params.species, params.baseline_temperature, params.threshold,
      params.diapause_daylength, params.diapause_temperature)
  else
    # Work out which species to return

    # Split the species name into bits based on separators '_', ',',',';',':' and ' '
    spNameBits = split(lowercase(speciesStr), ('_', ' ', ',', ';', ':'))


    # Find the index of the species in the list matching the speciesName
    species_idx = [all(occursin.(spNameBits, s)) for s in spList]


    if sum(species_idx) == 0
      # Can't find specie name in the parameter file
      @warn "Species " * speciesStr * " not found in species file" * speciesFile
      @info "!!!!!!! No parameters set for species " * speciesStr
      return parameters(species_name=missing,)

    elseif sum(species_idx) > 1
      @error "Multiple species found for " * speciesStr * " in file " * speciesFile *
             ". Provide a species string that is more specific."

    else
      # Use the values from the species parameter file
      return parameters(species_name=params.species[findfirst(species_idx)],
        base_temperature=params.baseline_temperature[findfirst(species_idx)],
        threshold=params.threshold[findfirst(species_idx)],
        diapause_daylength=params.diapause_daylength[findfirst(species_idx)],
        diapause_temperature=params.diapause_temperature[findfirst(species_idx)])
    end
  end

end


# ------------------------------------------------------------------------------------------




function read_meteo(meteoYear::Int64, meteoDirs::Vector, grid_thin::DataFrame, run_params::NamedTuple)
  # Import multiple years of daily min and max temperature for Republic of Ireland and 
  # Northern Ireland (and maybe UK)
  # Then creates the daily average temp for each eastings and northings of spatial locations
  # Locations are also given a unique ID
  #
  # Arguments:
  #   meteoYear     the starting year to be imported
  #   meteoDirs[1]  the directory containing the data for the Republic of Ireland 
  #                  (the daily maximum/minimum temps are assumed to be 
  #                   in folders maxtemp_grids and mintemp_grids)
  #   meteoDirs[2]  the directory containing the data for Northern Ireland 
  #                  (the daily maximum/minimum temps are assumed to be 
  #                   in folders NI_TX_daily_grid and NI_TN_daily_grid)
  #   grid_thin     the grid of locations to use
  #   run_params    parameters to define imported data. Use parameters are:
  #                     maxYears      the maximum number of years to be imported
  #                     lastMeteoYear the last year of meteo data available (default=2023)
  #                     country       string defining the country to import
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


  # An array containing the years to be imported
  # Make sure it is no greater than maxMeteoYear
  years = min.(collect(meteoYear:meteoYear+run_params.maxYears-1), run_params.lastMeteoYear)



  # ====================================================================
  # Import the weather data

  # ROI data
  if in(run_params.country, ["IE", "AllIreland"]) & !isnothing(meteoDirs[1])
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[1]))) > 0
      # Check for jld2 files and import them
      Tavg_IE, DOY_IE, ID_IE = read_JLD2_meteo(meteoDirs[1], years, grid_thin.ID, "IE")
    else
      # Otherwise import csv files
      Tavg_IE, DOY_IE, ID_IE = read_CSV_meteoIE(meteoDirs[1], grid_thin, years)
    end
  end

  # Northern Ireland or UK data
  if in(run_params.country, ["NI", "AllIreland"]) & !isnothing(meteoDirs[2])
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[2]))) > 0
      # Check for jld2 files and import them
      Tavg_NI, DOY_NI, ID_NI = read_JLD2_meteo(meteoDirs[2], years, grid_thin.ID, "NI")
    else
      # Otherwise import csv files
      Tavg_NI, DOY_NI, ID_NI = read_CSV_meteoNI(meteoDirs[2], grid_thin, years)
    end
  elseif in(run_params.country, ["UK", "EN", "SC", "WL"]) & !isnothing(meteoDirs[2])
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[2]))) > 0
      # Check for jld2 files and import them
      Tavg_NI, DOY_NI, ID_NI = read_JLD2_meteo(meteoDirs[2], years, grid_thin.ID, run_params.country)
    end
  end


  # =========================================================================================
  # Sort and clean the data

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




function read_meteo2(meteoYear::Int64, meteoDirs::Vector, grid_thin::DataFrame, run_params::NamedTuple)
  # Import multiple years of daily min and max temperature for Republic of Ireland and 
  # Northern Ireland (and maybe UK)
  #
  # Arguments:
  #   meteoYear     the starting year to be imported
  #   meteoDirs[1]  the directory containing the data for the Republic of Ireland 
  #                  (the daily maximum/minimum temps are assumed to be 
  #                   in folders maxtemp_grids and mintemp_grids)
  #   meteoDirs[2]  the directory containing the data for Northern Ireland 
  #                  (the daily maximum/minimum temps are assumed to be 
  #                   in folders NI_TX_daily_grid and NI_TN_daily_grid)
  #   grid_thin     the grid of locations to use
  #   run_params    parameters to define imported data. Use parameters are:
  #                     maxYears      the maximum number of years to be imported
  #                     lastMeteoYear the last year of meteo data available (default=2023)
  #                     country       string defining the country to import
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of minimum daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:    Matrix of maximum daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Third entry:     A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Fourth entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # An array containing the years to be imported
  # Make sure it is no greater than maxMeteoYear
  years = min.(collect(meteoYear:meteoYear+run_params.maxYears-1), run_params.lastMeteoYear)



  # ====================================================================
  # Import the weather data

  # ROI data
  if in(run_params.country, ["IE", "AllIreland"]) & !isnothing(meteoDirs[1])
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[1]))) > 0
      # Check for jld2 files and import them
      Tmin_IE, Tmax_IE, DOY_IE, ID_IE = read_JLD2_meteo2(meteoDirs[1], years, grid_thin.ID, "IE")
    else
      # Otherwise error
      error("No JLD2 meteo files found for IE")
    end
  end

  # Northern Ireland or UK data
  if in(run_params.country, ["NI", "AllIreland"]) & !isnothing(meteoDirs[2])
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[2]))) > 0
      # Check for jld2 files and import them
      Tmin_NI, Tmax_NI, DOY_NI, ID_NI = read_JLD2_meteo2(meteoDirs[2], years, grid_thin.ID, "NI")
    else
      # Otherwise error
      error("No JLD2 meteo files found for NI")
    end
  elseif in(run_params.country, ["UK", "EN", "SC", "WL"]) & !isnothing(meteoDirs[2])
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDirs[2]))) > 0
      # Check for jld2 files and import them
      Tmin_NI, Tmax_NI, DOY_NI, ID_NI = read_JLD2_meteo2(meteoDirs[2], years, grid_thin.ID, run_params.country)
    end
  end


  # =========================================================================================
  # Sort and clean the data

  if isnothing(meteoDirs[2])
    # Order meteo data by location ID in grid_thin
    sort_idx = sortperm(ID_IE)

    # Return sorted data
    return Tmin_IE[:, sort_idx], Tmax_IE[:, sort_idx], DOY_IE, ID_IE[sort_idx]

  elseif isnothing(meteoDirs[1])
    # Order meteo data by location ID in grid_thin
    sort_idx = sortperm(ID_NI)

    # Return sorted data
    return Tmin_NI[:, sort_idx], Tmax_NI[:, sort_idx], DOY_NI, ID_NI[sort_idx]

  elseif isnothing(meteoDirs[1]) & isnothing(meteoDirs[2])
    @error "No meteo directories have been given"

  else
    # Find duplicate locations and remove NI data
    NI_keep = [isdisjoint(ID_NI[i], ID_IE) for i in eachindex(ID_NI)]

    # Combine meteo data from NI and IE
    Tmin = hcat(Tmin_IE, Tmin_NI[:, NI_keep])
    Tmax = hcat(Tmax_IE, Tmax_NI[:, NI_keep])

    # Combine location data 
    ID = vcat(ID_IE, ID_NI[NI_keep])

    # Check days of year are the same for both data sets
    if (DOY_IE != DOY_NI)
      @error "IE and NI meteo files have different number of days!"
    end

    # Order meteo data by location ID in grid_thin
    sort_idx = sortperm(ID)

    # Return sorted data
    return Tmin[:, sort_idx], Tmax[:, sort_idx], DOY_IE, ID[sort_idx]
  end
end



# ------------------------------------------------------------------------------------------

function read_JLD2_meteo(meteoDir::String, years::Vector{Int64}, IDgrid::Vector{Int64},
  country::String)
  # Import multiple years of daily average temperature from a jld2 file 
  # for Republic of Ireland and create the daily average temp, and 
  # the ID, eastings and northings of spatial locations
  #
  # Arguments:
  #   meteoDir    the directory containing the meteo data in JLD2 format 
  #   years       years to be read
  #   IDgrid      ID of locations in grid_thin
  #   country     "IE" for Republic of Ireland, "NI" for Northern Ireland or "UK" for UK
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

  # CHeck to see if we're importing UK met office data or Met Eireann data
  if any(occursin.("meteoUK_Tavg_" * string(years[1]), readdir(meteoDir)))
    countryStr = "UK"
  else
    countryStr = country
  end

  for y in eachindex(years)
    # Get the correct filename

    meteoFile = filter(x -> occursin("meteo" * countryStr * "_Tavg_" * string(years[y]), x),
      readdir(meteoDir))

    # Import the data for Tavg
    f = jldopen(joinpath(meteoDir, meteoFile[1]), "r")
    TavgVec[y] = read(f, "Tavg")
    DOYVec[y] = read(f, "DOY")
    IDVec[y] = read(f, "ID")
    close(f)
  end

  # Check all ID's are equivalent if importing NI and IE data
  if length(IDVec) > 1
    if any(IDVec[1] .!= IDVec[2]) || any(IDVec[2] .!= IDVec[3])
      @error "Locations in different meteo files are not the same"
    end
  end

  # Find ID's that are in grid_thin
  ID = intersect(IDgrid, IDVec[1])

  # Put all the temperature data together into one Matrix
  Tavg = reduce(vcat, TavgVec)
  DOY = reduce(vcat, DOYVec)

  # Put location ID's (columns of Tavg) in order of ID
  # And leave out locations not in ID
  idx = [findfirst(IDVec[1] .== id) for id in ID]
  Tavg = Tavg[:, idx]
  return Tavg, DOY, ID
end

# ------------------------------------------------------------------------------------------


function read_JLD2_meteo2(meteoDir::String, years::Vector{Int64}, IDgrid::Vector{Int64},
  country::String)
  # Import multiple years of daily min and max temperature from a jld2 file 
  # and create the daily min/max temp, the ID, eastings and northings of spatial locations
  #
  # Arguments:
  #   meteoDir    the directory containing the meteo data in JLD2 format 
  #   years       years to be read
  #   IDgrid      ID of locations in grid_thin
  #   country     "IE" for Republic of Ireland, "NI" for Northern Ireland or "UK" for UK
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of minimum daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:    Matrix of maximum daily temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Third entry:     A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Fourth entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # Create empty array of arrays to hold average temps 
  TminVec = Vector{Array{Float32,2}}(undef, length(years))
  TmaxVec = Vector{Array{Float32,2}}(undef, length(years))
  DOYVec = Vector{Array{Int16,1}}(undef, length(years))
  IDVec = Vector{Array{Int64,1}}(undef, length(years))

  # Check to see if we're importing UK met office data or Met Eireann data
  if any(occursin.("meteoUK_Tminmax_" * string(years[1]), readdir(meteoDir)))
    countryStr = "UK"
  else
    countryStr = country
  end

  for y in eachindex(years)
    # Get the correct filename

    meteoFile = filter(x -> occursin("meteo" * countryStr * "_Tminmax_" * string(years[y]), x),
      readdir(meteoDir))

    # Import the data for Tmin and Tmax
    f = jldopen(joinpath(meteoDir, meteoFile[1]), "r")
    TminVec[y] = read(f, "Tmin")
    TmaxVec[y] = read(f, "Tmax")
    DOYVec[y] = read(f, "DOY")
    IDVec[y] = read(f, "ID")
    close(f)
  end

  # Check all ID's are equivalent if importing NI and IE data
  if length(IDVec) > 1
    if any(map(x -> !issetequal(IDVec[1], x), IDVec))
      @error "Locations in different meteo files are not the same"
      println(map(x -> setdiff(IDVec[1], x), IDVec))
    end
  end

  # Find ID's that are in grid_thin
  ID = intersect(IDgrid, IDVec[1])

  # Put all the temperature data together into one Matrix
  Tmin = reduce(vcat, TminVec)
  Tmax = reduce(vcat, TmaxVec)
  DOY = reduce(vcat, DOYVec)

  # Put location ID's (columns of Tavg) in order of ID
  # And leave out locations not in ID
  idx = [findfirst(IDVec[1] .== id) for id in ID]
  Tmin = Tmin[:, idx]
  Tmax = Tmax[:, idx]
  return Tmin, Tmax, DOY, ID
end

# ------------------------------------------------------------------------------------------



function read_JLD2_translate(meteoDir::String, rcp::String, period::String, IDgrid::Vector{Int64})
  # Import the TRANSLATE data for a given RCP and period from a jld2 file
  # and subset to the locations in IDgrid
  #
  # Arguments:
  #   meteoDir    the directory containing the meteo data in JLD2 format 
  #   rcp         the RCP (26, 45, 85)
  #   period      the 30 yr period (2021-2040, 2041-2070)
  #   IDgrid      ID of locations in grid_thin
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of average daily temperature (mean)
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:    Matrix of average daily temperature (standard deviation)
  #                      (Columns are spatial locations, rows are days of year)
  #     Third entry:     A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Fourth entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # Check RCP and period Arguments
  if isnothing(indexin(rcp, ["26", "45", "85"]))
    @error "RCP must be one of: '2.6', '4.5' or '8.5'"
  end

  if isnothing(indexin(period, ["2021-2055", "2041-2070"]))
    @error "Period must be one of: '2021-2055', '2041-2070'"
  end




  # Get the correct filename
  if !isnothing(meteoDir)
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDir))) > 0
      meteoFile = filter(x -> occursin("TRANSLATE_Tavg_rcp" * rcp * "_" * period, x),
        readdir(meteoDir))
    end
  end

  if length(meteoFile) == 0
    @error "No file found for RCP " * rcp * " and period " * period
  end

  # Import the TRANSLATE data
  f = jldopen(joinpath(meteoDir, meteoFile[1]), "r")
  Tavg_mean = read(f, "Tavg_mean")
  Tavg_sd = read(f, "Tavg_sd")
  DOY = read(f, "DOY")
  ID = read(f, "ID")
  close(f)

  # Sort IDgrid
  sort!(IDgrid)

  # Find ID's that are in grid_thin
  IDgrid2 = intersect(IDgrid, ID)

  # Put location ID's (columns of Tavg) in order of ID
  idx = [findfirst(ID .== id) for id in IDgrid2]

  return Tavg_mean[:,idx], Tavg_sd[:,idx], DOY, ID[idx]
end

# ------------------------------------------------------------------------------------------



function read_JLD2_translate2(meteoDir::String, rcp::String, period::String, IDgrid::Vector{Int64})
  # Import the max and min daily temperatures for the TRANSLATE data 
  # for a given RCP and period from a jld2 file
  # and subset to the locations in IDgrid
  #
  # Arguments:
  #   meteoDir    the directory containing the meteo data in JLD2 format 
  #   rcp         the RCP (26, 45, 85)
  #   period      the 30 yr period (2021-2040, 2041-2070)
  #   IDgrid      ID of locations in grid_thin
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of max daily temperature (mean)
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:    Matrix of max daily temperature (standard deviation)
  #                      (Columns are spatial locations, rows are days of year)
  #     Third entry:     Matrix of min daily temperature (mean)
  #                      (Columns are spatial locations, rows are days of year)
  #     Fourth entry:    Matrix of min daily temperature (standard deviation)
  #                      (Columns are spatial locations, rows are days of year)
  #     Fifth entry:     A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Sixth entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # Check RCP and period Arguments
  if isnothing(indexin(rcp, ["26", "45", "85"]))
    @error "RCP must be one of: '2.6', '4.5' or '8.5'"
  end

  if isnothing(indexin(period, ["2021-2055", "2041-2070"]))
    @error "Period must be one of: '2021-2055', '2041-2070'"
  end

  # Get the correct filename
  if !isnothing(meteoDir)
    if length(filter(x -> occursin(".jld2", x), readdir(meteoDir))) > 0
      meteoFile = filter(x -> occursin("TRANSLATE_Tmaxmin_rcp" * rcp * "_" * period, x),
        readdir(meteoDir))
    end
  end

  if length(meteoFile) == 0
    @error "No file found for RCP " * rcp * " and period " * period
  end

  # Import the TRANSLATE data
  f = jldopen(joinpath(meteoDir, meteoFile[1]), "r")
  Tmax_mean = read(f, "Tmax_mean")
  Tmax_sd = read(f, "Tmax_sd")
  Tmin_mean = read(f, "Tmin_mean")
  Tmin_sd = read(f, "Tmin_sd")
  DOY = read(f, "DOY")
  ID = read(f, "ID")
  close(f)

  # Sort IDgrid
  sort!(IDgrid)

  # Find ID's that are in grid_thin
  IDgrid2 = intersect(IDgrid, ID)

  # Put location ID's (columns of Tavg) in order of ID
  idx = [findfirst(ID .== id) for id in IDgrid2]

  return Tmax_mean[:, idx], Tmax_sd[:, idx], Tmin_mean[:, idx], Tmin_sd[:, idx], DOY, ID[idx]
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



function read_CSV_meteoIE2(meteoDir_IE::String, grid_thin::DataFrame, years)
  # Import multiple years of daily min and max temperature from a CSV file
  # for Republic of Ireland and place in a data frame with 
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
  #     First entry:     Matrix of daily min temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:     Matrix of daily max temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Third entry:   A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Fourth entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # Create empty array of arrays to hold daily max and min temps 
  TmaxVec = Vector{Array{Float32,2}}(undef, length(years) * 12)
  TminVec = Vector{Array{Float32,2}}(undef, length(years) * 12)

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
      TmaxVec[(y-1)*12+month] = permutedims(Array{Float32,2}(meteoTX[idx, :]))
      TminVec[(y-1)*12+month] = permutedims(Array{Float32,2}(meteoTN[idx, :]))
    end
  end

  # Put all the temperature data together into one Matrix
  Tmax = reduce(vcat, TmaxVec)
  Tmin = reduce(vcat, TminVec)

  # Calculate day of year (used for photoperiod calculations)
  DOM = size.(TmaxVec, 1)      # Day of month
  local DOY = reduce(vcat, [collect(1:sum(DOM[(1+12*(y-1)):12*y])) for y in 1:length(years)])  # Day of year
  DOY = convert.(Int16, DOY)

  # Return mean daily temp, day of year and location ID
  return Tmin, Tmax, DOY, grid_thin.ID[IDidx[idx]]
end

# ------------------------------------------------------------------------------------------



function read_netCDF_meteoUK(meteoDir::String, grid::DataFrame, years)
  # Import multiple years of daily min and max temperature from a netCDF file
  # for the UK and place in a data frame with the ID, eastings and northings of spatial locations
  #
  # Arguments:
  #   meteoDir   the directory containing the data for the UK 
  #               (the daily maximum/minimum temps are assumed to be 
  #                in folders maxtemp and mintemp)
  #   grid       the grid of locations to use
  #   years      years to be read (could be a vector or a scalar)
  #
  # Output:
  #   A list with three entries
  #     First entry:     Matrix of daily max temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:     Matrix of daily min temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Third entry:   A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Fourth entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  #
  # *************************************************************


  # Create empty array of arrays to hold daily max and min temps 
  TmaxVec = Vector{Array{Float32,2}}(undef, length(years) * 12)
  TminVec = Vector{Array{Float32,2}}(undef, length(years) * 12)

  # Read one meteo file to find locations information and calculate a filter
  ncFile = NCDataset(joinpath(meteoDir, "maxtemp", string(years[1]), "01", "01.nc"), "r")

  # Extract coordinates (remove 500m to put coordinate in bottom left of square)
  east_list = ncFile["projection_x_coordinate"][:] .- 500.0
  north_list = ncFile["projection_y_coordinate"][:] .- 500.0

  # Find indices of data to extract
  idx = [findfirst(grid.east[i] .== east_list) +
         length(east_list) * (findfirst(grid.north[i] .== north_list) - 1) for i in eachindex(grid.ID)]

  close(ncFile)

  for y in eachindex(years)
    for month in 1:12
      # Import data for month and year (y)

      # Get the correct filenames
      meteoFileTX = filter(x -> occursin(".nc", x), readdir(joinpath(meteoDir, "maxtemp", string(years[y]), string(month, pad=2))))
      meteoFileTN = filter(x -> occursin(".nc", x), readdir(joinpath(meteoDir, "mintemp", string(years[y]), string(month, pad=2))))

      # Import the data for min and max daily temperature
      for f in eachindex(meteoFileTX)
        # Import maximum daily temps for one month
        ncFileTX = NCDataset(joinpath(meteoDir, "maxtemp", string(years[y]), string(month, pad=2), meteoFileTX[f]), "r")
        if f == 1
          TmaxVec[(y-1)*12+month] = transpose(convert.(Float32, ncFileTX["daily_maxtemp"][idx]::Array{Float64,1}))
        else
          TmaxVec[(y-1)*12+month] = vcat(TmaxVec[(y-1)*12+month],
            transpose(convert.(Float32, ncFileTX["daily_maxtemp"][idx]::Array{Float64,1})))
        end
        close(ncFileTX)

        # Import minimum daily temps for one month
        ncFileTN = NCDataset(joinpath(meteoDir, "mintemp", string(years[y]), string(month, pad=2), meteoFileTN[f]), "r")
        if f == 1
          TminVec[(y-1)*12+month] = transpose(convert.(Float32, ncFileTN["daily_mintemp"][idx]::Array{Float64,1}))
        else
          TminVec[(y-1)*12+month] = vcat(TminVec[(y-1)*12+month], transpose(convert.(Float32, ncFileTN["daily_mintemp"][idx]::Array{Float64,1})))
        end
        close(ncFileTN)
      end
    end
  end

  # Put all the temperature data together into one Matrix
  Tmax = reduce(vcat, TmaxVec)
  Tmin = reduce(vcat, TminVec)

  # Calculate day of year (used for photoperiod calculations)
  DOM = size.(TmaxVec, 1)      # Day of month
  local DOY = reduce(vcat, [collect(1:sum(DOM[(1+12*(y-1)):12*y])) for y in 1:length(years)])  # Day of year
  DOY = convert.(Int16, DOY)

  # Return mean daily temp, day of year and location ID
  return Tmax, Tmin, DOY, grid.ID
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


function read_CSV_meteoNI2(meteoDir_NI::String, grid_thin::DataFrame, years)
  # Import multiple years of daily mean temperature for Northern Ireland
  # create the daily min and max temp, and the ID, eastings and northings of spatial locations
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
  #     First entry:     Matrix of daily minimum temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Second entry:    Matrix of daily maximum temperature for the period 
  #                      (Columns are spatial locations, rows are days of year)
  #     Third entry:     A vector of days of year (could span several years, length 
  #                      equals the number of rows of temp matrix)
  #     Fourth entry:    A vector of unique location ID's 
  #                      (length equals number of columns temp matrix)
  # *************************************************************

  # Create empty array of arrays to hold average temps 
  TminVec = Vector{Array{Float32,2}}(undef, length(years))
  TmaxVec = Vector{Array{Float32,2}}(undef, length(years))

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
    TminVec[y] = permutedims(Array{Float32,2}(meteoTN[idx, :]))
    TmaxVec[y] = permutedims(Array{Float32,2}(meteoTX[idx, :]))
  end

  # Put all the temperature data together into one Matrix
  Tmin = reduce(vcat, TminVec)
  Tmax = reduce(vcat, TmaxVec)


  # Calculate day of year (used for photoperiod calculations) 
  local DOY = reduce(vcat, [collect(1:size(TminVec[y], 1)) for y in 1:length(years)])
  DOY = convert.(Int16, DOY)  # Day of year

  # Return mean daily temp, day of year and location ID
  return Tmin, Tmax, DOY, grid_thin.ID[IDidx[idx]]
end



# ------------------------------------------------------------------------------------------



function read_grid(gridFilePath::String, thinFactor::Int, countryStr::String="IE")
  # Import the grid of locations and thin it using a thinning factor and the country of interest
  # 
  # Arguments:
  #   dataDir     the directory containing the grid data
  #   gridFile    the name of the file containing the grid data
  #   grid_thin   the thinning factor to use (default=1_)
  #   countryStr  one of "IE", "NI", "AllIreland" (default="IE")
  #
  # Output:
  #   A DataFrame containing the grid of locations
  # *************************************************************

  # Read in the grid data from a file
  grid = CSV.read(gridFilePath, DataFrame, missingstring="NA")

  # Sort locations in order of IDs
  grid = grid[sortperm(grid.ID), :]


  # Keep only locations for the required country
  if countryStr == "IE"
    subset!(grid, :country => c -> c .== "IE")
  elseif countryStr == "NI"
    subset!(grid, :country => c -> c .== "NI")
  end

  # Thin the locations  using the thinFactor
  thinInd = findall(mod.(grid.east, (thinFactor * 1e3)) .< 1e-8 .&&
                    mod.(grid.north, (thinFactor * 1e3)) .< 1e-8)


  return grid[thinInd, :]
end





# ------------------------------------------------------------------------------------------



function read_grid(run_params::NamedTuple)
  # Import the grid of locations and thin it using a thinning factor and the country of interest
  # 
  # Arguments:
  #   run_params   A tuple containing data on:
  #           gridFile     the file containing the grid data
  #           countyFile   the name of the file containing the county definitions
  #           thinFactor   the thinning factor to use (default=1_)
  #           country      one of "IE", "NI", "AllIreland" (default="IE")
  #
  # Output:
  #   A DataFrame containing the grid of locations
  # *************************************************************

  # Import county definitions
  county_defs = CSV.read(run_params.countyFile, DataFrame)

  # Rename county ID column
  rename!(county_defs, :Id => :countyID)

  # Read in the grid data from a file
  grid = CSV.read(run_params.gridFile, DataFrame,
    missingstring="NA",
    types=Dict(:ID => Int64, :east => Int32, :north => Int32,
      :hectad => String, :country => String, :county => String))

  # Sort locations in order of IDs
  grid = grid[sortperm(grid.ID), :]


  # Keep only locations for the required country
  if run_params.country == "UK"
    countryStr = ["NI", "EN", "SC", "WL"]
  elseif run_params.country == "AllIreland"
    countryStr = ["IE", "NI"]
  else
    countryStr = [run_params.country]
  end

  filter!(:country => c -> issubset([c], countryStr), grid)



  # Add in column with county ID
  leftjoin!(grid, county_defs, on=:county => :County)

  # Thin the locations  using the thinFactor
  thinInd = findall(mod.(grid.east, (run_params.thinFactor * 1e3)) .< 1e-8 .&&
                    mod.(grid.north, (run_params.thinFactor * 1e3)) .< 1e-8)

  return grid[thinInd, :]
end





# ------------------------------------------------------------------------------------------


function read_OPRAM_JLD2(inFiles::Union{String,Vector{String}}, years::Union{Int,Vector{Int}}, grid::DataFrame)
  # A function to import the OPRAM results from JLD2 files (as written by the run_model function)
  # and create a single DataFrame with the results
  #
  #  The function creates results for specified starting dates (usually the first of every month)
  #  The function checks to make sure all the required spatial locations are included. 
  #  Missing spatial locations are added as missing data
  #
  # Arguments:
  #   inFiles     the path to the results files
  #   years       the years to import (as a vector of integers)
  #   grid        the grid of locations to use (as a DataFrame)
  #
  # Output:
  #   A DataFrame with the results for the specified species and years
  #   Columns are:
  #        ID            ID for spatial location
  #        startDOY      starting day of year
  #        nGenerations  number of generations
  #        emergeDOY     day of year development is complete
  #        startMonth    the month development started
  # *************************************************************

  # # Check if the resultPath exists
  #   if !isdir(resultPath)
  #     @error "Result path does not exist: $(resultPath)"
  #   end

  if typeof(inFiles) == String
    # If a single file is given, convert it to a vector
    inFiles = [inFiles]
  end

  if typeof(years) == Int32
    # If a single file is given, convert it to a vector
    years = [years]
  end


  # Initialise output DataFrame
  df_1km = DataFrame(ID=Int32[],
    startDOY=Int32[],
    startDate=Date[],
    nGenerations=Vector{Union{Missing,Float64}}(missing, 0),
    emergeDOY=Vector{Union{Missing,Int32}}(missing, 0))

  for y in eachindex(inFiles)

    # Import results
    adult_emerge = load_object(inFiles[y])

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Create output for specific days of year
    dates = [Date(years[y], m, 01) for m in 1:12]         # The first of every month


    # If adult_emerge is a vector then it represents replicate runs
    if isa(adult_emerge, Vector)
      for r in eachindex(adult_emerge)
        @info "---- Generating output for specific starting dates for replicate " * string(r)
        # Create output for specific days of year for each replicate
        append!(df_1km, create_doy_results(dates, adult_emerge[r], grid))
      end
    else
      @info "---- Generating output for specific starting dates for year " * string(years[y])
      append!(df_1km, create_doy_results(dates, adult_emerge, grid))
    end
  end


  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Define starting month for start dates (and set as integer)
  idx = .!ismissing.(df_1km.startDate)
  df_1km.startMonth = Vector{Union{Missing,Int}}(missing, nrow(df_1km))
  df_1km.startMonth[idx] = Dates.month.(df_1km.startDate[idx])  # Use month rather than DOY to avoid leap year problems

  return df_1km

end





# ------------------------------------------------------------------------------------------




function save_OPRAM_1km_CSV(df_1km::DataFrame, species_name::String,
  yearStr::String, paths::NamedTuple, saveDates::Bool=false)
  # A function to save the final OPRAM results (at 1km scale) in a CSV file
  # The function is used in OPRAM_calculate_anomaly.jl
  #
  # Arguments:
  #   filename    the name of the file to save data to (including the path)
  #   df_1km      the DataFrame containing the results
  #   species_name the name of the species to save results for
  #   yearStr     the year string to use in the filename
  #   paths       a named tuple containing the paths to the output directories
  #
  # Output:
  #   None
  # *************************************************************


  # Find unique years and unique county names
  uniqueYears = unique(df_1km.year)
  uniqueCounty = sort(unique(df_1km.countyID))

  if length(uniqueYears) != 1
    @error "Data contains multiple years: $(uniqueYears)"
  end

  for c in uniqueCounty
    # @info "Writing CSV file for county " * string(c)
    # Extract data for the relevant county and year
    out_1km = select(subset(df_1km, :countyID => x1 -> x1 .== c),
      :ID, :east, :north, :startDate,
      [:ID, :east, :north, :startDate,
        :emergeDOY, :emergeDOY_median, :emergeDOY_anomaly,
        :nGenerations, :nGenerations_median, :nGenerations_anomaly])


    if saveDates
      # Add in Dates for emergence
      out_1km.emergeDate = Vector{Union{Missing,Date}}(missing, nrow(out_1km))
      idx = .!ismissing.(out_1km.startDate) .&& .!ismissing.(out_1km.emergeDOY)
      out_1km.emergeDate[idx] = Date.(year.(out_1km.startDate[idx]), 01, 01) .+ Day.(out_1km.emergeDOY[idx] .- 1)

      out_1km.emergeDate_30yr = Vector{Union{Missing,Date}}(missing, nrow(out_1km))
      idx = .!ismissing.(out_1km.startDate) .&& .!ismissing.(out_1km.emergeDOY_median)
      out_1km.emergeDate_30yr[idx] = Date.(year.(out_1km.startDate[idx]), 01, 01) .+ Day.(round.(Int, out_1km.emergeDOY_median[idx]) .- 1)
    end

    # Round number of generations
    out_1km.nGenerations = round.(out_1km.nGenerations, digits=2)
    out_1km.nGenerations_median = round.(out_1km.nGenerations_median, digits=2)
    out_1km.nGenerations_anomaly = round.(out_1km.nGenerations_anomaly, digits=2)

    # Round median and anomaly emergence DOY to 1 decimal place
    out_1km.emergeDOY_median = round.(out_1km.emergeDOY_median, digits=1)
    out_1km.emergeDOY_anomaly = round.(out_1km.emergeDOY_anomaly, digits=1)


    # Rename some columns
    rename!(out_1km,
      :east => "eastings",
      :north => "northings",
      :nGenerations => "nGen",
      :nGenerations_median => "nGen_30yr",
      :nGenerations_anomaly => "nGen_anomaly",
      :emergeDOY_median => "emergeDOY_30yr")

    # Create filename
    fileout1 = joinpath(paths.outDir, species_name, species_name * "_1_" * string(c) * "_" * yearStr * ".csv")

    CSV.write(fileout1, out_1km, missingstring="NA", dateformat="yyyy-mm-dd")

  end

  return nothing
end









# ------------------------------------------------------------------------------------------




function save_OPRAM_10km_CSV(out_10km::DataFrame, species_name::String,
  yearStr::String, paths::NamedTuple, saveDates::Bool=false)
  # A function to save the final OPRAM results (at 1km scale) in a CSV file
  # The function is used in OPRAM_calculate_anomaly.jl
  #
  # Arguments:
  #   out_10km      the DataFrame containing the results
  #   species_name the name of the species to save results for
  #   yearStr     the year string to use in the filename
  #   paths       a named tuple containing the paths to the output directories
  #   saveDates   whether to save emergence dates (true/false)
  #
  # Output:
  #   None
  # *************************************************************





  if saveDates
    # Add in Dates for emergence
    out_10km.emergeDate_min = Vector{Union{Missing,Date}}(missing, nrow(out_10km))
    idx = .!ismissing.(out_10km.startDate) .&& .!ismissing.(out_10km.emergeDOY_min)
    out_10km.emergeDate_min[idx] = Date.(year.(out_10km.startDate[idx]), 01, 01) .+ Day.(out_10km.emergeDOY_min[idx] .- 1)

    out_10km.emergeDate_median_min = Vector{Union{Missing,Date}}(missing, nrow(out_10km))
    idx = .!ismissing.(out_10km.startDate) .&& .!ismissing.(out_10km.emergeDOY_median_min)
    out_10km.emergeDate_median_min[idx] = Date.(year.(out_10km.startDate[idx]), 01, 01) .+ Day.(round.(Int, out_10km.emergeDOY_median_min[idx]) .- 1)

    # Rename some columns
    rename!(out_10km, :nGenerations_max => "nGen",
      :nGenerations_median_max => "nGen_30yr",
      :nGenerations_anomaly_median => "nGen_anomaly",
      :emergeDate_min => "emergeDate",
      :emergeDate_median_min => "emergeDate_30yr",
      :emergeDOY_min => "emergeDOY",
      :emergeDOY_median_min => "emergeDOY_30yr",
      :emergeDOY_anomaly_median => "emergeDOY_anomaly")

    # Select specific years
    select!(out_10km, [:hectad, :startDate,
      :emergeDOY, :emergeDOY_30yr, :emergeDOY_anomaly,
      :emergeDate, :emergeDate_30yr,
      :nGen, :nGen_30yr, :nGen_anomaly])

  else

    # Rename some columns
    rename!(out_10km, :nGenerations_max => "nGen",
      :nGenerations_median_max => "nGen_30yr",
      :nGenerations_anomaly_median => "nGen_anomaly",
      :emergeDOY_min => "emergeDOY",
      :emergeDOY_median_min => "emergeDOY_30yr",
      :emergeDOY_anomaly_median => "emergeDOY_anomaly")

    # Select specific years
    select!(out_10km, [:hectad, :startDate,
      :emergeDOY, :emergeDOY_30yr, :emergeDOY_anomaly,
      :nGen, :nGen_30yr, :nGen_anomaly])
  end


  # Round results for number of generations
  out_10km.nGen = round.(out_10km.nGen, digits=2)
  out_10km.nGen_30yr = round.(out_10km.nGen_30yr, digits=2)
  out_10km.nGen_anomaly = round.(out_10km.nGen_anomaly, digits=2)

  # Round emergence DOY to 1 decimal place
  out_10km.emergeDOY_30yr = round.(out_10km.emergeDOY_30yr, digits=1)
  out_10km.emergeDOY_anomaly = round.(out_10km.emergeDOY_anomaly, digits=1)


  # Create filename
  fileout3 = joinpath(paths.outDir, species_name, species_name * "_100_" * yearStr * ".csv")

  # Write the data
  CSV.write(fileout3, out_10km, missingstring="NA", dateformat="yyyy-mm-dd")

  return nothing
end



# ================ End of Function Definitions ======================
# ===================================================================
