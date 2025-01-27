# File to calculate day length (hours) across a set of locations 
# and days of year. The results are saved to a CSV file
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# Jan 2025
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Directory and file containing grid of locations
dataDir = "//Users//jon//git_repos//OPRAM//data"
gridFile = "FR_grid_locations.csv"

# Days of year for day length calculation
DOY = Vector{Int32}(1:366)

# Name of file for output
outFile = "daylength_FRANCE.csv"


# Load packages
using CSV;
using DataFrames;
using Plots;

# =========================================================
# ============= Defne functions ===========================

function daylength(latitude::Vector{Float64}, DOY::Vector{Int32})
  # Calculate day length in hours using the algorithm from the R package geosphere
  # This package uses the algorithm in 
  # Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and 
  # Robert M. Schoolfield, 1995. A model comparison for daylength as a function of 
  # latitude and day of the year. Ecological Modeling 80:87-95.
  #
  # Arguments:
  #   latitude   An array of latitudes  for which day length is to be calculated
  #   DOY        An array of days of year for which day length is to be calculated
  #
  # Output:
  #   A matrix giving the day length for each latitude and day of year 
  #             Number of rows is the number of unique latitudes, 
  #             Number of columns is the number of "days of year"
  #
  # ========================================================================================


  daylength = Array{Float64}(undef, (length(latitude), length(DOY)))

  for i in eachindex(DOY)
    P = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (convert(Float64, DOY[i]) - 186.0)))))

    a = (sind.(0.8333) .+ sind.(latitude) .* sin(P)) ./ (cosd.(latitude) .* cos(P))
    a = min.(max.(a, -1), 1)

    # Overwintering if daylength less than threshold     
    daylength[:, i] = 24.0 .- (24.0 / pi) .* acos.(a)
  end

  return daylength
end

# ------------------------------------------------------------------------------------------




# Load grid locations
grid = CSV.read(joinpath(dataDir, gridFile), DataFrame);

# Find unique latitudes
latitude_list = unique(grid.latitude)

# Calculate day length for all days of year at each latitude
day = daylength(latitude_list, DOY)





# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot daylength for a couple of latitudes
heatmap(day,
  cbar=true,
  xlabel="Day of Year",
  ylabel="Latitude")

# Visualise the data for the frst and last latitude
plot(day[1, :])
plot!(day[end, :],
  xlabel="Day of Year",
  ylabel="Day Length (hours)")







# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save the output as a CSV
out_df = DataFrame(latitude=repeat(latitude_list, inner=length(DOY)),
  DOY=repeat(DOY, outer=length(latitude_list)),
  daylength=reshape(day', prod(size(day))))

CSV.write(joinpath(dataDir, outFile), out_df)