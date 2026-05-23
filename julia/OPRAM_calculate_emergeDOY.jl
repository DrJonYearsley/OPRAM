# Use OPRAM model results to calculate dat of emergence with start date of 1st Jan

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#    Copyright (C) 2026  Jon Yearsley  (Jon.Yearsley@ucd.ie)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Load the required packages
using CSV
using DataFrames
using Statistics
using Dates
using JLD2
using Distributed
using SharedArrays
using TOML



include("OPRAM_io_functions.jl");
include("OPRAM_processing_functions.jl");


# ============================================================================================
# =============== Import parameter values =======================================================




# Model parameters are stored in a TOML file https://toml.io/en/
if length(ARGS) == 1
    nNodes, run_params, species_setup, paths = import_parameters(ARGS[1], false)

elseif length(ARGS) == 0 & isfile("parameters_UK_test.toml")
    nNodes, run_params, species_setup, paths = import_parameters("parameters_UK_test.toml", false)

else
    @error "No parameter file given"
end


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set species parameters from the parameter file (can be more than one species)
if species_setup.speciesStr[1] == "all"
    @info "Importing all species from the species file" * string(species_setup.speciesStr)
    species_params = import_species(species_setup.speciesFile, species_setup.speciesStr[1])
else
    @info "Importing species from the species file: " * string(species_setup.speciesStr)
    species_params = [import_species(species_setup.speciesFile, species_setup.speciesStr[s]) for s in eachindex(species_setup.speciesStr)]
end


# =============== End of parameter setup =====================================================
# ============================================================================================



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import location data and thin it using thinFactor
grid = read_grid(run_params)




for s in eachindex(species_params)
    # Find directory matching the species name in run_params
    regex = Regex(replace(lowercase(species_params[s].species_name),
        r"\s" => "\\w"))  # Replace spaces with reg expression
    speciesName = filter(x -> occursin(regex, x), readdir(paths.resultsDir))
    if length(speciesName) > 1
        @error "More than one species name found"
    elseif length(speciesName) == 0
        @error "Species not found"
    end

    @info "Importing data for $(speciesName[1])"


    # Create vector of files to import (one file for each year)
    inFiles = Vector{String}(undef, length(run_params.years))
    for y in eachindex(run_params.years)
        inFile = filter(x -> occursin(r"^" * speciesName[1] * "_" * run_params.country * "_" * string(run_params.years[y]) * "_1km_" * run_params.method * ".jld2", x),
            readdir(joinpath(paths.resultsDir, speciesName[1])))

        if length(inFile) > 1
            @error "More than one input file found"
        else
            inFiles[y] = joinpath(paths.resultsDir, speciesName[1], inFile[1])
        end

        # Import the data from the jld2 files
        df_1km = read_OPRAM_JLD2(inFiles[y], run_params.years, grid)

        # Retain only start date of 1st Jan
        filter!(row -> row.startDOY == 1, df_1km)

        # =========================================================
        # Write result to CSV file
        @info "Writing 1km results to CSV files"
        # Add eastings and northings into the 1km data frames (copy the data frames so we can round before saving)
        out_1km = rightjoin(grid[:, [:ID, :east, :north]], df_1km, on=[:ID])


        # Wrie CSV file
        fileout2 = joinpath(paths.outDir, speciesName[1], speciesName[1] * "_" * "01_Jan" * "_" *
                                                          string(run_params.years[y]) *
                                                          "_" * string(run_params.thinFactor) * "km.csv")
        CSV.write(fileout2, out_1km, missingstring="NA")
    end
end

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++