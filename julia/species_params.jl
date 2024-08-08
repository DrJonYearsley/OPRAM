# species_params.jl
#      Set parameters for the development model of various species
#
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================



#= Species : Agrilus anxius (bronze birch borer)
=#
agrilus_anxius = (base_temperature = 1.7f0,            # Degrees C
          threshold = 1004.0f0,                        # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C


#= Species Pseudips mexicanus 
=#
pseudips_mexicanus = (base_temperature = 8.5f0,        # Degrees C
          threshold = 889.2f0,                         # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C


#=  Species: Halyomorpha halys (brown marmorated stink bug)
# Faculative diapause 
 =#
halyomorpha_halys = (base_temperature = 12.9,       # Degrees C
          threshold = 625.0,                        # Degrees C
          diapause_photoperiod = 15.0,              # Hours
          diapause_temperature = 20.0)              # Degrees C


#=  Species: Leptinotarsa decemlineata (Colorado potato beetle)
# Obligate diapause 
=#
leptinotarsa_decemlineata = (base_temperature = 10.0,       # Degrees C
          threshold = 300.0,                                # Degrees C
          diapause_photoperiod = 12.0,                      # Hours
          diapause_temperature = missing)                   # Degrees C

