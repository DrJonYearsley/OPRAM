# species_params.jl
#      Set parameters for the development model of various species
#
# Species in this file are:
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
# 
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================



#= Species : Agrilus anxius (bronze birch borer)
Quiescence
=#
agrilus_anxius = (base_temperature = 1.7f0,            # Degrees C
          threshold = 1004.0f0,                        # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C


#= Species Pseudips mexicanus (Monterey Pine Engraver)
Quiescence
=#
pseudips_mexicanus = (base_temperature = 8.5f0,        # Degrees C
          threshold = 889.2f0,                         # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Spodoptera frugiperda (fall armyworm)
Quiescence
=#
spodoptera_frugiperda = (base_temperature = 13.01f0,   # Degrees C
          threshold = 391.02f0,                        # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Ips sexdentatus (six-toothed bark beetle)
Quiescence
=#
ips_sexdentatus = (base_temperature = 11.0f0,   # Degrees C
          threshold = 517.0f0,                        # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Ips duplicatus (double spined bark beetle)
Quiescence
=#
ips_duplicatus = (base_temperature = 11.1f0,           # Degrees C
          threshold = 217.3f0,                         # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Ips cembrae (large larch bark beetle)
Quiescence
=#
ips_cembrae = (base_temperature = 11.3f0,               # Degrees C
          threshold = 243.3f0,                          # Degrees C
          diapause_photoperiod = missing,               # Hours
          diapause_temperature = missing);              # Degrees C


# ============================================================


#=  Species: Leptinotarsa decemlineata (Colorado potato beetle)
# Obligate diapause 
=#
leptinotarsa_decemlineata = (base_temperature = 10.0,       # Degrees C
          threshold = 300.0,                                # Degrees C
          diapause_photoperiod = 12.0,                      # Hours
          diapause_temperature = missing)                   # Degrees C

#=   Oulema melanopus (cereal leaf beetle)
# Obligate diapause 
=#
oulema_melanopus = (base_temperature = 9.0,       # Degrees C
            threshold = 553.0,                                # Degrees C
            diapause_photoperiod = 12.0,                      # Hours
            diapause_temperature = missing)                   # Degrees C



#=  Species: Halyomorpha halys (brown marmorated stink bug)
# Faculative diapause 
 =#
 halyomorpha_halys = (base_temperature = 12.9,       # Degrees C
        threshold = 625.0,                        # Degrees C
        diapause_photoperiod = 15.0,              # Hours
        diapause_temperature = 20.0)              # Degrees C



#=   Ips typographus (European spruce bark beetle)
# Facultative diapause 
=#
ips_typographus = (base_temperature = 8.3,       # Degrees C
        threshold = 557.0,                             # Degrees C
        diapause_photoperiod = 14.5,                   # Hours
        diapause_temperature = 14.5)                   # Degrees C

