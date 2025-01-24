# species_params.jl
#      Set parameters for the development model of various species
#
# Species in this file are:

# Quiescence
#  :agrilus_anxius
#  :ips_cembrae
#  :ips_duplicatus
#  :ips_sexdentatus
#  :pseudips_mexicanus
#  :spodoptera_frugiperda

# Obligate diapause
#  :leptinotarsa_decemlineata
#  :oulema_melanopus_model1
#  :oulema_melanopus_model2

# Facultative
#  :halyomorpha_halys
#  :ips_typographus

# 
# Jon Yearsley (jon.yearsley@ucd.ie)
# 7th Aug 2024
#
# ====================================================================



#= Species : Agrilus anxius (bronze birch borer)
Quiescence
=#
agrilus_anxius = (name = "Agrilus anxius",
          base_temperature = 1.7f0,            # Degrees C
          threshold = 1004.0f0,                        # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C


#= Species Pseudips mexicanus (Monterey Pine Engraver)
Quiescence
=#
pseudips_mexicanus = (name = "Pseudips mexicanus",
          base_temperature = 8.5f0,        # Degrees C
          threshold = 889.2f0,                         # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Spodoptera frugiperda (fall armyworm)
Quiescence
=#
spodoptera_frugiperda = (name = "Spodoptera frugiperda",
          base_temperature = 13.01f0,   # Degrees C
          threshold = 391.02f0,                        # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Ips sexdentatus (six-toothed bark beetle)
Quiescence
=#
ips_sexdentatus = (name = "Ips sexdentatus",
          base_temperature = 11.0f0,   # Degrees C
          threshold = 517.0f0,                        # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Ips duplicatus (double spined bark beetle)
Quiescence
=#
ips_duplicatus = (name = "Ips duplicatus",
          base_temperature = 11.1f0,           # Degrees C
          threshold = 217.3f0,                         # Degrees C
          diapause_photoperiod = missing,              # Hours
          diapause_temperature = missing);             # Degrees C

#= Species : Ips cembrae (large larch bark beetle)
Quiescence
=#
ips_cembrae = (name = "Ips cembrae",
          base_temperature = 11.3f0,               # Degrees C
          threshold = 243.3f0,                          # Degrees C
          diapause_photoperiod = missing,               # Hours
          diapause_temperature = missing);              # Degrees C


# ============================================================


#=  Species: Leptinotarsa decemlineata (Colorado potato beetle)
# Obligate diapause 
=#
leptinotarsa_decemlineata = (name = "Leptinotarsa decemlineata",
          base_temperature = 10.0,       # Degrees C
          threshold = 300.0,                                # Degrees C
          diapause_photoperiod = 12.0,                      # Hours
          diapause_temperature = missing)                   # Degrees C

#=   Oulema melanopus (cereal leaf beetle)
# Obligate diapause 
=#

# This Paul's Model 2
oulema_melanopus_model2 = (name = "Oulema melanopus, model 1",
            base_temperature = 10.0,                   # Degrees C
            threshold = 491.7,                                # Degrees C
            diapause_photoperiod = 12.0,                      # Hours
            diapause_temperature = missing)                   # Degrees C

            # This Paul's Model 1
oulema_melanopus_model1 = (name = "Oulema melanopus, model 2",
            base_temperature = 9.0,       # Degrees C
            threshold = 553.0,                                # Degrees C
            diapause_photoperiod = 12.0,                      # Hours
            diapause_temperature = missing)                   # Degrees C


#=  Species: Halyomorpha halys (brown marmorated stink bug)
# Faculative diapause 
 =#
 halyomorpha_halys = (name = "Halyomorpha halys",
        base_temperature = 12.9,       # Degrees C
        threshold = 625.0,                        # Degrees C
        diapause_photoperiod = 15.0,              # Hours
        diapause_temperature = 20.0)              # Degrees C



#=   Ips typographus (European spruce bark beetle)
# Facultative diapause 
=#
# These values are the same as the RITY2 model (Ogris et al, 2019)
# https://doi.org/10.1016/j.ecolmodel.2019.108775
# Warning: 
#   The RITY2 model uses bark temperature NOT air temperature
#   RITY2 does not have a diapause temperature... but see
#       https://doi.org/10.1111/j.1439-0418.2006.01123.x

ips_typographus = (name = "Ips typographus",
        base_temperature = 8.3,             # Degrees C
        threshold = 557.0,                             # Degrees C
        diapause_photoperiod = 14.5,                   # Hours
        diapause_temperature = 14.5)                   # Degrees C

