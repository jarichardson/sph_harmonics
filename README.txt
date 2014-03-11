Readme File for sph_harmonic.py
parameter configuration file

Right now, this is the parameter set up, but with a bit more explanation.
More might come later.

#Available Models:
#gravity disturbance, "dg"
#     mGal
#     Difference between real gravity and "normal" gravity on spheriod
#gravitational acceleration, "g"
#     m s^-2
#     Acceleration due to gravity toward center of spheroid
#topography, "topo"
#     m
#     Difference between real elevation and spheroid elevation
#     (for unit conversion, change outside multiplier (out_coeff)
#        "multipliers" function)
#radial topography, "radialtopo"
#     m
#     Distance from spheroid center and surface
#     (for unit conversion, change outside multiplier (out_coeff)
#        "multipliers" function)
#To add models, change the functions "multipliers" and "model_formatting"

model_type = "topo"

#COEFFICIENTS TABLE SETUP
coeff_tbl_file = "lro_ltm01_pa_1080_sha_max50.tab" #input file
header_lines = 0
delim=" " #for whitespace delimited use " "
col_d = 0 #degree column (first column is 0, not 1)
col_o = 1 #order column
col_C = 2 #Cos coefficient column
col_S = 3 #Sin coefficient column

#MODEL SETUP
#to remove dynamic gravity (J2) for free-air, min deg should be 2.
#   gravity, min degree should likely be 2
min_degree = 0
max_degree = 50

west  = -180.0 #in lat,lon coords
east  = 178.0 #in degrees
south = -89  #global = w,-180; e,(180-degstep); s,-90; n,90
north = 89
degstep = 2 #for degree step based on max_degree of model (natural wavelength), set to 0.

ref_r =  1.738000000000000e06 #reference radius (meters)

#Gravitational Constant * Planet Mass (SI units)
#Needed for gravity models
GM = 4.902799806931690e12


#OUTPUT FILE
outfile = "out.llz"
#output format: tab delimited lat, long, value
