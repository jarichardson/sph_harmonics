Readme File for sph_harmonic.py, a spherical harmonics code
in python

Created by Jacob Richardson, jarichardson@mail.usf.edu
http://jarichardson.myweb.usf.edu
University of South Florida, 2013-4.

------------------------------------------------------------
Contents of the Readme File:
1. Running the code
2. Setting Parameters
	2.1. Models
		2.1.2. Weighting
	2.2. Coefficients Table
	2.3. Model Setup
	2.4. Physical Parameters
	2.5. Output file


------------------------------------------------------------

1. Running the code

The syntax for the code is 
./sph_harmonic.py <configuration file>


------------------------------------------------------------

2. Setting Parameters

Parameters are set in the configuration file. A sample 
configuration file is listed in this repository as
parameters.conf.

2.1. Models
Four models are currently available: gravity disturbance,
gravitational acceleration, topography, and radial
topography.

Gravity disturbance
   code: "dg"
   units: mGal
   Description: Difference between real gravity and "normal"
      gravity on spheroid.

Gravitational acceleration
   code: "g"
   units: m s^-2
   Description: Acceleration due to gravity toward center of
      spheroid.

Topography
   code: "topo"
   units: m
   Description: Difference between real elevation and 
      spheroid elevation. Currently assumes coefficients are
      in km, "multipliers" function converts to meters.   
      For unit conversion, change outside multiplier 
      (out_coeff) "multipliers" function.
      
Radial topography
   code: "radialtopo"
   units: m
   Distance between spheroid center and surface. Currently 
      assumes coefficients are in km, "multipliers" function 
      converts to meters. For unit conversion, change 
      outside multiplier (out_coeff) "multipliers" function)
      
To add models, change the functions: 
   "multipliers" and/or "model_formatting"

2.1.2. Weighting

All models can be unweighted (regular), weighted shallow or
weighted deep. Weighted "deep" favors long wavelengths which
should be compensated, based on physical parameters (see 
2.4). Weighted "shallow" favors short wavelengths which
should be supported by an elastic crust.


2.2. Coefficients Table

A four column coefficients table is required. Give the file
name in COEFFICIENT_TABLE_FILE.

The delimiter should be set as one character (e.g. ","). If
delimiter is white space, set DELIMITER = " ".


2.3. Model Setup

The minimum and maximum degrees are model dependent. For
instance, gravity anomaly models should generally start at
1 or 2. Max degree is either the largest degree in the
coefficients table or the largest degree needed for the
desired model.

The range is set with the WESN long/lat parameters. For a
global model, set SOUTH_LATITUDE = -90; NORTH_LATITUDE = 90;
WEST_LONGITUDE = -180. Set EAST_LONGITUDE to 180 minus the
degree step set by DEGREE_INTERVAL.

If DEGREE_INTERVAL = 0, it will be automatically set to
0.5 * Natural Wavelength, where the Natural Wavelength is
360/MAXIMUM_DEGREE.


2.4. Physical Parameters

All physical parameters should be in M-K-S Units.

Reference Radius is required for all models, but topography.
GM is required for gravity models
All other parameters are required for the current weighted
models (YOUNGS_MODULUS, POISSONS_RATIO, 
MANTLE_CRUST_DENSITY_CONTRAST, CRUSTAL_THICKNESS).


2.5. Output file

The output file is a headerless tab-delimited ASCII file, 
with the format:
longitude   latitude   value


