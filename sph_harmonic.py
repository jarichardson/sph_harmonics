#!/usr/bin/python
import sys, time
from numpy import *

#Evaluates fully normalized spherical harmonic coefficients from a table
#All variables needing to be changed are below (Lines 1-50), but use a configuration file.
#For Term Project, Potential Fields 2013, Jacob Richardson

#VARIABLE SET UP
global GM
global ref_r
global min_degree
global max_degree

#load configuration file
params = {}
try:
	execfile(sys.argv[1],params)
except: #if config file not correct, give syntax and exit
	sys.stdout.write("\nError: Cannot Load Parameters. Configuration file not given or does not exist.\nSyntax: ./sph_harmonic.py <configuration-file>\n\nExiting.\n")
	sys.exit()

model_type = params["MODEL_TYPE"]
model_weight = params["MODEL_WEIGHT"]

coeff_tbl_file = params["COEFFICIENT_TABLE_FILE"] #input file
header_lines = params["HEADER_LINES"]
delim= params["DELIMITER"]
col_d = params["DEGREE_COLUMN"]
col_o = params["ORDER_COLUMN"]
col_C = params["C_COEFFICIENT_COLUMN"]
col_S = params["S_COEFFICIENT_COLUMN"]

min_degree = params["MINIMUM_DEGREE"]
max_degree = params["MAXIMUM_DEGREE"]

west  = params["WEST_LONGITUDE"]
east  = params["EAST_LONGITUDE"]
south = params["SOUTH_LATITUDE"]
north = params["NORTH_LATITUDE"]
degstep = params["DEGREE_INTERVAL"]

ref_r =  params["REFERENCE_RADIUS"]
GM = params["GM"]

outfile = params["OUTFILE"]

#####################################
#FUNCTIONS
#####################################

def multipliers(code): #adds multipliers to a simple SH expansion
	in_coeff = ones(max_degree+1)
	out_coeff = 1
	
	if code=="dg":
		print "Gravity Disturbance (mGal)"
		for n in range(max_degree+1):
			in_coeff[n] = in_coeff[n]*n+1
		out_coeff = 1e5*GM/(ref_r**2) #1e5 converts to mGal (GM/R^2)
	elif code=="g": #g + dg
		print "Gravitational Acceleration (m s^-2)"
		for n in range(max_degree+1):
			in_coeff[n] = in_coeff[n]*n+1 
		out_coeff = GM/(ref_r**2) #(GM/R^2)
	elif code=="topo":
		print "Topography (meters)"
		out_coeff = 1000 #multiplies km model to meters
	elif code=="radialtopo":
		out_coeff = 1000
		print "Radial Topography (meters)"
	else:
		print "Simple Model (no multipliers/model type not identified)"
		
	return in_coeff,out_coeff
	
def model_formatting(v): #adds value to SH value
	if model_type=="radialtopo":
		v += ref_r #adds reference radius
	elif model_type=="g":
		v += (GM/ref_r**2) #adds average planetary gravity
	return v
	
def loadCoeffs(coeff_tbl_file,header_lines,col_d,col_o,col_C,col_S,delim):
	#Load data given user input columns at top.
	if delim==" ": #do not use delimeter parameter for white space
		degrees,orders,stokesC,stokesS = loadtxt(coeff_tbl_file,skiprows = header_lines,usecols=(col_d,col_o,col_C,col_S),unpack=True)
	else:
		degrees,orders,stokesC,stokesS = loadtxt(coeff_tbl_file,skiprows = header_lines,usecols=(col_d,col_o,col_C,col_S),delimiter=delim,unpack=True)

	#Set Stokes Coefficients with a loop given the raw data
	C = empty([max_degree+1,max_degree+1])
	S = empty([max_degree+1,max_degree+1])

	line = 0
	stopline = max_degree
	while degrees[line] <= stopline: #while the current degree is not bigger than max_degree
		d = degrees[line]
		o = orders[line]
		C[d][o] = stokesC[line] #set current line's C coefficient to C_n,m
		S[d][o] = stokesS[line]
		line += 1
		try:
			degrees[line]  #test for end of file
		except IndexError: #if end of file test if the last line was supposed to be the last line
			line -= 1
			stopline = 0
			if o<max_degree: #if order, m, is not max_degree, n cannot be max_degree either, so not enough coefficients.
				sys.stdout.write("\nCoefficient Table ends at n, m = %d, %d, not n, m = %d, %d!!" % (d,o,max_degree,max_degree))
				sys.stdout.write("\nError: Could not load desired coefficients (change max_degree).\n\n Exiting program :'(\n")
				sys.exit()
	return C,S
	
def nALF(mo,phi):
	#calculates the fully normalized associated legendre functions for a given phi
	#mo is max order, #phi is LATITUDE in radians (NOT Polar distance!)
	
	#            /1
	#           |
	#           | [NP(N,M;X)]^2 dX = 4 * pi    ,
	#           |
	#           /-1

	#   This program is based on a Fortran program by Robert L. Parker,
	#   Scripps Institution of Oceanography, Institute for Geophysics and 
	#   Planetary Physics, UCSD. February 1993.

	x = sin(phi)

	#error checking
	if mo < 0 or mo!=round(mo):
		print "n must be a positive integer."
		sys.exit(1)
	
	#some early constants
	tol = (finfo(x).tiny)**0.5 #calculate a sufficiently small number
	tstart = finfo(x).eps #what is the data type's resolution?
	
	normP = zeros([mo+1,mo+1])
	
	for n in range(mo+1):

		#compute list for normP(n,0:m;x)
		P = zeros(n+2) #0 to n including P[n] and an extra which will be removed at the end
	
		#The n=0 case, the only value returned is P(0,0)=1
		if n==0:
			P[0] = 1
		else:
			rootn = sqrt(range(0,(2*n+2)))
			s = sqrt(1-x**2) #equivalent of cos(phi)

			if x==-1: #calculate 2*cotangent(phi)
				twocot = inf
			elif x==1:
				twocot = -inf
			else:
				twocot = -2*x/s
	

	
			sn = (-s)**n
			if s>0 and abs(sn)<=tol:
				#if there will be underflow:
		
				v = 9.2-log(tol)/(n*s)
				w = 1/log(v)
				m1 = n*s*v*w*(1.0058+w*(3.819 - w*12.173))
				m1 = int(min(n, floor(m1)))
		
				P[m1] = sign(mod(m1+1,2)-0.5)*tstart
				if x<0:
					P[m1] = sign(mod(n+1,2)-0.5)*tstart
		
		
				mindex = range(m1)
				sumsq = tol
				for m in mindex[::-1]:
					P[m] = ((m+1)*twocot*P[m+1]-rootn[n+m+2]*rootn[n-m]*P[m+2])/ \
					(rootn[n+m+1]*rootn[n-m]); 
					sumsq = P[m]**2 + sumsq
				scale = 1/((2*sumsq - P[0]**2)**0.5)
				for m in range(m1+1):
					P[m] = scale*P[m]
	
	
			elif x!=1.0:
				#or if there is no underflow
		
				d = arange(2,2*n+1,2)
				c = prod(1.-1./d) #normalization constant
		
				P[n] = c**0.5*sn
				P[n-1] = P[n]*twocot*n/rootn[-1]
		
				mindex = range(n-1)
				for m in mindex[::-1]:
					P[m] =  (P[m+1]*twocot*(m+1) - P[m+2]*rootn[n+m+2]*rootn[n-m-1])/ \
					(rootn[n+m+1]*rootn[n-m]);
	
			elif abs(x)==1:
				#or if the phi is polar
		
				P[0] = 1
		
			#Now calculate fully normalized functions
			P[0] *= (2*n+1)**0.5
			P[1:-1] *= (2*(2*n+1))**0.5
	
			mindex = arange(1,n+1,2)
			for m in mindex[::1]:
				P[m] *= -1
	
		P = delete(P,n+1)	#delete the scrap
		normP[n][0:n+1] = P #append this n-row to the larger normalized P array

	return normP
	
def Potential(phi,lam,C,S,nP,in_coeff,out_coeff):
	#Evaluates the SH model at the given location
	#lat, long, C's, S's, norm. assoc. Legendre func.s, in_coefficients(with max_order entries), out_coefficient
	
	#convert lat,lon to radians here, never convert elsewhere
	phi = radians(phi)
	lam = radians(lam)
	
	#get cosine, sine multipliers
	Yc = [cos(m*lam) for m in range ((max_degree+1))]
	Ys = [sin(m*lam) for m in range ((max_degree+1))]

	pot = 0.0	#reset value
	for n in range(min_degree,(max_degree+1)): #for all n's
		mpot = 0.0 #reset order, m, linear combination sum
		for m in range((n+1)): #for all m's from 0 to n, including n
			sinprod = S[n][m]*Ys[m]             #S(n,m)*cos(m*longitude)
			cosprod = C[n][m]*Yc[m]             #C(n,m)*cos(m*longitude)
			mpot += (cosprod+sinprod)*nP[n][m]  #P(n,m)*(C*Y_c + S*Y_s)
		pot += in_coeff[n]*mpot                 #Inner multiplier (e.g. n+1 for dg)
		
	pot *= out_coeff                            #Outer multiplier (e.g. GM/R^2 for dg)
	return pot




##########################
#BEGIN SPHERICAL HARMONIC EVALUATION
##########################

startTime = time.time() #timer

#Identify model and add multipliers
sys.stdout.write("\nBeginning Spherical Harmonic Model: ")
(in_coeff,out_coeff) = multipliers(model_type)

#Load coefficient table
sys.stdout.write("\nLoading Coefficients from "+coeff_tbl_file+"...")
sys.stdout.flush()

(C,S) = loadCoeffs(coeff_tbl_file,header_lines,col_d,col_o,col_C,col_S,delim)

sys.stdout.write("Loaded!\n")
sys.stdout.flush()

#Check wavelength and degree step
nat_wavelength = 360./max_degree
nat_pixelwidth = nat_wavelength/2
if degstep==0:
	degstep=nat_pixelwidth
	sys.stdout.write("Degree step set to %.3f degrees (180/max_order)\n" % degstep)
elif nat_pixelwidth>degstep:
	sys.stdout.write("\nWarning: degree step (%.3f deg) < 0.5*natural wavelength (%.3f deg, 360/max_order) \n" % (degstep,nat_wavelength))
	sys.stdout.write("Warning: aliasing may occur.\n")
	sys.stdout.flush()
	
#Open output file	
oF = open(outfile,'w')

#Create and check range values
latrange = arange(south,(north+degstep),degstep) #create latitude values
lonrange = arange(west,(east+degstep),degstep)   #create longitude values

if latrange.size==0:
	sys.stdout.write("\nError: no latitude values created.\n")
	sys.stdout.write("check parameters south and north.\n\n Exiting program :'(\n")
	sys.exit()
if lonrange.size==0:
	sys.stdout.write("\nError: no longitude values created.\n")
	sys.stdout.write("check parameters east and west.\n\n Exiting program :'(\n")
	sys.exit()
	
sys.stdout.write("\nEvaluation Range: South "+str(south)+" -> North "+str(north)+"; West "+str(west)+" -> East "+str(east)+"\n")


#Run through latitude and longitude values to evaluate locations with SH model
for lat in latrange:
	#calculate fully normalized associated Legendre functions once for each latitude value
	normP = nALF(max_degree,radians(lat))
	
	for lon in lonrange:	
	
		sys.stdout.write("\rCurrently Evaluating: lat, %0.3f ; lon, %0.3f" % (lat,lon))
		sys.stdout.flush()
	
		Value = Potential(lat,lon,C,S,normP,in_coeff,out_coeff) #Evaluate using the SH expansion!
		Final = model_formatting(Value)                         #Add any constants (e.g. reference radius for radial topography)

		oF.write('%.6f\t%.6f\t%.5f\n' % (lon,lat,Final))        #Write value to tab-delimited file

oF.close() #close output file

sys.stdout.write("\rEvaluation Complete!                               \n\n")
sys.stdout.write("Solutions written to \""+outfile+"\"") #print out output file location
curTime = time.time()
elapsed = curTime-startTime
print("\n\nTime elapsed: %0.3f seconds" % elapsed) #print out elapsed time
