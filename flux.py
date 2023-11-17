import numpy as np

# ============================================ #

def Calc_MeanAnomaly(t_s):
	orbital_period_s = 365 * 86400

	n = (2*np.pi) / orbital_period_s
	tau = 0.0 # time when Earth is at its orbital periapsis point. 
	M = n*(t_s - tau)

	return M

# ============================================ #

def Calc_TrueAnomaly_Approximate(t_s, term_count=4):
	e = 0.01671

	# The series approximation for the true anomaly only converges if
	# the eccentricty of our elliptical orbit is less than the 'Laplace Limit'.
	LaplaceLimit = 0.6627
	assert(e < LaplaceLimit)

	# Utilise series approximation for true anomaly (f) [first 3 terms in series].
	# [ref: 'Celestial Mechanics', Moulton]
	M = Calc_MeanAnomaly(t_s)
	terms = [M, (2*e*np.sin(M)), (1.25*np.power(e,2)*np.sin(2*M)), (np.power(e,3)/12)*(13*np.sin(3*M)-3*np.sin(M))]
	
	return np.sum(terms[:term_count])

# ============================================ #

def Calc_TrueAnomaly_ODE(dt):
	P_s = 365 * 86400 # Period of Earth's orbit.
	e = 0.01671 # Earth's orbital eccentricity

	times = [0]
	pos = [0] 

	while(times[-1] < P_s):
		time, theta = times[-1], pos[-1]
		time_ROC = ( 2*np.pi * (1 + e * np.cos(theta))**2 ) / (P_s * np.power((1-e**2),1.5))
		pos.append(theta + time_ROC * dt) 
		times.append(time + dt)

	return times, pos	

# ============================================ #

def Calc_Declination(t_s, use_elliptical_orbit=True):
	# Define our orbital parameters.
	orbital_period_s = 365 * 86400
	e = 0.01671

	if(use_elliptical_orbit):
		# The 'ecliptic longitude' is the angle swept out by the Earth-Sun line,
		# relative to the periapsis point; ie: the 'true anomaly'.
		#
		# For a circular orbit, the true anomaly is simply (2*PI / period)*time.
		#
		# For an ellitpical orbit, we must solve Keplar's equation, relating
		# the mean anomaly (M) to the eccentricty (e) & eccentric anomaly (E);
		# it can be shown that we can obtain the true anomaly (f) from this.
		#
		# The solution to the equation requires numerical methods, or can be
		# approximated using a series expansion; we use the latter approach
		# for speed reasons.
		f = Calc_TrueAnomaly_Approximate(t_s)
		
		
		# Compute Earth-Sun declination angle, given the ecliptic longitude
		# for our desired orbit.
		sine_declination = -np.sin(np.radians(23.44)) * np.sin(f)
		declination = np.arcsin(sine_declination)

		return declination
	else:
		# Compute the mean anomaly for our orbit, at the given time.
		n = (2*np.pi) / orbital_period_s
		tau = 0.0 # time when Earth is at its orbital periapsis point. 
		M = n*(t_s - tau)

		sine_declination = -np.sin(np.radians(23.44)) * np.sin(M)
		declination = np.arcsin(sine_declination)

		return declination

# ============================================ #

def Calc_SolarEnergyAt(r):
	Luminosity = 3.846e26 # Sun's luminosity
	return Luminosity / (4 * np.pi * r**2)

# ============================================ #

def Calc_SolarRadiation(t_s):
	a = 149.6e9 # Earth-Sun Semi-major axis
	e = 0.01671 # Earth's orbital eccentricity

	# On an elliptical orbit, the distance between the Earth
	# & sun can be calculated via elliptical geometry, given 
	# that we know the true anomaly (f).
	f = Calc_TrueAnomaly_Approximate(t_s)

	# Here, f=0, corresponds to the Earth at the periapsis point.
	r = a*(1-np.power(e,2)) / (1 + e*np.cos(f))

	return Calc_SolarEnergyAt(r)

# ============================================ #

def ValidateSolarRadCalculation():
	# Given Earth's orbital parameters (a, e)
	# we can compute the periapsis and apsis
	# distances (relative to the sun), and
	# thus what our max and min solar flux
	# should be throughout the orbit.
	a = 149.6e9
	e = 0.0167

	periapsis = a*(1-e) # closest
	apsis = a*(1+e) # furthest

	exact_max_flux = Calc_SolarEnergyAt(periapsis)
	exact_min_flux = Calc_SolarEnergyAt(apsis)

	print("Theoretical Solar Flux values: ")
	print("Max flux: " + str(np.round(exact_max_flux,2)))
	print("Min flux: " + str(np.round(exact_min_flux,2)))

	fluxes = [Calc_SolarRadiation(86400 * t) for t in np.linspace(0, 365, 365)]
	model_max_flux = max(fluxes)
	mode_min_flux = min(fluxes)

	print("\nModel's Flux values: ")
	print("Max flux: " + str(np.round(model_max_flux,2)))
	print("Min flux: " + str(np.round(mode_min_flux,2)))


# ============================================ #

def Calc_DiurnalFlux(lat, t_s):
	Earth_e = 0.01671

	decl = Calc_Declination(t_s)
	q = Calc_SolarRadiation(t_s)

	flux = 0.0
	if(np.isclose(lat, np.pi/2)):
		flux = q*np.sin(decl) if (decl > 0) else 0.0
	elif(np.isclose(lat, -np.pi/2)):
		flux = q*np.sin(-decl) if (decl < 0) else 0.0
	else:
		H = 0.0
		tan_prod = np.tan(lat) * np.tan(decl)
		
		if(tan_prod < -1.0):
			return 0.0
		elif(tan_prod > 1.0):
			H = np.pi
		else:
			H = np.arccos(-tan_prod)

		term0 = H * np.sin(lat) * np.sin(decl)
		term1 = np.cos(lat) * np.cos(decl) * np.sin(H)
		flux = (q/np.pi) * (term0 + term1)

	return flux

# ============================================ #