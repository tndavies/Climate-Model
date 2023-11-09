import numpy as np

# ============================================ #

def Calc_MeanAnomaly(t_s):
	orbital_period_s = 365 * 86400

	n = (2*np.pi) / orbital_period_s
	tau = 0.0 # time when Earth is at its orbital periapsis point. 
	M = n*(t_s - tau)

	return M

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
		

		# The series approximation for the true anomaly only converges if
		# the eccentricty of our elliptical orbit is less than the 'Laplace Limit'.
		LaplaceLimit = 0.6627
		assert(e < LaplaceLimit)

		# Utilise series approximation for true anomaly (f) [first 3 terms in series].
		M = Calc_MeanAnomaly(t_s)
		f = M + (2*e*np.sin(M)) + (1.25*np.power(e,2)*np.sin(2*M)) + \
			(np.power(e,3)/12)*(13*np.sin(3*M)-3*np.sin(M))
		
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

def Calc_SolarRadiation(t_s):
	Luminosity = 3.846e26 # Sun's luminosity
	a = 149.6e9 # Earth-Sun Semi-major axis
	e = 0.01671 # Earth's orbital eccentricity

	# On an elliptical orbit, the distance between the Earth
	# & sun can be calculated given that we know the true anomaly (f).
	#
	# We can use the same approximation as for the declination calculation
	# to obtain (f), and in-fact we can get the distance (r) straight away
	# with a related series approximation [first 3 terms].
	M = Calc_MeanAnomaly(t_s)
	r_interim = 1 - (e*np.cos(M)) - (np.power(e,2)/2)*(np.cos(2*M)-1) \
					-(np.power(e,3)/8)*(3*np.cos(3*M) - 3*np.cos(M))  
	r = r_interim * a

	return Luminosity / (4 * np.pi * r**2)

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

def EllipticOrbit(SimTime_s, P_d, e):
	dt = 5000 if(e<0.8) else 50 # trial and error.
	P_s = P_d * 86400

	times = [0]
	pos = [0] 

	while(times[-1] < SimTime_s):
		time, theta = times[-1], pos[-1]
		time_ROC = ( 2*np.pi * (1 + e * np.cos(theta))**2 ) / (P_s * np.power((1-e**2),1.5))
		pos.append(theta + time_ROC * dt) 
		times.append(time + dt)

	return times, pos
	
# ============================================ #