import numpy as np

# ============================================ #

def Calc_Declination(t_s):
	obliquity, omega = np.radians(23.45), 1.98226e-7 # Earth orbital constants.
	sine_decl = -np.sin(obliquity) * np.sin(omega*t_s)

	return np.arcsin(sine_decl)

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

def Calc_SolarRadiation(t_s, e):
	Luminosity = 3.846e26 # sun's luminosity
	a = 149.6e9 # Earth_Sun_SMA
	#e = 0.01671 # Earth's orbital eccentricity

	# To find the angular position in the orbit,
	# at this time (t_s), we must solve Keplar's ODE,
	# like we do here ...
	# angular_pos = EllipticOrbit(t_s, 365, e)[-1][-1]
	# 
	# However, for Earth this is not required, as we
	# can see that as the eccentricity is close to zero,
	# the angular position is just given by (approximately)
	# ...
	angular_pos = ((2*np.pi) / (365 * 86400)) * t_s

	r = (a*(1-e**2)) / (1 + e * np.cos(angular_pos))

	return Luminosity / (4 * np.pi * r**2)


# ============================================ #

def Calc_DiurnalFlux(lat, t_s):
	Earth_e = 0.01671

	decl = Calc_Declination(t_s)
	q = Calc_SolarRadiation(t_s, Earth_e)

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