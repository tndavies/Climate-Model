import numpy as np

# ============================================ #

def Calc_Declination(t_s):
	obliquity, omega = np.radians(23.45), 1.98226e-7 # Earth orbital constants.
	sine_decl = -np.sin(obliquity) * np.cos(omega*t_s + np.pi/2)

	return np.arcsin(sine_decl)

# ============================================ #

def Calc_DiurnalFlux(lat, t_d):
	# earth-sun declination angle at this time
	decl = Calc_Declination(t_d * 86400)

	flux = 0.0
	q0 = 1360.0
	if(np.isclose(np.absolute(lat), np.pi/2)):
		# special case for the poles
		flux = (q0*np.sin(decl)) if (decl > 0) else 0.0
	else:
		cosine_H = -np.tan(lat) * np.tan(decl)
		H = np.arccos(cosine_H)

		term0 = H * np.sin(lat) * np.sin(decl)
		term1 = np.cos(lat) * np.cos(decl) * np.sin(H)
		flux = (q0/np.pi) * (term0 + term1)

	return flux

# ============================================ #

