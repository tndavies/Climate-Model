import numpy as np

# ============================================ #

def Calc_Declination(t_s):
	obliquity, omega = np.radians(23.45), 1.98226e-7 # Earth orbital constants.
	sine_decl = -np.sin(obliquity) * np.cos(omega*t_s + np.pi/2)

	return np.arcsin(sine_decl)

# ============================================ #

def Calc_DiurnalFlux(lat, t_d):
	decl = Calc_Declination(t_d * 86400)
	q0 = 1360.0

	flux = 0.0
	if(np.isclose(lat, np.pi/2)):
		flux = q0*np.sin(decl) if (decl > 0) else 0.0
	elif(np.isclose(lat, -np.pi/2)):
		flux = q0*np.sin(-decl) if (decl < 0) else 0.0
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
		flux = (q0/np.pi) * (term0 + term1)

	return flux

# ============================================ #

