from alive_progress import alive_bar
import numpy as np
import flux

# ============================================ #

def grad(lats, temps, idx):
	 # latitude we want the derivative evaluated at.
	lat = lats[idx]

	# the derivative is zero at the poles.
	if(np.isclose(np.absolute(lat), np.pi/2)):
		return 0.0

	# for any other latitude band, use central difference approximation.
	lat_above, lat_below = lats[idx + 1], lats[idx - 1]
	t_above, t_below = temps[idx + 1], temps[idx - 1]
	
	# we assume a uniform spacing in latitude coords!
	ss_above = lat_above - lat
	ss_below = lat - lat_below
	step = ss_above

	# central difference approximation
	return (t_above - t_below) / (2 * step)

# ============================================ #

def laplace(lats, temps, i):
	# latitude we want the derivative evaluated at.
	lat = lats[i]

	# north pole
	if(np.isclose(lat, np.pi/2)): 
		step = lat-lats[i-1]
		return -(temps[i] - temps[i-1]) / step**2

	# south pole
	elif(np.isclose(lat, -np.pi/2)): 
		step = lats[i+1]-lat
		return (temps[i+1] - temps[i]) / step**2

	# other latitudes
	else:
		return (temps[i+1] - 2*temps[i] + temps[i-1]) / (lats[i+1]-lat)**2

# ============================================ #

def Calculate_Albedo(T):
	return 0.525 - 0.245 * np.tanh(0.2*T-53.6)

# ============================================ #

def Calculate_IRCooling(T):
	SIGMA = 5.670374419e-8
	TauIR = 0.79 * np.power((T / 273), 3)
	return (SIGMA * np.power(T, 4)) / (1 + 0.75 * TauIR)

# ============================================ #

def Get_OceanFraction(lat):
	def IsIn(x, a, b): # returns true if 'x in [a, b)'.
		return (x <= a and x > b)

	lat = np.degrees(lat)

	# Table III of Earth's ocean fraction from WK97.
	if IsIn(lat, 90, 80):
		frac = 0.934
	elif IsIn(lat, 80, 70):
		frac = 0.713
	elif IsIn(lat, 70, 60):
		frac = 0.294
	elif IsIn(lat, 60, 50):
		frac = 0.428
	elif IsIn(lat, 50, 40):
		frac = 0.475
	elif IsIn(lat, 40, 30):
		frac = 0.572
	elif IsIn(lat, 30, 20):
		frac = 0.624
	elif IsIn(lat, 20, 10):
		frac = 0.736
	elif IsIn(lat, 10, 0):
		frac = 0.772
	elif IsIn(lat, 0, -10):
		frac = 0.764
	elif IsIn(lat, -10, -20):
		frac = 0.78
	elif IsIn(lat, -20, -30):
		frac = 0.769
	elif IsIn(lat, -30, -40):
		frac = 0.888
	elif IsIn(lat, -40, -50):
		frac = 0.97
	elif IsIn(lat, -50, -60):
		frac = 0.992
	elif IsIn(lat, -60, -70):
		frac = 0.896
	elif IsIn(lat, -70, -80):
		frac = 0.246
	else: # lat < -80
		frac = 0.0

	return frac

# ============================================ #
C_Land = 1e6 + 10.1e6
C_ml50 = 210e6
C_Ocean = C_ml50 + 10.1e6

def Calculate_HeatCapacity(lat, T):
	# heat capacity prescription from Vladilo1 et al (2013).
	C_Ice = (C_Land + 0.2*C_ml50) if(T>=263 and T<=273) else C_Land

	# fraction of the land/ocean that is ice covered.
	fIce = max(0.0, 1.0 - np.exp((T-273)/10))
	fOcean = Get_OceanFraction(lat)
	fLand = (1- fOcean)

	land_term = fLand*((1-fIce)*C_Land + fIce*C_Ice)
	ocean_term = fOcean*((1-fIce)*C_Ocean + fIce*C_Ice)
	return (land_term + ocean_term)

# ============================================ #

def Evaluate_DiffusionPDE(lats, temps, j, t):
	lat = lats[j]
	lat_temp = temps[j]

	sd1 = grad(lats, temps, j)
	sd2 = laplace(lats, temps, j)

	Albedo = 0.0
	Heat_Capacity = 0.0

	# Antartica special case
	if(lat >= -np.radians(90) and lat <= -np.radians(70)):
		Heat_Capacity = C_Land
		Albedo = 0.7
	else: # Not Antartica
		Albedo = Calculate_Albedo(lat_temp)
		Heat_Capacity = Calculate_HeatCapacity(lat, lat_temp)

	IR_Cooling = Calculate_IRCooling(lat_temp)
	Flux_In = flux.Calc_DiurnalFlux(lat, t)
	Diffusivity = 0.5394

	term0 = Flux_In * (1 - Albedo)
	term1 = Diffusivity * (sd2 - np.tan(lat) * sd1)
	term2 = IR_Cooling

	return (term0 + term1 - term2) / Heat_Capacity

# ============================================ #

def SimulateClimate(SimTime_yrs, iv=400, lat_step=6):
	lats = [np.radians(k) for k in np.arange(-90, 90+lat_step, lat_step)]
	TempFrames = [[iv for k in lats]]

	TIME_STEP = 86400 
	times_s = np.arange(TIME_STEP, SimTime_yrs * 365.25 * 86400, TIME_STEP)

	print("Simulating Earth's Climate ..")
	with alive_bar(times_s.size) as bar:
		for t in times_s:
			temps, tbuff = TempFrames[-1], []

			# evolve temperatures
			for j, temp in enumerate(temps):
				dTdt = Evaluate_DiffusionPDE(lats, temps, j, t)
				evolved_temp = temp + dTdt * TIME_STEP
				assert(evolved_temp >= 0.0)
				tbuff.append(evolved_temp)

			TempFrames.append(tbuff)
			bar()

	times_s = np.append(0, times_s)
	times_d = list(np.divide(times_s, 86400))

	return list(zip(times_d, TempFrames)), lats, (TIME_STEP / 86400)

# ============================================ #

def DumpSim(sim):
	for ds in sim:
		time = ds[0]
		temps = ds[1]

		print("At time=" + str(time) + " days")
		print(temps)
		print("")

# ============================================ #
