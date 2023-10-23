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

def Calculate_IRCooling(T):
	SIGMA = 5.670374419e-8
	TauIR = 0.79 * np.power((T / 273), 3)
	return (SIGMA * np.power(T, 4)) / (1 + 0.75 * TauIR)

# ============================================ #
def eval_pde(lats, temps, j, t):
	lat = lats[j]
	latitude_temp = temps[j]

	sd1 = grad(lats, temps, j)
	sd2 = laplace(lats, temps, j)

	Albedo = Calculate_Albedo(latitude_temp)
	IR_Cooling = Calculate_IRCooling(latitude_temp)
	Flux_In = flux.Calc_DiurnalFlux(lat, t)
	Diffusivity = 0.5394
	Heat_Capacity = 5.25e6

	term0 = Flux_In * (1 - Albedo)
	term1 = Diffusivity * (sd2 - np.tan(lat) * sd1)
	term2 = IR_Cooling

	return (term0 + term1 - term2) / Heat_Capacity

# ============================================ #

def EvolveGlobalTemperatures(lats, initial_temps, duration_d):
	DAY_SECS = 86400
	TIME_STEP = DAY_SECS 

	data = [(initial_temps, 0.0)]
	times = np.arange(0 + TIME_STEP, duration_d * DAY_SECS, TIME_STEP)

	print("Simulating Earth's Climate ..")
	with alive_bar(times.size) as bar:
		for t in times:
			temps = data[-1][0]
			tbuff = []

			# evolve temperatures
			for j, temp in enumerate(temps):
				dTdt = eval_pde(lats, temps, j, t)
				evolved_temp = temp + dTdt * TIME_STEP

				assert(evolved_temp >= 0.0)
				tbuff.append(evolved_temp)

			data.append((tbuff, t))
			bar()

	return data
# ============================================ #