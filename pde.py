from alive_progress import alive_bar
import numpy as np

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
	#assert(np.isclose(ss_above, ss_below)) 
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

def eval_pde(lats, temps, j, t):
	lat = lats[j]
	fos = grad(lats, temps, j)
	sos = laplace(lats, temps, j)

	return sos - np.tan(lat) * fos

# ============================================ #

def EvolveGlobalTemperatures(lats, initial_temps, duration_s):
	START_TIME, TIME_STEP = 0.0, 0.01
	
	data = [(initial_temps, START_TIME)]
	times = np.arange(START_TIME + TIME_STEP, duration_s, TIME_STEP)

	with alive_bar(times.size) as bar:
		for t in times:
			temps = data[-1][0]
			tbuff = []

			# evolve temperatures
			for j, temp in enumerate(temps):
				dTdt = eval_pde(lats, temps, j, t)
				evolved_temp = temp + dTdt * TIME_STEP
				tbuff.append(evolved_temp)

			data.append((tbuff, t))
			bar()


	return data
# ============================================ #