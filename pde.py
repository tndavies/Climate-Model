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
	assert(np.isclose(ss_above, ss_below)) 
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

# 'duration' & 'ts are in units of Earth days.
def EvolveGlobalTemperatures(lats, initial_temps, duration, tsm=1):
	DAY_CONST = 86400 # seconds in a day

	time = 0
	data = [(initial_temps, time)]
	tstep = tsm * DAY_CONST

	while (time < duration*DAY_CONST):
		temps = data[-1][0]
		tbuff = []

		# evolve temperatures
		for j, T in enumerate(temps):
			dTdt = eval_pde(lats, temps, j, time)
			tbuff.append(T + dTdt * tstep)

		# record calculated temperature profile
		record = (tbuff, time)
		data.append(record)

		# advance time
		time += tstep

	return data
# ============================================ #