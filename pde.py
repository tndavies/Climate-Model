import matplotlib.pyplot as plt
import numpy as np

# ============================================ #

def solve_pde(tprof, lats, lat_idx, t):
	heat_capacity = 30 	 	# C
	diffusivity = 10 		# D
	ir_flux = 230 			# I
	albedo = 0.3 			# A

	diurnal_flux = 0.0

	# approximate 1st order spatial derivative
	fo_spatial = 0.0
	if abs(lat) != 90.0:
		lat_diff = lats[lat_idx + 1] - lats[lat_idx - 1]
		fo_spatial = (tprof[lat_idx + 1] - tprof[lat_idx - 1]) / (2.0 * lat_diff)

	# approximate 2nd order spatial derivative
	so_spatial = 0.0

	return diurnal_flux*((1.0-albedo) / heat_capacity) - (diffusivity/heat_capacity)*np.tan(lat)*fo_spatial + (diffusivity/heat_capacity)*so_spatial - (ir_flux/heat_capacity)


# Initial temperature conditions for all latitude bands.
lats = np.arange(-90, 90, 10)
initial_conditions = [273.0 for k in lats]
data = [initial_conditions]

# timestep (1 day)
duration = 3e7 # evolve temps for a year.
dt = 86400.0
time = 0

while (time < duration):
	current_temp_profile = data[-1]
	data.append([])

	# evolve temperatures
	for idx, temp in enumerate(current_temp_profile):
		data[-1].append( temp + solve_pde(current_temp_profile, lats, idx, time) * dt )

	# advance time
	time += dt

# plot temperature profile of Earth across time.