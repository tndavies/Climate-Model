# ---------------------------------------------------------------- #
# 							Imports 							
# ---------------------------------------------------------------- #
from data import Earth_Geography
from data import Climate_Temperatures
from data import Climate_CO2_Concentrations

from alive_progress import alive_bar
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------- #
# 						Utility Functions 							
# ---------------------------------------------------------------- #
def days_to_seconds(x):		return x * 86400
def days_to_years(x): 		return x / 365
def years_to_days(x): 		return x * 365
def seconds_to_years(x): 	return x / (86400 * 365)
def to_kelvin(x): 			return x + 273.15
def to_celsius(x): 			return x - 273.15
def years_to_seconds(x):	return x * 365 * 86400

def Is_AntarcticLatitudeBand(lat: float):
	return lat >= -np.radians(90) and lat <= -np.radians(70)

def Correct_SeaLevelTemperature(T: float, altitude_gain: float):
	# @tidy: do we need to convert to celsius first here?
	Lapse_Rate = 6.5 # degrees celsius lost, per 1km above sea-level.
	Corrected_Temperature = to_celsius(T) - Lapse_Rate * (altitude_gain / 1000)
	return to_kelvin(Corrected_Temperature)

# ---------------------------------------------------------------- #
# 						Simulation Parameters 							
# ---------------------------------------------------------------- #

Latitude_Step = 10 	# [degrees]
Time_Step = days_to_seconds(1)

Starting_Year = 1961
Antarctica_Altitude = 2500			# [m]
Orbital_Obliquity = -0.40911 		# [rads]
Orbital_Period = 3.1558e7 			# [s]
Orbital_SemiMajorAxis = 149.6e9 	# [m]
Orbital_Eccentricity = 0.01671
Sun_Luminosity = 3.846e26 			# [W]
C_Transitioning_Ice = 5.31e7 		# [J/K]
C_Land = 1.11e7 					# [J/K]
C_Ocean = 2.201e8 					# [J/K]
Diffusivity = 0.5394				# [???]

# ---------------------------------------------------------------- #
# 					Diurnal Solar Radiation (S) 							
# ---------------------------------------------------------------- #

def calc_Declination(t: float):
	# Using the declination equation for a generically shaped orbit,
	# and applying to an elliptical orbit, we calculate the declination
	# angle.

	EarthSun_Angle = calc_EarthSunAngle(t)
	Sine_Declination = np.sin(Orbital_Obliquity) * np.sin(EarthSun_Angle)

	return np.arcsin(Sine_Declination)

def calc_EarthSunAngle(t: float):
	e = Orbital_Eccentricity # for concise-ness.

	# We begin by computing the mean anomaly angle for the
	# Earth-Sun orbit, at the desired time.
	m = ((2 * np.pi) / Orbital_Period)*t

	# The series approximation we use to calculate the angle between the two bodies
	# is only valid for orbital eccentricies less than the Laplace limit (0.6627).
	assert(e < 0.6627)

	# Evaulate the first three terms in the series approximation for the true anomaly angle,
	# derived in 'Celestial Mechanics' by Moulton.
	Body_Angle = m + (2*e*np.sin(m)) + (1.25*np.power(e,2)*np.sin(2*m)) + (np.power(e,3)/12)*(13*np.sin(3*m)-3*np.sin(m))

	return Body_Angle

def calc_RadiationIntensity(t: float):
	# Given the angle between the Sun and Earth (true anomaly angle), we compute
	# the distence between the two bodies via the polar equation for an ellipse,
	# and then compute how much of the sun's radition reaches that distance.

	EarthSun_Angle = calc_EarthSunAngle(t)
	EarthSun_Dist = Orbital_SemiMajorAxis * (1.0 - Orbital_Eccentricity**2.0) / (1.0 + Orbital_Eccentricity * np.cos(EarthSun_Angle))
	Radiation_Intensity = Sun_Luminosity / (4 * np.pi * EarthSun_Dist**2.0)

	return Radiation_Intensity

def calc_AverageRadiation(lat: float, t: float):
	# Computes the average solar radiation that hits Earth at a
	# specified point in time, where the average is taken over
	# a day.

	q = calc_RadiationIntensity(t)
	decl = calc_Declination(t)

	avg_radiation = 0.0
	if(np.isclose(lat, np.pi/2)):
		avg_radiation = q*np.sin(decl) if (decl > 0) else 0.0
	elif(np.isclose(lat, -np.pi/2)):
		avg_radiation = q*np.sin(-decl) if (decl < 0) else 0.0
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
		avg_radiation = (q/np.pi) * (term0 + term1)

	return avg_radiation

# ---------------------------------------------------------------- #
# 						Albedo Function (A) 							
# ---------------------------------------------------------------- #
def calculate_Albedo(T: float):
	# Computes the amount of incoming solar radiation 
	# that is reflected back by the Earth's surface;
	# accounting for the presence of snow/ice at
	# colder temperatures. (Spiegel, model 2)
	return 0.525 - 0.245 * np.tanh(0.2*T-53.6)

# ---------------------------------------------------------------- #
# 						IR Cooling Function (I) 							
# ---------------------------------------------------------------- #
def calculate_IRCooling(T, p):
	SIGMA = 5.670374419e-8
	TauIR = 0.79 * np.power((T / 273), 3)
	return (SIGMA * np.power(T, 4)) / (1 + 0.75 * TauIR)

def calculate_IRCooling2(T: float, t: float):
	Band_Emission = 5.67e-8 * (T**4)

	# Default CO2 concentration to the last historical value.
	LastDataYear = list(Climate_CO2_Concentrations)[-1]
	CO2_Concentration = Climate_CO2_Concentrations[LastDataYear]

	# If a CO2 concentration exists in the historical record, for
	# the period that we're in, then use this instead.
	Current_Period = int(Starting_Year + seconds_to_years(t))
	if Current_Period in Climate_CO2_Concentrations:
		CO2_Concentration = Climate_CO2_Concentrations[Current_Period]

	Atmospheric_Retension = ((T / 273)**3) * (CO2_Concentration/Climate_CO2_Concentrations[LastDataYear])**1.5

	return Band_Emission / (1 + Atmospheric_Retension)

# ---------------------------------------------------------------- #
# 						Heat Capacity (C) 							
# ---------------------------------------------------------------- #
def calc_IceFraction(T: float):
	# Computes the fraction of the land/ocean mass,
	# within a latitude band, that is covered in snow/ice.
	# (Vladilo, A.6)
	param = (T - 273.0) / 10.0
	return max(0, 1 - np.exp(param))

def get_OceanFraction(lat: float):
	# Performs a lookup into Earth's geography table to 
	# find the fraction of the given latitude band that is 
	# covered in ocean water. (table 3, WK97)

	lat_degrees = np.degrees(lat)
	assert(lat_degrees >= -90 and lat_degrees <= 90) # ensure lat-coord is valid. 

	for Band_Range in Earth_Geography:
		if(lat_degrees >= Band_Range[1] and lat_degrees <= Band_Range[0]):
			return Earth_Geography[Band_Range]

def calc_AverageHeatCapacity(lat: float, T: float):
	# Computes the heat capacity of a latitude band,
	# by averaging over all contributions from land,
	# ice and sea. (Vladilo, A.1)

	C_Ice = C_Transitioning_Ice if(T >= 263 and T <= 273) else C_Land

	Ocean_Frac = get_OceanFraction(lat)
	Ice_Frac = calc_IceFraction(T)
	Land_Frac = 1 - Ocean_Frac

	Land_Term = Land_Frac*((1-Ice_Frac)*C_Land + Ice_Frac*C_Ice)
	Ocean_Term = Ocean_Frac*((1-Ice_Frac)*C_Ocean + Ice_Frac*C_Ice)
	Average_HeatCapacity = Land_Term + Ocean_Term

	return Average_HeatCapacity

# ---------------------------------------------------------------- #
# 							1D-EBM PDE 							
# ---------------------------------------------------------------- #
def calc_Temp_sROC1(lats: float, temps: float, k: int):
	lat = lats[k]  # latitude we want the derivative evaluated at.
	
	if(np.isclose(np.absolute(lat), np.pi/2)):
		return 0.0 # the derivative is zero at the poles.

	# for any other latitude band, use central difference approximation.
	lat_above, lat_below = lats[k + 1], lats[k - 1]
	t_above, t_below = temps[k + 1], temps[k - 1]
	
	return (t_above - t_below) / (2 * (lat_above - lat))

def calc_Temp_sROC2(lats: float, temps: float, k: int):
	lat = lats[k] # latitude we want the derivative evaluated at.

	# north pole
	if(np.isclose(lat, np.pi/2)): 
		step = lat-lats[k-1]
		return -(temps[k] - temps[k-1]) / step**2

	# south pole
	elif(np.isclose(lat, -np.pi/2)): 
		step = lats[k+1]-lat
		return (temps[k+1] - temps[k]) / step**2

	# other latitudes
	else:
		return (temps[k+1] - 2*temps[k] + temps[k-1]) / (lats[k+1]-lat)**2

def calc_Temp_tROC(lats: list, temps: list, band_index: int, t: float):
	lat = lats[band_index]
	Band_Temp = temps[band_index]

	Temp_sROC1 = calc_Temp_sROC1(lats, temps, band_index)
	Temp_sROC2 = calc_Temp_sROC2(lats, temps, band_index)

	Albedo = calculate_Albedo(Band_Temp)
	Radiation_In = calc_AverageRadiation(lat, t)
	Heat_Capacity = calc_AverageHeatCapacity(lat, Band_Temp)
	IR_Cooling = calculate_IRCooling(Band_Temp, t)

	Temp_tROC = (Radiation_In * (1 - Albedo) + Diffusivity * \
		(Temp_sROC2 - np.tan(lat) * Temp_sROC1) - IR_Cooling) / Heat_Capacity

	return Temp_tROC

def SimClimate(sim_duration: float):
	# Simulates the climate, starting at the specified year,
	# for the specified number of years given by 'sim_duration';
	# returns each time-point, relative to starting years, along
	# with the temperature profile for each of these times.
	# ---------------------------------------------------

	Lat_Grid = []
	for Lat_Coord in np.arange(-90, 90 + Latitude_Step, Latitude_Step):
		Lat_Grid.append(np.radians(Lat_Coord))

	# Retrieve the historical global average temperature for
	# the specified starting year of the simulation.
	assert(Starting_Year in Climate_Temperatures)
	Starting_Temperature = Climate_Temperatures[Starting_Year]
	Starting_Temperature = to_kelvin(Starting_Temperature)

	Temperature_Sets = [[Starting_Temperature] * len(Lat_Grid)]
	Time_Points = [0] # [Days]

	# Metrics for the progress bar.
	Total_Iterations = sim_duration / Time_Step
	Num_Iterations_Done = 0

	with alive_bar(title="Climate Simulation", manual=True) as bar:
		# ---------------------------------------------------- 
		while(Time_Points[-1] < sim_duration):
			Curr_TemperatureProfile = Temperature_Sets[-1]
			Time = Time_Points[-1]

			# ----- 
			New_TemperatureProfile = []
			for k in range(len(Lat_Grid)):
				Temp_tROC =	calc_Temp_tROC(Lat_Grid, Curr_TemperatureProfile, k, Time)

				Band_Temp = Curr_TemperatureProfile[k]
				New_BandTemperature = Band_Temp + Temp_tROC * Time_Step
				assert(New_BandTemperature >= 0.0)

				# Correct sea-level temperature for Antarctica, so we get the
				# correct temperature at the surface of the continent. 
				if Is_AntarcticLatitudeBand(Lat_Grid[k]):
					Band_Temp = Correct_SeaLevelTemperature(Band_Temp, Antarctica_Altitude)

				New_TemperatureProfile.append(New_BandTemperature)
			# ----- 

			Time += Time_Step
			Time_Points.append(Time)
			Temperature_Sets.append(New_TemperatureProfile)

			# Update the progress bar.
			Num_Iterations_Done += 1
			Sim_Status = Num_Iterations_Done / Total_Iterations
			bar(Sim_Status)
		# ---------------------------------------------------- 

		# Remap time points in simultation to real world time,
		# which is relative to the starting year.
		for k in range(len(Time_Points)):
			Time_Points[k] = Starting_Year + seconds_to_years(Time_Points[k])

	return Time_Points, Temperature_Sets

# ---------------------------------------------------------------- #
# 							Main Code Path 							
# ---------------------------------------------------------------- #
times, temps_sets = SimClimate(years_to_seconds(45))

# foos = []
# for k in temps_sets:
# 	foos.append(k[0])

# plt.plot(times, foos)
# plt.show()

idxs = np.arange(0, len(temps_sets)-1, 365)
foo, temps = [], []
for k in idxs:
	foo.append(times[k])
	temps.append(np.mean(temps_sets[k])) # this needs to be weighted based on band area!

plt.figure()
plt.xlabel("Years")
plt.ylabel("Global Average Temperature [K]")
plt.plot(foo, temps, ".-", color=(0,0,0), label="model")

time_data = []
temps_data = []
for k in Climate_Temperatures:
	time_data.append(k)
	temps_data.append(to_kelvin(Climate_Temperatures[k]))

plt.plot(time_data, temps_data, "x--", color=(0.157, 0.89, 0.533),label="historical")
plt.grid()
plt.legend()
plt.show()