# ---------------------------------------------------------------- #
# 							Imports 							
# ---------------------------------------------------------------- #
from dataclasses import dataclass
from data import Earth_Geography
from data import Historic_Temperatures
from alive_progress import alive_bar
from scipy.optimize import minimize
from scipy.optimize import Bounds
import numpy as np

# ---------------------------------------------------------------- #
# 						Simulation Parameters 							
# ---------------------------------------------------------------- #

Latitude_Step = 10 	# [degrees]
Time_Step = 8600 	# [s]
Starting_Year = list(Historic_Temperatures)[0] # [wallclock]
HistoricCo2_Domain = (1958,2024)
Antarctica_Altitude = 2500			# [m]
Orbital_Obliquity = -0.40911 		# [rads]
Orbital_Period = 3.1558e7 			# [s]
Orbital_SemiMajorAxis = 149.6e9 	# [m]
Orbital_Eccentricity = 0.01671
Earth_Radius = 6371e3				# [m]
Sun_Luminosity = 3.846e26 			# [W]
C_Transitioning_Ice = 5.31e7 		# [J/K]
C_Land = 1.11e7 					# [J/K]
C_Ocean = 2.201e8 					# [J/K]
Diffusivity = 0.5394				# [???]
Prehistoric_Co2 = 290 				# [ppm]
Co2_Norm = 428						# [ppm]

# ---------------------------------------------------------------- #
# 						Misc' Functions 							
# ---------------------------------------------------------------- #
def days_to_seconds(x):		return x * 86400
def days_to_years(x): 		return x / 365
def years_to_days(x): 		return x * 365
def seconds_to_years(x): 	return x / (86400 * 365)
def to_kelvin(x): 			return x + 273.15
def to_celsius(x): 			return x - 273.15
def years_to_seconds(x):	return x * 365 * 86400
def to_wallclock(t):		return Starting_Year + seconds_to_years(t)
def from_wallclock(t):		return years_to_seconds(t - Starting_Year)
def GtC_to_ppm(x):			return x / 2.08

@dataclass
class Sim_Pack:
	lats: list
	times: list
	tdists: list

def Is_AntarcticLatitudeBand(lat: float):
	return lat >= -np.radians(90) and lat <= -np.radians(70)

def Correct_SeaLevelTemperature(T: float, altitude_gain: float):
	# @tidy: do we need to convert to celsius first here?
	Lapse_Rate = 6.5 # degrees celsius lost, per 1km above sea-level.
	Corrected_Temperature = to_celsius(T) - Lapse_Rate * (altitude_gain / 1000)
	return to_kelvin(Corrected_Temperature)

def calc_AntarcticAlbedo(sim: Sim_Pack):
	Avg_Albedos = []

	for X in sim.tdists:
		Antarctic_Temps = []
		for idx,l in enumerate(sim.lats):
			if(Is_AntarcticLatitudeBand(l)):
				Antarctic_Temps.append(X[idx])

		albedos = [calculate_Albedo(T) for T in Antarctic_Temps]
		Avg_Albedo = np.mean(albedos)

		Avg_Albedos.append(Avg_Albedo)

	return Avg_Albedos

def calc_AntarcticAverageTemperature(sim: Sim_Pack):
	Temperature_Averages = []

	for X in sim.tdists:
		Antarctic_Temps = []
		
		for idx,l in enumerate(sim.lats):
			if(Is_AntarcticLatitudeBand(l)):
				Antarctic_Temps.append(X[idx])

		Temperature_Averages.append(np.mean(Antarctic_Temps))

	return Temperature_Averages

# @todo: double check the logic of the bands works here!
def calc_SphericalAverage(sim: Sim_Pack, domain: tuple, func: callable):
	C = 2*np.pi*(Earth_Radius**2.0)
	idx0 = sim.lats.index(np.radians(domain[0]))
	idx1 = sim.lats.index(np.radians(domain[1]))
	Band_Areas = [C*np.absolute(np.cos(sim.lats[k])-np.cos(sim.lats[k+1])) \
		for k in range(idx0, idx1)]

	Variable_Averages = []
	for X in sim.tdists:
		Variable_List = []
		for k in range(idx0, idx1):
			Variable_List.append(func(X[k]))

		a = np.average(Variable_List, weights=Band_Areas)
		Variable_Averages.append(a)

	return Variable_Averages

# @todo: needs to be rewritten to work with new calc_SphericalAverage function. 
def calc_ComparisonMetric(lat_grid: list, times: list, temp_frames: list):
		# Computes the sum of the residuls between the 
		# historical record and the model.

		Sampled_Times = times[0::365]
		Sampled_Temps = temp_frames[0::365]

		Similarity_Metric = 0.0
		for k in range(len(Sampled_Times)):
			Historical_Temp = to_kelvin(Historic_Temperatures[Sampled_Times[k]])
			Global_Average_Temp = calc_GlobalAverageTemperature(lat_grid, Sampled_Temps[k])
			Similarity_Metric += (Global_Average_Temp - Historical_Temp)**2

		return Similarity_Metric

def find_OptimalCoolingFunction():
	def Objective_Func(model_params: list):
		(exp, coeff) = model_params
		Lat_Grid, Times, Temp_Sets = SimClimate(calc_HistoricPeriod(), exp, coeff)
		X = calc_ComparisonMetric(Lat_Grid, Times, Temp_Sets)
		print("exp={exp}, coeff={coeff}, X={x}".format(exp=exp, coeff=coeff, x=X))
		return X 

	# Co2 Exponent is unbounded (-inf, inf), but the retention factor
	# must be in the range (0, inf).
	result = minimize(Objective_Func, x0=[0,1], bounds=Bounds([-np.inf, 0.0], [np.inf, np.inf]))
	print(result)

def interp_HistoricCo2(t):
	Timestamp = to_wallclock(t)
	assert(Timestamp >= 1958 and Timestamp <= 2024) # valid range of interpolation to dataset
	return 1.31647498e-02*np.power(Timestamp,2) - 5.07940891e+01*Timestamp + 4.92987770e+04

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
	assert (e < 0.6627), "Orbital eccentricity is too large"

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
def calculate_IRCooling(T: float, t: float, emission_pathway: callable):
	Alpha, Beta = 0.6786585695589182, 0.05114966426233368 # picked for model to match climate record
	Present_Co2 = None

	Timestamp = to_wallclock(t)
	if(Timestamp < HistoricCo2_Domain[0]):
		Present_Co2 = Prehistoric_Co2
	elif(Timestamp >= HistoricCo2_Domain[0] and Timestamp <= HistoricCo2_Domain[1]):
		Present_Co2 = interp_HistoricCo2(t)
	elif(Timestamp >= HistoricCo2_Domain[1]):
		Present_Co2 = emission_pathway(t)

	Atmospheric_Absorption = Alpha * np.power((T/273),3) * np.power(Present_Co2 / Co2_Norm, Beta)
	Dissipated_Radiation = (5.67e-8 * np.power(T,4)) / (1 + Atmospheric_Absorption)

	return Dissipated_Radiation

# ---------------------------------------------------------------- #
# 						IPCC Emission Pathways 							
# ---------------------------------------------------------------- #
def A1FI_Pathway(t):
	Timestamp = to_wallclock(t)
	Co2_Emissions = -3.39e-03*np.power(Timestamp,2) + 14.2*Timestamp - 14790 
	return GtC_to_ppm(Co2_Emissions)

def B1_Pathway(t):
	Timestamp = to_wallclock(t)
	Co2_Emissions = -1.54*np.power(Timestamp,2) + 6.24*Timestamp - 6.32e+03
	return GtC_to_ppm(Co2_Emissions)

# ---------------------------------------------------------------- #
# 						Heat Capacity (C) 							
# ---------------------------------------------------------------- #
def calc_IceFraction(T: float):
	# Computes the fraction of the land/ocean mass,
	# within a latitude band, that is covered in snow/ice.
	# (Vladilo, A.6)
	param = (T - 273.0) / 10.0
	return np.maximum(0, 1 - np.exp(param))

def get_OceanFraction(lat: float):
	# Performs a lookup into Earth's geography table to 
	# find the fraction of the given latitude band that is 
	# covered in ocean water. (table 3, WK97)

	lat_degrees = np.degrees(lat)

	Latitude_IsValid = (lat_degrees >= -90 and lat_degrees <= 90)
	assert Latitude_IsValid, "Invalid Latitude Coordinate" 

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

def calc_Temp_tROC(lats: list, temps: list, band_index: int, t: float, emission_pathway: callable, enable_AltCorr: bool):
	lat = lats[band_index]
	Band_Temp = temps[band_index]

	if(enable_AltCorr and Is_AntarcticLatitudeBand(lat)):
		Band_Temp = Correct_SeaLevelTemperature(Band_Temp, 2500)

	Temp_sROC1 = calc_Temp_sROC1(lats, temps, band_index)
	Temp_sROC2 = calc_Temp_sROC2(lats, temps, band_index)

	IR_Cooling = calculate_IRCooling(Band_Temp, t, emission_pathway)
	Heat_Capacity = calc_AverageHeatCapacity(lat, Band_Temp)
	Radiation_In = calc_AverageRadiation(lat, t)
	Albedo = calculate_Albedo(Band_Temp)

	Temp_tROC = (Radiation_In * (1 - Albedo) + Diffusivity * \
		(Temp_sROC2 - np.tan(lat) * Temp_sROC1) - IR_Cooling) / Heat_Capacity

	return Temp_tROC

def SimClimate(emission_pathway: callable, End_Year: float, enable_AltCorr: bool = True):
	Simulation_Length = years_to_seconds(End_Year - Starting_Year)
	Lat_Grid = [np.radians(l) for l in np.arange(-90, 90 + Latitude_Step, Latitude_Step)]

	# Retrieve the historical global average temperature for
	# the specified starting year of the simulation.
	assert (Starting_Year in Historic_Temperatures), \
		"Starting year of simulation not present in historical record."
	Starting_Temperature = to_kelvin(Historic_Temperatures[Starting_Year])

	Temperature_Sets = [[Starting_Temperature] * len(Lat_Grid)]
	Time_Points = [0] # [Days]

	# Metrics for the progress bar.
	Total_Iterations = Simulation_Length / Time_Step
	Num_Iterations_Done = 0

	with alive_bar(title="Climate Simulation", manual=True) as bar:
		# ---------------------------------------------------- 
		while(Time_Points[-1] < Simulation_Length):
			Curr_TemperatureProfile = Temperature_Sets[-1]
			Time = Time_Points[-1]

			# ----- 
			New_TemperatureProfile = []
			for k in range(len(Lat_Grid)):
				Temp_tROC =	calc_Temp_tROC(Lat_Grid, Curr_TemperatureProfile, k, Time, emission_pathway, enable_AltCorr)
				New_BandTemperature = Curr_TemperatureProfile[k] + Temp_tROC * Time_Step
				New_TemperatureProfile.append(New_BandTemperature)
				assert (New_BandTemperature >= 0.0), "Unphysical band temperature."
			# ----- 

			Time += Time_Step
			Time_Points.append(Time)
			Temperature_Sets.append(New_TemperatureProfile)

			# Update the progress bar.
			Num_Iterations_Done += 1
			Sim_Status = Num_Iterations_Done / Total_Iterations
			bar(Sim_Status)
		# ---------------------------------------------------- 

	Times = [to_wallclock(t) for t in Time_Points]
	return Sim_Pack(Lat_Grid, Times, Temperature_Sets)
