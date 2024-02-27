# ---------------------------------------------------------------- #
# 							Imports 							
# ---------------------------------------------------------------- #
from alive_progress import config_handler
from numpy.polynomial import Polynomial
from alive_progress import alive_bar
from scipy.optimize import minimize
from scipy.optimize import Bounds
from dataclasses import dataclass
import numpy as np

from data import Historic_Temperatures
from data import Equilibrium_Config
from data import Earth_Geography
from data import Historic_Co2
from data import RCP85_Data
from data import RCP6_Data
from data import RCP45_Data
from data import RCP26_Data
from data import Serialise

config_handler.set_global(bar="classic2", spinner="classic") # Sets the style of the progress bar.

# ---------------------------------------------------------------- #
# 							Definitions 							
# ---------------------------------------------------------------- #
Antarctica_Altitude = 		2500		# [m]
Orbital_Obliquity = 		-0.40911 	# [rads]
Orbital_Period = 			3.1558e7 	# [s]
Orbital_SemiMajorAxis = 	149.6e9 	# [m]
Orbital_Eccentricity = 		0.01671
Earth_Radius = 				6371e3		# [m]
Sun_Luminosity = 			3.846e26 	# [W]
C_Transitioning_Ice = 		5.31e7 		# [J/K]
C_Land = 					1.11e7 		# [J/K]
C_Ocean = 					2.201e8 	# [J/K]
Diffusivity =				0.5394		# [???]
Co2_Norm = 					428			# [ppm]
Prerecord_Co2Level =		300			# [ppm]
Antarctic_Bounds = 			(-90,-70)	# [degrees]
Best_Alpha = 				0.6955093969474789
Best_Beta = 				0.05221459427741913
Stability_Duration = 		20 			# [years]

@dataclass
class Sim_Specification:
    Duration: float	# [years]
    # Optional (below)
    Initial_Year: float = list(Historic_Temperatures)[0]
    Altitude_Correction: bool = True
    Alpha: float = Best_Alpha
    Time_Step: float = 86400
    Beta: float = Best_Beta
    RCP: callable = None
    Lat_Step: int = 10

@dataclass
class Sim_Result:
	spec: Sim_Specification		# Copy of the Sim Spec used.
	times: list[float] 			# Array of timestamps.
	lats: list[float]			# Simulation latitude grid.
	tps: list[float]			# Temperature profile for each timestamp.

# ---------------------------------------------------------------- #
# 						Misc. Functions	 							
# ---------------------------------------------------------------- #
def Average(sim: Sim_Result, temp_map: callable, domain: tuple):
	C = 2*np.pi*(Earth_Radius**2.0)
	idx0 = sim.lats.index(np.radians(domain[0]))
	idx1 = sim.lats.index(np.radians(domain[1]))
	Band_Areas = [C*np.absolute(np.cos(sim.lats[k])-np.cos(sim.lats[k+1])) \
		for k in range(idx0, idx1)]

	Variable_Averages = []
	for X in sim.tps:
		Variable_List = []
		for k in range(idx0, idx1):
			Band_Temperature = (X[k] + X[k+1]) * 0.5
			Band_Variable = temp_map(Band_Temperature)
			Variable_List.append(Band_Variable)

		Avg = np.average(Variable_List, weights=Band_Areas)
		Variable_Averages.append(Avg)

	return Variable_Averages

def to_timestamp(initial_year: int, t: float):
    return initial_year + (t / 31536000)

# ---------------------------------------------------------------- #
# 						 Co2 Pathways 							
# ---------------------------------------------------------------- #
def Interpolate_Co2(co2_data: dict, poly_degree: int):
	Times, Concentrations = Serialise(co2_data)
	Interpolation = Polynomial.fit(Times, Concentrations, poly_degree)
	return Interpolation

RCP85 = Interpolate_Co2(RCP85_Data, 2)
RCP6 = Interpolate_Co2(RCP6_Data, 2)
RCP45 = Interpolate_Co2(RCP45_Data, 2)
RCP26 = Interpolate_Co2(RCP26_Data, 4)
RCP0 = Interpolate_Co2(Historic_Co2, 10)

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
	# We begin by computing the mean anomaly angle for the
	# Earth-Sun orbit, at the desired time.
	m = ((2 * np.pi) / Orbital_Period)*t

	# The series approximation we use to calculate the angle between the two bodies
	# is only valid for orbital eccentricies less than the Laplace limit (0.6627).
	e = Orbital_Eccentricity
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
def calc_Albedo(T: float):
	return 0.525 - 0.245 * np.tanh(0.2*T-53.6)

# ---------------------------------------------------------------- #
# 						IR Cooling Function (I) 							
# ---------------------------------------------------------------- #
def calculate_IRCooling(initial_year: int, t: float, rcp: callable, alpha: float, beta: float, temp: float):
	ts = to_timestamp(initial_year, t)
	
	Co2 = None
	if ts < list(Historic_Temperatures)[0]:
		Co2 = Prerecord_Co2Level
	elif ts >= list(Historic_Temperatures)[0] and ts <= list(Historic_Temperatures)[-1]:
		Co2 = RCP0(np.array(ts))
	elif ts > list(Historic_Temperatures)[-1]:
		Co2 = rcp(np.array(ts))

	Atmospheric_Absorption = alpha * np.power((temp/273),3) * np.power(Co2 / Co2_Norm, beta)
	Dissipated_Radiation = (5.67e-8 * np.power(temp,4)) / (1 + Atmospheric_Absorption)

	return Dissipated_Radiation

# ---------------------------------------------------------------- #
# 						Heat Capacity (C) 							
# ---------------------------------------------------------------- #
def Calc_IceFraction(T: float):
	param = (T - 273.0) / 10.0
	return np.maximum(0, 1 - np.exp(param))

def Get_OceanFraction(lat: float):
	lat_deg = np.degrees(lat)
	assert (lat_deg >= -90 and lat_deg <= 90), "Invalid Latitude Coordinate" 

	for Band_Range in Earth_Geography:
		if(lat_deg >= Band_Range[1] and lat_deg <= Band_Range[0]):
			return Earth_Geography[Band_Range]

def calc_AverageHeatCapacity(lat: float, T: float):
	# Computes the heat capacity of a latitude band,
	# by averaging over all contributions from land,
	# ice and sea. (Vladilo, A.1)
	C_Ice = C_Transitioning_Ice if(T >= 263 and T <= 273) else C_Land

	Ocean_Frac = Get_OceanFraction(lat)
	Ice_Frac = Calc_IceFraction(T)
	Land_Frac = 1 - Ocean_Frac

	Land_Term = Land_Frac*((1-Ice_Frac)*C_Land + Ice_Frac*C_Ice)
	Ocean_Term = Ocean_Frac*((1-Ice_Frac)*C_Ocean + Ice_Frac*C_Ice)
	Average_HeatCapacity = Land_Term + Ocean_Term

	return Average_HeatCapacity

# ---------------------------------------------------------------- #
# 							1D-EBM PDE 							
# ---------------------------------------------------------------- #
def Within_Antarctic(lat: float):
    return lat >= np.radians(Antarctic_Bounds[0]) and lat <= np.radians(Antarctic_Bounds[1])

def ApplyLapseRate(T: float, altitude_gain: float):
	# @tidy: do we need to convert to celsius first here?
	Lapse_Rate = 6.5 # degrees celsius lost, per 1km above sea-level.
	Corrected_Temperature = (T-273.15) - Lapse_Rate * (altitude_gain / 1000)
	return (Corrected_Temperature + 273.15)

def calc_Temp_sROC1(lats: list, temps: list, k: int):
	lat = lats[k]  # latitude we want the derivative evaluated at.
	
	if(np.isclose(np.absolute(lat), np.pi/2)):
		return 0.0 # the derivative is zero at the poles.

	# for any other latitude band, use central difference approximation.
	lat_above, lat_below = lats[k + 1], lats[k - 1]
	t_above, t_below = temps[k + 1], temps[k - 1]
	
	return (t_above - t_below) / (2 * (lat_above - lat))

def calc_Temp_sROC2(lats: list, temps: list, k: int):
	lat = lats[k] # latitude we want the derivative evaluated at.

	if(np.isclose(lat, np.pi/2)): # north pole
		step = lat-lats[k-1]
		return -(temps[k] - temps[k-1]) / step**2
	elif(np.isclose(lat, -np.pi/2)): # south pole
		step = lats[k+1]-lat
		return (temps[k+1] - temps[k]) / step**2
	else:
		return (temps[k+1] - 2*temps[k] + temps[k-1]) / (lats[k+1]-lat)**2

def calc_Temp_tROC(spec: Sim_Specification, lat_grid: list, tp: list, lat_idx: int, t: float):
	lat, temp = lat_grid[lat_idx], tp[lat_idx]

	if spec.Altitude_Correction and Within_Antarctic(lat):
		temp = ApplyLapseRate(temp, Antarctica_Altitude)

	sROC1 = calc_Temp_sROC1(lat_grid, tp, lat_idx)
	sROC2 = calc_Temp_sROC2(lat_grid, tp, lat_idx)

	IR_Cooling = calculate_IRCooling(spec.Initial_Year, t, spec.RCP, spec.Alpha, spec.Beta, temp)
	Heat_Capacity = calc_AverageHeatCapacity(lat, temp)
	Radiation_In = calc_AverageRadiation(lat, t)
	Albedo = calc_Albedo(temp)

	pde_lhs = Radiation_In*(1-Albedo) + Diffusivity*(sROC2-np.tan(lat)*sROC1)-IR_Cooling

	return pde_lhs / Heat_Capacity
	
def Simulate_Climate(spec: Sim_Specification) -> Sim_Result:
	assert spec.Duration >= 0, "Invalid Duration"
	Duration = spec.Duration * 31536000 # convert from years -> seconds.
    
	Time_Grid = np.arange(spec.Time_Step, Duration + spec.Time_Step, spec.Time_Step)
	Lat_Grid = [np.radians(l) for l in np.arange(-90, 90 + spec.Lat_Step, spec.Lat_Step)]

	assert len(Equilibrium_Config) == len(Lat_Grid), "Equilibrium-Config mismatched w/ latitude-grid size"
	Temperature_Profiles = [Equilibrium_Config]
	Time_Stamps = [spec.Initial_Year]

	def evolve_lat_temp(lat_idx: int, prior_tp: list[float], t: float):
		tROC =	calc_Temp_tROC(spec, Lat_Grid, prior_tp, lat_idx, t)
		new_temp = prior_tp[lat_idx] + tROC * spec.Time_Step
		assert new_temp >= 0.0, "Unphysical temperature"
		return new_temp

	with alive_bar(Time_Grid.size) as progress_bar:
		# ---------------------------------------------------------------------------- 
		for t in Time_Grid:
			TP_Prior = Temperature_Profiles[-1]
			TP_Now = [evolve_lat_temp(k,TP_Prior,t) for k in range(len(Lat_Grid))]
			Temperature_Profiles.append(TP_Now)
			Time_Stamps.append(to_timestamp(spec.Initial_Year, t))
			progress_bar()
		# ----------------------------------------------------------------------------

	return Sim_Result(spec, Time_Stamps, Lat_Grid, Temperature_Profiles)

# ---------------------------------------------------------------- #
# 						Model Calibration 							
# ---------------------------------------------------------------- #
def Optimize_Model():
	def compare(sim: Sim_Result):
		Gats = Average(sim, lambda x: x, (-90,90))
		Sampled_Timestamps = np.array(sim.times)[::365]
		Sampled_Gats = np.array(Gats)[::365]

		Similarity_Metric = 0.0
		for k in range(len(Sampled_Timestamps)):
			Timestamp, Gat = Sampled_Timestamps[k], Sampled_Gats[k]
			Lookup_Year = int(round(Timestamp))

			if Lookup_Year in Historic_Temperatures:
				Observed_Gat = Historic_Temperatures[Lookup_Year]
				Similarity_Metric += np.power((Gat - Observed_Gat),2.0)

		return Similarity_Metric

	def Objective_Func(params: list):
		ClimateRecord_Length = list(Historic_Temperatures)[-1] - list(Historic_Temperatures)[0]
		spec = Sim_Specification(ClimateRecord_Length, Alpha=params[0], Beta=params[1])
		result  = Simulate_Climate(spec)
		likeness = compare(result)

		Log_string = "Alpha={a}, Beta={b}, Likeness={l}".format(a=params[0],b=params[1],l=likeness)
		print(Log_string)
  
		return likeness

	# Note: Co2 Exponent (alpha) is unbounded (-inf, inf), but 
	# the retention factor (beta) must restricted to [0, inf).
	Param_Bounds = Bounds([-np.inf, 0.0], [np.inf, np.inf])
	Initial_Params = [Best_Alpha, Best_Beta]

	result = minimize(Objective_Func, x0=Initial_Params, bounds=Param_Bounds) 
	print(result)
 
def Find_Equilibrium():
	Spec = Sim_Specification(Stability_Duration, Initial_Year=list(Historic_Temperatures)[0]-Stability_Duration)
	Sim = Simulate_Climate(Spec)
	Gats = Average(Sim, lambda x: x, (-90,90))
	
	_, Observed_Gats = Serialise(Historic_Temperatures)
	Obs_Str = "First Gat in record: {x}".format(x=Observed_Gats[0])
	Gat_Str = "Equilibrium Gat: {x}".format(x=Gats[-1])

	print(Obs_Str)
	print(Gat_Str)
	print("Equilibrium Configuration:")
	print(Sim.tps[-1])