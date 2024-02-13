from climate import calculate_IRCooling
from climate import Sim_Specification
from climate import Simulate_Climate
from climate import Antarctic_Bounds
from climate import years_to_seconds
from climate import calc_Albedo
from climate import Best_Alpha
from climate import Best_Beta
from climate import Average
from climate import RCP85
from climate import RCP45
from climate import RCP26
from climate import RCP6
from climate import RCP0

from data import Historic_Temperatures
from data import Historic_Co2
from data import Serialise
from data import RCP85_Data
from data import RCP26_Data
from data import RCP45_Data
from data import RCP6_Data

from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt
import scienceplots
import numpy as np

# ---------------------------------------------------------------- #
# 						Utility Functions	 							
# ---------------------------------------------------------------- #
def Get_ClimateRecordLength():
    return list(Historic_Temperatures)[-1] - list(Historic_Temperatures)[0]


# ---------------------------------------------------------------- #
# 				Figure: Climate Observational Data	 							
# ---------------------------------------------------------------- #
def ClimateData_Plot():
	figure, (co2, gat) = plt.subplots(nrows=2,ncols=1,sharex=True)

	Times, Gats = Serialise(Historic_Temperatures)
	gat.plot(Times, Gats, "k-")

	Times, Concentrations = Serialise(Historic_Co2)
	co2.plot(Times, Concentrations, "k-")

	gat.set_ylabel("Temperature (k)")
	co2.set_ylabel("Concentration (ppm)")
	gat.set_xlabel("Year")
	
	plt.show()

# ---------------------------------------------------------------- #
# 		  			Figure: Model Calibration	 							
# ---------------------------------------------------------------- #
def ModelCalibration_Plot():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)
	
	Duration = years_to_seconds(Get_ClimateRecordLength())
	Sim = Simulate_Climate(Sim_Specification(Duration))

	Gats = Average(Sim, lambda x: x, (-90,90))
	axis.plot(np.array(Sim.times)[::365], np.array(Gats)[::365], "r-", label="Climate Model")

	Observation_Times, Observed_Gats = Serialise(Historic_Temperatures)
	axis.plot(Observation_Times, Observed_Gats, "k-", label="Observational Data")
	
	axis.set_ylabel("Temperature (k)")
	axis.set_xlabel("Year")
	axis.legend()

	plt.show()

# ---------------------------------------------------------------- #
# 				Figure: Projected Co2 Emissions	 							
# ---------------------------------------------------------------- #
def ProjectedEmissions_Plot():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	Observation_Times, Observed_Concentrations = Serialise(Historic_Co2)
	axis.plot(Observation_Times, Observed_Concentrations, "--", color="black", label="Observational Data")

	Times, Concentrations = Serialise(RCP85_Data)
	axis.plot(Times, Concentrations, color="firebrick", label="RCP 8.5")

	Times, Concentrations = Serialise(RCP6_Data)
	axis.plot(Times, Concentrations, color="darkorange", label="RCP 6")

	Times, Concentrations = Serialise(RCP45_Data)
	axis.plot(Times, Concentrations, color="royalblue", label="RCP 4.5")

	Times, Concentrations = Serialise(RCP26_Data)
	axis.plot(Times, Concentrations, color="seagreen", label="RCP 2.6")

	axis.set_ylabel("Concentration (ppm)")
	axis.set_xlabel("Year")
	axis.legend()

	plt.show()

# ---------------------------------------------------------------- #
# 			Figure: Antarctica Altitude Correction 	 							
# ---------------------------------------------------------------- #
def AntarcticaCorrection_Plot():
	figure, (Temp_axis, Alb_axis) = plt.subplots(nrows=2,ncols=1,sharex=True)

	def plot(correction_on: bool):
		Duration = years_to_seconds(Get_ClimateRecordLength())
		Sim = Simulate_Climate(Sim_Specification(Duration, Altitude_Correction=correction_on))

		Times = np.array(Sim.times)
		Temps = Average(Sim, lambda x: x, Antarctic_Bounds)
		Albedos = Average(Sim, calc_Albedo, Antarctic_Bounds)

		SamplesPerYear = 2
		N = int(365 / SamplesPerYear)
		Sampled_Times = np.array(Times)[::N]
		Sampled_Temps = np.array(Temps)[::N]
		Sampled_Albedos = np.array(Albedos)[::N]
  
		LineStyle = "r-" if correction_on else "k-"
		Label = "Altitude Correction: On" if correction_on else "Altitude Correction: Off"
		Temp_axis.plot(Sampled_Times, Sampled_Temps, LineStyle, label=Label)
		Alb_axis.plot(Sampled_Times, Sampled_Albedos, LineStyle, label=Label)

	Temp_axis.axhline(y=273, linestyle="-.", color='b', linewidth=3.0, alpha=0.5, label="Freezing point")
	plot(True)
	plot(False)

	Temp_axis.set_ylabel("Temperature (k)")
	Alb_axis.set_ylabel("Albedo")
	Alb_axis.set_xlabel("Year")
	Temp_axis.legend()

	plt.show()

# ---------------------------------------------------------------- #
# 					Figure: Co2 Pathway Fitting 	 							
# ---------------------------------------------------------------- #
def Co2PathwayInterpolation_Plot():
	fig = plt.figure()
	fig.supxlabel("Year")
	fig.supylabel("Concentration (ppm)")

	a = fig.add_subplot(321)
	b = fig.add_subplot(322)
	c = fig.add_subplot(323)
	d = fig.add_subplot(324)
	l = fig.add_subplot(313)

	def plot(Co2_Data: dict, Co2_Poly: Polynomial, Name: str, axis):
		Times, Concentrations = Serialise(Co2_Data)
		Sample_Times = np.linspace(Times[0], Times[-1], 100)
		Interpolated_Concentrations = Co2_Poly(Sample_Times)
		
		axis.plot(Times, Concentrations, ".", color="black", alpha=0.75, label=Name)
		axis.plot(Sample_Times, Interpolated_Concentrations, "--", color="red")
		axis.legend()

	plot(RCP85_Data, RCP85, "RCP 8.5", a)
	plot(RCP6_Data, RCP6, "RCP 6", b)
	plot(RCP45_Data, RCP45, "RCP 4.5", c)
	plot(RCP26_Data, RCP26, "RCP 2.6", d)
	plot(Historic_Co2, RCP0, "Historic Co2", l)

	plt.show()

# ---------------------------------------------------------------- #
# 					Figure: Temperature Forecasts 	 							
# ---------------------------------------------------------------- #
def TemperatureForecasts_Plot():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	Target_Year = 2100
	Duration = years_to_seconds(Target_Year - list(Historic_Temperatures)[0])
 
	Sim_RCP85 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP85))
	Sim_RCP6 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP6))
	Sim_RCP45 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP45))
	Sim_RCP26 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP26))

	def plot(sim, Name: str, col: str):
		gats = Average(sim, lambda x: x, (-90,90))
		Sample_Times = np.array(sim.times)[::365]
		Sample_Gats = np.array(gats)[::365]
		axis.plot(Sample_Times, Sample_Gats, label=Name, color=col)

	plot(Sim_RCP85, "RCP 8.5", "firebrick")
	plot(Sim_RCP6, "RCP 6", "darkorange")
	plot(Sim_RCP45, "RCP 4.5", "royalblue")
	plot(Sim_RCP26, "RCP 2.6", "seagreen")

	Times, Gats = Serialise(Historic_Temperatures)
	axis.plot(Times, Gats, "k--")	

	plt.ylabel("Temperature (k)")
	plt.xlabel("Year")
	plt.legend()
	plt.show()

# ---------------------------------------------------------------- #
# 					Figure: Model Stability 	 							
# ---------------------------------------------------------------- #
def ModelStability_Plot():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	for delta in np.linspace(0, 95, 3):
		yr0 = list(Historic_Temperatures)[0] - delta
		Duration = years_to_seconds(list(Historic_Temperatures)[-1] - yr0)
		spec = Sim_Specification(Duration, Initial_Year = yr0)
		Sim = Simulate_Climate(spec)

		Gats = Average(Sim, lambda x: x, (-90,90))
		axis.plot(np.array(Sim.times)[::365], np.array(Gats)[::365], "--")

	Observation_Times, Observed_Gats = Serialise(Historic_Temperatures)
	axis.plot(Observation_Times, Observed_Gats, "k-", label="Observational Data")

	plt.legend()
	plt.show()

# Heatmap over Earth map, for different projections to see how heat profiles vary between them.

# def CoolingModel_Plot():
# 	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

# 	SampleRate = 35
# 	Timestamps, _ = Serialise(Historic_Co2) 
# 	Sample_Timestamps = np.array(Timestamps)[::SampleRate]

# 	Temperatures = np.linspace(273,300)		

# 	for Timestamp in Sample_Timestamps:
# 		Coolings = calculate_IRCooling(Timestamp, 0.0, None, Best_Alpha, Best_Beta, Temperatures)
# 		Co2_Concentration = round(Historic_Co2[Timestamp], 2)
# 		axis.plot(Temperatures, Coolings, label=str(Co2_Concentration) + " ppm") 

# 	axis.legend()
# 	axis.set_ylabel("Intensity " + r"($W \; m^{-2}$)")
# 	axis.set_xlabel("Temperature (k)")
# 	plt.show()
 


# ---------------------------------------------------------------- #
# 						Main Code Path	 							
# ---------------------------------------------------------------- #
# plt.style.use('science')

# Thesis Ready:
# ClimateData_Plot()
# ModelCalibration_Plot()
# ProjectedEmissions_Plot()
# AntarcticaCorrection_Plot()
# Co2PathwayInterpolation_Plot()
TemperatureForecasts_Plot()

# In Progress:
# CoolingModel_Plot()
# ModelStability_Plot()

