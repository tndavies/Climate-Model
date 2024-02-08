import matplotlib.pyplot as plt
import numpy as np
import scienceplots

from climate import RCP_85
from climate import RCP_6
from climate import RCP_45
from climate import RCP_26
from climate import RCP_0

from climate import SimClimate
from climate import years_to_seconds
from climate import days_to_seconds
from climate import to_kelvin
from climate import calc_SphericalAverage
from climate import Starting_Year
from climate import Antarctic_Bounds
from climate import calc_Albedo

from data import Serialise
from data import Historic_Temperatures
from data import Historic_Co2
from data import RCP85_Data
from data import RCP6_Data
from data import RCP45_Data
from data import RCP26_Data

# ---------------------------------------------------------------- #
# 				Figure: Climate Observational Data	 							
# ---------------------------------------------------------------- #
def ClimateData_Plot():
	figure, (co2, gat) = plt.subplots(nrows=2,ncols=1,sharex=True)

	Times = [t for t in Historic_Temperatures]
	Gats = [Historic_Temperatures[t] for t in Times]
	gat.plot(Times, Gats, "k-")

	Times = [t for t in Historic_Co2]
	Concentrations = [Historic_Co2[t] for t in Times]
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
	
	Sim_Length = years_to_seconds(list(Historic_Temperatures)[-1] - Starting_Year)

	Sim = SimClimate(Sim_Length, None)
	Gats = calc_SphericalAverage(Sim, lambda Tband: Tband, (-90,90))
	axis.plot(np.array(Sim.times)[::365], np.array(Gats)[::365], "r-", label="Climate Model")

	Times = [t for t in Historic_Temperatures]
	Observed_Gats = [to_kelvin(Historic_Temperatures[t]) for t in Times]
	axis.plot(Times, Observed_Gats, "k-", label="Observational Data")
	
	axis.set_xlabel("Year")
	axis.set_ylabel("Temperature (k)")

	axis.legend()

	plt.show()

# ---------------------------------------------------------------- #
# 				Figure: Projected Co2 Emissions	 							
# ---------------------------------------------------------------- #
def ProjectedEmissions_Plot():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	Observation_Times = [t for t in Historic_Co2]
	Observed_Concentrations = [Historic_Co2[t] for t in Observation_Times]
	axis.plot(Observation_Times, Observed_Concentrations, "--", color="black", label="Observational Data")

	Times = [t for t in RCP85_Data]
	Concentrations = [RCP85_Data[t] for t in Times]
	axis.plot(Times, Concentrations, color="firebrick", label="RCP 8.5")

	Times = [t for t in RCP6_Data]
	Concentrations = [RCP6_Data[t] for t in Times]
	axis.plot(Times, Concentrations, color="darkorange", label="RCP 6")

	Times = [t for t in RCP45_Data]
	Concentrations = [RCP45_Data[t] for t in Times]
	axis.plot(Times, Concentrations, color="royalblue", label="RCP 4.5")

	Times = [t for t in RCP26_Data]
	Concentrations = [RCP26_Data[t] for t in Times]
	axis.plot(Times, Concentrations, color="seagreen", label="RCP 2.6")

	axis.set_xlabel("Year")
	axis.set_ylabel("Concentration (ppm)")

	axis.legend()

	plt.show()

# ---------------------------------------------------------------- #
# 			Figure: Antarctica Altitude Correction 	 							
# ---------------------------------------------------------------- #
def AntarcticaCorrection_Plot():
	figure, (Temp_axis, Alb_axis) = plt.subplots(nrows=2,ncols=1,sharex=True)

	Sim_Length = years_to_seconds(list(Historic_Temperatures)[-1] - Starting_Year)
	uSim = SimClimate(Sim_Length, None, enable_AltCorr=False)
	cSim = SimClimate(Sim_Length, None, enable_AltCorr=True)

	uTimes = np.array(uSim.times)
	uTemps = calc_SphericalAverage(uSim, lambda T: T, Antarctic_Bounds)
	uAlbedos = calc_SphericalAverage(uSim, calc_Albedo, Antarctic_Bounds)

	cTimes = np.array(cSim.times)
	cTemps = calc_SphericalAverage(cSim, lambda T: T, Antarctic_Bounds)
	cAlbedos = calc_SphericalAverage(cSim, calc_Albedo, Antarctic_Bounds)

	n = int(365 / 2)

	lbl = "Altitude Correction: Off"
	Sampled_Times = np.array(uTimes)[::n]
	Sampled_Temps = np.array(uTemps)[::n]
	Sampled_Albedos = np.array(uAlbedos)[::n]
	Temp_axis.plot(Sampled_Times, Sampled_Temps, "k-", label=lbl)
	Alb_axis.plot(Sampled_Times, Sampled_Albedos, "k-", label=lbl)

	lbl = "Altitude Correction: On"
	Sampled_Times = np.array(cTimes)[::n]
	Sampled_Temps = np.array(cTemps)[::n]
	Sampled_Albedos = np.array(cAlbedos)[::n]
	Temp_axis.plot(Sampled_Times, Sampled_Temps, "r-", label=lbl)
	Alb_axis.plot(Sampled_Times, Sampled_Albedos, "r-", label=lbl)

	Temp_axis.axhline(y=273, linestyle="-.", color='b', linewidth=3.0, alpha=0.5, label="Freezing point")

	Temp_axis.set_ylabel("Temperature (k)")
	Alb_axis.set_ylabel("Albedo")
	Alb_axis.set_xlabel("Year")
	
	Temp_axis.legend()

	plt.show()

# ---------------------------------------------------------------- #
# 					Figure: Co2 Pathway Fitting 	 							
# ---------------------------------------------------------------- #
def Co2PathwayFitting_Plot():
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

	plot(RCP85_Data, RCP_85, "RCP 8.5", a)
	plot(RCP6_Data, RCP_6, "RCP 6", b)
	plot(RCP45_Data, RCP_45, "RCP 4.5", c)
	plot(RCP26_Data, RCP_26, "RCP 2.6", d)
	plot(Historic_Co2, RCP_0, "Historic Co2", l)

	plt.show()


# ---------------------------------------------------------------- #
# 					Figure: Temperature Forecasts 	 							
# ---------------------------------------------------------------- #
def TemperatureForecasts_Plot():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	Sim_Length = years_to_seconds(2100 - Starting_Year)
	Sim_RCP85 = SimClimate(Sim_Length, RCP_85)
	Sim_RCP6 = SimClimate(Sim_Length, RCP_6)
	Sim_RCP45 = SimClimate(Sim_Length, RCP_45)
	Sim_RCP26 = SimClimate(Sim_Length, RCP_26)

	def plot(sim, Name: str, col: str):
		gats = calc_SphericalAverage(sim, lambda x: x, (-90,90))
		Sample_Times = np.array(sim.times)[::365]
		Sample_Gats = np.array(gats)[::365]
		axis.plot(Sample_Times, Sample_Gats, label=Name, color=col)

	plot(Sim_RCP85, "RCP 8.5", "firebrick")
	plot(Sim_RCP6, "RCP 6", "darkorange")
	plot(Sim_RCP45, "RCP 4.5", "royalblue")
	plot(Sim_RCP26, "RCP 2.6", "seagreen")

	Times, Gats = Serialise(Historic_Temperatures)
	Gats_Kelvin = [to_kelvin(gat) for gat in Gats]
	axis.plot(Times, Gats_Kelvin, "--", color="black")	

	plt.xlabel("Year")
	plt.ylabel("Temperature (k)")
	plt.legend()
	plt.show()

# ---------------------------------------------------------------- #
# 						Main Code Path	 							
# ---------------------------------------------------------------- #
# plt.style.use('science')

# ClimateData_Plot()
ModelCalibration_Plot()
# ProjectedEmissions_Plot()
# AntarcticaCorrection_Plot()
# Co2PathwayFitting_Plot()
# TemperatureForecasts_Plot()
