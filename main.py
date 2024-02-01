import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial
from cycler import cycler
import numpy as np
import scienceplots

from climate import SimClimate
from climate import years_to_seconds
from climate import days_to_seconds
from climate import A1FI_Pathway
from climate import B1_Pathway
from climate import to_kelvin
from climate import calc_SphericalAverage
from climate import Starting_Year
from climate import get_HistoricEmissions
from climate import Antarctic_Bounds
from climate import calc_Albedo

from data import Historic_Temperatures
from data import Historic_Co2
from data import A1FI_Projections
from data import A1T_Projections
from data import B1_Projections

# ---------------------------------------------------------------- #
# 				Figure: Historic Climate GAT/CO2 Data	 							
# ---------------------------------------------------------------- #
def GatCo2Data_Plot():
	figure, (co2, gat) = plt.subplots(nrows=2,ncols=1,sharex=True)

	Times = [t for t in Historic_Temperatures]
	Gats = [Historic_Temperatures[t] for t in Times]
	gat.plot(Times, Gats, "k-")

	Times = [t for t in Historic_Co2]
	Concentrations = [Historic_Co2[t] for t in Times]
	co2.plot(Times, Concentrations, "k--", label="Observational Data")
	co2.plot(Times, get_HistoricEmissions(np.array(Times)), "r-", label="Interpolation")

	gat.set_ylabel("Temperature (k)")
	co2.set_ylabel("Concentration (ppm)")
	gat.set_xlabel("Year")
	
	co2.legend()

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
# @todo: this plot needs some TLC.
def ProjectedEmissions_Plot():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	Observation_Times = [t for t in Historic_Co2]
	Observed_Concentrations = [Historic_Co2[t] for t in Observation_Times]
	axis.plot(Observation_Times, Observed_Concentrations, "k-", label="Observational Data")

	Times = [t for t in A1FI_Projections]
	Concentrations = [A1FI_Projections[t] for t in Times]
	axis.plot(Times, Concentrations, "r--", label="A1FI Pathway")

	Times = [t for t in A1T_Projections]
	Concentrations = [A1T_Projections[t] for t in Times]
	axis.plot(Times, Concentrations, "y--", label="A1T Pathway")

	Times = [t for t in B1_Projections]
	Concentrations = [B1_Projections[t] for t in Times]
	axis.plot(Times, Concentrations, "g--", label="B1 Pathway")

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

	n = 365

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
# 					Figure: Solar Intensity 	 							
# ---------------------------------------------------------------- #
def plot_SolarIntensities():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,figsize=Fig_Size,sharex=False)

	# ------------------------------------------
	Times = np.linspace(0, 365, 500)
	Latitudes = [90,0,-90]
	for l in Latitudes:
		LatBand_Intensities = []
		for t in Times:
			I = climate.calc_AverageRadiation(np.radians(l), days_to_seconds(t))
			LatBand_Intensities.append(I)

		Plot_Label = r"${lat}^\circ$".format(lat=str(l))
		axis.plot(Times, LatBand_Intensities, label=Plot_Label)
	# ------------------------------------------

	axis.set_xlabel("Day")
	axis.set_ylabel("Intensity " + r'(W$m^{-2}$)')
	
	plt.show()

# ---------------------------------------------------------------- #
# 					Figure: Albedo Model 	 							
# ---------------------------------------------------------------- #
def plot_AlbedoModel():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,figsize=Fig_Size,sharex=False)

	# ------------------------------------------
	temps = np.linspace(250,300, 100)
	As = climate.calculate_Albedo(temps)
	axis.plot(temps, As)
	# ------------------------------------------

	axis.set_xlabel("Temperature (k)")
	axis.set_ylabel("Albedo")
	
	plt.show()
	
# ---------------------------------------------------------------- #
# 					Figure: Ice Model 	 							
# ---------------------------------------------------------------- #
def plot_IceModel():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,figsize=Fig_Size,sharex=False)

	# ------------------------------------------
	temps = np.linspace(200,300,500)
	ice_frac = climate.calc_IceFraction(temps)
	axis.plot(temps, ice_frac)
	# ------------------------------------------

	axis.set_xlabel("Temperature (k)")
	axis.set_ylabel("Ice Fraction")
	
	plt.show()

# ---------------------------------------------------------------- #
# 						Main Code Path	 							
# ---------------------------------------------------------------- #
# plt.style.use('science')

# GatCo2Data_Plot()
# ModelCalibration_Plot()
# ProjectedEmissions_Plot()
AntarcticaCorrection_Plot()
