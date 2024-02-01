import matplotlib.pyplot as plt
from cycler import cycler
import numpy as np

from climate import SimClimate
from climate import years_to_seconds
from climate import days_to_seconds
from climate import A1FI_Pathway
from climate import B1_Pathway
from climate import to_kelvin
from climate import calc_SphericalAverage
from climate import Starting_Year
from data import Historic_Temperatures
from data import Historic_Co2Emissions
from climate import get_HistoricEmissions

import scienceplots

# plt.style.use('science')

# Default_LineStyles = cycler(linestyle=['-', '--', '-.', ':']) \
# 					+ cycler(color=['k','r','g','b'])
# plt.rc('axes', prop_cycle=Default_LineStyles)

Fig_Size = (8.8, 5.5)

# ---------------------------------------------------------------- #
# 				Figure: Historic Climate GAT/CO2 Data	 							
# ---------------------------------------------------------------- #
def plot_GATCO2():
	figure, (co2, gat) = plt.subplots(nrows=2,ncols=1,sharex=True)

	# ------------------------------------------
	times = [yr for yr in Historic_Temperatures]
	gats, co2_ppms = [], []
	for t in times:
		co2_ppms.append(Historic_Co2[t])
		gats.append(Historic_Temperatures[t])
	# ------------------------------------------

	gat.plot(times, gats)
	co2.plot(times, co2_ppms)

	gat.set_xlabel("Year")
	gat.set_ylabel("Temperature (k)")
	co2.set_ylabel("Concentration (ppm)")
	
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
# 			Figure: Antarctica Altitude Correction 	 							
# ---------------------------------------------------------------- #

def plot_AntarcticaCorrection():
	figure, (Temp_axis, Alb_axis) = plt.subplots(nrows=2,ncols=1,figsize=Fig_Size,sharex=True)

	Sim_End = 2022
	uSim = SimClimate(None, Sim_End, enable_AltCorr=False)
	cSim = SimClimate(None, Sim_End, enable_AltCorr=True)

	uTimes = np.array(uSim.times)
	uTemps = np.array(calc_AntarcticAverageTemperature(uSim))
	uICovs = np.array(calc_AntarcticAlbedo(uSim))
	
	cTimes = np.array(cSim.times)
	cTemps = np.array(calc_AntarcticAverageTemperature(cSim))
	cICovs = np.array(calc_AntarcticAlbedo(cSim))

	n = int(365 / 12)

	Temp_axis.plot(uTimes[::n], uTemps[::n], "k--", linewidth=1.0, label="No altitude correction")
	Alb_axis.plot(uTimes[::n], uICovs[::n], "k--", linewidth=1.0, label="No altitude correction")

	Temp_axis.plot(cTimes[::n], cTemps[::n], "r-", linewidth=1.2, label="With altitude correction")
	Alb_axis.plot(cTimes[::n], cICovs[::n], "r-", linewidth=1.2, label="With altitude correction")

	Temp_axis.axhline(y=273, linestyle=":", color='b', label="Freezing point")

	Temp_axis.set_ylabel("Temperature (k)")
	Alb_axis.set_ylabel("Albedo")
	Alb_axis.set_xlabel("Year")
	Temp_axis.legend()

	plt.show()

# ---------------------------------------------------------------- #
# 				Figure: Projected Co2 Emissions	 							
# ---------------------------------------------------------------- #
def plot_ProjectedCo2():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	obs_times, obs_co2 = [], []
	for t in Historic_Co2:
		emissions = Historic_Co2[t]
		obs_co2.append(emissions)
		obs_times.append(t)

	proj_times = np.linspace(obs_times[-1], 2100)

	A1FI_co2 = []
	Total_Emissions = obs_co2[-1]
	for t in proj_times:
		Total_Emissions += A1FI_Pathway(t)
		A1FI_co2.append(Total_Emissions)

	B1_co2 = []
	Total_Emissions = obs_co2[-1]
	for t in proj_times:
		Total_Emissions += B1_Pathway(t)
		B1_co2.append(Total_Emissions)


	axis.plot(obs_times, obs_co2, "k-", label="Observed")
	axis.plot(proj_times, A1FI_co2, "r--", label="A1FI")
	axis.plot(proj_times, B1_co2, "g--", label="B1")

	axis.set_xlabel("Year")
	axis.set_ylabel("Concentration (ppm)")
	plt.legend()
	plt.show()

# ---------------------------------------------------------------- #
# 		  			Figure: Model Calibration	 							
# ---------------------------------------------------------------- #
def plot_ModelCalibration():
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)
	
	Sim_Length = years_to_seconds(list(Historic_Temperatures)[-1] - Starting_Year)
	Sim = SimClimate(Sim_Length, None)

	Gats = calc_SphericalAverage(Sim, lambda Tband: Tband, (-90,90))
	axis.plot(np.array(Sim.times)[::365], np.array(Gats)[::365], "r-", linewidth=1.2, label="Climate Model")

	HR_Times = [t for t in Historic_Temperatures]
	HR_Gats = [to_kelvin(Historic_Temperatures[t]) for t in HR_Times]
	axis.plot(HR_Times, HR_Gats, "k-", linewidth=2.5, alpha=0.75, label="Historic Record")

	axis.set_xlabel("Year")
	axis.set_ylabel("Temperature (k)")

	plt.legend()
	plt.show()

# ---------------------------------------------------------------- #
# 		  			Figure: Co2 Data Interpolation	 							
# ---------------------------------------------------------------- #
def plot_InterpolatedCo2Data():
	times, co2 = [],[]
	for t in Historic_Co2Emissions:
		times.append(t)
		co2.append(Historic_Co2Emissions[t])

	plt.figure()
	plt.plot(times,co2, "k--", label="co2 observations")
	plt.plot(times, get_HistoricEmissions(np.array(times)), "r-", label="Interpolation")
	plt.legend()
	plt.show()

# ---------------------------------------------------------------- #
# 						Main Code Path	 							
# ---------------------------------------------------------------- #
# plot_GATCO2()
# plot_SolarIntensities()
# plot_AlbedoModel()
# plot_IceModel()
# plot_AntarcticaCorrection()
# plot_ProjectedCo2()
plot_ModelCalibration()

