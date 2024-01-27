import matplotlib.pyplot as plt
from cycler import cycler
import numpy as np

import climate
from climate import SimClimate
from climate import years_to_seconds
from climate import days_to_seconds
from climate import A1FI_Pathway
from climate import B1_Pathway

from data import Climate_Temperatures
from data import Climate_CO2_Concentrations

import scienceplots

plt.style.use('science')

Default_LineStyles = cycler(linestyle=['-', '--', '-.', ':']) \
					+ cycler(color=['k','r','g','b'])
plt.rc('axes', prop_cycle=Default_LineStyles)

Fig_Size = (8.8, 5.5)

# ---------------------------------------------------------------- #
# 				Figure: Historic Climate GAT/CO2 Data	 							
# ---------------------------------------------------------------- #
def plot_GATCO2():
	Plot_Name = "gat_co2"
	figure, (co2, gat) = plt.subplots(nrows=2,ncols=1,sharex=True)

	# ------------------------------------------
	times = [yr for yr in Climate_Temperatures]
	gats, co2_ppms = [], []
	for t in times:
		co2_ppms.append(Climate_CO2_Concentrations[t])
		gats.append(Climate_Temperatures[t])
	# ------------------------------------------

	gat.plot(times, gats)
	co2.plot(times, co2_ppms)

	gat.set_xlabel("Year")
	gat.set_ylabel("Temperature (k)")
	co2.set_ylabel("Concentration (ppm)")
	
	plt.show()
	# plt.savefig(fname="thesis/{name}.pdf".format(name=Plot_Name), bbox_inches="tight")

# ---------------------------------------------------------------- #
# 					Figure: Solar Intensity 	 							
# ---------------------------------------------------------------- #
def plot_SolarIntensities():
	Plot_Name = "solar_rad"
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
	# plt.savefig(fname="thesis/{name}.pdf".format(name=Plot_Name), bbox_inches="tight")

# ---------------------------------------------------------------- #
# 					Figure: Albedo Model 	 							
# ---------------------------------------------------------------- #
def plot_AlbedoModel():
	Plot_Name = "albedo"
	figure, (axis) = plt.subplots(nrows=1,ncols=1,figsize=Fig_Size,sharex=False)

	# ------------------------------------------
	temps = np.linspace(250,300, 100)
	As = climate.calculate_Albedo(temps)
	axis.plot(temps, As)
	# ------------------------------------------

	axis.set_xlabel("Temperature (k)")
	axis.set_ylabel("Albedo")
	
	plt.show()
	# plt.savefig(fname="thesis/{name}.pdf".format(name=Plot_Name), bbox_inches="tight")
	
# ---------------------------------------------------------------- #
# 					Figure: Ice Model 	 							
# ---------------------------------------------------------------- #
def plot_IceModel():
	Plot_Name = "ice"
	figure, (axis) = plt.subplots(nrows=1,ncols=1,figsize=Fig_Size,sharex=False)

	# ------------------------------------------
	temps = np.linspace(200,300,500)
	ice_frac = climate.calc_IceFraction(temps)
	axis.plot(temps, ice_frac)
	# ------------------------------------------

	axis.set_xlabel("Temperature (k)")
	axis.set_ylabel("Ice Fraction")
	
	plt.show()
	# plt.savefig(fname="thesis/{name}.pdf".format(name=Plot_Name), bbox_inches="tight")

# ---------------------------------------------------------------- #
# 			Figure: Antarctica Altitude Correction 	 							
# ---------------------------------------------------------------- #
def plot_AntarcticaCorrection():
	Plot_Name = "antrc_corr"
	figure, (axis) = plt.subplots(nrows=1,ncols=1,figsize=Fig_Size,sharex=False)

	# ------------------------------------------
	def calc_IceCoverages(enable_altitude_correction: bool):
		sim_pack = SimClimate(climate.calc_HistoricPeriod(), altitude_correction=enable_altitude_correction)

		Ice_Coverages = []
		for tdist in sim_pack.tdists:
			n=0
			IceFracs = []
			while climate.Is_AntarcticLatitudeBand(sim_pack.lats[n]):
				fIce = climate.calc_IceFraction(tdist[n])
				IceFracs.append(fIce)
				n += 1

			Ice_Coverages.append(np.mean(IceFracs)) 

		return sim_pack.times, Ice_Coverages
	# ------------------------------------------

	times, coverage = calc_IceCoverages(False)
	axis.plot(times, coverage)

	times, coverage = calc_IceCoverages(True)
	axis.plot(times, coverage)


	axis.set_xlabel("Year")
	axis.set_ylabel("Ice Coverage")
	
	plt.savefig(fname="thesis/{name}.pdf".format(name=Plot_Name), bbox_inches="tight")

# ---------------------------------------------------------------- #
# 				Figure: Projected Co2 Emissions	 							
# ---------------------------------------------------------------- #
def plot_ProjectedCo2():
	Plot_Name = "projected_co2"
	figure, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False)

	obs_times, obs_co2 = [], []
	for t in Climate_CO2_Concentrations:
		emissions = Climate_CO2_Concentrations[t]
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
# 						Main Code Path	 							
# ---------------------------------------------------------------- #
# plot_GATCO2()
# plot_SolarIntensities()
# plot_AlbedoModel()
# plot_IceModel()
# plot_AntarcticaCorrection()
plot_ProjectedCo2()