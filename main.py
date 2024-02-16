from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt
from climate import *
import scienceplots
import numpy as np

from data import Historic_Temperatures
from data import Historic_Co2
from data import Serialise
from data import RCP85_Data
from data import RCP26_Data
from data import RCP45_Data
from data import RCP6_Data

# ---------------------------------------------------------------- #
# 						Utility Functions	 							
# ---------------------------------------------------------------- #
def Get_ClimateRecordLength():
    return list(Historic_Temperatures)[-1] - list(Historic_Temperatures)[0]

def Save_Figure(fig, filename: str):
    path = ".\\figures\\" + filename + ".png"
    fig.savefig(path, dpi=300, bbox_inches="tight")

# ---------------------------------------------------------------- #
# 				Figure: Climate Observational Data	 							
# ---------------------------------------------------------------- #
def Fig_ClimateData():
	fig, (co2, gat) = plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(6.4, 4))

	Times, Gats = Serialise(Historic_Temperatures)
	gat.plot(Times, Gats, "k-")

	Times, Concentrations = Serialise(Historic_Co2)
	co2.plot(Times, Concentrations, "k-")

	gat.set_ylabel("Temperature (k)")
	co2.set_ylabel("Concentration (ppm)")
	gat.set_xlabel("Year")
	
	Save_Figure(fig, "climate_data")

# ---------------------------------------------------------------- #
# 		  			Figure: Model Calibration	 							
# ---------------------------------------------------------------- #
def Fig_ModelCalibration():
	fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False,figsize=(6.4, 4))
	
	Duration = years_to_seconds(Stability_Duration + Get_ClimateRecordLength())
	Sim = Simulate_Climate(Sim_Specification(Duration))

	Gats = Average(Sim, lambda x: x, (-90,90))
	axis.plot(np.array(Sim.times)[::365], np.array(Gats)[::365], "r-", label="Climate Model")

	Observation_Times, Observed_Gats = Serialise(Historic_Temperatures)
	axis.plot(Observation_Times, Observed_Gats, "k-", label="Observational Data")
	
	axis.set_ylabel("Temperature (k)")
	axis.set_xlabel("Year")
	axis.legend()

	Save_Figure(fig, "model_fit")

# ---------------------------------------------------------------- #
# 				Figure: Projected Co2 Emissions	 							
# ---------------------------------------------------------------- #
def Fig_Co2Projections():
	fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False,figsize=(6.4, 4))

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

	Save_Figure(fig, "rcp_curves")

# ---------------------------------------------------------------- #
# 			Figure: Antarctica Altitude Correction 	 							
# ---------------------------------------------------------------- #
def Fig_Antarctica():
	fig, (Temp_axis, Alb_axis) = plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(6.4, 4))

	# -----------------------------------------------------------------------------------------
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
  
		LineWidth = 1.0
		LineStyle = "r-" if correction_on else "k-"
		Label = "Altitude Correction: On" if correction_on else "Altitude Correction: Off"
		
		Temp_axis.plot(Sampled_Times, Sampled_Temps, LineStyle, linewidth=LineWidth, label=Label)
		Alb_axis.plot(Sampled_Times, Sampled_Albedos, LineStyle, linewidth=LineWidth, label=Label)
	# -----------------------------------------------------------------------------------------

	plot(True)
	plot(False)
	Temp_axis.axhline(y=273, linestyle="-.", color='b', linewidth=1.5, alpha=0.5, label="Freezing point")

	Temp_axis.set_ylabel("Temperature (k)")
	Alb_axis.set_ylabel("Albedo")
	Alb_axis.set_xlabel("Year")
	Temp_axis.legend()

	Save_Figure(fig, "antarctic")

# ---------------------------------------------------------------- #
# 					Figure: Co2 Pathway Fitting 	 							
# ---------------------------------------------------------------- #
def Fig_Co2Interpolations():
	fig = plt.figure(figsize=(6.4, 4))
	fig.supxlabel("Year")
	fig.supylabel("Concentration (ppm)")

	a = fig.add_subplot(321)
	b = fig.add_subplot(322)
	c = fig.add_subplot(323)
	d = fig.add_subplot(324)
	l = fig.add_subplot(313)

	# --------------------------------------------------------------------
	def plot(Co2_Data: dict, Co2_Poly: Polynomial, Name: str, axis):
		Times, Concentrations = Serialise(Co2_Data)
		Sample_Times = np.linspace(Times[0], Times[-1], 100)
		Interpolated_Concentrations = Co2_Poly(Sample_Times)
		
		axis.plot(Times, Concentrations, ".", color="black", alpha=0.75, label=Name)
		axis.plot(Sample_Times, Interpolated_Concentrations, "--", color="red")
		axis.legend()
	# --------------------------------------------------------------------

	plot(RCP85_Data, RCP85, "RCP 8.5", a)
	plot(RCP6_Data, RCP6, "RCP 6", b)
	plot(RCP45_Data, RCP45, "RCP 4.5", c)
	plot(RCP26_Data, RCP26, "RCP 2.6", d)
	plot(Historic_Co2, RCP0, "Historic Co2", l)

	Save_Figure(fig, "co2_interp")
 
# ---------------------------------------------------------------- #
# 					Figure: Forecasts 	 							
# ---------------------------------------------------------------- #
def Fig_Forecasts():
	# Simulations
	# -------------------------------------------------------------------
	Target_Year = 2100
	Duration = years_to_seconds(Target_Year - (list(Historic_Temperatures)[0] - Stability_Duration))

	Sim_RCP85 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP85))
	Sim_RCP6 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP6))
	Sim_RCP45 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP45))
	Sim_RCP26 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP26))
	# -------------------------------------------------------------------

	# Colour Map
	# -------------------------------------------------------------------
	RCP85_Col = "firebrick"
	RCP6_Col = "darkorange"
	RCP45_Col = "royalblue"
	RCP26_Col = "seagreen"
	# -------------------------------------------------------------------

	# Plot global temperature distributions
	# -------------------------------------------------------------------
	fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False,figsize=(6.4, 4))
  
	def plot_dist(sim: Sim_Result, name: str, col: str):
		axis.plot(sim.tps[-1], np.degrees(sim.lats), label=name, color=col)

	plot_dist(Sim_RCP85, "RCP 8.5", RCP85_Col)
	plot_dist(Sim_RCP6, "RCP 6", RCP6_Col)
	plot_dist(Sim_RCP45, "RCP 4.5", RCP45_Col)
	plot_dist(Sim_RCP26, "RCP 2.6", RCP26_Col)

	axis.set_xlabel("Temperature (k)")
	axis.set_ylabel("Latitude")
	axis.legend()

	Save_Figure(fig, "tdist_forecast")
 	# -------------------------------------------------------------------

	# Plot global average temperatures
	# -------------------------------------------------------------------
	fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False,figsize=(6.4, 4))
	
	def plot_gats(sim: Sim_Result, name: str, col: str):
		gats = Average(sim, lambda x: x, (-90,90))
		Sample_Times = np.array(sim.times)[::365]
		Sample_Gats = np.array(gats)[::365]
		axis.plot(Sample_Times, Sample_Gats, label=name, color=col)

	plot_gats(Sim_RCP85, "RCP 8.5", RCP85_Col)
	plot_gats(Sim_RCP6, "RCP 6", RCP6_Col)
	plot_gats(Sim_RCP45, "RCP 4.5", RCP45_Col)
	plot_gats(Sim_RCP26, "RCP 2.6", RCP26_Col)

	Times, Gats = Serialise(Historic_Temperatures)
	axis.plot(Times, Gats, "k--")	

	plt.ylabel("Temperature (k)")
	plt.xlabel("Year")
	plt.legend()

	Save_Figure(fig, "gat_forecast")
	# -------------------------------------------------------------------


# ---------------------------------------------------------------- #
# 						Main Code Path	 							
# ---------------------------------------------------------------- #
plt.style.use('science')

# Thesis Ready:
# Fig_ClimateData()
# Fig_ModelCalibration()
# Fig_Co2Projections()
# Fig_Antarctica()
# Fig_Co2Interpolations()
Fig_Forecasts()