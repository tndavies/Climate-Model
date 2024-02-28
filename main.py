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
	
	Sim = Simulate_Climate(Sim_Specification(Get_ClimateRecordLength()))
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

	# Times, Concentrations = Serialise(RCP6_Data)
	# axis.plot(Times, Concentrations, color="darkorange", label="RCP 6")

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
		Spec = Sim_Specification(Get_ClimateRecordLength(), Altitude_Correction=correction_on)
		Sim = Simulate_Climate(Spec)

		Times = np.array(Sim.times)
		Temps = Average(Sim, lambda x: x, Antarctic_Bounds)
		Albedos = Average(Sim, calc_Albedo, Antarctic_Bounds)

		SamplesPerYear = 2
		N = int(365 / SamplesPerYear)
		Sampled_Times = np.array(Times)[::N]
		Sampled_Temps = np.array(Temps)[::N]
		Sampled_Albedos = np.array(Albedos)[::N]
  
		LineWidth = 0.5
		LineStyle = "r-" if correction_on else "k-"
		Label = "Altitude Correction: On" if correction_on else "Altitude Correction: Off"
		
		Temp_axis.plot(Sampled_Times, Sampled_Temps, LineStyle, linewidth=LineWidth, label=Label)
		Alb_axis.plot(Sampled_Times, Sampled_Albedos, LineStyle, linewidth=LineWidth, label=Label)
	# -----------------------------------------------------------------------------------------

	plot(True)
	plot(False)
	Temp_axis.axhline(y=273, linestyle="-.", color='b', linewidth=1.2, alpha=0.75, label="Freezing point")

	Temp_axis.set_ylabel("Temperature (k)")
	Alb_axis.set_ylabel("Albedo")
	Alb_axis.set_xlabel("Year")
	# Temp_axis.legend()

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
# @think: Is this plotting the year 2099 average, or the year 2100 average? (temp dists plot)
# @todo: double check our averging code is working correctly.
def Fig_Forecasts():

	# Simulations
	Target_Year = 2100
	Duration = Target_Year - list(Historic_Temperatures)[0]
	Sim_RCP85 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP85))
	Sim_RCP45 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP45))
	Sim_RCP26 = Simulate_Climate(Sim_Specification(Duration, RCP=RCP26))
	Sim_Context = Sim_RCP85
 
	# Plot global average temperatures
  	# ----------------------------------------------------------------
	fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False,figsize=(6.4, 4))

	Time_Threshold = (list(Historic_Temperatures)[-1] - Sim_Context.spec.Initial_Year)*31536000
	N = round(Time_Threshold / Sim_Context.spec.Time_Step)

	def plot_gats(Sim: Sim_Result, Label: str, Colour: str):
		gats = np.array(Average(Sim, lambda x: x, (-90,90)))[N::365]
		times = np.array(Sim.times)[N::365]
		axis.plot(times, gats, "--", color=Colour, label=Label)

	gats = np.array(Average(Sim_Context, lambda x: x, (-90,90)))[:N:365]
	times = np.array(Sim_Context.times)[:N:365]
	axis.plot(times, gats, "k-")

	plot_gats(Sim_RCP85, "RCP 8.5", "firebrick")
	plot_gats(Sim_RCP45, "RCP 4.5", "royalblue")
	plot_gats(Sim_RCP26, "RCP 2.6", "seagreen")
  
	plt.ylabel("Temperature (k)")
	plt.xlabel("Year")
	plt.legend()

	Save_Figure(fig, "gat_forecast")

	# -------------------------------------------------------------------
 
	# Plot global temperature distributions
	# fig, (axis) = plt.subplots(nrows=1,ncols=1,sharex=False,figsize=(6.4, 4))
  
	# def plot_dist(sim: Sim_Result, name: str, col: str, line_pattern: str = "-"):
	# 	Average_Temps = []

	# 	Steps_Per_Year = int(3.154e7 / sim.spec.Time_Step) + 1
	# 	Temp_Dists = [sim.tps[-k] for k in range(1,Steps_Per_Year + 1)]
  
	# 	for k in range(len(sim.lats)):
	# 		buffer = [arr[k] for arr in Temp_Dists]
	# 		Average_Temps.append(np.mean(buffer))
  
	# 	axis.plot(np.degrees(sim.lats), Average_Temps, label=name, color=col, linestyle=line_pattern, linewidth=0.8)

	# plot_dist(Sim_Context, "2023", Context_Col, line_pattern="--")
	# plot_dist(Sim_RCP85, "RCP 8.5", RCP85_Col)
	# # plot_dist(Sim_RCP6, "RCP 6", RCP6_Col)
	# plot_dist(Sim_RCP45, "RCP 4.5", RCP45_Col)
	# plot_dist(Sim_RCP26, "RCP 2.6", RCP26_Col)

	# axis.set_xlabel("Latitude")
	# axis.set_ylabel("Temperature (k)")
	# axis.legend()

	# Save_Figure(fig, "tdist_forecast")
	
# ---------------------------------------------------------------- #
# 						Main Code Path	 							
# ---------------------------------------------------------------- #
plt.style.use('science')

# [Thesis Ready]:
# Fig_ClimateData()
# Fig_ModelCalibration()
# Fig_Co2Projections()
# Fig_Antarctica()
# Fig_Co2Interpolations()
Fig_Forecasts()
