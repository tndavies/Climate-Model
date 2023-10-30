from alive_progress import alive_bar
import matplotlib.pyplot as plt
from pde import grad, laplace
import numpy as np
import cmocean
import flux
import pde

# ======================================================== #
# 'dur' is in units of days.

def DeclinationPlot(dur):
	times_d = np.linspace(0, dur, 2*dur)
	times_s = np.array([t_d * 86400 for t_d in times_d])
	decls = flux.Calc_Declination(times_s)

	plt.figure()
	plt.title("Declination Angle vs. time", fontsize=20)
	plt.xlabel("time [days]", fontsize=20)
	plt.ylabel("declination [deg]", fontsize=20)

	plt.plot(times_d, np.degrees(decls))
	plt.grid()

	plt.show()


# ======================================================== #

def SolarFluxPlot():
	lats = [np.radians(k) for k in np.arange(-90,91,5)]
	times = np.arange(0, 365*86400, 86400)

	fig, ax = plt.subplots(1, 2)

	for L in lats:
		lat_str = str(np.round(np.degrees(L),1)) + r"$^\circ$"
		HemisphereID = 1 if (L >= -np.pi/2 and L <= 0) else 0

		fluxes = [flux.Calc_DiurnalFlux(L, day) for day in times]
		ax[HemisphereID].plot(times / 86400, fluxes, "--", label=lat_str)

	ax[0].set_title("Northen Hemisphere", fontsize=20)
	ax[0].legend()
	ax[0].set_xlabel("time [days]", fontsize=20)
	ax[0].set_ylabel("flux", fontsize=20)
	ax[0].grid()

	ax[1].set_title("Southern Hemisphere", fontsize=20)
	ax[1].set_xlabel("time [days]", fontsize=20)
	ax[1].set_ylabel("flux", fontsize=20)
	ax[1].legend()
	ax[1].grid()

	plt.show()


# ======================================================== #

def del_approximations(lats, ics):
	analytic_col = (0.58, 0.85, 0.48)
	LabelSize = 20
	thickness = 5.5

	lats_deg = np.degrees(lats)
	fig, ax = plt.subplots(3, 1)
	fig.suptitle("Numerical vs. Analytic Derivatives for Initial Temperature Function", fontsize=LabelSize)

	# ------------------------------------- #
	axis = ax[0] 
	axis.plot(lats_deg, ics, label="initial conditions", linewidth=thickness, color=analytic_col)
	axis.set_ylabel(r"$T(\lambda)$", fontsize=LabelSize)
	axis.set_xlabel(r"$\lambda$", fontsize=LabelSize)
	axis.legend()
	axis.grid()
	# ------------------------------------- #
	grads = [grad(lats, ics, k) for k in range(len(lats))]

	def analytical_grad(lat):
		return 273.0 * np.sin(2*lat)

	axis = ax[1] 
	axis.plot(lats, analytical_grad(lats), label="analytical", linewidth=thickness,color=analytic_col, alpha=0.5)
	axis.plot(lats, grads, "-", label="numerical", color=(0,0,0))
	axis.set_ylabel(r"$df/d\lambda$", fontsize=LabelSize)
	axis.set_xlabel(r"$\lambda$", fontsize=LabelSize)
	axis.legend()
	axis.grid()
	# ------------------------------------- #
	laplaces = [laplace(lats, ics, k) for k in range(len(lats))]

	def analytical_laplace(lat):
		return 546.0 * np.cos(2*lat)

	axis = ax[2] 
	axis.plot(lats, analytical_laplace(lats), label="analytical", linewidth=thickness,color=analytic_col, alpha=0.5)
	axis.plot(lats, laplaces, "-", label="numerical", color=(0,0,0))
	axis.set_ylabel(r"$d^2f/d\lambda^2$", fontsize=LabelSize)
	axis.set_xlabel(r"$\lambda$", fontsize=LabelSize)
	axis.legend()
	axis.grid()
	# ------------------------------------- #

	plt.show()

# ======================================================== #

def OceanProfile():
	lats = np.linspace(-np.pi/2, np.pi/2, 1000)
	fracs = [pde.Get_OceanFraction(lat) for lat in lats]

	plt.figure()
	plt.plot(lats, fracs)
	plt.title("Fraction of latitude band that is ocean")
	plt.ylabel("Ocean Fraction")
	plt.xlabel("latitude")
	plt.grid()
	plt.show()

# ======================================================== #

def TemporalHeatmap(times, tprofs):
	print("Generating temporal heatmap ..")

	def Celsius(temps_K):
		return [(t-273.15) for t in temps_K]

	heatmap = np.transpose([tp for tp in tprofs])
	t0, t1 = times[0] / 86400, times[-1] / 86400

	plt.figure()
	plt.imshow(heatmap, 
		origin="lower",
		extent=[t0,t1,-90,90], 
		interpolation="gaussian",
		cmap=cmocean.cm.thermal)

	plt.title("Global Temperature Simulation")
	plt.ylabel("latitude")
	plt.xlabel("days")
	plt.colorbar()
	plt.show()

# ======================================================== #

def CompareModel():
	def PHM(lat):
		return 302.3 - 45.3 * np.power(np.sin(lat),2.0)

	# simulate climate.
	LAT_STEP = 9
	SIM_TIME = 365*1 # simulate for 1 year.
	SIM_Lats = [np.radians(k) for k in np.arange(-90, 90+LAT_STEP, LAT_STEP)]
	IC_Temps = [400 for k in SIM_Lats] # global temperatures @ t=0.
	SIM_Data = pde.EvolveGlobalTemperatures(SIM_Lats, IC_Temps, SIM_TIME)

	SIM_AvgT = np.array(SIM_Data[0][0])
	for DataFrame in SIM_Data:
		SIM_AvgT = np.add(SIM_AvgT, np.array(DataFrame[0]))
	SIM_AvgT = np.divide(np.array(SIM_AvgT), SIM_TIME)

	# evaluate phenomological model.
	PHM_Lats = np.linspace(-np.pi/2, np.pi/2, 100)
	PHM_AvgT = PHM(PHM_Lats)

	# plot comparison.
	plt.figure()
	plt.title("Model Comparison")
	plt.ylabel("temperature")
	plt.xlabel("latitude")

	plt.plot(PHM_Lats, PHM_AvgT, color=(0,0,0), linewidth=5.0, alpha=0.25, label="phenomological model")
	plt.plot(SIM_Lats, SIM_AvgT, "x--", color=(1,0,0), linewidth=1.0, label="simulation")

	plt.grid()
	plt.legend()
	plt.show()

# ============================================ #
