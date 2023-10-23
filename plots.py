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

def GlobalTemperaturePlots(lats, data):
	print("Generating temperature visualisations ..")
	
	b0,b1 = np.degrees(np.min(lats)), np.degrees(np.max(lats))
	initial_min_temp = min(data[0][0])
	initial_max_temp = max(data[0][0])

	with alive_bar(len(data)) as bar:
		plt.figure()
		
		for j, ds in enumerate(data):
			temps, time = ds[0], ds[1]
			time_str = str(np.round(time, 2))

			plt.clf()
			plt.title("Global temperature @ t=" + time_str + " [s]")
			plt.xlabel("longitude [deg]")
			plt.ylabel("latitude [deg]")
			
			arr = np.transpose(np.array([temps for n in range(lats.size)]))
			plt.imshow(arr, origin="lower", 
					interpolation="gaussian",
					extent=[b0,b1,b0,b1],
					vmin=initial_min_temp,
					vmax=initial_max_temp,
					cmap=cmocean.cm.thermal)

			plt.colorbar()

			frame_name = "pde_video/" + str(j+1).zfill(4) + ".png"
			plt.savefig(frame_name, format="png", bbox_inches="tight", pad_inches=0.25)

			bar()
		
		plt.close()
# ============================================ #
def CompareModel():
	def PHM(lat):
		return 302.3 - 45.3 * np.power(np.sin(lat),2.0)

	# evaluate phenomological model for annually-averaged latitute temps.
	PHM_Lats = np.linspace(-np.pi/2, np.pi/2, 100)
	PHM_AvgT = PHM(PHM_Lats)

	# simulate climate to obtain annually-averaged latitude temps.
	SIM_TIME = 365 # simulate for 1 year.
	STEP = 9
	SIM_Lats = [np.radians(k) for k in np.arange(-90, 90+STEP, STEP)]
	IC_Temps = [400 for k in SIM_Lats] # global temperatures @ t=0.
	SIM_Data = pde.EvolveGlobalTemperatures(SIM_Lats, IC_Temps, SIM_TIME)

	SIM_AvgT = np.array(SIM_Data[0][0])
	for DataFrame in SIM_Data:
		SIM_AvgT = np.add(SIM_AvgT, np.array(DataFrame[0]))
	SIM_AvgT = np.divide(np.array(SIM_AvgT), SIM_TIME)

	# Plot comparison
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
