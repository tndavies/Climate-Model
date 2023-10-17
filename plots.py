from alive_progress import alive_bar
import matplotlib.pyplot as plt
from pde import grad, laplace
import numpy as np
import flux

# ======================================================== #
# 'dur' is in units of days.

def DeclinationPlot(dur):
	times_d = np.linspace(0, dur, 2*dur)
	times_s = np.array([t_d * 86400 for t_d in times_d])
	decls = flux.Calc_Declination(times_s)

	plt.figure()
	plt.title("Declination Angle vs. time")
	plt.xlabel("time [days]")
	plt.ylabel("declination [deg]")

	plt.plot(times_d, np.degrees(decls))
	plt.grid()

	plt.show()


# ======================================================== #

def SolarFluxPlot():
	lats = [np.radians(k) for k in np.arange(-90,91,15)]
	times = np.arange(0, 365)

	fig, ax = plt.subplots(1, 2)

	for L in lats:
		lat_str = str(np.round(np.degrees(L),1)) + r"$^\circ$"
		HemisphereID = 1 if (L >= -np.pi/2 and L <= 0) else 0

		fluxes = [flux.Calc_DiurnalFlux(L, day) for day in times]
		
		ax[HemisphereID].plot(times, fluxes, label=lat_str)

	ax[0].set_title("Northen Hemisphere")
	ax[0].legend()
	ax[0].set_xlabel("time [days]")
	ax[0].set_ylabel("flux")
	ax[0].grid()

	ax[1].set_title("Southern Hemisphere")
	ax[1].set_xlabel("time [days]")
	ax[1].set_ylabel("flux")
	ax[1].legend()
	ax[1].grid()

	plt.show()


# ======================================================== #

def del_approximations(lats, ics):
	analytic_col = (0.58, 0.85, 0.48)
	LabelSize = 14
	thickness = 5.5

	lats_deg = np.degrees(lats)
	fig, ax = plt.subplots(3, 1)

	# ------------------------------------- #
	axis = ax[0] 
	axis.plot(lats_deg, ics, linewidth=thickness, color=analytic_col)
	axis.set_ylabel(r"$T(\lambda)$", fontsize=LabelSize)
	axis.set_xlabel(r"$\lambda$", fontsize=LabelSize)
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
	axis.grid()
	# ------------------------------------- #

	plt.show()

# ======================================================== #

def GlobalTemperaturePlots(lats, data):
	print("Generating global temperature visualisations ..")
	
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
					interpolation="nearest",
					extent=[b0,b1,b0,b1],
					vmin=initial_min_temp,
					vmax=initial_max_temp)

			plt.colorbar()

			frame_name = "pde_video/" + str(j+1).zfill(4) + ".png"
			plt.savefig(frame_name, format="png", bbox_inches="tight", pad_inches=0.25)

			bar()
		
		plt.close()
# ============================================ #
