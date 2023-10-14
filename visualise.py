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

def SpatialDerivativePlots():
	# some functions f(x)
	def f0(x): return np.pi
	def f1(x): return np.power(x, 2.0)
	def f2(x): return np.sin(x)

	# their analytic derivatives f'(x)
	def f0_dash(x): return 0.0
	def f1_dash(x): return 2.0*x
	def f2_dash(x): return np.cos(x)

	# their analytic derivatives f''(x)
	def f0_dash2(x): return 0.0
	def f1_dash2(x): return 2.0
	def f2_dash2(x): return -np.sin(x)

	LabelSize = 18
	analytic_col = (0.58, 0.85, 0.48)
	labels = ["const.", r"$x^2$",r"$\sin(x)$"]

	functions = [f0,f1,f2]
	derivatives = [f0_dash,f1_dash,f2_dash]
	derivatives2 = [f0_dash2,f1_dash2,f2_dash2]

	fig, axes = plt.subplots(2, 3)
	fig.suptitle("Plots of function derivatives (analytic vs. numerical)", fontsize=25)

	PlotCol = 0
	xs = np.linspace(-np.pi/2,np.pi/2, 30)
	for j in range(len(functions)):
		f, dfdx, dfdx2 = functions[j], derivatives[j], derivatives2[j]
		f_of_x = [f(x) for x in xs]

		# first order derivative
		numerical_1 = [grad(xs, f_of_x, k) for k in range(len(xs))]
		analytical_1 = [dfdx(x) for x in xs]

		# second order derivative
		numerical_2 = [laplace(xs, f_of_x, k) for k in range(len(xs))]
		analytical_2 = [dfdx2(x) for x in xs]

		axis = axes[0, PlotCol]
		axis.grid()
		axis.set_title(r"$T(\lambda) = $" + labels[j], fontsize=LabelSize)
		axis.set_xlabel(r"$\lambda$", fontsize=LabelSize)
		axis.set_ylabel(r"$df/d\lambda$", fontsize=LabelSize)
		axis.plot(xs, analytical_1,label="Analytical", linewidth=4.5,color=analytic_col)
		axis.plot(xs, numerical_1, "x", label="Numerical", color=(0,0,0))
		axis.legend()

		axis = axes[1, PlotCol]
		axis.grid()
		axis.set_xlabel(r"$\lambda$", fontsize=LabelSize)
		axis.set_ylabel(r"$d^2f/d\lambda^2$", fontsize=LabelSize)
		axis.plot(xs, analytical_2,label="Analytical", linewidth=4.5,color=analytic_col)
		axis.plot(xs, numerical_2, "x", label="Numerical", color=(0,0,0))
		axis.legend()
		
		PlotCol += 1
		
	plt.show()

# ======================================================== #

def GlobalTemperaturePlots(lats, data):
	print("Generating global temperature visualisations ..")
	b0,b1 = np.degrees(np.min(lats)), np.degrees(np.max(lats))
	
	with alive_bar(len(data)) as bar:
		plt.figure()
		
		for j, ds in enumerate(data):
			time = (ds[1] / 86400) # convert from secs to days
			temps = ds[0]

			arr = np.transpose(np.array([temps for n in range(lats.size)]))
			frame_name = "pde_video/frame_" + str(j) + ".png"

			plt.clf()
			plt.title("Global temperature @ t={day} days".format(day=time))
			plt.xlabel("longitude [deg]")
			plt.ylabel("latitude [deg]")
			
			plt.imshow(arr, origin="lower", 
					interpolation="gaussian",
					extent=[b0,b1,b0,b1])
			plt.colorbar()

			plt.savefig(frame_name, format="png", bbox_inches="tight", pad_inches=0.25)

			bar()
		
		plt.close()
# ============================================ #
