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
	
	decls_ellip = [flux.Calc_Declination(t_s) for t_s in times_s]
	decls_circ = [flux.Calc_Declination(t_s, False) for t_s in times_s]

	plt.figure()
	plt.title("Declination Angle vs. time", fontsize=20)
	plt.xlabel("time [days]", fontsize=20)
	plt.ylabel("declination [deg]", fontsize=20)
	plt.plot(times_d, np.degrees(decls_ellip), label="elliptical orbit")
	plt.plot(times_d, np.degrees(decls_circ), label="circular orbit")
	plt.legend()
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

def TemporalHeatmap(sim, subset=0, sim_step=1):
	# @note: the 'subset' parameter is how many years relative
	# to the end of the simulation, we should include within
	# the plot.
	#
	# ie: subset=0, includes all the sim data
	# ie: subset=3, includes the last 3 years of sim data
	#
	# 'sim_step' is the timestep used when computing the
	# simulation, this is provided as a return value from
	# SimulateClimate (defaults to 1 day).

	print("Generating temporal heatmap ..")

	# Compute the subset index (sidx) from the subset parameter.
	tsteps_per_yr = 365 / sim_step
	sidx = int(subset * tsteps_per_yr)
	assert(sidx <= len(sim))

	# Gather the temperature array for each timestep
	# into a big array.
	temperature_frames = []
	for dframe in sim: 
		temperature_frames.append([T_k for T_k in dframe[1]])

	# Find the min/max temperatures across the entire sim
	temp_dump = []
	for k in temperature_frames:
		for T in k: temp_dump.append(T)
	Tmin, Tmax = min(temp_dump), max(temp_dump)

	# Extract the time scale of the simulation. 
	timescale = [(sim[-sidx])[0], (sim[-1])[0]]
	timescale = np.divide(timescale, 365)

	# Construct the heatmap plot.
	plt.figure()
	
	heatmap = temperature_frames[-sidx:]
	plt.imshow(np.transpose(heatmap), 
		origin="lower",
		extent=[timescale[0],timescale[1],-90, 90],
		vmin=Tmin, vmax=Tmax,
		interpolation="gaussian",
		cmap=cmocean.cm.thermal,
		aspect="auto")

	plt.title("Global Temperature Simulation")
	plt.ylabel("latitude")
	plt.xlabel("year")
	plt.colorbar()
	plt.show()

# ======================================================== #

def CompareModel(sim, lats, timestep=1):
	def PHM(lat):
		return 302.3 - 45.3 * np.power(np.sin(lat),2.0)

	LastNYrs = 5
	tsteps_per_yr = 365 / timestep
	SIM_Temps = [df[1] for df in sim]
	TempProfiles = SIM_Temps[int(-tsteps_per_yr*LastNYrs):]
	
	SIM_Avgs = np.array(TempProfiles[0])
	for k in range(1, len(TempProfiles)):
		SIM_Avgs = np.add(SIM_Avgs, TempProfiles[k])
	SIM_Avgs = np.divide(SIM_Avgs, len(TempProfiles))

	# evaluate phenomological model.
	PHM_Lats = np.linspace(-np.pi/2, np.pi/2, 100)
	PHM_AvgT = PHM(PHM_Lats)

	# plot comparison.
	plt.figure()
	plt.title("Model Comparison")
	plt.ylabel("temperature")
	plt.xlabel("latitude")

	plt.plot(lats, SIM_Avgs, "x--", color=(1,0,0), linewidth=1.0, label="simulation")
	plt.plot(PHM_Lats, PHM_AvgT, color=(0,0,0), linewidth=5.0, alpha=0.25, label="phenomological model")

	plt.grid()
	plt.legend()
	plt.show()

# ============================================ #

def OrbitalApproximation():

	plt.figure()
	plt.title("Elliptic Orbit Position")
	plt.xlabel("days", fontsize=20)
	plt.ylabel("orbital position", fontsize=20)
	
	for e in np.linspace(0, 0.99, 20):
		e_str = str(np.round(e, 2))
		times, angular_positions = flux.EllipticOrbit(365*86400, 365, e)
		plt.plot(times, angular_positions, "--", alpha=0.8, linewidth=0.8, color=(0.66, 0.34, 0.82))

	Earth_e = 0.01671
	times, angular_positions = flux.EllipticOrbit(365*86400, 365, Earth_e)
	plt.plot(times, angular_positions, label="(Earth) e="+str(Earth_e), linewidth=2, color=(1,0,0))
	
	plt.legend()
	plt.grid()
	plt.show()


# ============================================ #

def SolarFluxAlongOrbit():
	plt.figure()
	plt.title("Solar flux throughout various orbits")
	plt.xlabel("days", fontsize=20)
	plt.ylabel(r"$q(t)$", fontsize=20)

	# Calculate the actual flux Earth recived from the sun
	# throughout its orbit, first a circular orbit,
	# and then its true elliptical orbit.
	times = np.linspace(0, 365*86400, 1000)
	Flux_ActualOrbit = [flux.Calc_SolarRadiation(t) for t in times]

	plt.plot(np.divide(times, 86400), Flux_ActualOrbit, label="e=0.01671")

	plt.legend()
	plt.grid()
	plt.show()

# ============================================ #

def TrueAnomalyODE_Convergence():
	plt.figure()
	plt.title("Convergence of numerical integration for Earth's True Anomaly")
	plt.xlabel("days", fontsize=20)
	plt.ylabel(r"$f$", fontsize=20)

	for ss in np.arange(500, 3000, 100):
		times, pos = flux.Calc_TrueAnomaly_ODE(ss)
		plt.plot(times, pos, label=str(ss)+" steps")

	plt.legend()
	plt.grid()
	plt.show()

# ============================================ #

def TrueAnomalyResiduals():
	plt.figure()
	plt.title("True Anomaly: Numerical Integration vs. Series Approximation")
	plt.xlabel("days", fontsize=20)
	plt.ylabel("residuals rel. numeric result", fontsize=20)

	# Use 'TrueAnomalyODE_Convergence' plot to find a good step size
	# such that we get an accurate calculation here.
	ref_ts, ref_fs = flux.Calc_TrueAnomaly_ODE(dt=100)

	for tc in range(2,5):
		approx_fs = [flux.Calc_TrueAnomaly_Approximate(t, term_count=tc) for t in ref_ts]
		residuals = [(approx_fs[k] - ref_fs[k]) for k in range(len(ref_ts))]
		plt.plot(np.divide(ref_ts,86400), residuals, label=str(tc) + " terms")

	plt.grid()
	plt.legend()
	plt.show()


# ============================================ #