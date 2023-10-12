import matplotlib.pyplot as plt
import numpy as np

# ============================================ #
def grad(lats, temps, idx):
	 # latitude we want the derivative evaluated at.
	lat = lats[idx]

	# the derivative is zero at the poles.
	if(np.isclose(np.absolute(lat), np.pi/2)):
		return 0.0

	# for any other latitude band, use central difference approximation.
	lat_above, lat_below = lats[idx + 1], lats[idx - 1]
	t_above, t_below = temps[idx + 1], temps[idx - 1]
	
	# we assume a uniform spacing in latitude coords!
	ss_above = lat_above - lat
	ss_below = lat - lat_below
	assert(np.isclose(ss_above, ss_below)) 
	step = ss_above

	# central difference approximation
	return (t_above - t_below) / (2 * step)
# ============================================ #

# ============================================ #
def laplace(lats, temps, idx):
	# latitude we want the derivative evaluated at.
	lat = lats[idx]

	if(np.isclose(lat, np.pi/2)): # north pole
		return 0.0
	elif(np.isclose(lat, -np.pi/2)): # south pole
		return 0.0

	return 0.0

# ============================================ #


# ======================================================== #
# Plots showing grad approximation
# ======================================================== #
if(0):
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

	colours = [(1,0,0),(0,1,0), (0,0,1)]
	labels = ["const.", r"$x^2$",r"$\sin(x)$"]

	functions = [f0,f1,f2]
	derivatives = [f0_dash,f1_dash,f2_dash]
	derivatives2 = [f0_dash2,f1_dash2,f2_dash2]

	fig, axes = plt.subplots(3,2)
	fig.suptitle("Plots of function derivatives (analytic vs. numerical)")

	PlotRow = 0
	xs = np.linspace(-np.pi/2,np.pi/2, 60)
	for j in range(len(functions)):
		f, dfdx, dfdx2 = functions[j], derivatives[j], derivatives2[j]
		f_of_x = [f(x) for x in xs]

		# first order derivative
		numerical_1 = [grad(xs, f_of_x, k) for k in range(len(xs))]
		analytical_1 = [dfdx(x) for x in xs]

		# second order derivative
		numerical_2 = [laplace(xs, f_of_x, k) for k in range(len(xs))]
		analytical_2 = [dfdx2(x) for x in xs]

		axis = axes[PlotRow, 0]
		axis.grid()
		axis.set_title("f(x) = " + labels[j])
		axis.set_xlabel("x")
		axis.set_ylabel("df/dx(x)")
		axis.plot(xs, analytical_1,label="Analytic", color=colours[j], alpha=0.35)
		axis.plot(xs, numerical_1, "x", label="Numerical", color=colours[j])
		axis.legend()

		axis = axes[PlotRow, 1]
		axis.grid()
		axis.set_title("f(x) = " + labels[j])
		axis.set_xlabel("x")
		axis.set_ylabel(r"$d^2f/dx^2(x)$")
		axis.plot(xs, analytical_2,label="Analytic", color=colours[j], alpha=0.35)
		axis.plot(xs, numerical_2, "x", label="Numerical", color=colours[j])
		axis.legend()
		
		PlotRow += 1
		
	plt.show()
# ======================================================== #

# ============================================ #
if(0):
	# Initial temperature conditions for all latitude bands.
	lats = np.arange(-90, 91, 10)
	data = [[273.0 for k in lats]]

	# timestep (1 day)
	duration = 3e7 # evolve temps for a year.
	dt = 86400.0
	time = 0

	while (time < duration):
		current_temp_profile = data[-1]
		data.append([])

		# evolve temperatures
		for idx, T in enumerate(current_temp_profile):
			evolved_T = 0.0 # T + dt * solve_pde(current_temp_profile, lats, idx, time)
			data[-1].append( evolved_T )

		# advance time
		time += dt

	# plot temperature profile of Earth across time.
# ============================================ #
