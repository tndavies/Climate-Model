import numpy as np
import plots
import pde

# move to elliptical orbit
# put constants and flux into pde stuff
# compare sim temps to actual data.

# ============================================ #
lats = np.linspace(-np.pi/2,np.pi/2, 20)
ic = [200.0 for k in lats]

data = pde.EvolveGlobalTemperatures(lats, ic, 86400*3)
plots.GlobalTemperaturePlots(lats, data)

if(0):
	for ds in data:
		temps, time = ds[0], ds[1]

		print("t=" + str(time/86400) + " days")
		for T in temps:
			print(str(T) + " [k]")

		print("")
