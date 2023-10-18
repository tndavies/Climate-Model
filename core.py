import numpy as np
import plots
import pde

# fix flux stuff
# move to elliptical orbit
# put constants and flux into pde stuff
# compare sim temps to actual data.

# ============================================ #
lats = np.linspace(-np.pi/2,np.pi/2,20)
ic = [((np.sin(k)**2)*273.0) for k in lats]

#data = pde.EvolveGlobalTemperatures(lats, ic, 0.5)
#plots.GlobalTemperaturePlots(lats, data)

plots.SolarFluxPlot()
