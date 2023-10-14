import numpy as np
import visualise
import pde

# ============================================ #

lats = np.linspace(-np.pi/2,np.pi/2,60)
ic = [np.pi/2-abs(k) for k in lats]

#data = pde.EvolveGlobalTemperatures(lats, ic, 6)
#visualise.GenerateGlobalTemperaturePlots(lats, data)

visualise.SolarFluxPlot()