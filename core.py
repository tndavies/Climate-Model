import numpy as np
import plots
import pde

# ============================================ #
def foo(lat):
		return (np.sin(lat)**2)*273.0

lats = np.linspace(-np.pi/2,np.pi/2,20)
ic = [foo(k) for k in lats]

data = pde.EvolveGlobalTemperatures(lats, ic, 0.5)
plots.GlobalTemperaturePlots(lats, data)