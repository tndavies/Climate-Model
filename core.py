import numpy as np
import plots
import pde
import matplotlib.pyplot as plt
import flux

#LAT_STEP = 6
#lats = [np.radians(k) for k in np.arange(-90, 90+LAT_STEP, LAT_STEP)]
#ics = [400 for k in lats]
#times, tprofs = pde.EvolveGlobalTemperatures(lats, ics, 365*25)
#plots.TemporalHeatmap(times, tprofs)


plots.CompareModel()
