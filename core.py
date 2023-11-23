import numpy as np
import plots
import pde
import matplotlib.pyplot as plt
import flux

sim, lats, ts = pde.SimulateClimate(30)
#plots.CompareModel(sim, lats, timestep=ts)
#plots.TemporalHeatmap(sim, subset=0, sim_step=ts)

times = [df[0] for df in sim]
sp_temps = [df[1][0] for df in sim]

plt.figure()
plt.xlabel("years", fontsize=25)
plt.xlabel("temperature", fontsize=25)
plt.title("South Pole Temperature")
plt.plot(np.divide(times, 365), sp_temps)
plt.grid()
plt.show()