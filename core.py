import numpy as np
import plots
import pde
import matplotlib.pyplot as plt
import flux
import data

# yrs, temps = data.ImportTemps("Climate Temperature.csv")
yrs, co2 = data.ImportCO2("Climate CO2.csv")

plt.figure()
plt.xlabel("Year")
plt.ylabel("CO2 [ppm]")
plt.plot(yrs,co2, "o--")
plt.grid()
plt.show()

# sim, lats, ts = pde.SimulateClimate(15)
# plots.TemporalHeatmap(sim, sim_step=ts)
# plots.CompareModel(sim, lats, timestep=ts)

# temps = np.linspace(200,350, 1000)
# irs = [pde.Calculate_IRCooling(T) for T in temps]
# irs2 = [pde.foo_Calculate_IRCooling(T) for T in temps]
# plt.figure()
# plt.plot(temps, irs, label="greenhouse")
# plt.plot(temps, irs2, label="simple")
# plt.grid()
# plt.legend()
# plt.show()