import numpy as np
import plots
import pde
import matplotlib.pyplot as plt
import flux

# 1) What is the albedo function actually accounting for?
# 2) Add in ice fraction for land mass.
# 3) Fix C & A for Antartica.

sim, lats, ts = pde.SimulateClimate(60)
plots.CompareModel(sim, lats, timestep=ts)
# plots.TemporalHeatmap(sim, subset=10, sim_step=ts)