import numpy as np
import plots
import pde
import matplotlib.pyplot as plt
import flux

sim, ts = pde.SimulateClimate(60)
plots.TemporalHeatmap(sim, sim_step=ts, subset=10)