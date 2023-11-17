import numpy as np
import plots
import pde
import matplotlib.pyplot as plt
import flux

sim = pde.SimulateClimate(2)
plots.TemporalHeatmap(sim)