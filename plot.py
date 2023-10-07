from datapack import Datapack
import matplotlib.pyplot as plt
import numpy as np

### ==== Plotting the latitude datapacks === ###
Datapack_Count = 19
include_list = [i for i in range(1, Datapack_Count+1)]

heatmap = [[False] * n for y in range(0, n)]
n = int(np.ceil(np.sqrt(Datapack_Count)))
fig, axes = plt.subplots(n,n)

fig.suptitle("Plots of S(t) for various latitudes, in " + r"$10^\circ$" + "steps",
	fontsize=30)

axes[0, 0].text(0.5,0.65, r'$lat=-90^\circ$', fontsize=18, horizontalalignment='center',
     verticalalignment='center', transform=axes[0, 0].transAxes)

plot_x, plot_y = 0, 0
for pack_id in include_list:
	fname = "latitude_pack_" + str(pack_id)
	data = Datapack("datapacks/" + fname + ".dp")

	axis = axes[plot_y, plot_x]
	axis.plot(data.xs, data.ys, ".")

	heatmap[plot_y][plot_x] = True

	plot_x += 1
	if(plot_x > n - 1):
		plot_x = 0
		plot_y += 1

# remove unused axes from figure
for y in range(0, n):
	for x in range(0, n):
		if(not heatmap[y][x]):
			plt.delaxes(axes[y][x])

plt.show()

