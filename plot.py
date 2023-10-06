from datapack import Datapack
import matplotlib.pyplot as plt
import os
import re

### ==== Plotting the latitude datapacks === ###
include_list = [i for i in range(0, 19)]

Label_Size = 16

for i in include_list:
	fname = "latitude_pack_" + str(i)
	data = Datapack("datapacks/" + fname + ".dp")

	plt.figure()
	plt.title(data.desc, fontsize=Label_Size)
	plt.xlabel("t (s)", fontsize=Label_Size)
	plt.ylabel("daily averaged flux, S", fontsize=Label_Size)
	plt.plot(data.Xdata, data.Ydata, ".")
	plt.savefig("datapacks/" + fname + ".png", bbox_inches="tight")

