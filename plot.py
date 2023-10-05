from datapack import Datapack
import matplotlib.pyplot as plt

data = Datapack("datapacks/declination.dp")

plt.figure()
plt.title(data.desc)
plt.xlabel("time (s)")
plt.ylabel("declination")
plt.plot(data.Xdata, data.Ydata)
plt.show()
