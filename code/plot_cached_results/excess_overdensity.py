import numpy as np
import matplotlib.pyplot as plt
from toolkit import toolkit

toolkit.plt_use_tex()

error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:mint green", "pink"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple"]

overdensities = np.linspace(0.0, 0.2, 21)

data = np.load("./bias_data/bias_test_planck_const_range_planck.npy")

results = []

for i in range(len(overdensities)):
    mean = np.mean(data[i, :, 0, 0])
    # variation = np.mean(data[i, :, 0, 1] ** 2) / len(data[i])
    variation = np.var(data[i, :, 0, 0]) / len(data[i])
    results.append((mean, np.sqrt(variation)))

results = np.array(results).transpose()

plt.plot(overdensities, overdensities, color="black")
plt.errorbar(x=overdensities, y=results[0], yerr=results[1], label="Planck", capsize=3, capthick=1,
             color=line_colors[0], ecolor=error_bar_colors[0])
#plt.plot(overdensities, 0.01 * np.abs(results[0] - overdensities)/results[1])

data = np.load("./bias_data/bias_test_planck_const_range_act.npy")

results = []

for i in range(len(overdensities)):
    mean = np.mean(data[i, :, 0, 0])
    #variation = np.mean(data[i, :, 0, 1] ** 2) / len(data[i])
    variation = np.var(data[i, :, 0, 0]) / len(data[i])
    results.append((mean, np.sqrt(variation)))

results = np.array(results).transpose()

plt.plot(overdensities, overdensities, color="black")
plt.errorbar(x=overdensities, y=results[0], yerr=results[1], label="ACT", capsize=3, capthick=1,
             color=line_colors[1], ecolor=error_bar_colors[1])

plt.legend()

#plt.plot(overdensities, 0.01 * np.abs(results[0] - overdensities)/results[1])
#plt.plot(overdensities, np.ones(np.shape(overdensities)) * 0.045)
#plt.xlim(0, 0.05)
#plt.ylim(0, 0.05)

plt.title("Measuring the excess (const only) as a function of overdensity for 1000 \n realisations of data (per overdensity)")

plt.show()
print(np.shape(data))
print(results)
