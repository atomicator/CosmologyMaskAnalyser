import matplotlib.pyplot as plt
import numpy as np
from toolkit import toolkit

toolkit.plt_use_tex()

data = np.load("bias_data/normal_bias_act_400k.npy")

fig = plt.figure()
ax = fig.add_subplot(111)

print(data)

print(np.shape(data))

error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:mint green", "pink"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple"]

labels = ["Error on mean", "Deviation"]

mean = np.zeros(len(data[0]))
variance = np.zeros(len(data[0]))
weights = np.ones(len(data[0]))

print(data[0])

for data_set in data:
    mean += data_set[:, 0] * (1 / len(data))
    variance += (data_set[:, 1] ** 2) * (1 / len(data) ** 2)

mean = mean / weights
variance = variance / weights

print(mean, variance)

NSIDES = [1, 2, 4, 8, 16, 32, 64, 128]
x = [1/2, *NSIDES]

mean_error = np.zeros(np.shape(data)[1])
for i in range(np.shape(data)[1]):
    mean_error[i] = np.mean(data[:, i, 1])

ax.errorbar(x, mean, np.sqrt(variance), marker="+", ecolor=error_bar_colors[0],
            ls="none", color=line_colors[0], capsize=3, capthick=1, label=labels[0])
ax.errorbar(x, mean, mean_error, marker="+", ecolor=error_bar_colors[1],
            ls="none", color=line_colors[0], capsize=3, capthick=1, label=labels[1])
ax.plot(x, mean, color=line_colors[0])

ax.set_xscale("log", base=2)
ax.set_xlim(1/2 * np.sqrt(1/2), NSIDES[-1] * np.sqrt(2))

ax.plot([1/2, NSIDES[-1]], np.ones(2) * .1, color="k")
ax.set_xticks([0.5] + NSIDES, ["C"] + NSIDES)

#plt.ylim(-0.005, 0.005)

ax.legend()
ax.set_xlabel("NSIDE")
ax.set_ylabel("Excess")
ax.set_title(f"Fully random ACT data, 100 realisations")
plt.savefig("test.png")
plt.show()

# To do:
# Regenerate random catalogue
# Gen random catalogue with negative excess
# Add other catalogues to plot
# Create scatter from ACT data and pipeline
