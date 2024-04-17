import os
from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np


def gauss_monte_carlo():
    x = []
    for i in range(1000000):
        y = np.abs(np.random.normal(size=35))
        x.append((np.array(np.mean(y)), np.array(np.std(y))))
    x = np.array(x)
    print(x)
    print((np.mean(x[:, 0]), np.std(x[:, 0])))
    print((np.mean(x[:, 1]), np.std(x[:, 1])))


#gauss_monte_carlo()

color = ["xkcd:electric blue", "red"]
error_color = ["xkcd:aqua blue", "orange"]
toolkit.plt_use_tex()
fig = plt.figure()
ax = fig.add_subplot()

folders = ["./rotations/act/", "./rotations/planck/"][::1][0:1]
#folders = ["./rotations/planck/"]
labels = ["Point", "Galactic"]

NSIDES = np.array((1, 2, 4, 8, 16, 32, 64, 128))
x = np.append(np.array(1/2), NSIDES)

y = np.zeros((len(folders), len(x)))
y_std = np.zeros((len(folders), len(x)))

use_abs = False

for folder in folders:
    files = os.listdir(folder)
    sigma = np.zeros((len(files), len(x)))
    for file in files:
        data = np.load(folder + file)[0]
        if use_abs:
            sigma[files.index(file)] = np.abs(data[0] / data[1])
        else:
            sigma[files.index(file)] = data[0] / data[1]
    for i in range(len(y[0])):
        y[folders.index(folder)][i] = np.mean(sigma[:, i])
        y_std[folders.index(folder)][i] = np.std(sigma[:, i])
print(y)
print(y_std)

for i in range(len(folders)):
    plt.errorbar(x, y[i], y_std[i], label=labels[i], ecolor=error_color[i], color=color[i], capsize=3,
                 capthick=1,)

if use_abs:
    plt.plot((1/2, np.max(NSIDES)), (0.8, 0.8), color="k", linestyle="dashed")
    plt.plot((1/2, np.max(NSIDES)), (0.2, 0.2), color="k", linestyle="dotted")
    plt.plot((1/2, np.max(NSIDES)), (1.4, 1.4), color="k", linestyle="dotted")
else:
    plt.plot((1 / 2, np.max(NSIDES)), (1, 1), color="k", linestyle="dashed")
    plt.plot((1 / 2, np.max(NSIDES)), (0, 0), color="k", linestyle="dotted")
    plt.plot((1 / 2, np.max(NSIDES)), (-1, -1), color="k", linestyle="dotted")

plt.title("Normalised mean detection for the rotated ACT mask")
plt.xscale("log")
plt.legend()
ax.set_xticks(x, ["C"] + list(NSIDES))
plt.savefig("act_rot.pdf")
plt.show()
