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

folders = ["./rotations/act/", "./rotations/planck/"][::-1][0:1]
#folders = ["./rotations/planck/"]
labels = ["Point", "Galactic"]

NSIDES = np.array((1, 2, 4, 8, 16, 32, 64, 128))
x = np.append(np.array(1/2), NSIDES)

y = np.zeros((2, len(x)))
y_std = np.zeros((2, len(x)))

use_abs = False

for folder in folders:
    files = os.listdir(folder)
    sigma = np.zeros((2, len(files), len(x)))
    for file in files:
        data = np.load(folder + file)
        if use_abs:
            sigma[0][files.index(file)] = np.abs(data[0][0] / data[0][1])
            sigma[1][files.index(file)] = np.abs(data[1][0] / data[1][1])
        else:
            sigma[0][files.index(file)] = data[0][0] / data[0][1]
            sigma[1][files.index(file)] = data[1][0] / data[1][1]
    for i in range(len(y[0])):
        y[0][i] = np.mean(sigma[0][:, i])
        y_std[0][i] = np.std(sigma[0][:, i])
        y[1][i] = np.mean(sigma[1][:, i])
        y_std[1][i] = np.std(sigma[1][:, i])
print(y)
print(y_std)

plt.errorbar(x, y[0], y_std[0], label=labels[0], ecolor=error_color[0], color=color[0], capsize=3,
             capthick=1,)
#plt.errorbar(x, y[1], y_std[1], label=labels[1], ecolor=error_color[1], color=color[1], capsize=3,
#             capthick=1,)

if use_abs:
    plt.plot((1/2, np.max(NSIDES)), (0.8, 0.8), color="k", linestyle="dashed")
    plt.plot((1/2, np.max(NSIDES)), (0.2, 0.2), color="k", linestyle="dotted")
    plt.plot((1/2, np.max(NSIDES)), (1.4, 1.4), color="k", linestyle="dotted")
else:
    plt.plot((1 / 2, np.max(NSIDES)), (1, 1), color="k", linestyle="dashed")
    plt.plot((1 / 2, np.max(NSIDES)), (0, 0), color="k", linestyle="dotted")
    plt.plot((1 / 2, np.max(NSIDES)), (-1, -1), color="k", linestyle="dotted")

plt.title("Normalised mean detection for the rotated Planck masks")
plt.xscale("log")
#plt.legend()
plt.xticks(x, ["C"] + list(NSIDES))
plt.minorticks_off()
#ax.set_xticks(x, ["C"] + list(NSIDES))
plt.xlabel("NSIDE")
plt.ylabel("Mean (Excess / Error)")
plt.savefig("planck_rot.pdf")
plt.show()
