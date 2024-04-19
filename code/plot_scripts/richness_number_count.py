from toolkit import toolkit
import numpy as np
import matplotlib.pyplot as plt


def excess(log_r):
    return 1 + (log_r - (1/50) * log_r ** 2) / 100


def overdensity(excess):
    return excess

def log_richness(log_snr):
    alpha = 1
    beta = 2
    return alpha + beta * log_snr


def log_richness_dist(log_snr, log_r):
    return (np.exp(-(1/2) * np.square((log_r - log_richness(log_snr)) / scatter)) /
            (scatter * np.sqrt(2 * np.pi)))


scatter = 1
steps = 1000

x = np.linspace(4, 10, 1000)
y = np.zeros(np.shape(x))
bins = np.linspace(-100, 100, steps + 1)

for i in range(len(x)):
    y[i] = 0
    for j in range(steps):
        y[i] += ((bins[j+1] - bins[j]) * log_richness_dist(x[i], (bins[j+1] + bins[j]) / 2) *
                 overdensity(excess((bins[j+1] + bins[j]) / 2)))

plt.plot(x, y)
plt.plot(x, excess(x))
plt.show()
