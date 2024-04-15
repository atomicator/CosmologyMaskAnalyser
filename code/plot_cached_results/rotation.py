import os
from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np

folders = ["./rotations/act/", "./rotations/planck/"]

x = np.array((1/2, 1, 2, 4, 8, 16, 32, 64, 128))

y = np.zeros((len(folders), len(x)))
y_std = np.zeros((len(folders), len(x)))

for folder in folders:
    files = os.listdir(folder)
    sigma = np.zeros((len(files), len(x)))
    for file in files:
        data = np.load(folder + file)[0]
        sigma[files.index(file)] = np.abs(data[0] / data[1])
    for i in range(len(y[0])):
        y[folders.index(folder)][i] = np.mean(sigma[:, i])
        y_std[folders.index(folder)][i] = np.std(sigma[:, i])
print(y, y_std)
