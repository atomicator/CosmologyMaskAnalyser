import random
import numpy as np
import matplotlib.pyplot as plt


def bootstrapping(data, depth, iterations):
    mean_estimates = [np.mean(data)]
    n = len(data)
    new_data = data.copy()
    for i in range(depth):
        temp = random.choices(data, k=n)
        new_data = temp.copy()
        mean_estimates.append(np.mean(temp))
    plt.hist(mean_estimates)
    plt.show()
    return np.mean(mean_estimates), np.std(mean_estimates)

