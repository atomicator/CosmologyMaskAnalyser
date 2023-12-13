import numpy as np
import matplotlib.pyplot as plt
from toolkit import toolkit


toolkit.plt_use_tex()

data = np.loadtxt("../../sim_results/planck_point.csv", delimiter=",", dtype=np.complex_).transpose()
# TODO: Last line is broken - remove slice when fixed
min_row, max_row = (2, -1)
data = np.real(data)
new_data = data.copy()[min_row: max_row]

print(new_data)

for i in range(len(new_data)):
    new_data[i] -= np.mean(data[i + min_row])


max_value = np.max(new_data)
min_value = np.min(new_data)
value_range = (max_value - min_value) * 1.1

num_bins = 20

x = np.zeros((len(new_data), num_bins))
y = np.linspace(min_value, max_value, num_bins)

for i in range(len(new_data)):
    for data_point in new_data[i]:
        x[i][int(num_bins * (data_point - min_value) / value_range)] += 1

dp = 1 - int(np.floor(np.log10(np.std(data[0]))))

labels = [rf"Planck {np.mean(data[0]):.{dp}f} $\pm$ {np.std(data[0]):.{dp}f}",
          rf"SDSS {np.mean(data[1]):.{dp}f} $\pm$ {np.std(data[1]):.{dp}f}",
          rf"Planck only {np.mean(data[2]):.{dp}f} $\pm$ {np.std(data[2]):.{dp}f}",
          rf"SDSS only {np.mean(data[3]):.{dp}f} $\pm$ {np.std(data[3]):.{dp}f}",
          rf"Both {np.mean(data[4]):.{dp}f} $\pm$ {np.std(data[4]):.{dp}f}",
          rf"Neither {np.mean(data[5]):.{dp}f} $\pm$ {np.std(data[5]):.{dp}f}"][min_row:max_row]

for i in range(len(x)):
    plt.plot(y, x[i], linestyle="dotted", label=labels[i])

plt.legend()
plt.ylabel("Estimates")
plt.xlabel("Deviation from mean")
plt.title("The spread of results from the bootstrap simulations")
plt.savefig("../../graphs/bootstrap_results.png")
plt.show()

#TODO: generate random star sample in SDSS mask, then re run with actual data
