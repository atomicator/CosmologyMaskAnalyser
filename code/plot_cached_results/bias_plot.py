import matplotlib.pyplot as plt
import numpy as np
from toolkit import toolkit

toolkit.plt_use_tex()

step = 5
start = 4
stop = 5

files = ["bias_data_01/normal_bias_act_10k.npy", "bias_data/normal_bias_act_10k.npy", "bias_data/random_act_10k.npy",
         "bias_data/invert_bias_act_10k.npy", "bias_data_01/invert_bias_act_10k.npy"][start:stop:step]
#files = ["bias_data_01/normal_bias_act_400k.npy", "bias_data/normal_bias_act_400k.npy", "bias_data/random_act_400k.npy",
#         "bias_data/invert_bias_act_400k.npy", "bias_data_01/invert_bias_act_400k.npy"][start:stop:step]
#files = ["bias_data_01/normal_bias_planck_10k.npy", "bias_data/normal_bias_planck_10k.npy", "bias_data/random_planck_10k.npy",
#         "bias_data/invert_bias_planck_10k.npy", "bias_data_01/invert_bias_planck_10k.npy"][start:stop:step]
#files = ["bias_data_01/normal_bias_planck_400k.npy", "bias_data/normal_bias_planck_400k.npy", "bias_data/random_planck_400k.npy",
#         "bias_data/invert_bias_planck_400k.npy", "bias_data_01/invert_bias_planck_400k.npy"][start:stop:step]
alpha = [0.1, 0.05, 0, 1/19, 1/9][start:stop:step]
invert = [False, False, False, True, True][start:stop:step]

fig = plt.figure()
ax = fig.add_subplot(111)

error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:grass green", "pink"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple"]

labels = ["Deviation of sample", "Mean estimated error", "Error of mean"]

NSIDES = [1, 2, 4, 8, 16, 32, 64, 128]
x = [1/2, *NSIDES]

for file_num in range(len(files)):
    data = np.load(files[file_num])
    mean = np.zeros(len(data[0]))
    variance = np.zeros(len(data[0]))
    weights = np.ones(len(data[0]))

    print(data[0])

    for data_set in data:
        mean += data_set[:, 0] * (1 / len(data))
        variance += (data_set[:, 1] ** 2) * (1 / len(data) ** 2)

    print(np.shape(data))
    mean = mean / weights
    #variance = variance / weights

    #print(mean, variance)

    mean_error = np.zeros(np.shape(data)[1])
    variance = np.zeros(np.shape(data)[1])
    for i in range(np.shape(data)[1]):
        mean_error[i] = np.mean(data[:, i, 1])
        variance[i] = np.var(data[:, i, 0])

    ax.errorbar(x, mean, np.sqrt(variance), marker=None, ecolor=error_bar_colors[0],
                ls="none", color=line_colors[0], capsize=3, capthick=1, label=labels[0])
    ax.errorbar(x, mean, mean_error, marker=None, ecolor=error_bar_colors[1],
                ls="none", color=line_colors[0], capsize=3, capthick=1, label=labels[1])
    ax.errorbar(x, mean, np.sqrt(variance / np.shape(data)[0]), marker=None, ecolor=error_bar_colors[2],
                ls="none", color=line_colors[0], capsize=3, capthick=1, label=labels[2])
    ax.plot(x, mean, color=line_colors[0], marker="+")

    ax.set_xscale("log", base=2)
    ax.set_xlim(1/2 * np.sqrt(1/2), NSIDES[-1] * np.sqrt(2))
    #f_s = 0.01771962453092031  # ACT point
    f_s = 0.0138443493342983  # Planck point

    if not invert[file_num]:
        excess = alpha[file_num] * (1 - f_s) / (1 + alpha[file_num] * f_s)
    else:
        excess = - alpha[file_num] * (1 - f_s) / (1 + alpha[file_num] - alpha[file_num] * f_s)
    ax.plot([1/2, NSIDES[-1]], np.ones(2) * excess, color="k")


ax.set_xticks([0.5] + NSIDES, ["C"] + NSIDES)

#plt.ylim(-0.005, 0.005)

#ax.legend()
ax.set_xlabel("NSIDE")
ax.set_ylabel("Excess")
ax.set_title(f"A comparison of the errors for 100 realisations of 10,000 "
             f"clusters (Planck point)")
#ax.set_title(f"The ACT point mask, analysed with biased samples of 10,000"
#             f" clusters")
#plt.savefig("planck_random_10k.pdf")
plt.show()

print(np.mean(np.sqrt(variance) / mean_error))

# To do:
# Regenerate random catalogue
# Gen random catalogue with negative excess
# Add other catalogues to plot
# Create scatter from ACT data and pipeline
