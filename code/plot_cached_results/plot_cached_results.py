import matplotlib.pyplot as plt
import numpy as np
import toolkit.toolkit as toolkit

toolkit.plt_use_tex()
error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:mint green", "pink"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple"]
run_const = True

fig = plt.figure()
ax = fig.add_subplot(111)

NSIDES = [1, 2, 4, 8, 16, 32, 64, 128]
labels = ["1", "2", "3"]
title = "Title"
save_path = "test.png"
y_axis_label = "Excess"
files = ["./act_r/bin/act_r_bin_1.npy", "./act_r/bin/act_r_bin_2.npy", "./act_r/bin/act_r_bin_3.npy"]

for file in files:
    result_set = np.load(file)
    color = line_colors[files.index(file)]
    error_bar_color = error_bar_colors[files.index(file)]
    label = labels[files.index(file)]
    for i in range(len(result_set)):
        x = NSIDES.copy()
        if run_const:
            x = [0.5] + x
        print(result_set[i][1])
        print()
        ax.errorbar(x, result_set[i][0], result_set[i][1], marker="+", ecolor=error_bar_color,
                    ls="none", color=color, capsize=3, capthick=1, label=label)
        ax.plot(x, result_set[i][0], color=color)

ax.set_xscale("log", base=2)
ax.set_xlim(1/2 * np.sqrt(1/2), NSIDES[-1] * np.sqrt(2))

ax.plot([1/2, NSIDES[-1]], np.zeros(2), color="k")
ax.set_xticks([0.5] + NSIDES, ["C"] + NSIDES)

#plt.ylim(-0.005, 0.005)

ax.legend()
ax.set_xlabel("NSIDE")
ax.set_ylabel(y_axis_label)
ax.set_title(title)
plt.savefig(save_path)
plt.show()
