import matplotlib.pyplot as plt
import numpy as np
import toolkit.toolkit as toolkit

toolkit.plt_use_tex()
error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:mint green", "xkcd:light purple", "xkcd:burnt sienna", "xkcd: pink"]
#line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple", "xkcd:reddish brown", "xkcd:hot pink"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple", "xkcd:reddish brown", "xkcd:hot pink", "orange", "xkcd:aqua blue"]
#error_bar_colors = line_colors
run_const = True

fig = plt.figure()
ax = fig.add_subplot(111)

NSIDES = [1, 2, 4, 8, 16, 32, 64, 128]
save_path = "test.pdf"
#save_path = "final_excess_complete.pdf"
#title = "The excess of the ACT point mask, as a function of richness (bin)"
#save_path = "./act_bin_richness.png"
#files = ["./act_r/bin/act_r_bin_1.npy", "./act_r/bin/act_r_bin_2.npy", "./act_r/bin/act_r_bin_3.npy",
#         "./act_r/bin/act_r_bin_4.npy", "./act_r/min/act_r_min_4.npy"]

#title = "The excess of the Planck point mask, as a function of richness (bin)"
#save_path = "./planck_bin_richness.png"
#files = ["./planck_r/bin/planck_r_bin_1.npy", "./planck_r/bin/planck_r_bin_2.npy", "./planck_r/bin/planck_r_bin_3.npy",
#         "./planck_r/bin/planck_r_bin_4.npy", "./planck_r/min/planck_r_min_4.npy"]

#title = "The excess of the ACT point mask, as a function of minimum richness"
#save_path = "./act_min_richness.png"
#files = ["./act_r/min/act_excess.npy", "./act_r/min/act_r_min_1.npy", "./act_r/min/act_r_min_2.npy",
#         "./act_r/min/act_r_min_3.npy", "./act_r/min/act_r_min_4.npy"]

#files = ["./act_r/min/act_excess.npy", "./act_r/min/act_r_min_1.npy"]

#title = "The excess of the Planck point mask, as a function of minimum richness"
#save_path = "./planck_min_richness.png"
#files = ["./planck_r/bin/planck_excess.npy", "./planck_r/min/planck_r_min_1.npy", "./planck_r/min/planck_r_min_2.npy",
#         "./planck_r/min/planck_r_min_3.npy", "./planck_r/min/planck_r_min_4.npy"]

#title = "The excess of the ACT point mask, as a function of minimum redshift"
#save_path = "./act_min_z.png"
#save_path = "test.png"
#files = ["./act_z/min/act_excess.npy", "./act_z/min/act_z_min_1.npy", "./act_z/min/act_z_min_2.npy",
#         "./act_z/min/act_z_min_3.npy", "./act_z/min/act_z_min_4.npy"]
#files = ["./act_r/min/act_r_min_2.npy"]

#title = "The excess of the Planck point mask, as a function of minimum redshift"
#save_path = "./planck_min_z.png"
#files = ["./planck_z/bin/planck_excess.npy", "./planck_z/min/planck_z_min_1.npy", "./planck_z/min/planck_z_min_2.npy",
#         "./planck_z/min/planck_z_min_3.npy", "./planck_z/min/planck_z_min_4.npy"]

#labels = [r"$\lambda < 10$", r"$10 < \lambda < 20$", r"$20 < \lambda < 50$", r"$50 < \lambda < 100$", r"$100 < \lambda$"]
#labels = [r"All results", r"$10 < \lambda$", r"$20 < \lambda$", r"$50 < \lambda$", r"$100 < \lambda$"]
#labels = [r"$\lambda < 10$", r"$10 < \lambda < 20$", r"$20 < \lambda < 50$", r"$50 < \lambda < 100$", r"$100 < \lambda$"]
#labels = [r"All results", "$z > 0.319$", "$z > 0.319$", "$z > 0.420$", "$z > 0.489$", "$z > 0.500$"]
#labels = ["test"]

#files = ["./../../code/binned_results/planck_excess.npy"]

#files = ["./rotations/planck/sdss_planck_lon_shift_rot_10.npy", "./rotations/planck/sdss_planck_lon_shift_rot_50.npy", "./rotations/planck/sdss_planck_lon_shift_rot_100.npy", "./rotations/planck/sdss_planck_lon_shift_rot_150.npy"]
#labels = [r"$0^{\circ}$", r"$10^{\circ}$", r"$20^{\circ}$", r"$30^{\circ}$", r"$40^{\circ}$"]
#files = ["./act_r/min/act_excess.npy", "./rotations/sdss_act_lon_shift_rot_1.npy", "./rotations/sdss_act_lon_shift_rot_2.npy", "./rotations/sdss_act_lon_shift_rot_3.npy", "./rotations/sdss_act_lon_shift_rot_4.npy"]
#labels = [r"$0^{\circ}$", r"$1^{\circ}$", r"$2^{\circ}$", r"$3^{\circ}$", r"$4^{\circ}$"]
#title = "ACT Rotation"

#title = "ACT excess, as a function of ACT SNR ($s$)"
#files = ["./act_s/act_s_3.npy", "./act_s/act_s_4.npy", "./act_s/act_s_5.npy", "./act_s/act_s_6.npy", "./act_s/act_s_7.npy"]#, "./act_s/act_s_8.npy"]
#files = ["./act_s/act_s_4.npy", "./act_s/act_s_5.npy", "./act_s/act_s_6.npy", "./act_s/act_s_7.npy", "./act_s/act_s_8.npy"]

#title = "ACT excess, as a function of ACT SNR ($s$)"
#files = ["./act_s/act_s_12.npy", "./act_s/act_s_23.npy", "./act_s/act_s_34.npy", "./act_s/act_s_45.npy",
#         "./act_s/act_s_56.npy", "./act_s/act_s_67.npy", "./act_s/act_s_78.npy", "./act_s/act_s_8.npy"]
#files = ["./act_s/act_s_13.npy", "./act_s/act_s_35.npy", "./act_s/act_s_57.npy", "./act_s/act_s_7.npy"]
#files = ["./planck_s/planck_s_4.npy", "./planck_s/planck_s_5.npy", "./planck_s/planck_s_6.npy", "./planck_s/planck_s_7.npy", "./planck_s/planck_s_8.npy"]

#labels = ["$1 < s < 2$", "$2 < s < 3$", "$3 < s < 4$", "$4 < s < 5$", "$5 < s < 6$", "$6 < s < 7$", "$7 < s < 8$", "$8 < s$"]
#labels = ["$1 < s < 3$", "$3 < s < 5$", "$5 < s < 7$", "$7 < s$"]

#labels = ["Point", "Galactic"]
#files = ["./act_r/act_r_20+.npy"]
#files = ["./planck_r/planck_r_20+.npy"]
#title = "The measured excess for the two ACT masks"

#files = ["./act_r/act_r_5_10.npy", "./act_r/act_r_10_20.npy", "./act_r/act_r_20_40.npy", "./act_r/act_r_40_80.npy"]
#files = ["./planck_r/planck_r_5_10.npy", "./planck_r/planck_r_10_20.npy", "./planck_r/planck_r_20_40.npy", "./planck_r/planck_r_40_80.npy"]
#files = ["./act_z/act_z_0_15.npy", "./act_z/act_z_15_3.npy", "./act_z/act_z_3_45.npy", "./act_z/act_z_45_6.npy"]
files = ["./act_z/act_z_0_3.npy", "./act_z/act_z_3_45.npy", "./act_z/act_z_45_6.npy"]
#files = ["./planck_z/planck_z_0_3.npy", "./planck_z/planck_z_3_45.npy", "./planck_z/planck_z_45_6.npy"]
#labels = [r"$5 < \lambda < 10$", r"$10 < \lambda < 20$", r"$20 < \lambda < 40$", r"$40 < \lambda < 80$", r"$80 < \lambda$"]
#labels = [r"$0 < z < 0.15$", r"$0.15 < z < 0.3$", r"$0.3 < z < 0.45$", r"$0.45 < z < 0.6$"]
labels = [r"$0 < z < 0.3$", r"$0.3 < z < 0.45$", r"$0.45 < z < 0.6$"]
title = "The excess of the ACT point mask as a function of redshift"

point_only = False
y_axis_label = "Excess"

nside_eight_excess = []
nside_eight_error = []
nside_eight_snr = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5][::-1]
#nside_eight_snr = [2, 4, 6, 8][::-1]

n = 0
for file in files:
    result_set = np.load(file)
    print(result_set)
    for i in range(len(result_set)):
        x = NSIDES.copy()
        if run_const:
            x = [0.5] + x
        print(result_set[i][1])
        print()
        nside_eight_excess.append(result_set[i][0][x.index(8)])
        nside_eight_error.append(result_set[i][1][x.index(8)])
        ax.errorbar(x, result_set[i][0], result_set[i][1], marker="+", ecolor=error_bar_colors[n],
                    ls="none", color=line_colors[n], capsize=3, capthick=1, label=labels[n])
        ax.plot(x, result_set[i][0], color=line_colors[n])
        n += 1
        if point_only:
            break

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

print(len(nside_eight_excess), len(nside_eight_error), len(nside_eight_snr))

plt.clf()
exit()
plt.errorbar(nside_eight_snr, nside_eight_excess, nside_eight_error)
plt.show()
