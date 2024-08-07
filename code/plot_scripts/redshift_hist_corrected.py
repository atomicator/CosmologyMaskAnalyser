import numpy as np
import matplotlib.pyplot as plt
from toolkit import toolkit
import hmf
import copy

toolkit.plt_use_tex()


def load_richness_redshift(richness, redshift):
    if richness <= min_richness:
        return False
    else:
        sdss_richness.append(richness)
        sdss_redshift.append(redshift)
        return True


def planck_cat_selection_func(snr, pipe_det, cosmo, redshift, min_snr, pipe_det_val):
    if (snr > min_snr) and (pipe_det == pipe_det_val) and cosmo:
        sdss_redshift.append(redshift)
        return True


bins = np.linspace(0, 1, 51)

min_richness = 0

act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=0, table=True)
sdss_richness = []
sdss_redshift = []
act_cat.load_with_selection(load_richness_redshift, ("SNR", "redshift"), True)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = (np.sum(np.float_(np.bitwise_and(sdss_redshift > bins[i], sdss_redshift < bins[i + 1]))))# /
                      #(len(sdss_redshift) / 100))
#plt.stairs(sdss_height, bins, label=r"ACT" + f" $(N = {len(sdss_richness)})$", color="red")
#plt.stairs(sdss_height * (1 - 1 / np.sqrt(len(sdss_redshift))), bins, linestyle="dashed", color="orange")
#plt.stairs(sdss_height * (1 + 1 / np.sqrt(len(sdss_redshift))), bins, linestyle="dashed", color="orange")
f_s = 0.01771962453092031 # ACT

alpha = 0.13830263
beta = -8.24858781
poisson = 0.016

z = np.linspace(bins[0] + 1 / len(bins), bins[-1] - 1 / len(bins), len(bins) - 1)
excess = alpha * np.exp(beta * (z - 0.4))

plt.stairs(sdss_height * f_s * (1 + excess), bins, label=r"Calculated $N_{m}$", color="red")
plt.stairs(sdss_height * f_s, bins, label=r"Assumed $N_{m}$", color="orange", linestyle="solid")
plt.stairs(sdss_height * f_s * excess, bins, label=r"$N_{0} f_{s} E(z)$", color="xkcd:crimson", linestyle="solid")
plt.stairs(sdss_height, bins, label=r"$N_{0}/\sqrt{n}$", color="black", linestyle="solid")
#plt.stairs(sdss_height * (1 - 1 / np.sqrt(len(sdss_redshift))), bins, linestyle="dashed", color="orange")
#plt.stairs(sdss_height * (1 + 1 / np.sqrt(len(sdss_redshift))), bins, linestyle="dashed", color="orange")
print(np.sum(sdss_height * f_s * excess))

#plt.legend()
plt.title("The masked cluster count for the ACT point-source mask, as a function of redshift")
#plt.yscale("log")
plt.ylabel("Masked clusters per bin")
plt.xlabel("Redshift")
plt.savefig("redshift_histogram_corrected.pdf")
plt.show()
plt.clf()

plt.stairs(sdss_height * (1 + f_s * (1 + excess)), bins, label=r"Calculated $N_{T}$", linestyle="dashdot", color="red")
plt.stairs(sdss_height * (1 + f_s), bins, label=r"Assumed $N_{T}$", color="orange", linestyle="dashed")
plt.stairs(sdss_height * (f_s * excess), bins, label=r"$N_{T}$ bias", color="xkcd:crimson", linestyle="solid")
plt.stairs(sdss_height, bins, label=r"Detected clusters", color="black", linestyle="dotted")
plt.stairs(sdss_height / np.sqrt(len(sdss_redshift)), bins, label=r"Best-case error", color="black", linestyle="solid")
plt.yscale("log")
plt.legend()
plt.title("The total cluster count for the ACT point-source mask, as a function of redshift")
#plt.yscale("log")
plt.ylabel("Total clusters per bin")
plt.xlabel("Redshift")
plt.xlim(0, 1)
plt.savefig("redshift_number_count.pdf")
plt.show()
x = np.linspace(1 / (2 * len(bins)), 1 - 1 / (2 * len(bins)), len(bins) - 1)
plt.clf()
def correction(x):
    #r = (f_s - 1 - excess(x) * (f_s - 1)) / (f_s - 1 + excess(x))
    #return f_s * (r - 1)
    f_c = (1 + excess(x)) * f_s
    r = (f_c * (1 - f_s)) / (f_s * (1 - f_c))
    return f_s * (r - 1)


#def excess(x):
#    return alpha * np.exp(beta * (x - 0.4))


alpha = 0.138
beta = -8.2
f_s = 0.01771962453092031  # ACT
#plt.stairs(sdss_height * (1 - correction(x)), bins, linestyle="dotted", color="xkcd:crimson")
#plt.stairs(sdss_height * (1 + correction(x)), bins, linestyle="dotted", color="xkcd:crimson")


"""cat = toolkit.ClusterCatalogue("../../data/HFI_PCCS_SZ-union_R2.08.fits", hdu=1)
#cat.load_lon_lat()
min_snr = 6
pipe_det_value = 111
sdss_richness = []
sdss_redshift = []
cat.load_with_selection(planck_cat_selection_func, ["SNR", "PIPE_DET", "COSMO", "REDSHIFT"], True, min_snr, pipe_det_value)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = (np.sum(np.float_(np.bitwise_and(sdss_redshift > bins[i], sdss_redshift < bins[i + 1]))) / (len(sdss_redshift) / 100))
#plt.stairs(sdss_height, bins, label=r"Planck" + f" $(N = {len(sdss_redshift)})$", color="xkcd:electric blue")
#plt.stairs(sdss_height * (1 - 1 / np.sqrt(len(sdss_redshift))), bins, linestyle="dashed", color="xkcd:aqua blue")
#plt.stairs(sdss_height * (1 + 1 / np.sqrt(len(sdss_redshift))), bins, linestyle="dashed", color="xkcd:aqua blue")

#plt.stairs(sdss_height * (1 - 0.0028), bins, linestyle="dotted", color="xkcd:navy blue")
#plt.stairs(sdss_height * (1 + 0.0028), bins, linestyle="dotted", color="xkcd:navy blue")
"""
plt.legend()
plt.title("The number of clusters masked by the ACT point mask, as a function of redshift")
#plt.yscale("log")
plt.ylabel("Clusters per bin")
plt.xlabel("Redshift")
#plt.savefig("redshift_histogram_corrected.pdf")
#plt.show()
plt.clf()

sdss_redshift = np.array(sdss_redshift)

plt.plot(x, 100 * f_s * (1 + excess), color="red", label="Modelled masked clusters")
plt.plot(x, 100 * f_s * np.ones(np.shape(x)), color="orange", label="Assumed masked clusters")
plt.plot(x, 100 * f_s * excess, color="xkcd:crimson", label="Difference")
plt.plot(x, 100 * poisson * np.ones(np.shape(x)), color="black", label=r"Best-case Error estimate")
plt.stairs(100 * np.power(np.array((np.sum(sdss_redshift < 0.5), np.sum(sdss_redshift > 0.5))), -1/2),
           np.linspace(0, 1, 3), color="grey", linestyle="dashdot", label="Error estimate, 2 bins")
plt.stairs(100 * np.power(np.array((np.sum(sdss_redshift < 0.25),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.25, sdss_redshift < 0.5)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.5, sdss_redshift < 0.75)),
                                    np.sum(sdss_redshift > 0.75)
                                    )), -1/2), np.linspace(0, 1, 5), color="grey", linestyle="dashed",
           label="Error estimate, 4 bins")
plt.stairs(100 * np.power(np.array((np.sum(sdss_redshift < 0.125),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.125, sdss_redshift < 0.25)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.25, sdss_redshift < 0.375)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.375, sdss_redshift < 0.5)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.5, sdss_redshift < 0.625)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.625, sdss_redshift < 0.75)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.75, sdss_redshift < 0.875)),
                                    np.sum(sdss_redshift > 0.875)
                                    )), -1/2), np.linspace(0, 1, 9), color="grey", linestyle="dotted",
           label="Error estimate, 8 bins")
plt.legend()
plt.title("Fractional correction to the number-count ACT catalogue as a function of redshift")
plt.xlabel("$z$")
plt.ylabel(r"Correction to number count estimate $(\%)$")
plt.xlim(0, 1)
plt.savefig("fractional_correction_redshift.pdf")
plt.show()
