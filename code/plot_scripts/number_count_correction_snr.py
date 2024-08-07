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


bins = np.linspace(3.0, 20, 51)

min_richness = 0

act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=0, table=True)
sdss_richness = []
sdss_redshift = []
act_cat.load_with_selection(load_richness_redshift, ("SNR", "redshift"), True)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = (np.sum(np.float_(np.bitwise_and(sdss_richness > bins[i], sdss_richness < bins[i + 1]))))# /
                      #(len(sdss_richness) / 100))
#plt.stairs(sdss_height, bins, label=r"ACT" + f" $(N = {len(sdss_richness)})$", color="red")
#plt.stairs(sdss_height * (1 - 1 / np.sqrt(len(sdss_richness))), bins, linestyle="dashed", color="orange")
#plt.stairs(sdss_height * (1 + 1 / np.sqrt(len(sdss_richness))), bins, linestyle="dashed", color="orange")

x = np.linspace(bins[0] + 1 / (2 * len(bins)), bins[-1] - 1 / (2 * len(bins)), len(bins) - 1)


def correction(x):
    #r = (f_s - 1 - excess(x) * (f_s - 1)) / (f_s - 1 + excess(x))
    #return f_s * (r - 1)
    #f_c = (1 + excess(x)) * f_s
    #r = (f_c * (1 - f_s)) / (f_s * (1 - f_c))
    #return f_s * (r - 1)
    return f_s * (excess(x) + 1)


def excess(s):
    r = 10 ** 1.34 * s ** 0.57
    #return psi + phi * (r - 30)
    return (psi - 30 * phi) + phi * r * np.exp((sigma * np.log(10)**2) / 2)


psi = 0.116
phi = 0.0026
f_s = 0.01771962453092031  # ACT
poisson = 0.016
sigma = 0.17

sdss_redshift = np.array(sdss_redshift)

plt.stairs(sdss_height * f_s * (1 + excess(x)), bins, label=r"$N_{0} f_{s} (1 + E(s))$", color="red")
plt.stairs(sdss_height * f_s, bins, label=r"$N_{0} f_{s}$", color="orange", linestyle="solid")
plt.stairs(sdss_height * f_s * excess(x), bins, label=r"$N_{0} f_{s} E(s)$", color="xkcd:crimson", linestyle="solid")
plt.stairs(sdss_height * poisson, bins, label=r"$N_{0}/\sqrt{n}$", color="black", linestyle="solid")

plt.legend()
plt.title("The masked cluster count of for the ACT point-source mask, as a function of SNR")
#plt.yscale("log")
plt.ylabel("Cluster count")
plt.xlabel("SNR")
plt.savefig("richness_histogram_corrected.pdf")
plt.show()
plt.clf()

plt.stairs(sdss_height * (1 + f_s * (1 + excess(x))), bins, label=r"Calculated $N_{T}$", linestyle="dashdot", color="red")
plt.stairs(sdss_height * (1 + f_s), bins, label=r"Assumed $N_{T}$", color="orange", linestyle="dashed")
plt.stairs(sdss_height * (f_s * excess(x)), bins, label=r"$N_{T}$ bias", color="xkcd:crimson", linestyle="solid")
plt.stairs(sdss_height, bins, label=r"Detected clusters", color="black", linestyle="dotted")
plt.stairs(sdss_height / np.sqrt(len(sdss_redshift)), bins, label=r"Best-case error", color="black", linestyle="solid")
plt.yscale("log")
plt.legend()
plt.title("The total cluster count for the ACT point-source mask, as a function of SNR")
#plt.yscale("log")
plt.ylabel("Total clusters per bin")
plt.xlabel("SNR")
plt.xlim(bins[0], bins[-1])
plt.savefig("richness_number_count.pdf")
plt.show()
x = np.linspace(1 / (2 * len(bins)), 1 - 1 / (2 * len(bins)), len(bins) - 1)
plt.clf()

exit()
plt.plot(x, 100 * f_s * (1 + excess(x)), color="red", label="Modelled masked clusters")
plt.plot(x, 100 * f_s * np.ones(np.shape(x)), color="orange", label="Assumed masked clusters")
plt.plot(x, 100 * f_s * excess(x), color="xkcd:crimson", label="Difference")
plt.plot(x, 100 * poisson * np.ones(np.shape(x)), color="black", label=r"Standard error estimates")
plt.stairs(100 * np.power(np.array((np.sum(sdss_redshift < 0.5), np.sum(sdss_redshift > 0.5))), -1/2),
           np.linspace(x[0], x[-1], 3), color="grey", linestyle="dashdot", label="Error estimate, 2 bins")
plt.stairs(100 * np.power(np.array((np.sum(sdss_redshift < 0.25),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.25, sdss_redshift < 0.5)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.5, sdss_redshift < 0.75)),
                                    np.sum(sdss_redshift > 0.75)
                                    )), -1/2), np.linspace(x[0], x[-1], 5), color="grey", linestyle="dashed",
           label="Error estimate, 4 bins")
"""plt.stairs(100 * np.power(np.array((np.sum(sdss_redshift < 0.125),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.125, sdss_redshift < 0.25)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.25, sdss_redshift < 0.375)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.375, sdss_redshift < 0.5)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.5, sdss_redshift < 0.625)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.625, sdss_redshift < 0.75)),
                                    np.sum(np.bitwise_and(sdss_redshift > 0.75, sdss_redshift < 0.875)),
                                    np.sum(sdss_redshift > 0.875)
                                    )), -1/2), np.linspace(x[0], x[-1], 9), color="grey", linestyle="dotted",
           label="Error estimate, 8 bins")"""
plt.legend()
plt.title("Fractional correction to the number-count ACT catalogue as a function of SNR")
plt.xlabel("SNR")
plt.ylabel(r"Correction to number count estimate $(\%)$")
plt.savefig(r"fractional_correction_snr.pdf")
plt.show()
