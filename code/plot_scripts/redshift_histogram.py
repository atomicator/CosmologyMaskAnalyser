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

sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = (np.sum(np.float_(np.bitwise_and(sdss_redshift > bins[i], sdss_redshift < bins[i + 1]))) /
                      (len(sdss_redshift) / 100))
plt.stairs(sdss_height, bins, label="SDSS, complete" + f" $(N = {len(sdss_richness)})$")


min_richness = 20

sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = (np.sum(np.float_(np.bitwise_and(sdss_redshift > bins[i], sdss_redshift < bins[i + 1]))) /
                      (len(sdss_redshift) / 100))
plt.stairs(sdss_height, bins, label=r"SDSS, $\lambda > 20$" + f" $(N = {len(sdss_richness)})$")

min_richness = 0

act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=0, table=True)
sdss_richness = []
sdss_redshift = []
act_cat.load_with_selection(load_richness_redshift, ("SNR", "redshift"), True)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = (np.sum(np.float_(np.bitwise_and(sdss_redshift > bins[i], sdss_redshift < bins[i + 1]))) /
                      (len(sdss_redshift) / 100))
plt.stairs(sdss_height, bins, label=r"ACT" + f" $(N = {len(sdss_richness)})$")

cat = toolkit.StarCatalogue("../../data/HFI_PCCS_SZ-union_R2.08.fits", hdu=1)
#cat.load_lon_lat()
min_snr = 6
pipe_det_value = 111
sdss_richness = []
sdss_redshift = []
cat.load_with_selection(planck_cat_selection_func, ["SNR", "PIPE_DET", "COSMO", "REDSHIFT"], True, min_snr, pipe_det_value)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = (np.sum(np.float_(np.bitwise_and(sdss_redshift > bins[i], sdss_redshift < bins[i + 1]))) / (len(sdss_redshift) / 100))
plt.stairs(sdss_height, bins, label=r"Planck" + f" $(N = {len(sdss_redshift)})$")

plt.legend()
plt.title("The cluster count of the catalogues, as a function of redshift")
#plt.yscale("log")
plt.ylabel("Percentage")
plt.xlabel("Redshift")
plt.savefig("redshift_histogram.pdf")
plt.show()
