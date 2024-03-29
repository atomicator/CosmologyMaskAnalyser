import numpy as np
import matplotlib.pyplot as plt
import scipy
from toolkit import toolkit


def load_act_richness(flag, cluster_richness, cluster_snr, cluster_redshift, y_error, y):
    if flag == 0.0:
        return False
    else:
        richness.append(cluster_richness)
        snr.append(cluster_snr)
        redshift.append(cluster_redshift)
        y_err.append(y_error)
        y_c.append(y)


def linear_one_d(x, alpha, beta):
    return alpha + x * beta


def linear_two_d(x, alpha, beta, gamma):
    return alpha + x[0] * beta + x[1] * gamma


def load_richness_redshift(richness, redshift):
    if richness <= min_richness:
        return False
    else:
        sdss_richness.append(richness)
        sdss_redshift.append(redshift)
        return True

min_richness = 0
act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1, table=True)
sdss_richness = []
sdss_redshift = []
act_cat.load_with_selection(load_richness_redshift, ("fixed_y_c", "fixed_SNR"), True)
sdss_richness = np.array(sdss_richness) * 1e-4
sdss_redshift = np.array(sdss_redshift)
act_noise = np.mean(sdss_richness / sdss_redshift)

act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1, table=True)
richness = []
snr = []
redshift = []
y_err = []
y_c = []

act_cat.load_with_selection(load_act_richness, ("RM", "RM_LAMBDA", "SNR", "redshift", "err_y_c", "y_c"), True)
snr = np.array(snr)
richness = np.array(richness)
redshift = np.array(redshift)
y_err = np.array(y_err)
y_c = np.array(y_c)
s_err = y_err * (snr / y_c)
log_s_err = s_err / snr

log_richness = np.log(richness)
log_snr = np.log(snr)
log_z = np.log(1+redshift)

popt_one, pcov_one = scipy.optimize.curve_fit(linear_one_d, log_richness, log_snr, sigma=log_s_err)

print(popt_one, np.sqrt(np.diag(pcov_one)))
print(f"Chi^2: {np.sum(np.square((linear_one_d(log_richness, *popt_one) - log_snr) / log_s_err))}")
print(f"Chi^2 / N: {np.sum(np.square((linear_one_d(log_richness, *popt_one) - log_snr) / log_s_err)) / len(snr)}")

popt_two, pcov_two = scipy.optimize.curve_fit(linear_two_d, np.array((log_richness, log_z)), log_snr, sigma=log_s_err)

print(popt_two, np.sqrt(np.diag(pcov_two)))
print(f"Chi^2: {np.sum(np.square((linear_two_d((log_richness, log_z), *popt_two) - log_snr) / log_s_err))}")
print(f"Chi^2 / N: {np.sum(np.square((linear_two_d((log_richness, log_z), *popt_two) - log_snr) / log_s_err)) / len(snr)}")

plt.scatter(log_snr, linear_one_d(log_richness, *popt_one), marker="+", color="r")
plt.plot((1.5, 3.5), (1.5, 3.5))

bins = np.linspace(0, 25, 51)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
#plt.stairs(sdss_height, bins, label="ACT, complete")
print(np.sum(sdss_height))
sdss_height = np.zeros(len(bins) - 1)
snr = snr[snr > 5]
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
#plt.stairs(sdss_height, bins, label="ACT, SNR > 5")
print(np.sum(sdss_height))

min_richness = 0
sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_richness = np.array(sdss_richness)
snr = np.exp(linear_one_d(np.log(sdss_richness), *popt_one))

print(snr)
print(np.min(snr), np.max(snr))

sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
#plt.stairs(sdss_height, bins, label="SDSS, complete")

min_richness = np.exp((np.log(5) - popt_one[0]) / popt_one[1])
print(min_richness)
sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_richness = np.array(sdss_richness)
snr = np.exp(linear_one_d(np.log(sdss_richness), *popt_one))
print(f"N: {len(snr)}")
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
#plt.stairs(sdss_height, bins, label=r"SDSS, $\lambda > \lambda_{min}$")

plt.legend()
#plt.title("The cluster count of the catalogues, as a function of ACT SNR")
#plt.yscale("log")
#plt.ylabel("Percentage")
#plt.xlabel("SNR")

plt.xlabel("ln SNR")
plt.ylabel("ln f(richness)")
plt.title("Measured vs mapped SNR")

plt.show()
