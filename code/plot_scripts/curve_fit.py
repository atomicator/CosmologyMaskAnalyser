import numpy as np
import matplotlib.pyplot as plt
import scipy
from toolkit import toolkit
import matplotlib

toolkit.plt_use_tex()


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


def inverse_linear_one_d(x, alpha, beta):
    return (1 / beta) * (x - alpha)


def linear_two_d(x, alpha, beta, gamma):
    return alpha + x[0] * beta + x[1] * gamma


def power_law_two(x, alpha, beta, gamma, delta):
    return alpha * x ** beta + gamma * x ** delta


def load_richness_redshift(richness, redshift):
    if richness <= min_richness:
        return False
    else:
        sdss_richness.append(richness)
        sdss_redshift.append(redshift)
        return True


def exponential(x, alpha, beta):
    return alpha * np.exp(x * beta)


def inverse_exponential(x, alpha, beta):
    return beta * np.log(x / alpha)


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

log_richness = np.log10(richness)
log_snr = np.log10(snr)

log_richness, log_snr = log_snr, log_richness

log_z = np.log(1+redshift)

popt_one, pcov_one = scipy.optimize.curve_fit(linear_one_d, log_richness, log_snr, sigma=log_s_err)

print(popt_one, np.sqrt(np.diag(pcov_one)))
print(f"Chi^2: {np.sum(np.square((linear_one_d(log_richness, *popt_one) - log_snr) / log_s_err))}")
print(f"Chi^2 / N: {np.sum(np.square((linear_one_d(log_richness, *popt_one) - log_snr) / log_s_err)) / len(snr)}")

popt_two, pcov_two = scipy.optimize.curve_fit(linear_two_d, np.array((log_richness, log_z)), log_snr, sigma=log_s_err)

print(popt_two, np.sqrt(np.diag(pcov_two)))
print(f"Chi^2: {np.sum(np.square((linear_two_d((log_richness, log_z), *popt_two) - log_snr) / log_s_err))}")
print(f"Chi^2 / N: {np.sum(np.square((linear_two_d((log_richness, log_z), *popt_two) - log_snr) / log_s_err)) / len(snr)}")

popt_three, pcov_three = scipy.optimize.curve_fit(power_law_two, log_richness, log_snr, sigma=log_s_err)

print(popt_three, np.sqrt(np.diag(pcov_three)))
print(f"Chi^2: {np.sum(np.square((power_law_two(log_richness, *popt_three) - log_snr) / log_s_err))}")
print(f"Chi^2 / N: {np.sum(np.square((power_law_two(log_richness, *popt_three) - log_snr) / log_s_err)) / len(snr)}")

popt_four, pcov_four = scipy.optimize.curve_fit(exponential, log_richness, log_snr, sigma=log_s_err)

print(popt_four, np.sqrt(np.diag(pcov_four)))
print(f"Chi^2: {np.sum(np.square((exponential(log_richness, *popt_four) - log_snr) / log_s_err))}")
print(f"Chi^2 / N: {np.sum(np.square((exponential(log_richness, *popt_four) - log_snr) / log_s_err)) / len(snr)}")

sorted_richness = np.sort(log_richness)
sorted_snr = np.zeros(np.shape(snr))
log_richness = list(log_richness)
for i in range(len(snr)):
    sorted_snr[i] = log_snr[log_richness.index(sorted_richness[i])]
print(len(log_richness))
variance_bins = 20
sigma = np.zeros(np.shape(snr))
deviation = []
for i in range(variance_bins):
    bin_min, bin_max = int(i * len(snr) / variance_bins), int((i + 1) * len(snr) / variance_bins)
    deviation.append(np.std(sorted_snr[bin_min:bin_max]))
    sigma[bin_min:bin_max] = np.ones(np.shape(sigma[bin_min:bin_max])) * deviation[-1]

print("test")

popt_one, pcov_one = scipy.optimize.curve_fit(linear_one_d, sorted_richness, sorted_snr, sigma=sigma)

print(popt_one, np.sqrt(np.diag(pcov_one)))
print(f"Chi^2: {np.sum(np.square((linear_one_d(sorted_richness, *popt_one) - sorted_snr) / sigma))}")
print(f"Chi^2 / N: {np.sum(np.square((linear_one_d(sorted_richness, *popt_one) - sorted_snr) / sigma)) / len(snr)}")

#popt_three, pcov_three = scipy.optimize.curve_fit(power_law_two, log_richness, sorted_snr, sigma=sigma)

#print(popt_three, np.sqrt(np.diag(pcov_three)))
#print(f"Chi^2: {np.sum(np.square((power_law_two(sorted_richness, *popt_three) - sorted_snr) / sigma))}")
#print(f"Chi^2 / N: {np.sum(np.square((power_law_two(sorted_richness, *popt_three) - sorted_snr) / sigma)) / len(snr)}")

popt_four, pcov_four = scipy.optimize.curve_fit(exponential, sorted_richness, sorted_snr, sigma=sigma)

lin_reg = scipy.stats.linregress(log_richness, log_snr)

print(lin_reg.slope, lin_reg.intercept, lin_reg.rvalue, lin_reg.pvalue,
      lin_reg.stderr, lin_reg.intercept_stderr)

print(popt_four, np.sqrt(np.diag(pcov_four)))
print(f"Chi^2: {np.sum(np.square((exponential(sorted_richness, *popt_four) - sorted_snr) / sigma))}")
print(f"Chi^2 / N: {np.sum(np.square((exponential(sorted_richness, *popt_four) - sorted_snr) / sigma)) / len(snr)}")

print(f"N: {len(sorted_richness)}")
plt.scatter(np.power(10, log_richness), np.power(10, log_snr), color="r", marker="+")
func = linear_one_d
inverse_func = inverse_linear_one_d
#func = exponential
#inverse_func = inverse_exponential
pcov = np.matrix(((lin_reg.intercept_stderr ** 2, 0.0), (0.0, lin_reg.stderr ** 2)))
params = (lin_reg.intercept, lin_reg.slope)
#pcov = pcov_four
#params = popt_four
rotated = True

if not rotated:
    x = np.linspace(1.5, 3, 1000)
else:
    x = np.linspace(0.55, 1.8, 1000)

middle = np.power(10, func(x, *params))
lower = np.power(10, func(x, *(params - np.sqrt(np.diag(pcov)))))
upper = np.power(10, func(x, *(params + np.sqrt(np.diag(pcov)))))
x = np.power(10, x)

plt.plot(x, middle, color="b")
plt.plot(x, lower, color="b", linestyle="dashed")
plt.plot(x, upper, color="b", linestyle="dashed")

if not rotated:
    plt.xlabel("Richness")
    plt.ylabel("SNR")
else:
    plt.ylabel("Richness")
    plt.xlabel("SNR")
plt.title("Power law fit to SNR - richness data")
plt.xscale("log")
plt.yscale("log")
x_ticks = [4, 6, 8, 10, 20, 40, 60]
plt.xticks(x_ticks, labels=x_ticks)
y_ticks = [20, 40, 60, 80, 100, 200]
plt.yticks(y_ticks, labels=y_ticks)
plt.show()

print(np.mean(sigma), np.std(sigma))
print(np.min(richness))
plt.savefig("snr_richness.pdf")
#exit()

plt.clf()
#plt.scatter(richness, snr, color="r", marker="+")
#x = np.exp(x)
#plt.plot(x, np.exp(middle), color="b")
#plt.plot(x, np.exp(lower), color="b", linestyle="dashed")
#plt.plot(x, np.exp(upper), color="b", linestyle="dashed")

bins = np.linspace(0, 30, 31)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
plt.stairs(sdss_height, bins, label=rf"ACT, complete $(N = {len(snr)})$")
print(np.sum(sdss_height))
sdss_height = np.zeros(len(bins) - 1)
snr = snr[snr > 5]
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
plt.stairs(sdss_height, bins, label=rf"ACT, SNR $>$ 5 $(N = {len(snr)})$")
print(np.sum(sdss_height))

min_richness = 0
sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_richness = np.array(sdss_richness)
snr = np.power(10, inverse_func(np.log10(sdss_richness), *params))

print(snr)
print(np.min(snr), np.max(snr))

sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
print(sdss_height)
plt.stairs(sdss_height, bins, label=rf"SDSS, complete $(N = {len(snr)})$")

#min_richness = np.exp((np.log(5) - popt_one[0]) / popt_one[1])
min_snr = 4
bin_search_min = 0
bin_search_max = 1000
bin_search_iterations = 100
for i in range(bin_search_iterations):
    if inverse_func((bin_search_max + bin_search_min) / 2, *params) > np.log10(min_snr):
        bin_search_max = (bin_search_min + bin_search_max) / 2
    else:
        bin_search_min = (bin_search_min + bin_search_max) / 2
min_richness = 10 ** ((bin_search_max + bin_search_min) / 2)
print(f"min richness: {min_richness}")
sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_richness = np.array(sdss_richness)
snr = np.power(10, inverse_func(np.log10(sdss_richness), *params))
print(f"N: {len(snr)}")
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
plt.stairs(sdss_height, bins, label=fr"SDSS, $f(\lambda) > {min_snr}$ " + f"$(N = {len(snr)})$")

min_richness = 20
print(f"min richness: {min_richness}")
sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_richness = np.array(sdss_richness)
snr = np.power(10, inverse_func(np.log10(sdss_richness), *params))
print(f"N: {len(snr)}")
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1])))
sdss_height = 100 * sdss_height / np.sum(sdss_height)
plt.stairs(sdss_height, bins, label=fr"SDSS, $\lambda > {min_richness}$ " + f"$(N = {len(snr)})$")

plt.legend()
plt.title("The cluster count of the catalogues, as a function of ACT SNR")
#plt.yscale("log")
plt.ylabel("Percentage")
plt.xlabel("SNR")

#plt.xlabel("ln SNR")
#plt.ylabel("ln f(richness)")
#plt.title("Measured vs mapped SNR")

#plt.xlim(0, 200)
#plt.ylim(0, 40)

plt.savefig("cluster_count_snr.pdf")
plt.show()
