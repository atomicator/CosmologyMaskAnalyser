import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt


def linear_weighting(f_s, _f_c, n, **_kwargs):
    variance = f_s * (1 - f_s) / n
    weights = np.ones(len(n)) / len(n)
    mean = np.sum(f_s * weights) * 100
    mean_error = np.sqrt(np.sum(variance * weights ** 2)) * 100
    return mean, mean_error


def inverse_variance_weighting(f_s, _f_c, n, **_kwargs):
    variance = f_s * (1 - f_s) / n
    weights = (1 / variance) / (np.sum(1 / variance))
    mean = np.sum(f_s * weights) * 100
    mean_error = np.sqrt(np.sum(variance * weights ** 2)) * 100
    return mean, mean_error


def density_weighting(f_s, _f_c, n, **_kwargs):
    print(f_s, _f_c, n)
    variance = f_s * (1 - f_s) / n
    weights = n / (np.sum(n))
    print(np.sum(weights))
    mean = np.sum(f_s * weights) * 100
    mean_error = np.sqrt(np.sum(variance * weights ** 2)) * 100
    return mean, mean_error


def excess_measurement(f_s, f_c, n, skip_n_filter=False, minimum_n=0, **_kwargs):
    print(f"Mean N: {np.mean(n)}")
    if not skip_n_filter:
        bin_filter = np.bitwise_not(np.bitwise_or(f_s == 0, f_s == 1))
        bin_filter = np.bitwise_and(bin_filter, n > minimum_n)
        n = n[bin_filter]
        f_c = f_c[bin_filter]
        f_s = f_s[bin_filter]
    print(np.min(f_s), np.max(f_s))
    excess = (f_c - f_s) / f_s
    print(f"Mean N: {np.mean(n)}, {np.sum(n)}")
    variance = (1 - f_s) / (f_s * n)
    weight = (1 / variance) / np.sum(1 / variance)
    #weight = (np.sin(np.pi * f_s) / variance) / np.sum(np.sin(np.pi * f_s) / variance)
    mean = np.sum(excess * weight)
    mean_error = np.sqrt(np.sum(variance * weight ** 2))
    if len(n) == 0:
        mean, mean_error = np.nan, np.nan
    return mean, mean_error


def excess_measurement_frac(f_s, f_c, n, **_kwargs):
    excess = (f_c - f_s) / f_s
    variance = (1 - f_s) / (f_s * n)
    weight = (1 / variance) / np.sum(1 / variance)
    mean = np.sum(excess * weight)
    mean_error = np.sqrt(np.sum(variance * weight ** 2))
    cluster_masked_frac = (mean + 1) * np.sum(n * f_s) / np.sum(n) * 100
    cluster_masked_frac_error = mean_error * np.sum(n * f_s) / np.sum(n) * 100
    return cluster_masked_frac, cluster_masked_frac_error


def _line(x, a, **_kwargs):
    return x * a


def regression_weighting(f_s, f_c, n, **_kwargs):
    variance = f_c * (1 - f_c) / n
    """plt.errorbar(f_s, f_c, np.sqrt(variance), marker="+", ls="none")
    plt.plot([0, 1], [0, 1])
    plt.show()"""
    #variance = (1 - f_c) / n
    # Some sets are skewed towards zero - finding the unmasked fraction can help this
    f_s = 1 - f_s
    f_c = 1 - f_c
    print(np.min(n), np.max(n))
    weight = 1 / variance
    popt, pcov = scipy.optimize.curve_fit(_line, f_s, f_c, sigma=np.sqrt(1 / weight))
    grad = popt[0]
    grad_error = np.sqrt(pcov[0][0])
    print(f"Grad: {grad} +/-  {grad_error}")
    # Original version didn't have the n dependence, caused weird errors later on
    final = np.array([grad, grad_error]) * np.sum((1 - f_s) * n) * 100 / np.sum(n)
    plt.errorbar(f_s, f_c, np.sqrt(variance), marker="+", ls="none")
    plt.plot([0, 1], [0, 1])
    plt.plot([0, 1], [0, grad])
    plt.show()
    return final


def _skewness_function(x, a, **_kwargs):
    return (a + 1) * x / (a * x + 1)


def skewness_match(f_s, f_c, n, **_kwargs):
    minimum_n = 0
    f_s = f_s[n > minimum_n]
    f_c = f_c[n > minimum_n]
    n = n[n > minimum_n]
    #variance = f_c * (1 - f_c) / n
    variance = 1 / n
    popt, pcov = scipy.optimize.curve_fit(_skewness_function, f_s, f_c, sigma=np.sqrt(variance))
    final = np.array([popt[0], np.sqrt(pcov[0][0])])
    return final


def scatter(f_s, f_c, _n, **_kwargs):
    plt.plot((0, 0.6), (0, 0.6))
    plt.scatter(f_c, f_s, marker="+")
    plt.xlim(0, 0.3)
    plt.ylim(0, 0.3)
    plt.show()


def ratio(_f_s, f_c, n, **_kwargs):
    weight = n / np.sum(n)
    mean = np.sum(f_c * weight) * 100
    variance = np.sum(1/n * weight ** 2) * 100
    return np.array((mean, np.sqrt(variance)))
