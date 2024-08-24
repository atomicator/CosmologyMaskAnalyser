# This file contains various weight settings - the functions are passed as arguments to the methods in toolkit.
# Note: we decided that excess_measurement was the best option

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
    #print((f_s, f_c, n))
    #print(np.min(f_s), np.max(f_s))
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
    variance = np.sum(1 / n * weight ** 2) * 100
    return np.array((mean, np.sqrt(variance)))


def overdensity(f_s, f_c, n, **_kwargs):
    try:
        valid = np.bitwise_and(np.bitwise_and(n > 1, f_c > 0), f_c < 1)
        n = n[valid]
        f_s = f_s[valid]
        f_c = f_c[valid]
    except IndexError:  # catches the case where n is not an array
        n = np.array((n,))
        f_s = np.array((f_s,))
        f_c = np.array((f_c,))
    plt.scatter(f_s, f_c, 1 + np.int_(n / 10), "r", "+")
    plt.xlabel(r"$f_{s}$")
    plt.ylabel(r"$f_{c}$")
    plt.plot([np.min((np.min(f_s), np.min(f_c))), np.max((np.max(f_s), np.max(f_c)))],
             [np.min((np.min(f_s), np.min(f_c))), np.max((np.max(f_s), np.max(f_c)))], linestyle='--')
    plt.show()
    f_c_error = np.sqrt(f_s * (1 - f_s) / n)
    #  f_c_error = np.sqrt(1 / n)
    #if len(n) > 1:
    #    print("check used")
    #    valid = np.bitwise_and((f_c + f_c_error) < 1, (f_c - f_c_error) > 0)
    #    f_s = f_s[valid]
    #    f_c = f_c[valid]
    #    f_c_error = f_c_error[valid]
    alpha = (f_c - f_s) / (f_s * (1 - f_c))
    print(f"alpha min / max: {np.min(alpha)} / {np.max(alpha)}")
    # alpha_variance = np.abs((1 / (f_s * (1 - f_s)) + (f_s - f_s) / (f_s * (1 - f_s) ** 2)) ** 2 * f_s * (1 - f_s) / n)
    alpha_variance = np.square((((f_s + f_c_error - f_s) / (f_s * (1 - f_s - f_c_error))) -
                                ((f_s - f_c_error - f_s) / (f_s * (1 - f_s + f_c_error)))) / 2)
    weight = (1 / alpha_variance) / np.sum(1 / alpha_variance)
    mean = np.sum(alpha * weight)
    mean_error = np.sqrt(np.sum(alpha_variance * weight ** 2))
    return np.array((mean, mean_error))


def overdensity_manual(f_s, f_c, n, **kwargs):
    f_c_error = np.sqrt(f_s * (1 - f_s) / n)
    valid = f_c_error > 0
    f_s = f_s[valid]
    f_c = f_c[valid]
    f_c_error = f_c_error[valid]
    error_tolerance = 1e-5
    chi_square_increase_error = np.sqrt(len(f_c))
    # This is assuming golden section search will work - assume only one minimum exists

    def calc_chi_square(a):
        return np.sum(((f_c - (((1 + a) * f_s) / (1 + a * f_s))) / f_c_error) ** 2)

    x1, x2, x3 = -10, 0, 10
    f1, f2, f3 = calc_chi_square(x1), calc_chi_square(x2), calc_chi_square(x3)
    while True:
        x4, x5 = (x1 + x2) / 2, (x2 + x3) / 2
        f4, f5 = calc_chi_square(x4), calc_chi_square(x5)
        if f4 < f2:
            x2, f2 = x4, f4
            x3, f3 = x2, f2
        elif f5 < f2:
            x1, f1 = x2, f2
            x2, f2 = x5, f5
        else:
            x1, f1 = x4, f4
            x3, f3 = x5, f5
        if x3 - x1 < error_tolerance:
            break
    alpha = x2
    chi_square_min = f2
    chi_square_target = chi_square_min + chi_square_increase_error
    max_error = 1
    lower1 = alpha - max_error
    upper1 = alpha
    lower2 = alpha
    upper2 = alpha + max_error
    x = np.linspace(-1, 1, 1000)
    y = []
    for data_point in x:
        y.append(calc_chi_square(data_point))
    plt.plot(x, y)
    plt.plot((-1, 1), (chi_square_target, chi_square_target), linestyle='--')
    plt.show()
    #(fl1, fl2, fu1, fu2) = (calc_chi_square(lower1), calc_chi_square(lower2), calc_chi_square(upper1),
    #                        calc_chi_square(upper2))
    #while True:
    #    m1 = (lower1 + upper1) / 2
    #    m2 = (lower2 + upper2) / 2
    #    fm1, fm2 = calc_chi_square(m1), calc_chi_square(m2)
    #    if fm1 < chi_square_target:
    #        upper1, fu1 = m1, fm1
    #    else:
    #        lower1, fl1 = m1, fm1

    #    if fm2 < chi_square_target:
    #        upper2, fu2 = m2, fm2
    #    else:
    #        lower2, fl2 = m2, fm2

    #    if (upper1 - lower1) < error_tolerance and (upper2 - lower2) < error_tolerance:
    #        break

    i = 0
    while True:
        i += 1
        midpoint1 = (lower1 + upper1) / 2
        midpoint2 = (lower2 + upper2) / 2
        chi_square1 = calc_chi_square(midpoint1)
        chi_square2 = calc_chi_square(midpoint2)
        if chi_square1 < chi_square_target:
            upper1 = midpoint1
        else:
            lower1 = midpoint1
        if chi_square2 < chi_square_target:
            upper2 = midpoint2
        else:
            lower2 = midpoint2
        #if np.abs(upper1 - lower1) + np.abs(upper2 - lower2) < error_tolerance:
        if i > 10:
            break
    print(alpha, midpoint1, midpoint2)
    return np.array((alpha, np.abs(midpoint1 - midpoint2) / 2))
