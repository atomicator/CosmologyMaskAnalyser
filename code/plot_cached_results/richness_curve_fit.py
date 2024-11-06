from toolkit import toolkit
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

toolkit.plt_use_tex()


def exponential(x, alpha, beta):
    return alpha * np.exp(x * beta)


def linear(x, alpha, beta):
    return alpha + (x - 30) * beta


def power_law(x, alpha, beta):
    return alpha * (1 + x) ** beta


files = ["./act_r/act_r_5_10.npy", "./act_r/act_r_10_20.npy", "./act_r/act_r_20_40.npy", "./act_r/act_r_40_80.npy",
         "./act_r/act_r_80+.npy"]

NSIDES = [1 / 2, 1, 2, 4, 8, 16, 32, 64, 128]
NSIDE_used = 8

richness = []
min_richness = 0

y = np.zeros(len(files))
y_std = np.zeros(len(files))

for i in range(len(files)):
    data = np.load(files[i])[0]
    y[i] = data[0][NSIDES.index(NSIDE_used)]
    y_std[i] = data[1][NSIDES.index(NSIDE_used)]


def load_func(r):
    if r > min_richness:
        richness.append(r)


cat = toolkit.load_catalogue("sdss")
cat.load_with_selection(load_func, ("LAMBDA_CHISQ",), lon_lat=True)
richness = np.array(richness)

x = np.array((np.mean(richness[np.bitwise_and(richness > 5, richness < 10)]),
              np.mean(richness[np.bitwise_and(richness > 10, richness < 20)]),
              np.mean(richness[np.bitwise_and(richness > 20, richness < 40)]),
              np.mean(richness[np.bitwise_and(richness > 40, richness < 80)]),
              np.mean(richness[richness > 80])))

print(x)
func = linear

popt, pcov = scipy.optimize.curve_fit(func, x, y, sigma=y_std, absolute_sigma=True)
#popt, pcov = scipy.optimize.curve_fit(func1, np.log(1+x), np.log(1 + y), sigma=0.01*np.ones(len(x)), absolute_sigma=True)
popt = np.array(popt)
print(np.sqrt(np.diag(pcov)))
perr = np.array(np.sqrt(np.diag(pcov)))# / (np.sqrt(2)) # Divide due to combining two equivalent errors
#perr[1] = - perr[1]
#perr[1] = 0
print(popt, perr)
print(f"chi-squared: {np.sum(np.square((y - func(x, *popt)) / y_std)) / (len(x) - len(popt))}")
plt.errorbar(x, y, y_std, color="xkcd:electric blue", ecolor="xkcd:aqua blue", capsize=3, capthick=1)
plt.plot((x[0], x[-1]), (0, 0), color="k")
base = np.linspace(x[0], x[-1], 1000)
plt.plot(base, func(base, *popt), color="red")
plt.plot(base, func(base, *(popt + perr)), color="orange")
plt.plot(base, func(base, *(popt - perr)), color="orange")
plt.xscale("log")
#plt.yscale("log")
#plt.ylim(-0.1, 0.5)
plt.title("The excess of the ACT point-source mask as a function of richness")
plt.ylabel("Excess")
plt.xlabel("Richness")
plt.savefig("richness_excess_fit.pdf")
plt.show()
