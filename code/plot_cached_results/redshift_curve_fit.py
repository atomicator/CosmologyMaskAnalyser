from toolkit import toolkit
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

toolkit.plt_use_tex()


def exponential(x, alpha, beta):
    return alpha * np.exp(x * beta)


def linear(x, alpha, beta):
    return alpha + x * beta


def power_law(x, alpha, beta):
    return 10 ** alpha * (1 + x) ** beta

#files = ["./act_z_fine/act_z_1.npy", "./act_z_fine/act_z_2.npy", "./act_z_fine/act_z_3.npy", "./act_z_fine/act_z_4.npy",
#         "./act_z_fine/act_z_5.npy", "./act_z_fine/act_z_6.npy", "./act_z_fine/act_z_7.npy", "./act_z_fine/act_z_8.npy",
#         "./act_z_fine/act_z_9.npy", "./act_z_fine/act_z_10.npy", "./act_z_fine/act_z_11.npy", "./act_z_fine/act_z_12.npy"]
files = ["./act_z_fine/act_z_0.npy", "./act_z_fine/act_z_7.npy", "./act_z_fine/act_z_8.npy",
         "./act_z_fine/act_z_9.npy", "./act_z_fine/act_z_10.npy", "./act_z_fine/act_z_11.npy", "./act_z_fine/act_z_12.npy"]

NSIDES = [1 / 2, 1, 2, 4, 8, 16, 32, 64, 128]
NSIDE_used = 8

redshift = []
min_richness = 20

y = np.zeros(len(files))
y_std = np.zeros(len(files))

for i in range(len(files)):
    data = np.load(files[i])[0]
    y[i] = data[0][NSIDES.index(NSIDE_used)]
    y_std[i] = data[1][NSIDES.index(NSIDE_used)]


def load_func(richness, z):
    if richness > min_richness:
        redshift.append(z)


cat = toolkit.load_catalogue("sdss")
cat.load_with_selection(load_func, ("LAMBDA_CHISQ", "ZRED"), lon_lat=True)
redshift = np.array(redshift)

#x = np.array((np.mean(redshift[redshift < 0.05]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.05, redshift < 0.1)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.1, redshift < 0.15)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.15, redshift < 0.2)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.2, redshift < 0.25)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.25, redshift < 0.3)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.3, redshift < 0.35)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.35, redshift < 0.4)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.4, redshift < 0.45)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.45, redshift < 0.5)]),
#              np.mean(redshift[np.bitwise_and(redshift > 0.5, redshift < 0.55)]),
#              np.mean(redshift[redshift > 0.55])))
x = np.array((np.mean(redshift[redshift < 0.3]),
              np.mean(redshift[np.bitwise_and(redshift > 0.3, redshift < 0.35)]),
              np.mean(redshift[np.bitwise_and(redshift > 0.35, redshift < 0.4)]),
              np.mean(redshift[np.bitwise_and(redshift > 0.4, redshift < 0.45)]),
              np.mean(redshift[np.bitwise_and(redshift > 0.45, redshift < 0.5)]),
              np.mean(redshift[np.bitwise_and(redshift > 0.5, redshift < 0.55)]),
              np.mean(redshift[redshift > 0.55])))


func = power_law

popt, pcov = scipy.optimize.curve_fit(func, x, y, sigma=y_std, absolute_sigma=True)
#popt, pcov = scipy.optimize.curve_fit(func1, np.log(1+x), np.log(1 + y), sigma=0.01*np.ones(len(x)), absolute_sigma=True)
popt = np.array(popt)
perr = np.array(np.sqrt(np.diag(pcov))) / np.sqrt(2)
#perr[1] = - perr[1]
print(popt, perr)
plt.errorbar(x, y, y_std, color="xkcd:electric blue", ecolor="xkcd:aqua blue", capsize=3, capthick=1)
plt.plot((x[0], x[-1]), (0, 0), color="k")
base = np.linspace(x[0], x[-1], 1000)
plt.plot(base, func(base, *popt), color="red")
plt.plot(base, func(base, *(popt + perr)), color="orange")
plt.plot(base, func(base, *(popt - perr)), color="orange")
#plt.yscale("log")
#plt.xscale("log")
#plt.xscale("log")
#plt.ylim(-0.5, 0.5)
plt.show()

print(f"Chi squared: {np.sum((((func(x, *popt) - y) / y_std) ** 2 / (len(files) - len(popt))))}")
exit()
step = 0.1
i = step
iterations = 100
initial_min, initial_max = 0, 0.6
while i < 1:
    min, max = initial_min, initial_max
    for j in range(iterations):
        s = np.sum(redshift > (min + max) / 2) / len(redshift)
        if s > i:
            min = (min + max) / 2
        else:
            max = (min + max) / 2
    print(f"{i}: {(min + max) / 2}")
    i += step
print(redshift)
