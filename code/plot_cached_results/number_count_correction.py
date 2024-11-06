import numpy as np
from toolkit import toolkit
import matplotlib.pyplot as plt

toolkit.plt_use_tex()
f_s = 0.01771962453092031 # ACT


def excess(x):
    return alpha * np.exp(beta * (x - 0.4))


def correction_old(x):
    #r = (f_s - 1 - excess(x) * (f_s - 1)) / (f_s - 1 + excess(x))
    #return f_s * (r - 1)
    f_c = (1 + excess(x)) * f_s
    r = (f_c * (1 - f_s)) / (f_s * (1 - f_c))
    return f_s * (r - 1)


def correction(x):
    return f_s * excess(x)

z = np.linspace(0, 1, 1000)

alpha = 0.13830263
beta = -8.24858781
plt.plot(z, correction(z) * 100, color="red")
alpha = 0.13830263 + 0.04718889
beta = -8.24858781 + 2.11872194
data1 = correction(z)
alpha = 0.13830263 - 0.04718889
beta = -8.24858781 + 2.11872194
data2 = correction(z)
alpha = 0.13830263 + 0.04718889
beta = -8.24858781 - 2.11872194
data3 = correction(z)
alpha = 0.13830263 - 0.04718889
beta = -8.24858781 - 2.11872194
data4 = correction(z)
plt.plot(z, np.amin((data1, data2, data3, data4), axis=0) * 100, color="orange")
plt.plot(z, np.amax((data1, data2, data3, data4), axis=0) * 100, color="orange")
plt.xlabel(r"Redshift")
plt.ylabel(r"Correction (\%)")
plt.yscale("log")
plt.title("Number count correction for the redshift model")
#plt.savefig("redshift_number_count.pdf")
#plt.show()
plt.clf()


def excess(x):
    #return alpha * np.exp(beta * x)
    return alpha + (x - 30) * beta


def log_richness_dist(log_snr, log_r):
    return (np.exp(-(1/2) * np.square((log_r - log_richness(log_snr)) / scatter)) /
            (scatter * np.sqrt(2 * np.pi)))


def log_richness(log_snr):
    return 1.34 + 0.57 * log_snr


def num_int(snr):
    log_snr = np.log10(snr)
    y = np.zeros(np.shape(snr))
    bins = np.linspace(-10, 50, steps + 1)
    for i in range(len(snr)):
        y[i] = 0
        for j in range(steps):
            log_r = (bins[j+1] + bins[j]) / 2
            y[i] += (bins[j+1] - bins[j]) * log_richness_dist(log_snr[i], log_r) * correction(np.exp(log_r * np.log(10)))
    return y


steps = 1000
snr = np.linspace(0.5, 10, 1000)
scatter = 0.17
r = 10 ** 1.34 * snr ** 0.57
alpha = 0.11569115
beta = 0.002644
#plt.plot(snr, correction(r) * 100, color="red", linestyle="solid", label="Direct approximation")
plt.plot(snr, ((alpha - 30 * beta) + beta * r * np.exp((np.log(10) * scatter) ** 2 / 2)) * f_s * 100,
         color="red", label="Scatter approximation")
#plt.plot(snr, num_int(snr) * 100, color="blue", linestyle="solid", label="Numerical integral")
#plt.plot(snr, ((alpha - 30 * beta) + beta * np.exp(scatter ** 2 / 2) * r) * f_s * 100, linestyle="dashed")
alpha = 0.11569115 + 0.01870632
beta = 0.002644 + 0.00098928
plt.plot(snr, ((alpha - 30 * beta) + beta * r * np.exp((np.log(10) * scatter) ** 2 / 2)) * f_s * 100,
         color="orange", label="Scatter approximation", linestyle="solid")
#plt.plot(snr, num_int(snr) * 100, color="blue", linestyle="dashed")
alpha = 0.11569115 - 0.01870632
beta = 0.002644 - 0.00098928
plt.plot(snr, ((alpha - 30 * beta) + beta * r * np.exp((np.log(10) * scatter) ** 2 / 2)) * f_s * 100,
         color="orange", label="Scatter approximation", linestyle="solid")
#plt.plot(snr, num_int(snr) * 100, color="blue", linestyle="dashed")

plt.xlabel(r"SNR")
plt.ylabel(r"Number Count Correction (\%)")
plt.title("The number count correction for the ACT point-source mask as a function of richness")
#plt.legend()
plt.savefig("richness_number_count.pdf")
plt.show()
plt.clf()

plt.plot(excess(z), excess(z) * f_s)
plt.plot(excess(z), correction(z))
#plt.show()
