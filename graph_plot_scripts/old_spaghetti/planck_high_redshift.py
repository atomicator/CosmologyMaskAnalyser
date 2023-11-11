import numpy as np
import matplotlib.pyplot as plt
import pixell
import healpy as hp
import astropy.cosmology
import astropy.units
import pandas as pd

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})


# DIST_OPTION (the one considered optimal) is variable 78
wmap_map_I = hp.read_cl("../../data/raw_planck_data/HFI_PCCS_HZ_R2.00.fits")[35:59:4]
data = wmap_map_I.astype(np.float_).transpose()

averages = []

# Average the data
'''
a = len(data[0])
print(a)
print(data)
for array in data:
    averages += [sum(array) / a]

# Bin the data
num_bins = 30
min_bin = 0
max_bin = 6

bin_width = (max_bin - min_bin) / num_bins
x = np.linspace(min_bin, max_bin, num_bins + 1) #+ (bin_width / 2)
bins = np.zeros(num_bins)
n = len(averages)

for point in averages:
    bins[int((point - min_bin) / bin_width)] += 1

print(bins)
'''

plot_data = []

# Bin the data
num_bins = 30
min_bin = 0
max_bin = 8

bin_width = (max_bin - min_bin) / num_bins
x = np.linspace(min_bin, max_bin, num_bins + 1) #+ (bin_width / 2)
bins = np.zeros(num_bins)

for data_set in data.transpose():
    for data_point in data_set:
        bins[int((data_point - min_bin) / bin_width)] += 1
    print("test")
    plot_data.append(list(bins))
    bins = np.zeros(num_bins)

plt.title("Distribution of Redshifts")
plt.xlabel(r"$z$")
plt.ylabel("Number of Stars")

print(len(plot_data))

for data_set in plot_data[::-1]:
    plt.stairs(data_set, x, fill=False)
    print("test")

#plt.stairs(bins * 2, x, fill=True)
#plt.stairs(bins, x, fill=True)
plt.xlim(min_bin, max_bin)

plt.savefig("./graphs/redshift_histogram.pdf", format="pdf")

plt.show()
