import numpy
import numpy as np
import matplotlib.pyplot as plt
import pixell
import pixell.enmap as enmap
import pixell.enplot
import healpy as hp
import astropy.cosmology
import astropy.units
import pandas as pd
import matplotlib.projections

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})

imap = enmap.read_fits("../../data/ACT_mask.fits", hdu=1)
new_data = numpy.array(imap)

masked_area = 0
total_area = 0

lat = len(new_data)
lon = len(new_data[0])

print(lat, lon)

for i in range(lat):
    s = np.sin(np.pi * (i + 0.5) / lat)
    total_area += s * lon
    masked_area += s * (lon - np.sum(new_data[i]))

masked_fraction = masked_area / total_area
print(f"Masked fraction: {masked_fraction}")

print(np.sum(new_data) / (lat * lon))

# Generate a random set of cartesian co-ords and check if masked
N = int(1e6)

masked_points = 0
total_points = 0
n = 0

while n < N:
    # Constrain so that the points must lie within a unit sphere - this would skew the reaults nif it wasn't done
    x = - 1 + 2 * np.random.random()
    y = - 1 + 2 * np.random.random()
    z = - 1 + 2 * np.random.random()
    if (x ** 2 + y ** 2 + z ** 2) <= 1:
        theta = np.arctan2(np.sqrt(x ** 2 + y ** 2), z)
        phi = np.arctan2(y, x)
        #print(theta / np.pi, phi/np.pi)
        masked_points += 1 - imap[int((theta / np.pi) * lat)][int(((phi / (2 * np.pi)) + 1/2) * lon)]
        total_points += 1
        n += 1

simulation_masked_fraction = masked_points / total_points

print(f"Fraction of {N} points that where masked: {simulation_masked_fraction}")
print(f"Difference: {np.abs(simulation_masked_fraction - masked_fraction)}")
# Assuming a binomial distribution on the theoretical result
print(f"Predicted deviation: {np.sqrt(N * masked_fraction * (1 - masked_fraction)) / N}")
