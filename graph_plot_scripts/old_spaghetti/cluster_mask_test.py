import numpy as np
import matplotlib.pyplot as plt
import pixell
import healpy as hp
import astropy
import pandas as pd
import toolkit

wmap_map_I, header = hp.read_map("../../data/planck_galactic_mask.fits", h=True, nest=False)
print(header)
new_header = {}
for array in header:
    new_header[array[0]] = array[1]
NAXIS1 = new_header["NAXIS1"]
NAXIS2 = new_header["NAXIS2"]
NPIX = new_header["NPIX"]
NSIDE = new_header["NSIDE"]

# Pixel lookup test
i = hp.vec2pix(NSIDE, 1, 1, 1, nest=False)
print(i)
print(wmap_map_I[i])

# From healpix website: pixels cover the same surface area, meaning:

n = sum(wmap_map_I)
fraction_masked = (NPIX - n) / NPIX

print(f"Fraction masked: {fraction_masked}")

# Generate a random set of cartesian co-ords and check if masked
N = int(1e6)

masked_points = 0
total_points = 0
n = 0

while n < N:
    '''# Constrain so that the points must lie within a unit sphere - this would skew the results if it wasn't done
    x = - 1 + 2 * np.random.random()
    y = - 1 + 2 * np.random.random()
    z = - 1 + 2 * np.random.random()
    if (x ** 2 + y ** 2 + z ** 2) <= 1:
        pix = hp.vec2pix(NSIDE, x, y, z, nest=False)
        masked_points += 1 - wmap_map_I[pix]
        total_points += 1
        n += 1'''
    ra, dec = toolkit.gen_random_coord()
    pix = hp.ang2pix(NSIDE, dec, ra, nest=False)
    masked_points += 1 - wmap_map_I[pix]
    total_points += 1
    n += 1

simulation_masked_fraction = masked_points / total_points
print(f"Number of points:           {N}")
print(f"Fraction that where masked: {simulation_masked_fraction}")
print(f"Difference:                 {np.abs(simulation_masked_fraction - fraction_masked)}")
# Assuming a binomial distribution on the theoretical result
print(f"Predicted deviation:        {np.sqrt(N * fraction_masked * (1 - fraction_masked)) / N}")
