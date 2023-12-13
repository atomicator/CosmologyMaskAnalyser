import matplotlib.pyplot as plt
import numpy as np
from toolkit import toolkit

mask = toolkit.load_mask("planck_galactic")
cat = toolkit.load_catalogue("sdss")
cat.load_lon_lat()

comparison_mask = toolkit.load_mask("comparison_sdss_planck_galactic")

"""
fig = plt.figure()
ax = fig.add_subplot(111)
mask.set_fig_ax(fig, ax)
mask.plot(cbar=True, show=True)
plt.show()
"""
f = np.sum(mask.lookup_point(*cat.lon_lat.transpose())) / len(cat.lon_lat)
print(f"actual fraction of sample: {f, 1 - f}")

NSIDES = [1] # , 2, 4, 8, 16, 32, 64, 128, 256]

for n in NSIDES:
    bin = toolkit.HealpixBinMap(n)
    bin.bin_catalogue(cat)
    bin.load_catalogue(cat)
    bin.calc_masked_fraction(mask)
