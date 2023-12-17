from toolkit import toolkit
import numpy as np
import matplotlib.pyplot as plt

mask1 = toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_1.fits")
mask2 = toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_2.fits")
mask3 = toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_3.fits")
mask4 = toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_4.fits")

data = mask1.map + mask2.map + mask3.map + mask4.map

fig = plt.figure()
ax = fig.add_subplot(111)

mask1.set_fig_ax(fig, ax)
mask1.plot(show=True, cbar=True, cmap="rainbow")

# [masked by both , masked by planck only, masked by SDSS only, Allowed by both]

mask = toolkit.load_mask("comparison_sdss_planck_galactic")

print(mask.map)
