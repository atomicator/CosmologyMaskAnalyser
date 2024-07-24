import numpy as np
import matplotlib.pyplot as plt
from toolkit import toolkit, toolkit_new
import gc
import healpy as hp
import pixell.enmap

toolkit.plt_use_tex()

print(*())

data = toolkit.plt_use_tex

if data:
    print(True)
else:
    print(False)

#planck_mask = toolkit.load_mask("planck_modified_total")
#planck_mask = toolkit_new.PixellMask("../../data/ACT_mask.fits", mask_using_latlon=False, hdu=1)
#print(planck_mask.calc_exact_unmasked_fraction())
#sdss_mask = toolkit.load_mask("sdss_mask")
#print(planck_mask.calc_exact_unmasked_fraction())
#planck_mask = toolkit_new.HealpyMask("../../data/HFI_PCCS_SZ-selfunc-MMF3-cosmolog_R2.08.fits", hdu=1)
#planck_mask = toolkit_new.load_mask("planck_galactic")
#imap = pixell.enmap.read_map("../../data/ACT_mask.fits", hdu=1)
#print(imap.shape)
#print(imap.wcs)
#cat = toolkit.ClusterCatalogue("../../data/sdss_catalogue.fits", hdu=1)
"""
rotator = hp.Rotator(coord=["G", "C"])
#sdss_mask.map = rotator.rotate_map_pixel(sdss_mask.map)

#print(len(cat.data))

resolution = (1000, 2000)

fig = plt.figure()
ax = fig.add_subplot(111)

#cat.set_fig_ax(fig, ax)
#cat.set_cat_mask(sdss_mask)
#cat.load_lon_lat()
#cat.load_ra_dec()
#cat.plot_heatmap(128, cmap="rainbow", resolution=resolution, cbar=True, cbar_label="Clusters per "
#                                                                                     "square degree", clear=False,
#                 show=False)
gc.collect()

planck_mask.set_fig_ax(fig, ax)
planck_mask.plot(show=False, cbar=False, cmap=toolkit.bw_heatmap, clear=False, alpha=0.45)

plt.xlabel("Longitude")
plt.ylabel("Latitude")

#plt.xlabel("Right Ascension")
#plt.ylabel("Declination")

plt.title(r"The SDSS clusters and the point mask")

#plt.savefig("../../graphs/PlanckMaskSDSSData.png")
#plt.savefig("../../graphs/ACTMaskSDSSData.png")

plt.show()
"""