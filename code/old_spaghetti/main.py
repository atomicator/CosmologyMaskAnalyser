import numpy as np
import matplotlib.pyplot as plt
import pixell
import healpy as hp
import astropy
import pandas as pd

NSIDE = 32
print(
    "Approximate resolution at NSIDE {} is {:.2} deg".format(
        NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
    )
)

NPIX = hp.nside2npix(NSIDE)
print(NPIX)

wmap_map_I = hp.read_map("../../data/planck_galactic_mask.fits")

hp.mollview(
    wmap_map_I,
    xsize=int(1e4),
    #coord=["G", "E"],
    #title="Histogram equalized Ecliptic",
    #unit="mK",
    #norm="hist"
)
hp.graticule()

plt.savefig("./graphs/mask_plot.png", format="png", dpi=1.5e3)

plt.show()
