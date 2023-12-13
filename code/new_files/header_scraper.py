from toolkit import toolkit
import numpy as np
# point - cosmology intersection
# galactic - MMF3 cosmology
# survey - MMF3 survey
import astropy
import healpy as hp

mask = astropy.io.fits.open("../../data/planck_galactic_mask.fits")

toolkit.get_header_info(mask)

mask = astropy.io.fits.open("../../data/planck_point_mask.fits")

print(mask)
