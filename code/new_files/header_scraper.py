from toolkit import toolkit
import numpy as np
# point - cosmology intersection
# galactic - MMF3 cosmology
# survey - MMF3 survey
import astropy
import healpy as hp

#toolkit.load_catalogue("sdss")

mask = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1)

#mask = astropy.io.fits.open("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1)
#mask = astropy.io.fits.open("../../data/ACT_mask.fits", hdu=1)
#mask = astropy.table.Table.read("../../data/ACT_mask.fits")
#mask = astropy.table.Table.read("../binned_results/test.fits")
#mask = astropy.table.Table.read("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1)
#toolkit.get_header_info(mask)
#print(mask.info)
#data = mask.field("RA---CAR")

#print(np.min(data), np.max(data))

#mask = astropy.io.fits.open("../binned_results/test.fits")

#cat = toolkit.StarCatalogue("../binned_results/test.fits", hdu=1)
#cat.load_lon_lat()

#print(np.shape(cat.lon_lat))
