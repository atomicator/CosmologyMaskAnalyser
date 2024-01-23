from toolkit import toolkit
import numpy as np
# point - cosmology intersection
# galactic - MMF3 cosmology
# survey - MMF3 survey
import astropy
import healpy as hp

toolkit.load_catalogue("sdss")

#mask = astropy.io.fits.open("../binned_results/test.fits", hdu=1)
#mask = astropy.io.fits.open("../../data/sdss_catalogue.fits")
mask = astropy.table.Table.read("../binned_results/test.fits")

#toolkit.get_header_info(mask)
print(mask.info)
print(mask.field("GLON"))

mask = astropy.io.fits.open("../binned_results/test.fits")

cat = toolkit.StarCatalogue("../binned_results/test.fits", hdu=1)
cat.load_lon_lat()

print(np.shape(cat.lon_lat))
