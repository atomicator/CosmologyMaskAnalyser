import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from toolkit import toolkit

mask = toolkit.load_mask("sdss_mask")
sample_cat = toolkit.load_catalogue("sdss")
sample_cat.load_lon_lat()

target = len(sample_cat.lon_lat) * 10
nside = 32

points = toolkit.gen_random_coords(target, mask)[::-1].transpose()
print(points)
random_cat = toolkit.StarCatalogue()
random_cat.lon_lat = points
data = toolkit.match_distribution(random_cat, sample_cat, nside)

print(data)
print(np.shape(data))

"""col1 = fits.Column(name="GLON", format="D")
col2 = fits.Column(name="GLAT", format="D")
hdu = fits.BinTableHDU.from_columns([col1, col2])
hdu.data = np.float64(data.transpose())
print(hdu.columns)
hdu.writeto("test.fits")"""
a = Table(data, names=["GLON", "GLAT"], dtype=[data.dtype, data.dtype])
a.write(f"random_sdss_nside_{nside}.fits", overwrite=True)
