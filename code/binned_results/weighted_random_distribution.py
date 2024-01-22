import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from toolkit import toolkit
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--mask", default="sdss_mask")
parser.add_argument("--sample_catalogue", default="sdss")
parser.add_argument("--nside", type=int, default=8)
parser.add_argument("--save_path")
parser.add_argument("--target", type=int)
parser.add_argument("--raise_path", type=int, default=2)

args = parser.parse_args()

mask = toolkit.load_mask(args.mask, raise_dir=args.raise_path)
sample_cat = toolkit.load_catalogue(args.sample_catalogue, raise_dir=args.raise_path)
sample_cat.load_lon_lat()

if args.target:
    target = args.target
else:
    target = len(sample_cat.lon_lat) * 10

nside = args.nside

points = toolkit.gen_random_coords(target, mask)[::-1].transpose()
print(points)
random_cat = toolkit.StarCatalogue()
random_cat.lon_lat = points
data = toolkit.match_distribution(random_cat, sample_cat, nside)

print(data)
print(np.shape(data))

"""col1 = fits.Column(name="GLON", format="D")
col2 = fits.Column(name="GLAT", format="D")
hdu = fits.BinTableHDU.from_columns([col1, col2])({args.catalogue})
hdu.data = np.float64(data.transpose())
print(hdu.columns)
hdu.writeto("test.fits")"""
a = Table(data, names=["GLON", "GLAT"], dtype=[data.dtype, data.dtype])

if args.save_path:
    a.write(args.save_path, overwrite=True)
else:
    a.write(f"random_sdss_nside_{nside}.fits", overwrite=True)