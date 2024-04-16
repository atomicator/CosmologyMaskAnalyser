from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np
import astropy
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--min_z", type=float, default=0.0)
parser.add_argument("--min_r", type=float, default=20.0)
parser.add_argument("--max_z", type=float, default=20.0)
parser.add_argument("--max_r", type=float, default=10000.0)
args = parser.parse_args()


def filter(redshift, richness):
    global args
    return (args.min_z < redshift < args.max_z) and (args.min_r < richness < args.max_r)


class SquareMaskFilter:
    def __init__(self, lat_range=(-90, 90), lon_range=(0, 360), lonlat=True):
        self.lat_range = lat_range
        self.lon_range = lon_range
        self.lonlat = lonlat

    def lookup_point(self, lon, lat):
        if not self.lonlat:
            c = astropy.coordinates.SkyCoord(lon, lat, unit="deg", frame="galactic")  # convert to galactic co ords
            lon = c.icrs.ra.degree
            lat = c.icrs.dec.degree
        lon[lon < 0] = lon[lon < 0] + 360
        return np.int_(np.bitwise_not(np.bitwise_and(np.bitwise_and(self.lon_range[0] < lon, self.lon_range[1] > lon),
                                      np.bitwise_and(self.lat_range[0] < lat, self.lat_range[1] > lat))))


toolkit.plt_use_tex()

fig = plt.figure()
ax = fig.add_subplot(111)

#galactic_mask = toolkit.load_mask("planck_galactic")
#act_mask = toolkit.load_mask("act")
#act_mask = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=True)
#mask = toolkit.CombinationMask(galactic_mask, act_mask, invert=True)

#mask.set_fig_ax(fig, ax)
#mask.plot(cbar=False, label="Galactic", cmap="bwr", show=False, clear=False)

#galactic_mask.map = 1 - galactic_mask.map
#mask = toolkit.CombinationMask(galactic_mask, act_mask, invert=True)

# Galactic mask
"""act_mask = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
act_graph_filter = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
act_graph_filter.map = np.load("../../data/act_galactic_mask_array.npy")
test = SquareMaskFilter((-10, 20), (40, 100), lonlat=False)
final_filter = toolkit.CombinationMask(act_mask, test, use_and=False, invert=False)
galactic_mask = toolkit.CombinationMask(act_graph_filter, final_filter, use_and=True, invert=False)"""
#galactic_mask = toolkit.load_mask("act_galactic")

#plot_mask = galactic_mask
#plot_mask.set_fig_ax(fig, ax)
#plot_mask.plot(cbar=False, label="Galactic", cmap="bwr", show=False, clear=False)
lon_shift = -50
# Point mask
"""galactic_mask.invert = True
point_mask = toolkit.CombinationMask(act_mask, galactic_mask, use_and=False)"""
point_mask = toolkit.load_mask("sdss_mask", lon_shift=lon_shift)

plot_mask = point_mask
plot_mask.set_fig_ax(fig, ax)
plot_mask.plot(cbar=False, cmap="bwr", show=False, clear=False)
#print(plot_mask.lon_shift)
"""
point_mask = toolkit.load_mask("act_point")
plot_mask = point_mask
plot_mask.set_fig_ax(fig, ax)
plot_mask.plot(cbar=False, label="2", cmap="bwr_r", show=False, clear=False)"""
#print(plot_mask.lon_shift)

plt.ylabel("Latitude")
plt.xlabel("Longitude")
plt.title(r"The \enquote{Galactic} and \enquote{Point} components of the ACT mask")

cat = toolkit.load_catalogue("sdss")
cat.load_with_selection(filter, ["ZRED", "LAMBDA_CHISQ"], lon_lat=True)
cat.lon_lat[:, 0] -= lon_shift
cat.lon_lat[cat.lon_lat[:, 0] < 0, 0] += 360
cat.lon_lat[cat.lon_lat[:, 0] > 360, 0] -= 360

plt.scatter(*cat.lon_lat.transpose(), s=1)

#plt.legend()
#plt.ylabel("Declination")
#plt.xlabel("Right Ascension")
#plt.title(r"The ACT mask with the fill filter")

#plt.savefig("act_mask_final_lat_lon.png")
plt.show()
