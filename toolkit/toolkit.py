import numpy as np
import matplotlib.pyplot as plt
import pixell
import pixell.enmap
import healpy as hp
import astropy.cosmology
import astropy
import matplotlib


def plt_use_tex():
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif"})
    return None


class __Mask(object):  # overwite this with pixell and healpy specific versions
    # Do not use this class outside of this file
    def lookup_point(self, lon, lat):
        pass

    def set_fig_ax(self, fig, ax):
        self.fig = fig
        self.ax = ax

    def gen_random_point(self):
        x = self.lookup_point(*gen_random_coord())
        return x

    def calc_exact_unmasked_fraction(self):
        pass

    def simulate_unmasked_fraction(self, num_points):
        unmasked_points = 0
        for i in range(int(num_points)):
            unmasked_points += self.gen_random_point()
        return unmasked_points / num_points


class PixellMask(__Mask):
    def __init__(self, path, step=1, mask_using_latlon=True, **kwargs):
        # Load data from fits file
        self.map = pixell.enmap.read_fits(path, **kwargs)[::step, ::step]
        self.x = len(self.map)
        self.y = len(self.map[0])
        self.lat_min = -90
        self.lat_max = 90
        self.lon_min = 0
        self.lon_max = 360
        self.fig = None
        self.ax = None
        self.mask_using_latlon = mask_using_latlon

    def update_angle_range(self, lat_min=-90, lat_max=90, lon_min=0, lon_max=360):
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.lon_min = lon_min
        self.lon_max = lon_max

    def lookup_point(self, lat, lon):
        # print(lon, lat)
        point = self.map[int(((90 + lat) / 180) * self.x)][int((lon / 360) * self.y)]
        return point

    def plot_on_ax(self, ax, alpha, cmap="plasma"):
        lat = np.linspace(self.lat_min, self.lat_max, self.x + 1)
        lon = np.linspace(self.lon_min, self.lon_max, self.y + 1)
        lon, lat = np.meshgrid(lon, lat)
        c = ax.pcolormesh(lon, lat, self.map, alpha=alpha, cmap=cmap)
        return c

    def old_plot(self, fig, ax, save_path=None, cmap="plasma", title=None, show=True, cbar=True, clear=True, **kwargs):
        if clear:
            plt.clf()
        lat = np.linspace(self.lat_min, self.lat_max, self.x + 1)
        lon = np.linspace(self.lon_min, self.lon_max, self.y + 1)
        lon, lat = np.meshgrid(lon, lat)
        c = ax.pcolormesh(lon, lat, self.map, cmap=cmap, **kwargs)
        print("Test")
        if cbar:
            fig.colorbar(c, orientation="horizontal", cmap=cmap)
        if title:
            plt.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            fig.show()

    def plot_old(self, ax, save_path=None, cmap="plasma", title=None, show=True, cbar=True, clear=True, **kwargs):
        if clear:
            plt.clf()
        lat = np.linspace(self.lat_min, self.lat_max, self.x + 1)
        lon = np.linspace(self.lon_min, self.lon_max, self.y + 1)
        lon, lat = np.meshgrid(lon, lat)
        c = self.ax.pcolormesh(lon, lat, self.map, cmap=cmap, **kwargs)
        print("Test")
        if cbar:
            self.fig.colorbar(c, orientation="horizontal", cmap=cmap)
        if title:
            ax.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        return c

    def plot(self, clear=True, cbar=None, title=None, save_path=None, show=False, cmap="plasma", alpha=1):
        if clear:
            self.ax.cla()
        c = self.plot_on_ax(self.ax, alpha, cmap=cmap)
        if cbar:
            self.fig.colorbar(c, orientation="horizontal", cmap=cmap)
        if title:
            self.ax.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()

    def calc_exact_unmasked_fraction(self):
        total_area = 0
        unmasked_area = 0
        for i in range(self.x):
            s = np.sin(np.pi * (i + 0.5) / self.x)
            total_area += s * self.y
            unmasked_area += s * np.sum(self.map[i])
        return unmasked_area / total_area


class HealpyMask(__Mask):
    def __init__(self, path, nest=False, mask_using_latlon=True, **kwargs):  # Load data from fits file
        self.map, header = hp.read_map(path, h=True, nest=nest, **kwargs)
        self.header = dict(header)
        self.NSIDE = int(self.header["NSIDE"])
        self.NPIX = 12 * self.NSIDE ** 2
        self.nest = nest
        self.mask_using_latlon = mask_using_latlon

    def lookup_point(self, lon, lat, correction_applied=False):
        if self.mask_using_latlon or correction_applied:
            pix = hp.ang2pix(self.NSIDE, lon, lat, lonlat=True, nest=self.nest)
        else:
            c = astropy.coordinates.SkyCoord(lon, lat, unit="deg", frame="galactic")  # convert to galactic co ords
            ra = c.icrs.ra.degree
            dec = c.icrs.dec.degree
            # print(ra, dec)
            pix = hp.ang2pix(self.NSIDE, ra, dec, lonlat=True, nest=self.nest)
        return self.map[pix]

    def plot_quick(self, save_path=None, xsize=1e4, title="", graticule=False, clear=True, show=True):
        if clear:
            plt.clf()
        hp.mollview(
            self.map,
            xsize=xsize,
            title=title
        )
        if graticule:
            hp.graticule()
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()

    def calc_exact_unmasked_fraction(self):
        unmasked_points = np.sum(self.map)
        unmasked_fraction = unmasked_points / self.NPIX
        return unmasked_fraction

    def plot_on_ax(self, alpha, resolution, cmap="plasma"):
        # x is lat, y is lon
        print("Converting units")
        x, y = resolution
        print(resolution)
        pixel_data = np.zeros((x, y))
        lon = np.linspace(-90, 90, x)
        lat = np.linspace(0, 360, y)
        lat, lon = np.meshgrid(lat, lon)
        c = astropy.coordinates.SkyCoord(lat, lon, unit="deg", frame="galactic")  # convert to galactic co ords
        ra = c.icrs.ra.degree
        dec = c.icrs.dec.degree
        print("converted units")
        for i in range(x):
            for j in range(y):
                # pixel_data[i, j] = self.map[hp.ang2pix(self.NSIDE, (j / y) * 360, (i / x) * 180 - 90, lonlat=True)]
                # print(ra[j], dec[i])
                pixel_data[i, j] = self.lookup_point(ra[i, j], dec[i, j], correction_applied=True)
        lon = np.linspace(0, 360, y + 1)
        lat = np.linspace(- 90, 90, x + 1)
        lon, lat = np.meshgrid(lon, lat)
        print(x, y)
        print(pixel_data)
        c = self.ax.pcolormesh(lon, lat, pixel_data, alpha=alpha, cmap=cmap)
        return c

    def plot(self, resolution=(1000, 2000), clear=True, cbar=None, title=None, save_path=None, show=False,
             cmap="plasma", alpha=1):
        if clear:
            self.ax.cla()
        c = self.plot_on_ax(alpha, resolution, cmap=cmap)
        if cbar:
            self.fig.colorbar(c, orientation="horizontal", cmap=cmap)
        if title:
            self.ax.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()

    def compare(self, mask):
        data = self.map + 1j * mask.map
        return data


class NotMaskError(TypeError):
    """Raise for errors relating to calling mask related functions on a star catalogue"""


class NoFigAx(TypeError):
    """Raise if fig and ax are not defined for a plot"""


class StarCatalogue(object):
    def __init__(self, path, hdu=0, **kwargs):
        self.fits = astropy.io.fits.open(path, **kwargs)
        self.h = dict(self.fits[hdu].header)
        self.data = self.fits[hdu].data

        self.lon_lat = None
        self.heatmap = None
        self.plot_data = None
        self.NSIDE = None
        self.fig = None
        self.ax = None
        self.cat_mask = None

    def set_cat_mask(self, mask):
        self.cat_mask = mask

    def set_fig_ax(self, fig, ax):
        self.fig = fig
        self.ax = ax

    def load_ra_dec(self):
        dec = self.data.field('DEC')
        ra = self.data.field('RA')
        self.lon_lat = np.array([ra, dec]).transpose()

    def load_lon_lat(self):
        try:
            lon = self.data.field('GLON')
            lat = self.data.field('GLAT')
        except KeyError:
            dec = self.data.field('DEC')
            ra = self.data.field('RA')
            c = astropy.coordinates.SkyCoord(ra, dec, unit="deg")  # convert to galactic co ords
            lon = c.galactic.l.degree
            lat = c.galactic.b.degree
        self.lon_lat = np.array([lon, lat]).transpose()

    def load_with_selection(self, selection_function, requested_fields, lon_lat, *args, **kwargs):
        if lon_lat:
            lon = self.data.field('GLON')
            lat = self.data.field('GLAT')
        else:
            dec = self.data.field('DEC')
            ra = self.data.field('RA')
            c = astropy.coordinates.SkyCoord(ra, dec, unit="deg")  # convert to galactic co ords
            lon = c.galactic.l.degree
            lat = c.galactic.b.degree

        requested_data = []
        for field in requested_fields:
            requested_data += [self.data.field(field)]
        requested_data = np.array(requested_data).transpose()

        new_data = []
        for i in range(len(requested_data)):
            if selection_function(*requested_data[i], *args, **kwargs):
                new_data += [[lon[i], lat[i]]]
        self.lon_lat = np.array(new_data)

    def update_lon_lat(self, new_lon_lat):
        self.lon_lat = new_lon_lat

    def plot_heatmap_quick(self, nside=128, show=True, **kwargs):
        self.gen_heatmap_data(nside)
        hp.mollview(self.heatmap, **kwargs)
        if show:
            plt.show()

    def plot_on_ax(self, ax, alpha, resolution, cmap="plasma", vmax=None):
        x, y = resolution
        print(np.sum(self.heatmap))
        pixel_data = np.zeros((x, y))
        for i in range(x):
            print(i)
            for j in range(y):
                pixel_data[i, j] = self.heatmap[hp.ang2pix(self.NSIDE, (j / y) * 360, (i / x) * 180 - 90, lonlat=True)]
        print(f"max: {np.max(pixel_data)}")
        lat = np.linspace(- 90, 90, x + 1)
        lon = np.linspace(0, 360, y + 1)
        lon, lat = np.meshgrid(lon, lat)
        # TODO: Make NaN filter mask based
        if not self.cat_mask:
            pixel_data[pixel_data == 0] = np.nan
        else:
            for i in range(x):
                print(i)
                for j in range(y):
                    if self.cat_mask.lookup_point((j / y) * 360, (i / x) * 180 - 90) == 0:
                        pixel_data[i, j] = np.nan
        # c = ax.pcolormesh(lon, lat, pixel_data, alpha=alpha, cmap=cmap,
        #                  norm=matplotlib.colors.LogNorm(vmin=None, vmax=None))
        c = ax.pcolormesh(lon, lat, pixel_data, vmin=1, alpha=alpha, cmap=cmap)
        return c

    def gen_heatmap_data(self, nside=2048):
        self.NSIDE = nside
        self.heatmap = np.zeros(hp.nside2npix(self.NSIDE))
        for cluster in self.lon_lat:
            self.heatmap[hp.ang2pix(self.NSIDE, *cluster[::1], lonlat=True)] += 1
        resolution = hp.nside2resol(self.NSIDE, arcmin=True) / 60
        print(f"Resolution: {resolution} pixels")
        print(f"Pixels in heatmap: {len(self.heatmap)}")
        print(max(f"Max value in heatmap: {self.heatmap}"))
        print(f"Sum: {np.sum(self.heatmap)}")
        self.heatmap = self.heatmap / (resolution ** 2)

    def plot_heatmap(self, nside, resolution=(1000, 2000), clear=True, cbar=None, title=None, save_path=None,
                     show=False, cmap="plasma", alpha=1, cbar_label="", vmax=None):
        self.gen_heatmap_data(nside)
        if clear:
            self.ax.cla()
        c = self.plot_on_ax(self.ax, alpha, resolution, cmap=cmap, vmax=vmax)
        if cbar:
            cbar = self.fig.colorbar(c, orientation="horizontal", location="bottom", cmap=cmap, label=cbar_label)
        if title:
            self.ax.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()

    def plot_scatter(self, **kwargs):
        self.ax.scatter(*self.lon_lat.transpose(), **kwargs)


class HealpyPlot(object):
    def __init__(self, healpy_object):
        self.object = healpy_object
        self.heatmap = None
        self.NSIDE = None

    def gen_heatmap_from_catalogue(self, nside=128):
        # Create a blank healpix plot and add the data points to it
        self.heatmap = np.zeros(hp.nside2npix(self.NSIDE))
        self.NSIDE = nside
        for cluster in self.object.lon_lat:
            self.heatmap[hp.ang2pix(self.NSIDE, *cluster[::-1], lonlat=True)] += 1
        # Scale based on resolution (normalise to value per square degree)
        resolution = hp.nside2resol(self.NSIDE, arcmin=True) / 60
        print(resolution)
        self.heatmap = self.heatmap / (resolution ** 2)

    def gen_heatmap_from_mask(self):
        self.heatmap = self.object.map
        self.NSIDE = self.object.NSIDE

    def plot_on_ax(self, ax, alpha, resolution, cmap="plasma"):
        x, y = resolution
        pixel_data = np.zeros((x, y))
        for i in range(x):
            for j in range(y):
                pixel_data[i, j] = self.heatmap[hp.ang2pix(self.NSIDE, (i / x) * 180 - 90, (j / y) * 2 * 360,
                                                           lonlat=True)]
        lat = np.linspace(- 90, 90, x + 1)
        lon = np.linspace(0, 360, y + 1)
        lon, lat = np.meshgrid(lon, lat)
        c = ax.pcolormesh(lon, lat, pixel_data, alpha=alpha, cmap=cmap)
        return c


def gen_random_coord():
    finished = False
    (lon, lat) = (0, 0)
    while not finished:
        x = - 1 + 2 * np.random.random()
        y = - 1 + 2 * np.random.random()
        z = - 1 + 2 * np.random.random()
        if (x ** 2 + y ** 2 + z ** 2) <= 1:
            finished = True
            lat, lon = hp.vec2ang(np.array((x, y, z)), lonlat=True)
    return lon, lat


def get_file_info(path, **kwargs):
    fits_image = astropy.io.fits.open(path, **kwargs)
    fits_image.info()
    return fits_image


def get_header_info(hdu_list):
    for hdu in hdu_list:
        print(hdu.header)


def load_mask(mask, raise_dir=2):
    value = None
    if mask == "planck_galactic":
        value = HealpyMask("../" * raise_dir + "data/planck_galactic_mask.fits")
    if mask == "planck_point":
        value = HealpyMask("../" * raise_dir + "data/planck_point_mask.fits")
    if mask == "sdss_mask":
        value = HealpyMask("../" * raise_dir + "data/redmapper_dr8_public_v6.3_zmask.fits", mask_using_latlon=False,
                           hdu=1, partial=True)
        value.map[value.map > 0.4] = 1.0
        value.map[value.map < 0.3] = 0
        # value.map = (value.map - 1) * -1
    return value


def bootstrap(data, samples):
    n = len(data)
    mean_estimates = []
    for i in range(samples):
        print(f"Bootstrap iteration: {i}")
        new_data = np.random.choice(data, n)
        estimate = np.sum(new_data)
        print(estimate)
        mean_estimates.append(estimate)
    return mean_estimates


def fraction_masked_pair(mask1, mask2, n=int(1e3)):
    if isinstance((mask1, mask2), HealpyMask) and mask2.NPIX == mask1.NPIX and mask1.using_latlon and mask2.using_latlon:
        data = mask1.compare(mask2)
        k = np.max((mask1.NPIX, mask2.NPIX))
        results = [np.real(np.sum(data) / k), np.imag(np.sum(data) / k), np.sum(data == 1 + 1j) / k,
                   np.sum(data == 1) / k, np.sum(data == 1j) / k,
                   np.sum(data == 0) / k]
    else:
        x = np.linspace(-np.pi, np.pi, n)
        y = np.linspace(-np.pi/2, np.pi/2, n)
        x, y = np.meshgrid(x, y)
        data1 = mask1.lookup_point(x * 180 / np.pi, y * 180 / np.pi)
        data2 = mask2.lookup_point(x * 180 / np.pi, y * 180 / np.pi)
        data = np.array(data1 + 1j * data2)
        sums = np.array([0, 0, 0, 0, 0])
        for i in range(len(y)):
            s = np.sin(np.pi * (i + 0.5) / len(y))
            sums[0] += len(x) * s
            sums[1] += np.sum(data[i] == 1) * s
            sums[2] += np.sum(data[i] == 1j) * s
            sums[3] += np.sum(data[i] == 1+1j) * s
            sums[4] += np.sum(data[i] == 0) * s
        results = sums[1:] / sums[0]
    return results


def gen_random_coords(n, mask=None):
    if mask:
        n = n / mask.calc_exact_unmasked_fraction()  # correct for masked fraction
    n = int(n * 6 / np.pi)  # correct for volume generation
    x = 2 * np.random.rand(n) - 1
    y = 2 * np.random.rand(n) - 1
    z = 2 * np.random.rand(n) - 1
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # filter out points outside of unit sphere
    x = x[r < 1]
    y = y[r < 1]
    z = z[r < 1]
    r = r[r < 1]
    theta = np.arcsin(z / r) * 180 / np.pi
    phi = np.arctan2(x, y) * 180 / np.pi
    data = np.array((theta, phi))
    # if a mask is being used, filter out masked points
    if mask:
        data = data[:, np.bool_(mask.lookup_point(phi, theta))]
    return data


#  Define a black and white heatmap
__cdict_bw = [(0, 0, 0), (1, 1, 1)]
bw_heatmap = matplotlib.colors.LinearSegmentedColormap.from_list("black and white", __cdict_bw, N=2)
rb_heatmap = matplotlib.colors.LinearSegmentedColormap.from_list("rb", [(1, 0, 0), (0, 0, 1)])
