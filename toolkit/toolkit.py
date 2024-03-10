import numpy as np
import matplotlib.pyplot as plt
import pixell
import pixell.enmap
import healpy as hp
import astropy.cosmology
import astropy
import matplotlib
import copy
import toolkit.filters as filters


def plt_use_tex():
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        'text.latex.preamble': r'\usepackage{csquotes}'
    })
    return None


class __Mask(object):  # overwite this with pixell and healpy specific versions
    # Do not use this class outside of this file
    fig = None
    ax = None

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
    def __init__(self, path, step=1, mask_using_latlon=True, invert=True, **kwargs):
        # Load data from fits file
        self.imap = pixell.enmap.read_fits(path, **kwargs)
        self.map = np.array(self.imap)
        print(np.shape(self.map))
        self.x = len(self.map)
        self.y = len(self.map[0])
        self.lat_min = -90
        self.lat_max = 90
        self.lon_min = 00
        self.lon_max = 360
        self.fig = None
        self.ax = None
        self.mask_using_latlon = mask_using_latlon
        self.invert = invert

    def update_angle_range(self, lat_min=-90, lat_max=90, lon_min=0, lon_max=360):
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.lon_min = lon_min
        self.lon_max = lon_max

    def lookup_point(self, lon, lat):
        if not self.mask_using_latlon:
            c = astropy.coordinates.SkyCoord(lon, lat, unit="deg", frame="galactic")  # convert to galactic co ords
            lon = c.icrs.ra.degree
            lat = c.icrs.dec.degree
        # point = self.map[np.int_(((90 + lat) / 180) * self.x), np.int_(((lon - 180) / 360) * self.y)]
        pix = np.int_(pixell.enmap.sky2pix(self.imap.shape, self.imap.wcs, np.array((lat, lon)) * np.pi / 180))
        point = np.zeros(np.shape(pix)[1:])
        # points_in_range = pix[
        #    np.bitwise_and(np.bitwise_and(0 < np.array(pix[0]), np.array(pix[0]) < np.shape(self.map)[0]),
        #                   np.bitwise_and(0 < np.array(pix[1]), np.array(pix[1]) < np.shape(self.map)[1]))]
        points_in_range = np.bitwise_and(np.bool_(np.int_(0 < np.array(pix[0])) * np.int_(np.array(pix[0]) < np.shape(self.map)[0] - 1)),
                                             np.bool_(np.int_(0 < np.array(pix[1])) * np.int_(np.array(pix[1]) < np.shape(self.map)[1] - 1)))
        point[points_in_range] = self.map[pix[0, points_in_range], pix[1, points_in_range]]
        #point[points_in_range] = self.map
        if self.invert:
            point = 1 - point
        return point

    def plot_on_ax(self, ax, alpha, resolution, cmap="plasma", **kwargs):
        print(resolution)
        lat = np.linspace(self.lat_min, self.lat_max, resolution[0] + 1)[1:-2]
        lon = np.linspace(self.lon_min, self.lon_max, resolution[1] + 1)[1:-2]
        lon, lat = np.meshgrid(lon, lat)
        data = self.lookup_point(lon, lat)
        data[data == 1] = np.nan
        c = ax.pcolormesh(lon, lat, data, alpha=alpha, cmap=cmap, **kwargs)
        return c

    def old_plot(self, fig, ax, save_path=None, cmap="plasma", title=None, show=True, cbar=True, clear=True, **kwargs):
        if clear:
            plt.clf()
        lat = np.linspace(self.lat_min, self.lat_max, self.x + 1)
        lon = np.linspace(self.lon_min, self.lon_max, self.y + 1)
        lon, lat = np.meshgrid(lon, lat)
        c = ax.pcolormesh(lon, lat, self.map, cmap=cmap, **kwargs)
        if cbar:
            fig.colorbar(c, orientation="horizontal", cmap=cmap)
        if title:
            plt.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            fig.show()

    def plot_old(self, ax, save_path=None, cmap="plasma", title=None, show=True, alpha=0, cbar=True, clear=True, **kwargs):
        if clear:
            plt.clf()
        lat = np.linspace(self.lat_min, self.lat_max, self.x + 1)
        lon = np.linspace(self.lon_min, self.lon_max, self.y + 1)
        lon, lat = np.meshgrid(lon, lat)
        c = self.ax.pcolormesh(lon, lat, self.map, cmap=cmap, **kwargs)
        if cbar:
            self.fig.colorbar(c, orientation="horizontal", cmap=cmap)
        if title:
            ax.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()
        return c

    def plot(self, resolution=(1000, 2000), clear=True, cbar=None, title=None, save_path=None, show=False,
             cmap="plasma", alpha=1, **kwargs):
        if clear:
            self.ax.cla()
        print(resolution)
        c = self.plot_on_ax(self.ax, alpha, resolution=resolution, cmap=cmap, **kwargs)
        if cbar:
            self.fig.colorbar(c, orientation="horizontal", cmap=cmap)
        if title:
            self.ax.set_title(title)
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
        # ra = c.icrs.ra.degree
        # dec = c.icrs.dec.degree
        print("converted units")
        """for i in range(x):
            for j in range(y):
                # pixel_data[i, j] = self.map[hp.ang2pix(self.NSIDE, (j / y) * 360, (i / x) * 180 - 90, lonlat=True)]
                # print(ra[j], dec[i])
                pixel_data[i, j] = self.lookup_point(lat[i, j], lat[i, j], correction_applied=False)
        """
        pixel_data = self.lookup_point(lat, lon, correction_applied=False)
        pixel_data[pixel_data == 1] = np.nan
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

    def clone(self):
        fig, ax = self.fig, self.ax
        new = copy.deepcopy(self)
        new.fig, new.ax = fig, ax
        return new


class CombinationMask(PixellMask):
    def __init__(self, mask1, mask2, invert=False, use_and=True):
        self.mask1 = mask1
        self.mask2 = mask2
        self.lat_min = -90
        self.lat_max = 90
        self.lon_min = 00
        self.lon_max = 360
        self.fig = None
        self.ax = None
        self.mask_using_latlon = True
        self.invert = invert
        self.use_and = use_and

    def lookup_point(self, lon, lat):
        if self.use_and:
            data = np.float_(np.bitwise_and(np.bool_(self.mask1.lookup_point(lon, lat)),
                                                           np.bool_(self.mask2.lookup_point(lon, lat))))
        else:
            data = np.float_(np.bitwise_or(np.bool_(self.mask1.lookup_point(lon, lat)),
                                            np.bool_(self.mask2.lookup_point(lon, lat))))
        #data = 1 - self.mask1.lookup_point(lon, lat) * self.mask2.lookup_point(lon, lat)
        if self.invert:
            data = 1 - data
        return data

    def calc_exact_unmasked_fraction(self):
        pixels = np.int_(np.linspace(0, hp.nside2npix(32) - 1, hp.nside2npix(32)))
        num_masked = np.sum(self.lookup_point(*hp.pix2ang(32, pixels, lonlat=True)))
        unmasked_fraction = num_masked / len(pixels)
        return unmasked_fraction


class NotMaskError(TypeError):
    """Raise for errors relating to calling mask related functions on a star catalogue"""


class NoFigAx(TypeError):
    """Raise if fig and ax are not defined for a plot"""


class StarCatalogue(object):
    def __init__(self, path=None, hdu=0, table=False, **kwargs):
        if path:
            if not table:
                self.fits = astropy.io.fits.open(path, **kwargs)
                self.h = dict(self.fits[hdu].header)
                self.data = self.fits[hdu].data
            if table:
                self.fits = astropy.table.Table.read(path)
                self.data = self.fits
                self.h = None
        else:
            self.fits = None
            self.h = None
            self.data = None

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
        try:
            dec = self.data.field('DEC')
            ra = self.data.field('RA')
        except KeyError:
            dec = self.data.field('decDeg')
            ra = self.data.field('RADeg')
        self.lon_lat = np.array([ra, dec]).transpose()

    def load_lon_lat(self):
        try:
            lon = self.data.field('GLON')
            lat = self.data.field('GLAT')
        except KeyError:
            try:
                dec = self.data.field('DEC')
                ra = self.data.field('RA')
            except KeyError:
                dec = self.data.field('decDeg')
                ra = self.data.field('RADeg')
            c = astropy.coordinates.SkyCoord(ra, dec, unit="deg")  # convert to galactic co ords
            lon = c.galactic.l.degree
            lat = c.galactic.b.degree
        self.lon_lat = np.array([lon, lat]).transpose()

    def load_with_selection(self, selection_function, requested_fields, lon_lat, *args, **kwargs):
        if lon_lat:
            try:
                lon = self.data.field('GLON')
                lat = self.data.field('GLAT')
            except KeyError:
                dec = self.data.field('DEC')
                ra = self.data.field('RA')
                c = astropy.coordinates.SkyCoord(ra, dec, unit="deg")  # convert to galactic co ords
                lon = c.galactic.l.degree
                lat = c.galactic.b.degree
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
        """x, y = resolution
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
        return c"""
        x = np.linspace(-90, 90, resolution[0])
        y = np.linspace(0, 360, resolution[1])
        x, y = np.meshgrid(x, y)
        pixel_data = self.heatmap[hp.ang2pix(self.NSIDE, y, x, lonlat=True)]
        if not self.cat_mask:
            pixel_data[pixel_data == 0] = np.nan
        else:
            pixel_data[self.cat_mask.lookup_point(y, x) == 0] = np.nan
        lat = np.linspace(-90, 90, resolution[0] + 1)
        lon = np.linspace(0, 360, resolution[1] + 1)
        lat, lon = np.meshgrid(lat, lon)
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


class _BinMap(object):
    map = []
    weighted_map = []
    mask_fraction_map = []
    NSIDE = 0
    binned_sample = [[]]

    def load_catalogue(self, catalogue):
        if self.NSIDE <= 8:
            data = self.lookup_pix(*catalogue.lon_lat.transpose())
            for i in range(len(self.map)):
                self.map[i] = np.sum(data == i)
        else:
            for cluster in catalogue.lon_lat:
                self.map[hp.ang2pix(self.NSIDE, *cluster, lonlat=True)] += 1

    def lookup_point(self, lon, lat):
        pass

    def lookup_pix(self, lon, lat):
        pass

    def lookup_weighted_point(self, lon, lat):
        pass

    def calc_weighted_map(self, mask, resolution=(1000, 1000)):
        pass

    def bin_catalogue(self, catalogue):
        self.binned_sample = []
        for map_bin in range(hp.nside2npix(self.NSIDE)):
            self.binned_sample.append([])
        for cluster in catalogue.lon_lat:
            self.binned_sample[self.lookup_pix(*cluster)].append(cluster)
        for map_bin in range(hp.nside2npix(self.NSIDE)):
            self.binned_sample[map_bin] = np.array(self.binned_sample[map_bin])

    def calc_masked_fraction(self, mask, area_masked_fraction, weight_function, const=False):
        # Divide the sample into relevant bins
        pass
        # Run a switch-case on v, to run the method for the relevant weight set.

    def divide_sample(self, mask, area_masked_fraction, const=False, filter_set="all", a=5):
        if not const:
            sample_masked_fraction = np.zeros(hp.nside2npix(self.NSIDE))
            for pixel in range(hp.nside2npix(self.NSIDE)):
                try:
                    sample_masked_fraction[pixel] = 1 - np.sum(
                        mask.lookup_point(*self.binned_sample[pixel].transpose())) / len(
                        self.binned_sample[pixel])
                except TypeError:
                    sample_masked_fraction[pixel] = np.NAN
        else:
            sample_masked_fraction = np.array(
                1 - np.sum(mask.lookup_point(*self.binned_sample[0].transpose())) / len(self.binned_sample[0]))

        sdss_allowed_region = area_masked_fraction[1] + area_masked_fraction[3]
        planck_masked_only = area_masked_fraction[1]

        nan_filter = np.bitwise_or(np.isnan(sample_masked_fraction), np.isinf(sample_masked_fraction))
        self.map = np.array(self.map)

        if filter_set == "all":
            not_masked = np.bitwise_or(np.bitwise_or(nan_filter, sample_masked_fraction == 0), planck_masked_only == 0)
            masked = np.bitwise_and(
                np.bitwise_or(sample_masked_fraction == 1, sdss_allowed_region == planck_masked_only),
                np.bitwise_not(not_masked))
            mixed_bins = np.bitwise_not(np.bitwise_or(masked, not_masked))

            """nan_filter = np.bitwise_or(np.isnan(sample_masked_fraction), np.isinf(sample_masked_fraction))
            masked = np.bitwise_and(np.bitwise_or(sample_masked_fraction == 1, sdss_allowed_region == planck_masked_only),
                                    np.bitwise_not(nan_filter))
            not_masked = np.bitwise_and(np.bitwise_or(np.bitwise_or(nan_filter, sample_masked_fraction == 0),
                                                      sdss_allowed_region == 0), np.bitwise_not(masked))
            mixed_bins = np.bitwise_not(np.bitwise_or(masked, not_masked))"""

            """nan_filter = np.bitwise_not(np.bitwise_or(np.isnan(sample_masked_fraction), np.isinf(sample_masked_fraction)))
            masked = sample_masked_fraction == 1
            not_masked = sample_masked_fraction == 0
            mixed_bins = (np.bitwise_not(np.bitwise_or(masked, np.bitwise_not(nan_filter))) * np.bitwise_not(not_masked) *
                          np.bitwise_not(np.bitwise_or(planck_masked_only == 0, planck_masked_only == sdss_allowed_region)))"""

        elif filter_set == "n_only":
            not_masked = np.bitwise_or(np.bitwise_or(nan_filter, self.map == 0), sdss_allowed_region == 0)
            masked = np.bitwise_and(not_masked, np.bitwise_not(not_masked))
            mixed_bins = np.bitwise_and(np.bitwise_not(not_masked), np.bitwise_not(masked))
        elif filter_set == "overkill":
            with np.errstate(divide="ignore", invalid="ignore"):
                not_masked = np.bitwise_or(np.bitwise_or(nan_filter, sample_masked_fraction == 0),
                                           planck_masked_only == 0)
                masked = np.bitwise_and(
                    np.bitwise_or(sample_masked_fraction == 1, sdss_allowed_region == planck_masked_only),
                    np.bitwise_not(not_masked))
                mixed_bins = np.bitwise_not(np.bitwise_or(masked, not_masked))
                masked = np.bitwise_or(masked, np.bitwise_and(mixed_bins, 1 - (
                            planck_masked_only / sdss_allowed_region) < a / self.map))
                not_masked = np.bitwise_or(not_masked, np.bitwise_and(mixed_bins, (
                            planck_masked_only / sdss_allowed_region) < a / self.map))
                mixed_bins = np.bitwise_not(np.bitwise_or(masked, not_masked))
                masked = np.bitwise_or(masked, np.bitwise_and(mixed_bins, 1 - sample_masked_fraction < a / self.map))
                not_masked = np.bitwise_or(not_masked,
                                           np.bitwise_and(mixed_bins, sample_masked_fraction < a / self.map))
                mixed_bins = np.bitwise_not(np.bitwise_or(masked, not_masked))
        elif filter_set == "excess":
            print("test")
            not_masked = np.bitwise_or(np.bitwise_or(nan_filter, self.map == 0), sdss_allowed_region == 0)
            filtered = np.bitwise_or(np.bitwise_or(sdss_allowed_region == 0, planck_masked_only == 0), sdss_allowed_region <= planck_masked_only)
            masked = np.bitwise_and(not_masked, filtered)
            mixed_bins = np.bitwise_and(np.bitwise_not(not_masked), np.bitwise_not(masked))
        else:
            raise ValueError("filter set not recognised")

        if not const:
            fully_masked_cluster_fraction = np.sum(self.map[masked]) / np.sum(self.map)
            not_masked_cluster_fraction = np.sum(self.map[not_masked]) / np.sum(self.map)
            fully_masked_sky_fraction = np.sum(sdss_allowed_region[masked]) / np.sum(sdss_allowed_region)
            not_masked_sky_fraction = np.sum(sdss_allowed_region[not_masked]) / np.sum(sdss_allowed_region)
            sky_fraction = planck_masked_only[mixed_bins] / sdss_allowed_region[mixed_bins]
            cluster_fraction = sample_masked_fraction[mixed_bins]
            n = self.map[mixed_bins]
        else:
            fully_masked_cluster_fraction = 0
            not_masked_cluster_fraction = 0
            fully_masked_sky_fraction = 0
            not_masked_sky_fraction = 0
            sky_fraction = np.array((planck_masked_only / sdss_allowed_region))
            cluster_fraction = np.array((sample_masked_fraction))
            n = self.map

        return [[fully_masked_cluster_fraction, not_masked_cluster_fraction, fully_masked_sky_fraction,
                 not_masked_sky_fraction], sky_fraction, cluster_fraction, n]


class BinaryMap(_BinMap):
    def __init__(self):
        self.map = np.array((0.0, 0.0))
        self.weighted_map = np.array((0.0, 0.0))

    def lookup_point(self, _lon, lat):
        return self.map[np.int_(lat > 0)]

    def lookup_weighted_point(self, _lon, lat):
        return self.weighted_map[np.int_(lat > 0)]

    def lookup_pix(self, _lon, lat):
        return np.int_(lat > 0)

    def calc_weighted_map(self, mask, resolution=(1000, 1000)):
        bin = self.lookup_pix(*hp.pix2ang(mask.NSIDE, np.arange(hp.nside2npix(mask.NSIDE)), lonlat=True))
        print(np.min(bin), np.max(bin))
        self.mask_fraction_map = np.array([np.sum(mask.map[bin == 0]) / len(mask.map[bin == 0]),
                                           np.sum(mask.map[bin == 1]) / len(mask.map[bin == 1])])
        print(f"Mask fraction map: {self.mask_fraction_map}")
        with np.errstate(divide="ignore"):
            inverted_fraction_map = np.power(self.mask_fraction_map,
                                             -1)  # spits out inf/NaN values, ignore RuntimeWarning
        inverted_fraction_map[np.bitwise_or(np.isnan(inverted_fraction_map), np.isinf(inverted_fraction_map))] = 0
        self.weighted_map = self.map * inverted_fraction_map


class ConstantMap(_BinMap):
    def __init__(self):
        self.map = np.ones(1)
        self.weighted_map = np.ones(1)

    def lookup_point(self, _lon, lat):
        return np.ones(len(lat))

    def lookup_weighted_point(self, _lon, lat):
        return np.ones(len(lat))

    def lookup_pix(self, _lon, lat):
        try:
            return np.zeros(len(lat))
        except TypeError:
            return 0

    def bin_catalogue(self, catalogue):
        self.binned_sample = [[]]
        for cluster in catalogue.lon_lat:
            self.binned_sample[self.lookup_pix(*cluster)].append(cluster)
        for map_bin in range(len(self.binned_sample)):
            self.binned_sample[map_bin] = np.array(self.binned_sample[map_bin])


class HealpixBinMap(_BinMap):
    def __init__(self, NSIDE):
        self.NSIDE = NSIDE
        self.map = np.zeros(hp.nside2npix(NSIDE))
        self.weighted_map = np.zeros(hp.nside2npix(NSIDE))

    def lookup_point(self, lon, lat):
        return self.map[hp.ang2pix(self.NSIDE, lon, lat, lonlat=True)]

    def lookup_weighted_point(self, lon, lat):
        return self.weighted_map[hp.ang2pix(self.NSIDE, lon, lat, lonlat=True)]

    def lookup_pix(self, lon, lat):
        return hp.ang2pix(self.NSIDE, lon, lat, lonlat=True)

    def calc_weighted_map(self, mask, resolution=(1000, 1000)):
        self.mask_fraction_map = hp.pixelfunc.ud_grade(mask.map, self.NSIDE)
        # with np.seterr(all="ignore"):
        inverted_fraction_map = np.power(self.mask_fraction_map, -1)  # spits out inf/NaN values, ignore RuntimeWarning
        inverted_fraction_map[np.bitwise_or(np.isnan(inverted_fraction_map), np.isinf(inverted_fraction_map))] = 0
        self.weighted_map = self.map * inverted_fraction_map


def gen_mask_comparison_map(mask1, mask2, NSIDE=512, NSIDE_internal=2048, name="", res=int(1e4)):
    """pix = np.linspace(0, NSIDE_internal - 1, NSIDE_internal)
    print(mask1.NSIDE, mask2.NSIDE)
    sum = mask1.map + 1j * mask2.map
    data1 = hp.ud_grade(sum == 0.0, NSIDE)
    data2 = hp.ud_grade(sum == 1.0, NSIDE)
    data3 = hp.ud_grade(sum == 1j, NSIDE)
    data4 = hp.ud_grade(sum == 1 + 1j, NSIDE)
    results = np.float_(np.array([data1, data2, data3, data4]))"""
    """x = np.linspace(-180, 180, 2 * res)
    y = np.linspace(-90, 90, res)
    x, y = np.meshgrid(x, y)
    data1 = mask1.lookup_point(x, y)
    data2 = mask2.lookup_point(x, y)
    data1[data1 > 0] = 1.0
    data2[data2 > 0] = 1.0
    data = data1 + 1j * data2
    map = np.zeros((5, hp.nside2npix(NSIDE)))
    bins = hp.ang2pix(NSIDE, x, y, lonlat=True)
    for pixel in range(len(map[0])):
        points = bins == pixel
        print(pixel / len(map[0]))
        reduced_data = data[points]
        temp = y[points] * np.pi / 180
        map[0, pixel] += np.sum(np.cos(temp))
        map[1, pixel] += np.sum(np.cos(temp[reduced_data == 0.0]))
        map[2, pixel] += np.sum(np.cos(temp[reduced_data == 1.0]))
        map[3, pixel] += np.sum(np.cos(temp[reduced_data == 1j]))
        map[4, pixel] += np.sum(np.cos(temp[reduced_data == 1.0 + 1j]))
    print(map)
    results = map[1:] / (map[0] + 1e-100)
    print(np.max(results[0]), np.max(results[1]), np.max(results[2]), np.max(results[3]))"""
    pix = np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1, hp.nside2npix(NSIDE_internal)))
    points = hp.pix2ang(NSIDE_internal, pix, lonlat=True)
    sum = mask1.lookup_point(*points) + 1j * mask2.lookup_point(*points)
    mask1_masked = mask1.lookup_point(*points) == 0.0
    mask2_masked = mask2.lookup_point(*points) == 0.0
    print(np.shape(sum))
    print(sum)
    """data1 = hp.ud_grade(sum == 0.0, NSIDE)
    data2 = hp.ud_grade(sum == 1.0, NSIDE)
    data3 = hp.ud_grade(sum == 1j, NSIDE)
    data4 = hp.ud_grade(sum == 1 + 1j, NSIDE)"""
    data1 = hp.ud_grade(np.float_(np.bitwise_and(mask1_masked, mask2_masked)), NSIDE)
    data2 = hp.ud_grade(np.float_(np.bitwise_and(np.bitwise_not(mask1_masked), mask2_masked)), NSIDE)
    data3 = hp.ud_grade(np.float_(np.bitwise_and(mask1_masked, np.bitwise_not(mask2_masked))), NSIDE)
    data4 = hp.ud_grade(np.float_(np.bitwise_and(np.bitwise_not(mask1_masked), np.bitwise_not(mask2_masked))), NSIDE)
    results = np.float_(np.array([data1, data2, data3, data4]))
    hp.fitsfunc.write_map(f"./{name}_{NSIDE}_1.fits", results[0], overwrite=True)
    hp.fitsfunc.write_map(f"./{name}_{NSIDE}_2.fits", results[1], overwrite=True)
    hp.fitsfunc.write_map(f"./{name}_{NSIDE}_3.fits", results[2], overwrite=True)
    hp.fitsfunc.write_map(f"./{name}_{NSIDE}_4.fits", results[3], overwrite=True)
    print(results)
    print(np.min(data1 + data2 + data3 + data4), np.max(data1 + data2 + data3 + data4))


def run_const(data_set, mask, filter_set, a, cat, weight_function, convert_to_mask_frac):
    data = np.array((
        np.mean(data_set[0]),
        np.mean(data_set[1]),
        np.mean(data_set[2]),
        np.mean(data_set[3])
    ))
    binmap = ConstantMap()
    binmap.bin_catalogue(cat)
    binmap.load_catalogue(cat)
    output = binmap.divide_sample(mask, data, True, filter_set, a)
    mixed = weight_function(*output[1:], skip_n_filter=True)
    print(f"Mixed C: {mixed[0]} +/- {mixed[1]}")
    if convert_to_mask_frac:
        final = np.array(
            [output[0][0] + (1 - output[0][0] - output[0][1]) * mixed[0], (1 - output[0][0] - output[0][1]) * mixed[1]])
    else:
        final = np.array(mixed)
    print(f"Final C: {final[0]} +/- {final[1]}")
    return final


def run_nside(n, data_set, mask, filter_set, a, cat, weight_function, convert_to_mask_frac):
    try:
        data = np.array((
            hp.ud_grade(data_set[0], n),
            hp.ud_grade(data_set[1], n),
            hp.ud_grade(data_set[2], n),
            hp.ud_grade(data_set[3], n)
        ))
        binmap = HealpixBinMap(n)
        binmap.bin_catalogue(cat)
        binmap.load_catalogue(cat)
        output = binmap.divide_sample(mask, data, False, filter_set, a)
        mixed = weight_function(*output[1:])
        print(f"Mixed {n}: {mixed[0]} +/- {mixed[1]}")
        print(output[0])
        if convert_to_mask_frac:
            final = np.array([output[0][0] * 100 + (1 - output[0][0] - output[0][1]) * mixed[0],
                              (1 - output[0][0] - output[0][1]) * mixed[1]])
        else:
            final = np.array(mixed)
        print(f"Final {n}: {final[0]} +/- {final[1]}")
    except NotMaskError:
        final = np.array([np.NaN, np.NaN])
    #results.append(final)
    return final


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


def load_mask(mask, raise_dir=2, nside=8, invert=False):
    value = None
    if mask == "planck_galactic":
        # value = HealpyMask("../" * raise_dir + "data/HFI_PCCS_SZ-selfunc-union-survey_R2.08.fits", mask_using_latlon=True, hdu=1, partial=False)
        value = HealpyMask("../" * raise_dir + "data/planck_galactic_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "planck_point" or mask == "planck_modified_total":
        value = HealpyMask("../" * raise_dir + "data/HFI_PCCS_SZ-selfunc-inter-cosmo_2.02.fits", partial=False,
                           mask_using_latlon=True)
    elif mask == "planck_survey":
        value = HealpyMask("../" * raise_dir + "data/planck_survey_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "sdss_mask":
        value = HealpyMask("../" * raise_dir + "data/redmapper_dr8_public_v6.3_zmask.fits", mask_using_latlon=False,
                           hdu=1, partial=True)
        value.map[value.map > 0.4] = 1.0
        value.map[value.map < 0.3] = 0
        rotator = hp.Rotator(coord=["C", "G"])
        value.map = rotator.rotate_map_pixel(value.map)
        value.mask_using_latlon = True
        value.map[value.map > 0.5] = 1.0
        value.map[value.map < 0.5] = 0
        # value.map = (value.map - 1) * -1
    elif mask == "planck_modified_galactic":
        point_mask = load_mask("planck_point", raise_dir=raise_dir)
        galactic_mask = load_mask("planck_galactic", raise_dir=raise_dir)
        compound_map = load_mask("planck_point", raise_dir=raise_dir)
        value = load_mask("planck_point", raise_dir=raise_dir)
        compound_map.map = 2 * point_mask.map + galactic_mask.map
        value.map = np.float_(np.bitwise_not(compound_map.map == 0))
    elif mask == "planck_modified_point":
        point_mask = load_mask("planck_point", raise_dir=raise_dir)
        galactic_mask = load_mask("planck_galactic", raise_dir=raise_dir)
        compound_map = load_mask("planck_point", raise_dir=raise_dir)
        value = load_mask("planck_point", raise_dir=raise_dir)
        compound_map.map = 2 * point_mask.map + galactic_mask.map
        value.map = np.float_(np.bitwise_not(compound_map.map == 1))
    elif mask == "comparison_sdss_planck_galactic":
        # masked by
        planck_only = HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_2.fits")
        neither = HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_4.fits")
        data = planck_only.map / (planck_only.map + neither.map + 1e-30)
        data = hp.pixelfunc.ud_grade(data, nside)
        planck_only.map = data
        planck_only.NSIDE = nside
        planck_only.NPIX = hp.nside2npix(nside)
        return planck_only
    elif mask.lower() == "act":
        value = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, mask_using_latlon=False, invert=invert)
    elif mask == "sdss_planck_point_only":
        mask1 = load_mask("sdss_mask", raise_dir)
        mask2 = load_mask("planck_modified_point", raise_dir)
        mask2.map = 1 - mask2.map
        value = CombinationMask(mask1, mask2)
    elif mask == "planck_galactic_test":
        value = HealpyMask("../" * raise_dir + "data/planck_galactic_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "planck_point_test":
        value = HealpyMask("../" * raise_dir + "data/planck_point_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "act_galactic":
        act_mask = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
        act_graph_filter = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
        act_graph_filter.map = np.load("../" * raise_dir + "data/act_galactic_mask_array.npy")
        test = filters.SquareMaskFilter((-10, 20), (40, 100), lonlat=False)
        final_filter = CombinationMask(act_mask, test, use_and=False, invert=False)
        value = CombinationMask(act_graph_filter, final_filter, use_and=True, invert=False)
    elif mask == "act_point":
        act_mask = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
        galactic_mask = load_mask("act_galactic", raise_dir, nside, invert)
        galactic_mask.invert = True
        value = CombinationMask(act_mask, galactic_mask, use_and=False)
    else:
        raise ValueError(f"{mask} is not a recognised mask")
    return value


def load_catalogue(cat, raise_dir=2, **kwargs):
    if cat == "sdss":
        data = StarCatalogue(raise_dir * "../" + "data/sdss_catalogue.fits", hdu=1)
        data.load_lon_lat()
    else:
        raise ValueError(f"{cat} is not a recognised catalogue")
    return data


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


def fraction_masked_pair(mask1, mask2, n=int(1e3), ram_limited=False, weight_map=None):
    if (isinstance((mask1, mask2),
                   HealpyMask) and mask2.NPIX == mask1.NPIX and mask1.using_latlon and mask2.using_latlon
            and not weight_map):
        data = mask1.compare(mask2)
        k = np.max((mask1.NPIX, mask2.NPIX))
        results = [np.real(np.sum(data) / k), np.imag(np.sum(data) / k), np.sum(data == 1 + 1j) / k,
                   np.sum(data == 1) / k, np.sum(data == 1j) / k,
                   np.sum(data == 0) / k]
    elif not ram_limited and not weight_map:
        x = np.linspace(-np.pi, np.pi, n)
        y = np.linspace(-np.pi / 2, np.pi / 2, n)
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
            sums[3] += np.sum(data[i] == 1 + 1j) * s
            sums[4] += np.sum(data[i] == 0) * s
        results = sums[1:] / sums[0]
    else:
        if weight_map is None:
            weight_map = ConstantMap()
        x = np.linspace(-np.pi, np.pi, n)
        y = np.linspace(-np.pi / 2, np.pi / 2, n)
        sums = np.array([0, 0, 0, 0, 0])
        for i in range(len(y)):
            theta = np.pi * (i + 0.5) / len(y)
            s = np.sin(theta)
            data1 = mask1.lookup_point(x * 180 / np.pi, (theta * 180 / np.pi) - 90)
            data2 = mask2.lookup_point(x * 180 / np.pi, (theta * 180 / np.pi) - 90)
            weight = weight_map.lookup_weighted_point(x * 180 / np.pi, ((theta * 180 / np.pi) - 90) * np.ones(len(x)))
            data = np.array(data1 + 1j * data2)
            # data = np.array(data1 + 1j * data2)
            sums[0] += np.sum(np.ones(n) * weight) * s
            sums[1] += np.sum(np.float_(data == 1) * weight) * s
            sums[2] += np.sum(np.float_(data == 1j) * weight) * s
            sums[3] += np.sum(np.float_(data == 1 + 1j) * weight) * s
            sums[4] += np.sum(np.float_(data == 0) * weight) * s
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


def match_distribution(random_coords, sample, nside):
    binmap_sample = HealpixBinMap(nside)
    binmap_sample.bin_catalogue(sample)
    binmap_random = HealpixBinMap(nside)
    binmap_random.bin_catalogue(random_coords)
    ratio = np.zeros(hp.nside2npix(nside))

    for bin_num in range(hp.nside2npix(nside)):
        try:
            ratio[bin_num] = len(binmap_sample.binned_sample[bin_num]) / len(binmap_random.binned_sample[bin_num])
        except ZeroDivisionError:
            ratio[bin_num] = 0

    new_ratio = 1 / np.max(ratio)
    final_data_set = []
    first = True
    for bin_num in range(hp.nside2npix(nside)):
        if ratio[bin_num] != 0:
            if first:
                final_data_set = binmap_random.binned_sample[bin_num][
                                 :int(new_ratio * len(binmap_sample.binned_sample[bin_num]))]
                first = False
            else:
                final_data_set = np.append(final_data_set, binmap_random.binned_sample[bin_num][
                                                           :int(new_ratio * len(binmap_sample.binned_sample[bin_num]))],
                                           axis=0)
    print(final_data_set)
    return final_data_set


def ra_dec_to_lon_lat(ra, dec, reverse=False):
    if not reverse:
        c = astropy.coordinates.SkyCoord(ra=ra, dec=dec, unit="deg")  # convert to galactic co ords
        lon = c.galactic.l.degree
        lat = c.galactic.b.degree
    else:
        c = astropy.coordinates.SkyCoord(lon=ra, lat=dec)  # convert to galactic co ords
        lon = c.icrs.ra.degree
        lat = c.icrs.dec.degree
    return lon, lat


#  Define a black and white heatmap
__cdict_bw = [(0, 0, 0), (1, 1, 1)]
bw_heatmap = matplotlib.colors.LinearSegmentedColormap.from_list("black and white", __cdict_bw, N=2)
rb_heatmap = matplotlib.colors.LinearSegmentedColormap.from_list("rb", [(1, 0, 0), (0, 0, 1)])

# To Do
"""Minimum richness, larger data set, binning data (z plus richness), then start ACT"""
