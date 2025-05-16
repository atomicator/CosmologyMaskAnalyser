# This defines all the classes used by the rest of the project
import multiprocessing
import warnings
from idlelib.autocomplete import ATTRS
from multiprocessing.pool import ThreadPool
from os import error
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pixell
import pixell.enmap
import healpy as hp
import astropy.cosmology
import astropy.coordinates
import astropy.io.fits
import astropy.table
import astropy.wcs.utils
import copy
from astropy import units as u


def plt_use_tex():  # updates variables in the matplotlib library so that it uses LaTeX
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        'text.latex.preamble': r'\usepackage{csquotes}'
                               r'\usepackage{amsmath}'
    })
    return None


class __Mask(object):  # A template class, used to define general methods. Separate inherited classes are needed for
    # each type of mask used.

    fig = None  # The references to a matplotlib plot
    ax = None

    lon_shift = 0  # The shift when querying a co-ordinate
    mask_using_latlon = True  # Stores if the mask is using latlon or RA-DEC
    invert = False  # Stores if the output should be inverted

    def __init__(self, lon_shift=0, mask_using_latlon=True, invert=False):
        self.lon_shift = lon_shift
        self.mask_using_latlon = mask_using_latlon
        self.invert = invert

    def lookup_point(self, lon, lat):  # A wrapper around _lookup_point, which handles the mask file itself
        if self.lon_shift:  # If the co-ords are shifted
            #lon = lon.clone()
            lon = np.copy(lon)
            lon += self.lon_shift
            lon[lon > 360] -= 360
            lon[lon < 0] += 360
        if not self.mask_using_latlon:  # If the co-ords need to be converted to RA-DEC
            c = astropy.coordinates.SkyCoord(lon, lat, unit="deg", frame="galactic")  # convert to galactic co ords
            lon = c.icrs.ra.degree
            lat = c.icrs.dec.degree
        points = self._lookup_point(lon, lat)  # _lookup_point handles the methods directly related to the mask file
        if self.invert:  # If the mask needs inverting
            points = 1 - points
        return points

    def _lookup_point(self, lon, lat):  # Needs to be overwritten when inherited - querying a point depends on the
        # format of the mask
        raise NotMaskError("__Mask._lookup_point() called")

    def set_fig_ax(self, fig, ax):  # updates the matplotlib environment
        self.fig = fig
        self.ax = ax

    def calc_exact_unmasked_fraction(self):  # Needs to be overwritten when inherited
        raise NotMaskError("__Mask.calc_exact_unmasked_fraction() called")

    def plot(self, resolution=(1000, 2000), clear=True, cbar=None, title=None, save_path=None, show=False,
             cmap="plasma", alpha=1, nan_filter=True, lon_range=(0, 360), lat_range=(-90, 90), use_lon_lat=True):  # Resolution is the resolution of the plot. Clear determines if
        # the plot environment needs clearing before it gets plotted. Cbar adds a colourbar to the plot. Title adds a
        # title to the plot. Nan_filter determines if the non-masked region should be transparent. Cmap is the colourmap
        # to be used for the plot.
        if clear:
            self.ax.cla()  # Clear the plot environment
        x, y = resolution
        lat = np.linspace(lat_range[0], lat_range[1], x)  # Create a 2D grid for a heatmap
        lon = np.linspace(lon_range[0], lon_range[1], y)
        lon, lat = np.meshgrid(lon, lat)
        if not use_lon_lat:
            #temp = astropy.coordinates.SkyCoord(l=lon*u.degree, b=lat*u.degree, frame='galactic')
            #pixel_data = self.lookup_point(temp.icrs.ra.degree, temp.icrs.dec.degree)
            temp = astropy.coordinates.SkyCoord(ra=lon*u.degree, dec=lat*u.degree, frame='icrs')
            pixel_data = self.lookup_point(temp.galactic.l.degree, temp.galactic.b.degree)
        else:
            pixel_data = self.lookup_point(lon, lat)
        #plt.imshow(pixel_data)
        #plt.title("pixel_data")
        #plt.colorbar()
        #plt.show()
        if nan_filter:
            pixel_data[pixel_data == 1] = np.nan  # Convert the non-masked regions to be transparent
        #plt.imshow(pixel_data)
        #plt.title("pixel_data_2")
        #plt.colorbar()
        #plt.show()
        lon = np.linspace(lon_range[0], lon_range[1], y + 1)  # Create a co-ordinate grid for a heatmap
        lat = np.linspace(lat_range[0], lat_range[1], x + 1)
        lon, lat = np.meshgrid(lon, lat)
        c = self.ax.pcolormesh(lon, lat, pixel_data, alpha=alpha, cmap=cmap)  # Plot a colourmap
        if cbar:
            self.fig.colorbar(c, orientation="horizontal", cmap=cmap)  # add a colourbar to the plot
        if title:
            self.ax.title(title)  # add a title to the plot
        if save_path:
            plt.savefig(save_path)  # save the plot
        if show:
            plt.show()  # show the plot

    def clone(self):  # clones the mask
        fig, ax = self.fig, self.ax
        new = copy.deepcopy(self)
        new.fig, new.ax = fig, ax
        return new


class PixellMask(__Mask):  # Defines the methods for the PixellMask classes
    def __init__(self, path, mask_using_latlon=True, invert=False, lon_shift=0, init_val=0, suppress_warnings=False, **kwargs):
        # Load data from fits file
        super().__init__(lon_shift, mask_using_latlon, invert)  # Calls the constructor of the base class
        if suppress_warnings:
            with warnings.catch_warnings():
                warnings.filterwarnings(action="ignore", category=RuntimeWarning)
                self.imap = pixell.enmap.read_map(path, **kwargs)  # Loads the file
        else:
            self.imap = pixell.enmap.read_map(path, **kwargs)  # Loads the file
        self.map = np.array(self.imap)  # A numpy version of the map
        self.x = len(self.map)  # The dimensions of the map - useful to have them saved
        self.y = len(self.map[0])
        self.lat_min = -90  # The range used for plotting the graphs - updated by a different function
        self.lat_max = 90
        self.lon_min = 00
        self.lon_max = 360
        self.init_val = init_val
        #print(self.imap.wcs)

    def _lookup_point(self, lon, lat):
        #print(np.min(lon), np.max(lon))
        #print(np.min(lat), np.max(lat))
        #print(self.imap.wcs)
        #wcs = astropy.wcs.WCS(self.imap.wcs)
        # Note: co-ord are transformed by wrapper so they match the mask
        if not self.mask_using_latlon:
        #if True:
            c = astropy.coordinates.SkyCoord(ra=lon * u.degree, dec=lat * u.degree)
        else:
            c = astropy.coordinates.SkyCoord(lon * u.degree, lat * u.degree)
        pix = np.int_(self.imap.wcs.world_to_pixel(c))[::-1]
        # Line below works for ACT mask, not for SPT
        #pix = np.int_(pixell.enmap.sky2pix(self.imap.shape, self.imap.wcs, np.array((lat, lon)) * np.pi / 180))
        #print(pix, lon, lat)
        # convert lat lon co-ords to pixel co-ords
        point = np.ones(np.shape(pix)[1:]) * self.init_val  # initialise the array as masked
        points_in_range = np.bitwise_and(
            np.bool_(np.int_(0 < np.array(pix[0])) * np.int_(np.array(pix[0]) < np.shape(self.map)[0] - 1)),
            np.bool_(np.int_(0 < np.array(pix[1])) * np.int_(np.array(pix[1]) < np.shape(self.map)[1] - 1)))
        #plt.imshow(points_in_range)
        #plt.title("Points in range")
        #plt.show()
        #print(f"Points in range: {np.sum(points_in_range)} / {points_in_range.shape}")
        #print(f"Default: {point[0][0]}")
        #print(f"Pix range: {pix[0][0]}, {pix[0][-1]}, {pix[-1][0]}, {pix[-1][-1]}")
        # Check which co-ords are within the range detailed in the mask
        #print(f"{np.sum(point)}")
        point[points_in_range] = self.map[pix[0, points_in_range], pix[1, points_in_range]]  # update the points in the
        #print(f"{np.sum(point)}")
        # region covered by the mask
        #plt.imshow(point)
        #plt.title("Point")
        #plt.show()
        return point

    def calc_exact_unmasked_fraction(self):
        total_area = 0
        unmasked_area = 0
        boundary = pixell.enmap.pix2sky(self.imap.shape, self.imap.wcs, [np.array((0, self.x)),
                                                                         np.array((0, self.y))])  # The range covered
        # by the map in sky co-ords
        dec_range = boundary[0]
        for i in range(self.x):
            s = np.cos(dec_range[0] + (i / self.x) * (dec_range[1] - dec_range[0]))
            total_area += s * self.y
            unmasked_area += s * np.sum(self.map[i])
        # fraction of covered sky masked * fraction of sky covered by map
        return (unmasked_area / total_area) * (np.sin(dec_range[1]) - np.sin(dec_range[0])) / 2


class HealpyMask(__Mask):
    def __init__(self, path, nest=False, mask_using_latlon=True, lon_shift=None, invert=False, **kwargs):  # path is the
        # path to the file. Nest is if the mask is using a nested co-ord system. Mask_using_latlon stores if the mask is
        # using LAT / LON or RA / DEC. Lon_shift stores if the points should be shifted.
        super().__init__(lon_shift, mask_using_latlon, invert)
        self.map, header = hp.read_map(path, h=True, nest=nest, **kwargs)  # load the file
        self.header = dict(header)  # convert the header to a dictionary, so it can be read
        self.NSIDE = int(self.header["NSIDE"])
        self.NPIX = 12 * self.NSIDE ** 2
        self.nest = nest

    def _lookup_point(self, lon, lat):
        pix = hp.ang2pix(self.NSIDE, lon, lat, lonlat=True, nest=self.nest)
        return self.map[pix]

    def calc_exact_unmasked_fraction(self):
        unmasked_points = np.sum(self.map)
        unmasked_fraction = unmasked_points / self.NPIX
        return unmasked_fraction


class CombinationMask(__Mask):  # An object that behaves like a mask - it combines two masks with a logic operator
    # Currently supports AND, NAND, OR, NOR
    def __init__(self, mask1, mask2, invert=False, use_and=True, lon_shift=0):
        super().__init__(invert=invert, lon_shift=lon_shift, mask_using_latlon=True)
        self.mask1 = mask1  # The masks to be compared
        self.mask2 = mask2
        self.lat_min = -90
        self.lat_max = 90
        self.lon_min = 00
        self.lon_max = 360
        self.fig = None
        self.ax = None
        self.mask_using_latlon = True
        self.invert = invert  # Bool - converts AND to NAND etc.
        self.use_and = use_and  # Bool - stores if AND or OR should be used

    def _lookup_point(self, lon, lat):
        if self.use_and:
            data = np.float64(np.bitwise_and(np.bool_(self.mask1.lookup_point(lon, lat)),
                                             np.bool_(self.mask2.lookup_point(lon, lat))))
        else:
            data = np.float64(np.bitwise_or(np.bool_(self.mask1.lookup_point(lon, lat)),
                                            np.bool_(self.mask2.lookup_point(lon, lat))))
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


class ClusterCatalogue(object):
    def __init__(self, path=None, hdu=0, table=False, **kwargs):  # Path is the path to the file. HDU is the hdu of the
        # file to load. Table stores if the data is stored in an ordinary fits file.
        if path:
            if not table:
                self.fits = astropy.io.fits.open(path, **kwargs)  # load the data
                self.h = dict(self.fits[hdu].header)
                self.data = self.fits[hdu].data
            if table:
                self.fits = astropy.table.Table.read(path)  # load the table
                self.data = self.fits
                self.h = None
        else:
            self.fits = None
            self.h = None
            self.data = None

        self.lon_lat = None  # stores the lon / lat co-ords of the
        self.heatmap = None  # stores data of the heatmap
        self.NSIDE = None  # stores the NSIDE used to create the heatmap
        self.fig = None  # stores references to the matplotlib environment
        self.ax = None
        self.cat_mask = None  # The mask used for the catalogue - this is used for the NaN filter

    def set_cat_mask(self, mask):  # updates the catalogue mask
        self.cat_mask = mask

    def set_fig_ax(self, fig, ax):  # updates the matplotlib environment
        self.fig = fig
        self.ax = ax

    def load_data(self, lon_lat=True, selection_function=None, requested_fields=(), *args, **kwargs):
        # List of co-ord column names, in the format [lon_name, lat_name, "lon_lat"] or [dec_name, ra_name, "RA_DEC"]
        column_list = [["GLON", "GLAT", "lon_lat"], ["RA", "DEC", "RA_DEC"], ["decDeg", "RADeg", "RA_DEC"],]
        loaded = False
        data_format = None
        data = None
        for column_set in column_list:  # for each column
            try:
                data = [self.data.field(column_set[0]), self.data.field(column_set[1])]
                data_format = column_set[2]
                loaded = True
                break
            except KeyError:  # if the column isn't in the table
                pass
        if not loaded:  # if none of the keys matched
            raise KeyError("None of the keys match the data")

        if lon_lat and data_format == "RA_DEC":
            c = astropy.coordinates.SkyCoord(ra=data[0], dec=data[1], unit="deg")  # convert ICRS to galactic co ords
            data = [c.galactic.l.degree, c.galactic.b.degree]
        elif not lon_lat and data_format == "lon_lat":
            c = astropy.coordinates.SkyCoord(lon=data[0], lat=data[1], unit="deg", frame="galactic")  # convert galactic
            # to ICRS co ords
            data = [c.icrs.dec.degree, c.icrs.ra.degree]
        data = np.array(data).transpose()
        if selection_function:  # if a filter is being applied to the data
            requested_data = []
            for field in requested_fields:  # for each field to load from the data
                requested_data += [self.data.field(field)]
            requested_data = np.array(requested_data).transpose()
            new_data = []
            for i in range(len(requested_data)):  # for each cluster
                if selection_function(*requested_data[i], *args, **kwargs):  # if it passes the selection function
                    new_data += [data[i]]  # add the data to the updated list
            data = new_data
        self.lon_lat = np.array(data)

    def plot_on_ax(self, ax, alpha, resolution, cmap="plasma"):
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
        self.heatmap = self.heatmap / (resolution ** 2)

    def plot_heatmap(self, nside, resolution=(1000, 2000), clear=True, cbar=None, title=None, save_path=None,
                     show=False, cmap="plasma", alpha=1, cbar_label=""):
        self.gen_heatmap_data(nside)
        if clear:
            self.ax.cla()
        c = self.plot_on_ax(self.ax, alpha, resolution, cmap=cmap)
        if cbar:
            self.fig.colorbar(c, orientation="horizontal", location="bottom", cmap=cmap, label=cbar_label)
        if title:
            self.ax.title(title)
        if save_path:
            plt.savefig(save_path)
        if show:
            plt.show()

    def plot_scatter(self, **kwargs):
        self.ax.scatter(*self.lon_lat.transpose(), **kwargs)


def gen_random_coords(n, mask=None):
    np.random.seed()
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


class _BinMap(object):
    n = np.array([])
    binned_sample = [[]]
    sky_masked_fraction = []
    cluster_masked_fraction = []
    mask = None
    final_filter = None

    def lookup_pix(self, lon, lat):
        raise NotImplementedError

    def set_mask(self, mask):
        self.mask = mask

    def bin_catalogue(self, catalogue):
        print("Looking up clusters")
        masked = self.mask.lookup_point(*catalogue.lon_lat.transpose())
        masked_count_pixels = np.zeros(np.shape(self.n))
        print("Binning Catalogue")
        for cluster_number in range(len(catalogue.lon_lat)):
            cluster = catalogue.lon_lat[cluster_number]
            pix_location = self.lookup_pix(*cluster)
            self.binned_sample[pix_location].append(cluster)
            masked_count_pixels[pix_location] += (1 - masked[cluster_number])
        print("Sorting Catalogue")
        for i in range(len(self.binned_sample)):
            self.binned_sample[i] = np.array(self.binned_sample[i])
            self.n[i] = len(self.binned_sample[i])
        has_clusters_filter = self.n > 0
        self.cluster_masked_fraction[has_clusters_filter] = masked_count_pixels[has_clusters_filter] / self.n[has_clusters_filter]
        self.cluster_masked_fraction[np.bitwise_not(has_clusters_filter)] = np.nan
        print("Sorted catalogue")
        #self.calc_masked_fraction()

    def calc_masked_fraction(self):
        print("Calculating masked fraction")
        for pixel in range(len(self.binned_sample)):
            #try:
            if len(self.binned_sample[pixel]) > 0:
                temp = np.array(self.binned_sample[pixel])
                self.cluster_masked_fraction[pixel] = 1 - np.sum(self.mask.lookup_point(*temp.transpose())) \
                                                    / len(self.binned_sample[pixel])
            #except TypeError:  # runs if there's no clusters in the pixel
            else:
                self.cluster_masked_fraction[pixel] = np.nan
        self.cluster_masked_fraction = np.array(self.cluster_masked_fraction)

    def divide_sample(self, mask, area_masked_fraction, filter_fully_masked=True, filter_empty=True):
        #  Cat allowed region - the fraction of (total) sky in which the clusters could have been detected in.
        #  Survey masked only - the fraction of the (total) sky in which the clusters could have spawned and is masked
        #                       by the survey mask
        cat_allowed_region = area_masked_fraction[1] + area_masked_fraction[3]
        survey_masked_only = area_masked_fraction[1]
        with np.errstate(divide='ignore', invalid='ignore'):  # supress warnings for dividing by zero - these pixels are
            # filtered below
            self.sky_masked_fraction = np.array(survey_masked_only / cat_allowed_region)

        # Filter pixels that have no clusters, or are fully / not masked
        nan_filter = np.bitwise_or(np.isnan(self.sky_masked_fraction), np.isinf(self.sky_masked_fraction))
        if filter_empty:
            nan_filter = np.bitwise_or(nan_filter, self.n == 0)
        else:
            self.cluster_masked_fraction[self.n == 0] = 0
        fully_masked_filter = self.sky_masked_fraction == 1.
        not_masked_filter = self.sky_masked_fraction == 0.

        # Calculate the sky / cluster fractions masked by the filters
        fully_masked_cluster_fraction = np.sum(self.n[fully_masked_filter]) / np.sum(self.n)
        not_masked_cluster_fraction = np.sum(self.n[not_masked_filter]) / np.sum(self.n)
        fully_masked_sky_fraction = np.sum(cat_allowed_region[fully_masked_filter]) / np.sum(cat_allowed_region)
        not_masked_sky_fraction = np.sum(cat_allowed_region[not_masked_filter]) / np.sum(cat_allowed_region)

        # Sky fractions ignore the no clusters filter - add this as a separate value?
        if filter_fully_masked:
            self.final_filter = np.bitwise_not(np.bitwise_or(nan_filter, np.bitwise_or(fully_masked_filter, not_masked_filter)))
        else:
            self.final_filter = np.bitwise_not(nan_filter)
        sky_fraction = self.sky_masked_fraction[self.final_filter]
        cluster_fraction = self.cluster_masked_fraction[self.final_filter]
        n = self.n[self.final_filter]

        #print("test")
        #print(f"{len(self.n)}: {n[survey_masked_only[final_filter] < 0.05]}")

        return [[fully_masked_cluster_fraction, not_masked_cluster_fraction, fully_masked_sky_fraction,
                 not_masked_sky_fraction], sky_fraction, cluster_fraction, n]


class BinaryBinMap(_BinMap):
    n = np.zeros(2)
    cluster_masked_fraction = np.zeros(2)
    binned_sample = [[], []]

    def lookup_pix(self, _lon, lat):
        value = None
        if np.shape(_lon) != np.shape(lat):
            raise ValueError
        elif np.shape(_lon) == ():
            if lat > 0:
                value = 0
            else:
                value = 1
        else:
            value = np.int_(lat < 0)
        return value


class ConstantBinMap(_BinMap):
    def __init__(self):
        self.n = np.zeros(1)
        self.cluster_masked_fraction = np.zeros(1)
        self.binned_sample = [[],]

    def lookup_pix(self, _lon, _lat):
        if np.shape(_lon) != np.shape(_lat):
            raise ValueError
        elif np.shape(_lon) == ():
            return 0
        return np.array(np.zeros(np.shape(_lon)))


class HealpixBinMap(_BinMap):
    def __init__(self, nside):
        self.nside = nside
        self.n = np.zeros(hp.nside2npix(self.nside))
        self.cluster_masked_fraction = np.zeros(hp.nside2npix(self.nside))
        self.binned_sample = [[] for i in range(hp.nside2npix(self.nside))]

    def lookup_pix(self, lon, lat):
        value = None
        if np.shape(lon) != np.shape(lat):
            raise ValueError
        else:
            value = hp.ang2pix(self.nside, lon, lat, lonlat=True)
        return value


def gen_mask_comparison_map(load_func, args, kwargs, mask1_name, mask2_name, NSIDE=512, NSIDE_internal=2048, name="", write=True,
                            copy=False, num_thread=1, skip=0):
    #print("Initialising pix array")
    #pix = np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1, hp.nside2npix(NSIDE_internal)))
    #print("Loading masks")
    #mask1 = load_func(mask1_name, *args, **kwargs)
    #mask2 = load_func(mask2_name, *args, **kwargs)
    #print("Creating point array")
    #points = hp.pix2ang(NSIDE_internal, pix, lonlat=True)
    #del pix
    #points = hp.pix2ang(NSIDE_internal, np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1,
    #                                                        hp.nside2npix(NSIDE_internal))), lonlat=True)
    #print("Allocating memory for mask1")
    #mask1_masked = np.int_(np.zeros(pix.size))
    #print("Allocating memory for mask2")
    #mask2_masked = np.int_(np.zeros(pix.size))
    print("Allocating memory for data")
    data = np.float32(np.zeros(hp.nside2npix(NSIDE_internal)))

    steps = 100 * num_thread
    divisions = np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1, steps + 1))
    count = skip * steps

    def numpy_func(_mask1_masked, _mask2_masked):
        raise AttributeError

    def func_one(mask1_masked, mask2_masked):
        return np.bitwise_and(mask1_masked, mask2_masked)

    def func_two(mask1_masked, mask2_masked):
        return np.bitwise_and(np.bitwise_not(mask1_masked), mask2_masked)

    def func_three(mask1_masked, mask2_masked):
        return np.bitwise_and(mask1_masked, np.bitwise_not(mask2_masked))

    def func_four(mask1_masked, mask2_masked):
        return np.bitwise_and(np.bitwise_not(mask1_masked), np.bitwise_not(mask2_masked))

    funcs = [func_one, func_two, func_three, func_four]

    def scope_func(i):
        nonlocal count
        #nonlocal divisions
        #nonlocal mask1
        #nonlocal mask2
        print(f"{25 * (count / steps)}%")
        count += 1
        points = hp.pix2ang(NSIDE_internal, pix[divisions[i]:divisions[i + 1]], lonlat=True)
        mask1_masked = mask1.lookup_point(*points) == 0.0
        mask2_masked = mask2.lookup_point(*points) == 0.0
        #return [i, np.int_(numpy_func(mask1_masked, mask2_masked))]
        return_func((i, np.int_(numpy_func(mask1_masked, mask2_masked))))
        return None
        # segfault (somehow)

    lock = multiprocessing.Lock()

    def return_func(result):
        nonlocal data
        nonlocal lock
        lock.acquire()
        i = result[0]
        #print(result)
        #print(result[0])
        #print(result[1])
        data[divisions[i]:divisions[i + 1]] = result[1]
        lock.release()

    for j in range(skip, 4):
        pool = ThreadPool(num_thread)
        print("Initialising pix array")
        pix = np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1, hp.nside2npix(NSIDE_internal)))
        print("Loading masks")
        mask1 = load_func(mask1_name, *args, **kwargs)
        mask2 = load_func(mask2_name, *args, **kwargs)
        print("Querying masks")
        numpy_func = funcs[j]
        #for result in pool.map(scope_func, range(steps)):
        #    i = result[0]
        #    data[divisions[i]:divisions[i + 1]] = result[1]
        pool.map_async(scope_func, range(steps))#, callback=return_func)
        pool.close()
        pool.join()
        print("Rescaling")
        pix = None
        mask1 = None
        mask2 = None
        temp = hp.ud_grade(data, NSIDE)
        print("Writing")
        if write:
            print(f"Writing to: \"./{name}_{NSIDE}_{j+1}.fits\"")
            hp.fitsfunc.write_map(f"./{name}_{NSIDE}_{j+1}.fits", temp, overwrite=True)
        print("Finished writing")


    """print("Initialising pix array")
    pix = np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1, hp.nside2npix(NSIDE_internal)))
    print("Loading masks")
    mask1 = load_func(mask1_name, *args, **kwargs)
    mask2 = load_func(mask2_name, *args, **kwargs)
    
    numpy_func = func_two
    for result in pool.map(scope_func, range(steps)):
        i = result[0]
        data[divisions[i]:divisions[i + 1]] = result[1]
    print("Rescaling")
    pix = None
    mask1 = None
    mask2 = None
    temp = hp.ud_grade(data, NSIDE)
    print("Writing")
    if write:
        print(f"Writing to: \"./{name}_{NSIDE}_2.fits\"")
        hp.fitsfunc.write_map(f"./{name}_{NSIDE}_2.fits", temp, overwrite=True)
    print("Finished writing")
    print("Initialising pix array")
    pix = np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1, hp.nside2npix(NSIDE_internal)))
    print("Loading masks")
    mask1 = load_func(mask1_name, *args, **kwargs)
    mask2 = load_func(mask2_name, *args, **kwargs)

    numpy_func = func_three
    for result in pool.map(scope_func, range(steps)):
        i = result[0]
        data[divisions[i]:divisions[i + 1]] = result[1]
    print("Rescaling")
    pix = None
    mask1 = None
    mask2 = None
    temp = hp.ud_grade(data, NSIDE)
    print("Writing")
    if write:
        print(f"Writing to: \"./{name}_{NSIDE}_3.fits\"")
        hp.fitsfunc.write_map(f"./{name}_{NSIDE}_3.fits", temp, overwrite=True)
    print("Finished writing")
    print("Initialising pix array")
    pix = np.int_(np.linspace(0, hp.nside2npix(NSIDE_internal) - 1, hp.nside2npix(NSIDE_internal)))
    print("Loading masks")
    mask1 = load_func(mask1_name, *args, **kwargs)
    mask2 = load_func(mask2_name, *args, **kwargs)

    numpy_func = func_four
    for result in pool.map(scope_func, range(steps)):
        i = result[0]
        #print(f"{75 + i/(steps * 0.04)}%")
        data[divisions[i]:divisions[i + 1]] = result[1]
    print("Rescaling")
    pix = None
    mask1 = None
    mask2 = None
    temp = hp.ud_grade(data, NSIDE)
    print("Writing")
    if write:
        print(f"Writing to: \"./{name}_{NSIDE}_4.fits\"")
        hp.fitsfunc.write_map(f"./{name}_{NSIDE}_4.fits", temp, overwrite=True)
    print("Finished writing"))"""
    return None


def run_const(data_set, mask, cat, weight_function, convert_to_mask_frac):
    data = np.array((
        (np.mean(data_set[0]),),
        (np.mean(data_set[1]),),
        (np.mean(data_set[2]),),
        (np.mean(data_set[3]),),
    ))
    binmap = ConstantBinMap()
    binmap.set_mask(mask)
    binmap.bin_catalogue(cat)
    #binmap.load_catalogue(cat)
    #output = binmap.divide_sample(mask, data, True, filter_set, a)
    output = binmap.divide_sample(mask, data)
    mixed = weight_function(*output[1:], skip_n_filter=True)
    print(f"Mixed C: {mixed[0]} +/- {mixed[1]}")
    if convert_to_mask_frac:
        final = np.array(
            [output[0][0] + (1 - output[0][0] - output[0][1]) * mixed[0], (1 - output[0][0] - output[0][1]) * mixed[1]])
    else:
        final = np.array(mixed)
    print(f"Final C: {final[0]} +/- {final[1]}")
    return final


def run_nside(n, data_set, mask, cat, weight_function, convert_to_mask_frac):
    try:
        data = np.array((
            hp.ud_grade(data_set[0], n),
            hp.ud_grade(data_set[1], n),
            hp.ud_grade(data_set[2], n),
            hp.ud_grade(data_set[3], n)
        ))
        binmap = HealpixBinMap(n)
        binmap.set_mask(mask)
        binmap.bin_catalogue(cat)
        #binmap.load_catalogue(cat)
        #output = binmap.divide_sample(mask, data, False, filter_set, a)
        output = binmap.divide_sample(mask, data)
        mixed = weight_function(*output[1:])
        print(f"Mixed {n}: {mixed[0]} +/- {mixed[1]}")
        print(output[0])
        if convert_to_mask_frac:
            final = np.array([output[0][0] * 100 + (1 - output[0][0] - output[0][1]) * mixed[0],
                              (1 - output[0][0] - output[0][1]) * mixed[1]])
        else:
            final = np.array(mixed)
        print(f"Final {n}: {final[0]} +/- {final[1]}")
    except ValueError:
        final = np.array([np.nan, np.nan])
    return final

#  Define a black and white heatmap
__cdict_bw = [(0, 0, 0), (1, 1, 1)]
bw_heatmap = matplotlib.colors.LinearSegmentedColormap.from_list("black and white", __cdict_bw, N=2)
rb_heatmap = matplotlib.colors.LinearSegmentedColormap.from_list("rb", [(1, 0, 0), (0, 0, 1)])
