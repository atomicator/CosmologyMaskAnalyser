import numpy as np
import astropy

class SquareMaskFilter:
    def __init__(self, lat_range=(-90, 90), lon_range=(0, 360), lonlat=True, lon_shift=None):
        self.lat_range = lat_range
        self.lon_range = lon_range
        self.lonlat = lonlat
        self.lon_shift = lon_shift

    def lookup_point(self, lon, lat):
        if self.lon_shift:
            lon += self.lon_shift
            lon[lon > 180] -= 360
            lon[lon < -180] += 360
        if not self.lonlat:  # convert to RA DEC
            c = astropy.coordinates.SkyCoord(lon, lat, unit="deg", frame="galactic")  # convert to galactic co ords
            lon = c.icrs.ra.degree
            lat = c.icrs.dec.degree
        lon[lon < 0] = lon[lon < 0] + 360
        return np.int_(np.bitwise_not(np.bitwise_and(np.bitwise_and(self.lon_range[0] < lon, self.lon_range[1] > lon),
                                      np.bitwise_and(self.lat_range[0] < lat, self.lat_range[1] > lat))))