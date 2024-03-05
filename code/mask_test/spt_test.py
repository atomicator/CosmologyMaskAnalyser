import copy
from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np
import pixell.enmap
import pixell.utils
import os
import astropy


class SPTMask(toolkit.PixellMask):
    def __init__(self, file_path):
        self.imap = pixell.enmap.read_fits(file_path)
        self.map = np.array(self.imap)
        self.corner = np.rad2deg(self.imap.box())
        self.mask_using_latlon = True
        self.lat_min = -90
        self.lon_min = -180
        self.lat_max = 90
        self.lon_max = 180

    def lookup_point(self, lon, lat):
        if not self.mask_using_latlon:
            c = astropy.coordinates.SkyCoord(lon, lat, unit="deg", frame="galactic")  # convert to galactic co ords
            lon = c.icrs.ra.degree
            lat = c.icrs.dec.degree
        lon = copy.copy(lon)
        lat = copy.copy(lat)
        #lon[lon > 180] = lon[lon > 180] - 360
        min_lat = np.min(self.corner[:, 0])
        min_lon = np.min(self.corner[:, 1])
        max_lat = np.max(self.corner[:, 0])
        max_lon = np.max(self.corner[:, 1])
        #print(self.corner)
        #print((min_lat, max_lon),
        #      (max_lat, min_lon))
        points = np.zeros(np.shape(lon))
        lon[lon > 180] -= 360
        valid_points = np.bitwise_and(np.bitwise_and(lon > min_lon, lon < max_lon),
                                      np.bitwise_and(lat > min_lat, lat < max_lat))
        #pix = np.array((np.int_(((-lat[valid_points] + max_lat) * len(self.map) / (max_lat - min_lat))),
        #                np.int_(((-lon[valid_points] + max_lon) * len(self.map[0]) / (max_lon - min_lon)))
        #                ))
        pix = np.int_(pixell.enmap.sky2pix(self.imap.shape, self.imap.wcs, np.array((lat[valid_points], lon[valid_points])) * np.pi / 180))
        # pix = np.array((np.int_(((lat[valid_points] - min_lat) * len(self.map) / (max_lat - min_lat))),
        #                np.int_(((lon[valid_points] - min_lon) * len(self.map[0]) / (max_lon - min_lon)))
        #                ))
        #pix = np.int_(pixell.enmap.sky2pix(self.imap.shape, self.imap.wcs, np.array((lon[valid_points], lat[valid_points])) * np.pi / 180))
        #lon[lon > 0] += 360
        try:
            print(np.min(pix[0]), np.max(pix[0]), np.min(pix[1]), np.max(pix[1]))
        except ValueError:
            print("Zero length array")
        print(np.shape(valid_points), np.shape(pix), np.shape(lat))
        points[valid_points] = self.map[pix[0], pix[1]]
        return points


#path = "../../data/pixel_galactic_masks"
path = "../../data/pixel_masks_new"
files = os.listdir(path)

fig = plt.figure(figsize=(8, 4), dpi=1000)
ax = fig.add_subplot(111)

masks = []

for file in files:
    mask = SPTMask(path + "/" + file)
    masks.append(mask)

final_mask = toolkit.CombinationMask(masks.pop(0), masks.pop(0), use_and=False)
while len(masks) > 0:
    final_mask = toolkit.CombinationMask(final_mask, masks.pop(0), use_and=False)
#final_mask = mask
#mask2 = SPTMask("../../data/pixel_masks/pixel_masks/ra1hdec-60_pixel_mask_ps.fits.gz")
#mask1 = SPTMask("../../data/pixel_masks/pixel_masks/ra21hdec-42.5_pixel_mask_ps.fits.gz")
#mask = toolkit.CombinationMask(mask1, mask2, use_and=False, invert=False)
#final_mask.update_angle_range(-70, -15, -80, 150)
final_mask.update_angle_range(-90, 90, 0, 360)
final_mask.set_fig_ax(fig, ax)
final_mask.plot(resolution=(4000, 2000), cbar=False)
#plt.show()

cat = toolkit.StarCatalogue("../../data/sptecs_catalog_oct919.fits", hdu=1)
cat.load_ra_dec()
cat1 = toolkit.StarCatalogue("../../data/2500d_cluster_sample_Bocquet19.fits", hdu=1)
cat1.load_ra_dec()
cat.lon_lat = np.append(cat.lon_lat, cat1.lon_lat, axis=0)
#plt.scatter(*(cat.lon_lat.transpose()), marker="+")
#plt.ylim(-90, 90)
#plt.imshow(mask.map)
#plt.show()
print("done")

data = final_mask.lookup_point(*cat.lon_lat.transpose())
print(np.sum(data), len(data))

plt.scatter(*(cat.lon_lat[data == 1].transpose()), marker="+", label="Not masked")
plt.scatter(*(cat.lon_lat[data == 0].transpose()), marker="+", label="Masked")
plt.legend()
plt.ylim(-90, 90)
plt.xlim(0, 360)
plt.savefig("test.png")
plt.show()
