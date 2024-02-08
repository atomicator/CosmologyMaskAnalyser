import numpy
import numpy as np
import matplotlib.pyplot as plt
import pixell
import pixell.enmap as enmap
import pixell.enplot
import healpy as hp
import astropy.cosmology
import astropy.units
import matplotlib.projections

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})

#dec_min = 170
#dec_max = 180
#ra_min = 170
#ra_max = 180

dec_min_test = -90
dec_max_test = 90
ra_min_test = -180
ra_max_test = 180

#box = np.array([[ra_min, ra_max], [dec_min, dec_max]]) * pixell.utils.degree

#print(box)

imap = enmap.read_fits("../../data/ACT_mask.fits", hdu=1)
#imap = enmap.read_map("./data/ACT_mask.fits", hdu=1)
#data = astropy.io.fits.open("./data/planck_galactic_mask.fits", hdu=2)[2]

#print(data)

#print(imap.shape, imap.dtype)
print(imap)
x = len(imap)
y = len(imap[0])
print(x, y)

new_data = imap[::6, ::6]
print(new_data)

new_data = numpy.array(new_data)
#fig, ax = plt.subplots(figsize=(1, 1))
fig = plt.figure()
ax = fig.add_subplot(111, projection="mollweide")
#ax = fig.add_subplot(111)
#lat = np.linspace(-np.pi/2, np.pi/2, len(new_data) + 1)
#lon = np.linspace(-np.pi/2, np.pi/2, len(new_data[0]) + 1)
lat = np.linspace(np.pi * (dec_min_test / 180), np.pi * (dec_max_test / 180), len(new_data) + 1)
lon = np.linspace(np.pi * (ra_min_test / 180), np.pi * (ra_max_test / 180), len(new_data[0]) + 1)
Lon, Lat = np.meshgrid(lon, lat)
#ax = fig.add_subplot(111)
#c = ax.imshow(new_data, extent=(- np.pi / 2, np.pi / 2, -np.pi, np.pi))
#c = ax.pcolormesh(new_data, extent=(- np.pi / 2, np.pi / 2, -np.pi, np.pi))
print(Lon.shape, Lat.shape, new_data.shape)
#c = ax.pcolormesh(Lat, Lon, new_data)
c = ax.pcolormesh(Lon, Lat, new_data)
#ax.colorbar(orientation="horizontal")
fig.colorbar(c, orientation="horizontal")
plt.title("Part of the ACT mask")
#plt.savefig("./graphs/act_mask/mollweide.png")
fig.show()


#temp = pixell.enplot.plot(new_data)
#pixell.enplot.show(pixell.enplot.plot(new_data))
