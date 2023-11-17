from toolkit import toolkit
import matplotlib.pyplot as plt

fig = plt.figure()
# ax = fig.add_subplot(111, projection="mollweide")
ax = fig.add_subplot(111)

# load the catalogue
cat = toolkit.StarCatalogue("../../data/sdss_catalogue.fits", hdu=1)

# load the masks
#mask = toolkit.PixellMask("../../data/ACT_mask.fits", step=10, hdu=1)
#mask = toolkit.HealpyMask("../../data/planck_galactic_mask.fits")
mask = toolkit.HealpyMask("../../data/planck_point_mask.fits")


sdss_mask = toolkit.HealpyMask("../../data/redmapper_dr8_public_v6.3_zmask.fits", mask_using_latlon=False, hdu=1, partial=True)

# Make the values of the masked / not masked regions more normal - use the convention 1 = masked
sdss_mask.map[sdss_mask.map > 0.4] = 1.0
sdss_mask.map[sdss_mask.map < 0.3] = 0

sdss_mask.set_fig_ax(fig, ax)
sdss_mask.plot()
plt.title("SDSS mask in RA-DEC")
#plt.savefig("../../graphs/SDSS_mask.png", dpi=1e3)
plt.show()
#sdss_mask = toolkit.StarCatalogue("../../data/redmapper_dr8_public_v6.3_zmask.fits", hdu=1)

#print(sdss_mask.data.field("HPIX"))

'''
cat.set_fig_ax(fig, ax)
cat.set_cat_mask(sdss_mask)
#cat.load_ra_dec()
cat.load_lon_lat()
cat.plot_heatmap(256, cmap="rainbow", resolution=(1000, 2000), cbar=True, cbar_label="Clusters per "
                                                                                     "square degree")

mask.set_fig_ax(fig, ax)
mask.plot(alpha=0.35, cbar=False, cmap=toolkit.bw_heatmap, clear=False)
plt.xlabel("Galactic Longitude")
plt.ylabel("Galactic Latitude")
plt.title("SDSS clusters and the Planck point mask")
plt.savefig("../../graphs/PlanckPointMaskSDSSData.png", dpi=1e3)
plt.show()
'''