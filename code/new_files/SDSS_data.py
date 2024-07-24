from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
#ax = fig.add_subplot(111, projection="mollweide")
ax = fig.add_subplot(111)

# load the catalogue
#cat = toolkit.ClusterCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1)
cat = toolkit.StarCatalogue("../binned_results/test.fits", hdu=1, table=True)
cat.load_lon_lat()
data = np.load("../../data/random_catalogue_400k_inigo.npy")
print(data)
cat.lon_lat = np.array([data[0] * 180 / np.pi, 90 - data[1] * 180/np.pi]).transpose()
print(cat.lon_lat)

# load the masks
#mask = toolkit.PixellMask("../../data/ACT_mask.fits", step=10, hdu=1, mask_using_latlon=False)

#mask = toolkit.load_mask("sdss_planck_point_only")
#galactic_mask = toolkit.load_mask("planck_galactic")
#galactic_mask.map = 1 - galactic_mask.map
#planck = toolkit.load_mask("planck_modified_point")
#act_mask = toolkit.load_mask("act")
#mask = toolkit.CombinationMask(galactic_mask, act_mask, invert=True)
#mask = toolkit.load_mask("planck_modified_point")
#mask = sdss_mask
#mask = toolkit.HealpyMask("../../data/planck_galactic_mask.fits")
#mask = toolkit.HealpyMask("../../data/planck_point_mask.fits")
#mask = toolkit.load_mask("planck_modified_galactic")
#mask = toolkit.load_mask("planck_modified_galactic")
#mask = toolkit.load_mask("act")

#sdss_mask = toolkit.HealpyMask("../../data/redmapper_dr8_public_v6.3_zmask.fits", mask_using_latlon=False, hdu=1, partial=True)

# Make the values of the masked / not masked regions more normal - use the convention 1 = masked
#sdss_mask.map[sdss_mask.map > 0.4] = 1.0
#sdss_mask.map[sdss_mask.map < 0.3] = 0

#sdss_mask.set_fig_ax(fig, ax)
#sdss_mask.plot(cbar=True)
#mask.set_fig_ax(fig, ax)
#mask.plot(cbar=False, label="Galactic", cmap="bwr", show=False)
#plt.title("SDSS mask in GLAT-GLON")
#plt.savefig("../../graphs/SDSS_mask.png", dpi=1e3)
#plt.show()
#sdss_mask = toolkit.ClusterCatalogue("../../data/redmapper_dr8_public_v6.3_zmask.fits", hdu=1)

#print(sdss_mask.data.field("HPIX"))


cat.set_fig_ax(fig, ax)
#cat.set_cat_mask(m)
#cat.load_ra_dec()
#cat.load_lon_lat()
cat.plot_heatmap(32, cmap="rainbow", resolution=(1000, 2000), cbar=True, cbar_label="Clusters per "
                                                                                     "square degree")
#mask.set_fig_ax(fig, ax)
#mask.plot(alpha=0.35, cbar=True, cmap="plasma", clear=False)
#plt.xlabel("Galactic Longitude")
#plt.ylabel("Galactic Latitude")
#plt.title("SDSS clusters and the Planck point mask")
#plt.savefig("../../graphs/PlanckPointMaskSDSSData.png", dpi=1e3)
plt.show()

#print(mask.lookup_point(0, 0))
