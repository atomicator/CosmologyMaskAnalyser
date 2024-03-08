from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np

toolkit.plt_use_tex()
fig = plt.figure()
ax = fig.add_subplot(111)

cat = toolkit.load_catalogue("sdss")
#cat = toolkit.StarCatalogue("../binned_results/act_bias2.fits", table=True)
cat.load_lon_lat()
cat.set_fig_ax(fig, ax)

#act_mask = toolkit.load_mask("act_point")
#act_mask = toolkit.load_mask("sdss_mask")
#act_mask = toolkit .load_mask("act")
#act_mask.set_fig_ax(fig, ax)

mask = toolkit.load_mask("planck_modified_point")
mask.set_fig_ax(fig, ax)
mask.plot(cmap="bwr_r", clear=False)

mask = toolkit.load_mask("planck_modified_galactic")
mask.set_fig_ax(fig, ax)
mask.plot(cmap="bwr", clear=False)

#cat.plot_heatmap(128, cmap="rainbow", resolution=(1000, 2000), cbar=True, cbar_label="Clusters per "
#                                                                                    "square degree")
#cat.plot_scatter(marker="+")
#act_mask.plot(clear=False, alpha=0.35, cmap=toolkit.bw_heatmap)
plt.ylabel("Latitude")
plt.xlabel("Longitude")
plt.title(r"A plot of the Point and Galactic components of the Planck mask")
plt.savefig("planck_point_galactic.png", dpi=1000)
plt.show()
