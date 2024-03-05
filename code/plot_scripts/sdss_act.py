from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np

toolkit.plt_use_tex()
fig = plt.figure()
ax = fig.add_subplot(111)

#cat = toolkit.load_catalogue("sdss")
cat = toolkit.StarCatalogue("../binned_results/act_bias2.fits", table=True)
cat.load_lon_lat()
cat.set_fig_ax(fig, ax)

#act_mask = toolkit.load_mask("act_point")
#act_mask = toolkit.load_mask("sdss_mask")
act_mask = toolkit .load_mask("act")
act_mask.set_fig_ax(fig, ax)

#cat.plot_heatmap(32, cmap="rainbow", resolution=(1000, 2000), cbar=True, cbar_label="Clusters per "
#                                                                                    "square degree")
cat.plot_scatter(marker="+")
act_mask.plot(clear=False, alpha=0.35, cmap=toolkit.bw_heatmap)
plt.ylabel("Latitude")
plt.xlabel("Longitude")
plt.title(r"The SDSS clusters and the ACT mask")
plt.savefig("SDSS_ACT.png")
plt.show()
