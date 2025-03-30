from toolkit import toolkit, data
import matplotlib.pyplot as plt
import numpy as np

toolkit.plt_use_tex()
fig = plt.figure()
ax = fig.add_subplot(111)


def planck_cat_selection_func(snr, pipe_det, cosmo, min_snr, pipe_det_val):
    return (snr > min_snr) and (pipe_det == pipe_det_val) and cosmo


#cat = toolkit.load_catalogue("sdss")
#cat = toolkit.StarCatalogue("../binned_results/sdss_random_400k.fits", table=True)
#cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1)
cat = toolkit.ClusterCatalogue("../../data/HFI_PCCS_SZ-union_R2.08.fits", hdu=1)
#cat.load_lon_lat()
min_snr = 6
pipe_det_value = 111
cat.load_data(True, planck_cat_selection_func, ["SNR", "PIPE_DET", "COSMO"], min_snr, pipe_det_value)
cat.set_fig_ax(fig, ax)

#act_mask = toolkit.load_mask("act_point")
#act_mask = toolkit.load_mask("sdss_mask")
#act_mask = toolkit .load_mask("act")
#act_mask.set_fig_ax(fig, ax)

#mask = toolkit.load_mask("planck_modified_point")
#mask.set_fig_ax(fig, ax)
#mask.plot(cmap="bwr_r", clear=False)

#mask = toolkit.load_mask("planck_modified_galactic")
#act_mask = toolkit.load_mask("act")
act_mask = data.load_mask("planck_modified_total")
act_mask.set_fig_ax(fig, ax)
#mask.plot(cmap="bwr", clear=False)

#cat.plot_heatmap(128, cmap="rainbow", resolution=(1000, 2000), cbar=True, cbar_label="Clusters per "
#                                                                                    "square degree")

cat.plot_scatter(marker=".", s=1)
act_mask.plot(clear=False, alpha=0.35, cmap=toolkit.bw_heatmap)
plt.ylabel("Latitude")
plt.xlabel("Longitude")
plt.title(r"The Planck catalogue and mask")
plt.savefig("planck_cat_mask.png", dpi=1000)
plt.show()
