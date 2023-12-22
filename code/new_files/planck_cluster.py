from toolkit import toolkit
import matplotlib.pyplot as plt


def planck_cat_selection_func(snr, pipe_det, cosmo, min_snr, pipe_det_val):
    return (snr > min_snr) and (pipe_det == pipe_det_val) and cosmo


fig = plt.figure()
#ax = fig.add_subplot(111, projection="mollweide")
ax = fig.add_subplot(111)

# load the catalogue
cat = toolkit.StarCatalogue("../../data/HFI_PCCS_SZ-union_R2.08.fits", hdu=1)

# Filter the data
min_snr = 6
pipe_det_value = 111

cat.load_with_selection(planck_cat_selection_func, ["SNR", "PIPE_DET", "COSMO"], True, min_snr, pipe_det_value)

print(len(cat.lon_lat))

# load the mask
#mask = toolkit.PixellMask("../../data/ACT_mask.fits", step=2, hdu=1)
#mask = toolkit.HealpyMask("../../data/planck_galactic_mask.fits")
mask = toolkit.load_mask("planck_modified_galactic")
mask2 = toolkit.load_mask("planck_galactic")
mask.map = mask.map + 2 * mask2.map

#cat.set_fig_ax(fig, ax)
#cat.plot_scatter(marker="+", color="red", label=r"$\sigma$ > 6")


mask.set_fig_ax(fig, ax)
mask.plot(alpha=1, cbar=True, cmap="rainbow", clear=False)
plt.xlabel("Galactic Longitude")
plt.ylabel("Galactic Latitude")
plt.title("A heatmap of the stars in the Planck cluster against the ACT mask")
#plt.savefig("../../graphs/ACTMaskPlanckData.png", dpi=1e3)
plt.show()
