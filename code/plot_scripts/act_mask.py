from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np

toolkit.plt_use_tex()

fig = plt.figure()
ax = fig.add_subplot(111)

#galactic_mask = toolkit.load_mask("planck_galactic")
#act_mask = toolkit.load_mask("act")
act_mask = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1, invert=False)
#mask = toolkit.CombinationMask(galactic_mask, act_mask, invert=True)

#mask.set_fig_ax(fig, ax)
#mask.plot(cbar=False, label="Galactic", cmap="bwr", show=False, clear=False)

#galactic_mask.map = 1 - galactic_mask.map
#mask = toolkit.CombinationMask(galactic_mask, act_mask, invert=True)
act_mask.map = np.ones(np.shape(act_mask.map))
act_mask.set_fig_ax(fig, ax)
act_mask.plot(cbar=False, label="Galactic", cmap="bwr_r", show=False, clear=False)

plt.ylabel("Latitude")
plt.xlabel("Longitude")
plt.title(r"The \enquote{Galactic} and \enquote{Point} components of the ACT mask")

plt.savefig("act_mask_lat_lon.png")
plt.show()
