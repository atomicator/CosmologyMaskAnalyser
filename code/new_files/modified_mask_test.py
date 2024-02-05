import numpy as np
from toolkit import toolkit
import matplotlib.pyplot as plt

toolkit.plt_use_tex()

fig = plt.figure()
ax = fig.add_subplot(111)

resolution = (1000, 2000)

galactic = toolkit.load_mask("planck_modified_galactic")
#point = toolkit.load_mask("planck_point")
point = toolkit.load_mask("planck_modified_point")

galactic.set_fig_ax(fig, ax)
point.set_fig_ax(fig, ax)

#temp.map = 1 - np.float_(data == 0)
galactic.plot(show=False, cbar=False, cmap="bwr", resolution=resolution, clear=False)

#temp.map = 1 - np.float_(data == 2)
point.plot(show=False, cbar=False, cmap="bwr_r", resolution=resolution, clear=False)

plt.title("A comparison of the two masks")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.savefig("../../graphs/comparison.png", dpi=1e3)

plt.show()
