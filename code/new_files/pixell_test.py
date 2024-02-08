from toolkit import toolkit
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

mask = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1, step=1, mask_using_latlon=False)
mask.set_fig_ax(fig, ax)
mask.plot(title="Toolkit plot test", show=True, cbar=True, resolution=(10000, 5000))
#print(mask.calc_exact_unmasked_fraction())
#print(mask.simulate_unmasked_fraction(1e6))
