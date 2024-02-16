from toolkit import toolkit
import numpy as np
import matplotlib.pyplot as plt

mmf1 = toolkit.HealpyMask("../../data/HFI_PCCS_SZ-selfunc-MMF1-cosmolog_R2.08.fits")
mmf3 = toolkit.HealpyMask("../../data/HFI_PCCS_SZ-selfunc-MMF3-cosmolog_R2.08.fits")
ps = toolkit.HealpyMask("../../data/HFI_PCCS_SZ-selfunc-PwS-cosmolog_R2.08.fits")
total = toolkit.HealpyMask("../../data/HFI_PCCS_SZ-selfunc-inter-cosmo_2.02.fits")
galactic = mmf3.clone()
point = mmf3.clone()

fig = plt.figure()
ax = fig.add_subplot(111)

error = galactic.clone()

mmf1.set_fig_ax(fig, ax)
mmf3.set_fig_ax(fig, ax)
ps.set_fig_ax(fig, ax)
total.set_fig_ax(fig, ax)
point.set_fig_ax(fig, ax)
galactic.set_fig_ax(fig, ax)

galactic.map = np.float_((mmf3.map + mmf1.map + ps.map + total.map) == 0)
error.set_fig_ax(fig, ax)

print(np.min(mmf3.map), np.max(mmf3.map))
point.map = np.float_((mmf3.map * 2 + total.map) != 2)
error.map = np.float_((mmf3.map * 2 + total.map) != 1)
galactic.map = np.float_((mmf3.map * 2 + total.map) != 0)

#galactic.plot(cmap="bwr", clear=False)
point.plot(cmap="rainbow_r", clear=False)
error.plot(cmap="bwr_r", clear=False)

plt.show()
