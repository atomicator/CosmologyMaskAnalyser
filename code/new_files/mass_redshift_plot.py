from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np


def planck_cat_selection_func(snr, pipe_det, cosmo, redshift, mass):
    global min_snr, max_snr, pipe_det_val, m, z
    x = False
    if (max_snr > snr > min_snr) and (pipe_det == pipe_det_val) and cosmo:
        z.append(redshift)
        m.append(mass)
        x = True
    return x


max_snr = 1e10
min_snr = 10
pipe_det_val = 111
# max snr is 48.9
cat = toolkit.StarCatalogue("../../data/HFI_PCCS_SZ-union_R2.08.fits", hdu=1)

labels = [r"$\sigma > 10$", r"$\sigma > 9$", r"$\sigma > 8$", r"$\sigma > 7$", r"$\sigma > 6$"]

min_z = 1e100
max_z = 0

while min_snr >= 6:
    z = []
    m = []
    cat.load_with_selection(planck_cat_selection_func, ["SNR", "PIPE_DET", "COSMO", "REDSHIFT", "MSZ"], True)
    plt.scatter(z, m, marker="+", s=10, label=labels[10 - min_snr])
    max_snr = min_snr
    min_snr -= 1
    if np.max(z) > max_z:
        max_z = np.max(z)
    if 0 < np.min(z) < min_z:
        min_z = np.min(z)

print(min_z, max_z)

plt.xlabel("Redshift")
plt.ylabel(r"Mass $\left( 10^{14} M_{\odot} \right)$")
plt.title(r"Cluster Mass of the Planck detections as a function of redshift and signal to noise $(\sigma)$")
plt.xlim(0, 1)
plt.legend()
plt.savefig("../../graphs/RedshiftSNRPlanck.png", dpi=1e3)
plt.show()
