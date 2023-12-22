import matplotlib.pyplot as plt
from toolkit import toolkit
import numpy as np
import healpy as hp

toolkit.plt_use_tex()

mask = toolkit.load_mask("planck_galactic")
cat = toolkit.load_catalogue("sdss")
cat.load_lon_lat()

fig = plt.figure()
ax = fig.add_subplot(111)

f = (1 - np.sum(mask.lookup_point(*cat.lon_lat.transpose())) / len(cat.lon_lat)) * 100
sample_error = 0.017

print(f)
NSIDES = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
#NSIDES = [512]

results = []
for n in NSIDES:
    binmap = toolkit.HealpixBinMap(n)
    binmap.bin_catalogue(cat)
    binmap.load_catalogue(cat)
    data = np.array((
        hp.pixelfunc.ud_grade(toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_galactic_512_1.fits").map,
                              n),
        hp.pixelfunc.ud_grade(toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_galactic_512_2.fits").map,
                              n),
        hp.pixelfunc.ud_grade(toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_galactic_512_3.fits").map,
                              n),
        hp.pixelfunc.ud_grade(toolkit.HealpyMask("../../data/cached_results/sdss_mask_planck_galactic_512_4.fits").map,
                              n)
    ))
    results.append(binmap.calc_masked_fraction_new(mask, data))

results = np.array(results).transpose() * 100

#ax.plot(NSIDES, np.abs(f - results[0]), label="Galactic")
#ax.errorbar(NSIDES, np.abs(f - results[0]), results[1], fmt="none", capsize=1, capthick=1, ecolor="k")

ax.plot(NSIDES, results[2], label="Galactic")
ax.errorbar(NSIDES, results[2], results[3], fmt="none", capsize=1, capthick=1, ecolor="k")

#ax.plot(NSIDES, np.ones(len(NSIDES)) * sample_error, label="Sample error", linestyle="dashed", color="black")
ax.plot(NSIDES, np.ones(len(NSIDES)) * 100, label="Sample error", linestyle="dashed", color="black")

ax.set_xscale("log", base=2)
ax.legend()
#plt.ylim(0, 0.08)
ax.set_xticks(NSIDES, NSIDES)

ax.set_xlabel("NSIDE")
ax.set_ylabel(r"Absolute masked fraction difference $(\%)$")
ax.set_title("The effects of weighting the masked fraction using a binning algorithm to\nestimate the ratio")
plt.show()

print(f, results)
