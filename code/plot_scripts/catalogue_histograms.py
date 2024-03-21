import numpy as np
import matplotlib.pyplot as plt
from toolkit import toolkit

toolkit.plt_use_tex()


def test(*args):
    sdss_data.append(np.max(args))
    return True


bins = np.linspace(-500, 500, 51)

sdss_cat = toolkit.load_catalogue("sdss")
sdss_data = []
sdss_cat.load_with_selection(test, ("LAMBDA_CHISQ",), True)
sdss_data = np.array(sdss_data)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(sdss_data > bins[i], sdss_data < bins[i + 1]))) / (len(sdss_data) / 100)
plt.stairs(sdss_height, bins, label="SDSS")

sdss_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1, table=True)
sdss_data = []
sdss_cat.load_with_selection(test, ("RM_LAMBDA", "RMDESY3_LAMBDA_CHISQ", "CAMIRA_N_mem"), True)
sdss_data = np.array(sdss_data)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(sdss_data > bins[i], sdss_data < bins[i + 1]))) / (len(sdss_data) / 100)
plt.stairs(sdss_height, bins, label="ACT")
print(np.sum(sdss_height))
print(np.min(sdss_data), np.max(sdss_data))
print(np.sum(sdss_data < 0))

plt.title("The cluster count of the catalogues, as a function of redshift")
#plt.yscale("log")
plt.ylabel("Percentage")
plt.xlabel("Richness")

plt.show()
