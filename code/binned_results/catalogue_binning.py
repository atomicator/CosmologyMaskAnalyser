from toolkit import toolkit
import matplotlib.pyplot as plt
import numpy as np

print("Loading catalogue")

cat = toolkit.load_catalogue("sdss", raise_dir=2)

print("Binning catalogue")

bin = toolkit.HealpixBinMap(8)
#bin = toolkit.BinaryMap()
bin.load_catalogue(cat)

#print(bin.map)

print("Loading mask")

map1 = toolkit.load_mask("planck_point")
map2 = toolkit.load_mask("sdss_mask")

print("Weighting bins")

bin.calc_weighted_map(map2)

print(bin.weighted_map)

print("Calculating fractions")

#old_data = toolkit.fraction_masked_pair(map1, map2, weight_map=None, n=10000, ram_limited=True)
data = toolkit.fraction_masked_pair(map1, map2, weight_map=bin, n=10000, ram_limited=True)

#print(old_data)
#print(old_data[1] / (old_data[2] + old_data[1]))

print(data)
print(100 * data[1] / (data[2] + data[1]))

map = toolkit.load_mask("planck_point")
map.map = bin.weighted_map
map.NSIDE = bin.NSIDE

fig = plt.figure()
ax = fig.add_subplot(111)

map.set_fig_ax(fig, ax)

map.plot(cbar=True, show=True)

#data = map.map[np.bitwise_and(map.map > 1200, map.map < 3000)]
data = map.map[bin.mask_fraction_map > 0.9]
mean = np.mean(data)
std = np.std(data)

plt.hist(map.map[map.map > 0], bins=40)
plt.title(f"{mean} $\pm$ {std}")
plt.show()
print(np.sqrt(mean))
