from toolkit import toolkit
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("--mask_one", choices=["sdss_mask", "planck_point", "planck_galactic"],
                    help="The first mask to use, all points will be allowed through this mask", default="sdss_mask")
parser.add_argument("--mask_two", choices=["sdss_mask", "planck_point", "planck_galactic"],
                    help="The second mask to use", default="planck_galactic")
parser.add_argument("--path_raise", type=int, default=2,
                    help="The number of times ../ should be added to file paths - 2 on laptop, 0 on cluster")
parser.add_argument("--bootstrap_iterations", type=int, default=5)
parser.add_argument("--save_path", help="The path to save the output")
parser.add_argument("--num_clusters", type=int, default=500000,
                    help="The number of clusters (roughly) to generate")  # Laptop runs out of memory for 50 mil,
# shouldn't be a problem on cluster

args = parser.parse_args()
mask_names = [args.mask_one, args.mask_two]

masks = [toolkit.load_mask(mask_names[0], args.path_raise), toolkit.load_mask(mask_names[1], args.path_raise)]

masks[0] = toolkit.load_mask("planck_galactic")

print("Number 1: using the mask directly")
print(np.sum(masks[0].map) / len(masks[0].map))

print("Number 2: sampling the map")
n = 5000
x = np.linspace(-180, 180, n)
y = np.linspace(-90, 90, n)
x, y = np.meshgrid(x, y)
data = masks[0].lookup_point(x, y)
plt.pcolormesh(x, y, data, cmap="rainbow")
plt.show()
print(data.shape)
print(np.sum(data[0]))
total = 0
unmasked = 0
for i in range(len(y)):
    s = np.sin(np.pi * (i + 0.5) / len(y))
    total += n * s
    unmasked += np.sum(data[i]) * s
print(unmasked, total)
print(unmasked / total)
