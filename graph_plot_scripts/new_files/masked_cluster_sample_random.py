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

print(mask_names)
print(f"Bootstrap iterations: {args.bootstrap_iterations}")
print(f"Number of simulated clusters: {args.num_clusters}")

masks = [toolkit.load_mask(mask_names[0], args.path_raise), toolkit.load_mask(mask_names[1], args.path_raise)]
mask_frac = toolkit.fraction_masked_pair(masks[0], masks[1])
unmasked_fraction_exact = mask_frac[2] / (mask_frac[0] + mask_frac[2])
print(f"The exact unmasked fraction is {unmasked_fraction_exact}")

# Generate a random sample that is allowed through mask[0]
co_ords = toolkit.gen_random_coords(args.num_clusters, masks[0])
data = masks[1].lookup_point(*co_ords[::-1])

unmasked_fraction = np.sum(data) / len(data)

print(f"The unmasked fraction is {unmasked_fraction}")

mean_estimates = np.array(toolkit.bootstrap(data, args.bootstrap_iterations)) / len(data)

print(f"Final results: {mean_estimates[0]} +/- {np.std(mean_estimates)}")

if args.save_path:
    np.savetxt("../" * args.path_raise + args.save_path, mean_estimates, delimiter=",")

