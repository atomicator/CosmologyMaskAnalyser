from toolkit import toolkit
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("--mask", choices=["sdss_mask", "planck_point", "planck_galactic"],
                    help="The mask to use", default="planck_galactic")
parser.add_argument("--path_raise", type=int, default=2,
                    help="The number of times ../ should be added to file paths - 2 on laptop, 0 on cluster")
parser.add_argument("--bootstrap_iterations", type=int, default=5)
parser.add_argument("--save_path", help="The path to save the output")
parser.add_argument("--catalogue", help="The cluster map to use", default="data/sdss_catalogue.fits")

args = parser.parse_args()

print(args.mask)
print(args.catalogue)
print(f"Bootstrap iterations: {args.bootstrap_iterations}")

mask = toolkit.load_mask(args.mask, args.path_raise)
cat = toolkit.StarCatalogue("../" * args.path_raise + args.catalogue, hdu=1)

cat.load_lon_lat()

print(len(cat.lon_lat))

data = mask.lookup_point(*cat.lon_lat.transpose())

print(min(data))

print(np.sum(data) / len(data))

mean_estimates = np.array(toolkit.bootstrap(data, args.bootstrap_iterations)) / len(data)

print(f"Final results: {mean_estimates[0]} +/- {np.std(mean_estimates)}")

if args.save_path:
    np.savetxt("../" * args.path_raise + args.save_path, mean_estimates, delimiter=",")
