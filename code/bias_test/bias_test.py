from toolkit import toolkit, weights
import argparse
import multiprocessing.pool
import numpy as np


parser = argparse.ArgumentParser()

parser.add_argument("--catalogue", default="sdss")
parser.add_argument("--save_path", default="test.npy")
parser.add_argument("--raise_dir", type=int, default=2)
parser.add_argument("--nside", type=int, default=2)
parser.add_argument("--weight_function", default="excess")
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("--iterations", type=int, default=1)
parser.add_argument("--target", type=int, default=1000)
parser.add_argument("--data_mask", default="sdss_act")
parser.add_argument("--overdensity", type=float, default=0.05)
args = parser.parse_args()

if args.catalogue == "sdss":
    random_mask = toolkit.load_mask("sdss_mask", raise_dir=args.raise_dir)
else:
    raise ValueError

if args.data_mask == "sdss_act":
    point_mask = toolkit.load_mask("act_point", raise_dir=args.raise_dir)
    temp1 = toolkit.load_mask("sdss_mask", raise_dir=args.raise_dir)
    temp1.map = np.int_(temp1.map)
    temp2 = toolkit.load_mask("act_point", raise_dir=args.raise_dir)
    temp1.map = 1 - temp1.map
    overdensity_mask = toolkit.CombinationMask(temp1, temp2, invert=True, use_and=False)
    sky_mask_frac = 0.02024328631080069
    data_set = np.float_(np.array((
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_1.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_2.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_3.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_4.fits").map
    )))
    filter_set = "n_only"
else:
    raise ValueError

pool = multiprocessing.pool.ThreadPool(processes=args.threads)
threads = []


def func():
    random_points = toolkit.gen_random_coords(args.target, random_mask).transpose()[::-1]
    bias_points = toolkit.gen_random_coords(len(random_points) * args.overdensity * sky_mask_frac * 5, overdensity_mask).transpose()[::-1]
    cat = toolkit.StarCatalogue()
    cat.lon_lat = np.append(random_points, bias_points[:int(len(random_points) * args.overdensity * sky_mask_frac)], axis=0)
    return toolkit.run_nside(n=args.nside, data_set=data_set, mask=point_mask, filter_set=filter_set, a=0, cat=cat,
                             weight_function=weights.excess_measurement, convert_to_mask_frac=False)



for i in range(args.iterations):
    threads.append(pool.apply_async(func))

#pool.join()
print(pool.close())
print(pool.join())

results = []
for thread in threads:
    results.append(thread.get())
print(results)
results = np.array(results)
np.save(args.save_path, results)
