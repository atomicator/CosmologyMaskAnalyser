import matplotlib.pyplot as plt
from toolkit import toolkit, weights
import argparse
import multiprocessing.pool
import numpy as np
import healpy as hp
import random


parser = argparse.ArgumentParser()

parser.add_argument("--catalogue", default="sdss")
parser.add_argument("--save_path", default="test.npy")
parser.add_argument("--raise_dir", type=int, default=2)
parser.add_argument("--nside", type=int, default=2)
parser.add_argument("--weight_function", default="excess")
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("--iterations", type=int, default=1)
parser.add_argument("--target", type=int, default=400000)
parser.add_argument("--data_mask", default="sdss_act")
parser.add_argument("--overdensity", type=float, default=0.05)
args = parser.parse_args()


def test_function():
    global data_set
    sky_mask_frac = 0.015567657065637092
    random_points = toolkit.gen_random_coords(args.target, random_mask)[::-1].transpose()[::-1]
    bias_points = toolkit.gen_random_coords(len(random_points) * args.overdensity * sky_mask_frac * 5, overdensity_mask)[::-1].transpose()[::-1]
    #print(np.shape(random_points), np.shape(bias_points))
    #print(np.min(random_points[:, 0]), np.max(random_points[:, 0]))
    #print(np.min(random_points[:, 1]), np.max(random_points[:, 1]))
    #print(np.min(bias_points[0]), np.max(bias_points[0]))
    #print(np.min(bias_points[1]), np.max(bias_points[1]))
    cat = toolkit.StarCatalogue()
    cat.lon_lat = np.append(random_points, bias_points[:int(len(random_points) * args.overdensity * sky_mask_frac)], axis=0)
    #print(np.shape(cat.lon_lat))
    #print(np.sum(random_mask.lookup_point(*cat.lon_lat.transpose()[::1])) / len(cat.lon_lat))
    #print(np.sum(point_mask.lookup_point(*cat.lon_lat.transpose()[::1])) / len(cat.lon_lat))
    #print(np.sum(point_mask.lookup_point(*bias_points.transpose()[::1])) / len(bias_points))
    #print(np.min(cat.lon_lat[0]), np.max(cat.lon_lat[0]))
    #print(np.min(cat.lon_lat[1]), np.max(cat.lon_lat[1]))
    #cat.lon_lat = cat.lon_lat.transpose()
    try:
        data = np.array((
            hp.ud_grade(data_set[0], args.nside),
            hp.ud_grade(data_set[1], args.nside),
            hp.ud_grade(data_set[2], args.nside),
            hp.ud_grade(data_set[3], args.nside)
        ))
        binmap = toolkit.HealpixBinMap(args.nside)
        binmap.bin_catalogue(cat)
        binmap.load_catalogue(cat)
        output = binmap.divide_sample(point_mask, data, False, filter_set, 0)
        mixed = weights.excess_measurement(*output[1:])
        print(f"Mixed {args.nside}: {mixed[0]} +/- {mixed[1]}")
        print(output[0])
        final = np.array(mixed)
        print(f"Final {args.nside}: {final[0]} +/- {final[1]}")
    except ValueError:
        final = np.array([np.NaN, np.NaN])
    # results.append(final)
    return final


if args.catalogue == "sdss":
    random_mask = toolkit.load_mask("sdss_mask", raise_dir=args.raise_dir)
else:
    raise ValueError

if args.data_mask == "sdss_act":
    point_mask = toolkit.load_mask("act_point", raise_dir=args.raise_dir)
    temp1 = toolkit.load_mask("sdss_mask", raise_dir=args.raise_dir)
    temp1.map = np.int_(temp1.map)
    #temp2 = toolkit.load_mask("act_point", raise_dir=args.raise_dir)
    temp1.map = 1 - temp1.map
    overdensity_mask = toolkit.CombinationMask(temp1, point_mask, invert=True, use_and=False)
    sky_mask_frac = 0.009859934289099422
    data_set = np.float_(np.array((
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_1.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_2.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_3.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_4.fits").map
    )))
    filter_set = "n_only"
else:
    raise ValueError

pool = multiprocessing.pool.Pool(processes=args.threads, initializer=random.seed)
threads = []


for i in range(args.iterations):
    threads.append(pool.apply_async(test_function))

#pool.join()
#print(pool.close())
#print(pool.join())

results = []
for thread in threads:
    results.append(thread.get())
print(results)
results = np.array(results)
np.save(args.save_path, results)
