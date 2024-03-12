import matplotlib.pyplot as plt
from toolkit import toolkit, weights
import argparse
import multiprocessing.pool
import numpy as np
import healpy as hp
import random
import multiprocessing.managers


parser = argparse.ArgumentParser()

parser.add_argument("--catalogue", default="sdss")
parser.add_argument("--save_path", default="test.npy")
parser.add_argument("--raise_dir", type=int, default=2)
parser.add_argument("--weight_function", default="excess")
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("--iterations", type=int, default=1)
parser.add_argument("--target", type=int, default=400000)
parser.add_argument("--data_mask", default="sdss_act")
parser.add_argument("--overdensity", type=float, default=0.05)
args = parser.parse_args()

NSIDES = [1, 2, 4, 8, 16, 32, 64, 128]

multiprocessing.managers.BaseManager.register("HealpyMask", toolkit.HealpyMask)
multiprocessing.managers.BaseManager.register("PixellMask", toolkit.PixellMask)
multiprocessing.managers.BaseManager.register("CombinationMask", toolkit.CombinationMask)


def test_function():
    global data_set
    random_points = toolkit.gen_random_coords(args.target, random_mask)[::-1].transpose()
    #bias_points = toolkit.gen_random_coords(len(random_points) * args.overdensity * sky_mask_frac * 5, overdensity_mask)[::-1].transpose()
    bias_points = toolkit.gen_random_coords(args.target * args.overdensity, random_mask)[::-1].transpose()
    bias_points = bias_points[point_mask.lookup_point(*bias_points.transpose()) == 0]
    cat = toolkit.StarCatalogue()
    cat.lon_lat = np.append(random_points, bias_points[:int(len(random_points) * args.overdensity * sky_mask_frac)], axis=0)
    temp = []
    data = np.array((
        np.mean(data_set[0]),
        np.mean(data_set[1]),
        np.mean(data_set[2]),
        np.mean(data_set[3])
    ))
    binmap = toolkit.ConstantMap()
    binmap.bin_catalogue(cat)
    binmap.load_catalogue(cat)
    output = binmap.divide_sample(point_mask, data, True, filter_set, 0)
    mixed = weights.excess_measurement(*output[1:], skip_n_filter=True)
    print(f"Mixed C: {mixed[0]} +/- {mixed[1]}")
    final = np.array(mixed)
    print(f"Final C: {final[0]} +/- {final[1]}")
    temp.append(final)
    return final
    for n in NSIDES:
        try:
            data = np.array((
                hp.ud_grade(data_set[0], n),
                hp.ud_grade(data_set[1], n),
                hp.ud_grade(data_set[2], n),
                hp.ud_grade(data_set[3], n)
            ))
            binmap = toolkit.HealpixBinMap(n)
            binmap.bin_catalogue(cat)
            binmap.load_catalogue(cat)
            output = binmap.divide_sample(point_mask, data, False, filter_set, 0)
            mixed = weights.excess_measurement(*output[1:])
            print(f"Mixed {n}: {mixed[0]} +/- {mixed[1]}")
            print(output[0])
            final = np.array(mixed)
            print(f"Final {n}: {final[0]} +/- {final[1]}")
        except ValueError:
            final = np.array([np.NaN, np.NaN])
        temp.append(final)
    return temp


if args.catalogue == "sdss":
    random_mask = toolkit.load_mask("sdss_mask", raise_dir=args.raise_dir)
else:
    raise ValueError

if args.data_mask == "sdss_act":
    point_mask = multiprocessing.managers.BaseManager.toolkit.load_mask("act_point", raise_dir=args.raise_dir)
    temp1 = multiprocessing.managers.BaseManager.toolkit.load_mask("sdss_mask", raise_dir=args.raise_dir)
    temp1.map = np.int_(temp1.map)
    #temp2 = toolkit.load_mask("act_point", raise_dir=args.raise_dir)
    temp1.map = 1 - temp1.map
    overdensity_mask = multiprocessing.managers.BaseManager.toolkit.CombinationMask(temp1, point_mask, invert=True, use_and=False)
    sky_mask_frac = 0.009859934289099422
    data_set = np.float_(np.array((
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_1.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_2.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_3.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_4.fits").map
    )))
    filter_set = "n_only"
elif args.data_mask == "sdss_planck":
    point_mask = toolkit.load_mask("planck_modified_point", raise_dir=args.raise_dir)
    temp1 = toolkit.load_mask("sdss_mask", raise_dir=args.raise_dir)
    temp1.map = np.int_(temp1.map)
    # temp2 = toolkit.load_mask("act_point", raise_dir=args.raise_dir)
    temp1.map = 1 - temp1.map
    overdensity_mask = toolkit.CombinationMask(temp1, point_mask, invert=True, use_and=False)
    sky_mask_frac = 0.0138443493342983
    sdss_mask = toolkit.load_mask("sdss_mask", args.raise_dir)
    temp = point_mask.map + 1j * sdss_mask.map
    data_set = np.float_(np.array((
        temp == 0,
        temp == 1j,
        temp == 1,
        temp == 1 + 1j
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
