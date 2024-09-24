from toolkit import toolkit, data
import argparse
import multiprocessing.pool
import numpy as np
import healpy as hp
import random

parser = argparse.ArgumentParser()

parser.add_argument("--catalogue", default="sdss")
parser.add_argument("--save_path", default="test.npy")
parser.add_argument("--raise_dir", type=int, default=2)
parser.add_argument("--weight_function", default="overdensity")
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("--iterations", type=int, default=1)
parser.add_argument("--target", type=int, default=400000)
parser.add_argument("--data_mask", default="sdss_act")
parser.add_argument("--overdensity", type=float, default=0.05)
parser.add_argument("--const_only", default=False, type=lambda x: (str(x).lower() == 'true'))
parser.add_argument("--invert_bias", default=False, type=lambda x: (str(x).lower() == 'true'))
args = parser.parse_args()

NSIDES = [0, 1, 2, 4, 8, 16, 32, 64, 256]

density_min = 0
density_max = 2 * args.target
density_steps = 1000

overdensity_min = -1
overdensity_max = 1
overdensity_steps = 1000


def test_function(const_only=args.const_only, overdensity=args.overdensity):
    global data_set
    random_points = toolkit.gen_random_coords(args.target, random_mask)[::-1].transpose()
    #bias_points = toolkit.gen_random_coords(len(random_points) * args.overdensity * sky_mask_frac * 5, overdensity_mask)[::-1].transpose()
    bias_points = toolkit.gen_random_coords(args.target * overdensity, random_mask)[::-1].transpose()
    if not args.invert_bias:
        bias_points = bias_points[point_mask.lookup_point(*bias_points.transpose()) == 0]
    else:
        bias_points = bias_points[point_mask.lookup_point(*bias_points.transpose()) != 0]
    print(len(bias_points))
    cat = toolkit.ClusterCatalogue()
    #cat.lon_lat = np.append(random_points, bias_points[:int(len(random_points) * overdensity * sky_mask_frac)], axis=0)
    cat.lon_lat = np.append(random_points, bias_points, axis=0)
    temp = []

    for NSIDE in NSIDES:
        if NSIDE != 0:
            data_array = np.array((
                hp.ud_grade(data_set[0], NSIDE),
                hp.ud_grade(data_set[1], NSIDE),
                hp.ud_grade(data_set[2], NSIDE),
                hp.ud_grade(data_set[3], NSIDE)
            ))
            binmap = toolkit.HealpixBinMap(NSIDE)
        else:
            data_array = np.array((
                (np.mean(data_set[0]),),
                (np.mean(data_set[1]),),
                (np.mean(data_set[2]),),
                (np.mean(data_set[3]),),
            ))
            binmap = toolkit.ConstantBinMap()
        binmap.set_mask(point_mask)
        binmap.bin_catalogue(cat)
        output = binmap.divide_sample(point_mask, data_array, filter_fully_masked=False, filter_empty=False)

        sky_masked_fraction = output[1]
        cluster_masked_fraction = output[2]
        n = output[3]
        sky_surveyed_fraction = (data_array[1] + data_array[3])[binmap.final_filter]

        masked_clusters = np.int_(np.round(n * cluster_masked_fraction))
        unmasked_clusters = np.int_(np.round(n * (1 - cluster_masked_fraction)))
        results = np.zeros((overdensity_steps, density_steps))

        #pool = multiprocessing.pool.Pool()
        thread_objects = []
        for i in range(overdensity_steps):
            thread_objects.append([])
            for j in range(density_steps):
                #thread_objects[-1].append(pool.apply_async(func, args=(density_min + ((j + 0.5) / density_steps) *
                #                                                       (density_max - density_min), overdensity_min +
                #                                                       ((i + 0.5) / overdensity_steps) *
                #                                                       (overdensity_max - overdensity_min), i, j)))
                results[i, j] = func(density_min + ((j + 0.5) / density_steps) *
                                     (density_max - density_min), overdensity_min +
                                     ((i + 0.5) / overdensity_steps) *
                                     (overdensity_max - overdensity_min), i, j, NSIDE, sky_masked_fraction,
                                                  sky_surveyed_fraction, masked_clusters,
                                                  unmasked_clusters)
        #for i in range(overdensity_steps):
        #    for j in range(density_steps):
        #        results[i][j] = thread_objects[i][j].get()

        results = results - np.max(results)
        results = np.exp(results)
        results = ((results / np.sum(results)) / ((overdensity_max - overdensity_min) * (density_max - density_min) /
                                                  (overdensity_steps * density_steps)))
        x = np.linspace(overdensity_min, overdensity_max, overdensity_steps)
        y = np.sum(results, axis=1)
        y = y * (density_max - density_min) / density_steps
        y_cum = np.cumsum(y / np.sum(y)) * 100  # convert to percentiles
        search_percentiles = [16, 50, 84]
        x_vals = []

        for percentile in search_percentiles:
            lower = 0
            upper = len(y_cum) - 1
            while True:
                # print(lower, upper)
                if y_cum[int((upper + lower) / 2)] < percentile:
                    lower = int((upper + lower) / 2)
                else:
                    upper = int((upper + lower) / 2)
                if (upper - lower) in (0, 1):
                    x_vals.append(x[lower])
                    break

        temp.append([x_vals[1], (x_vals[2] - x_vals[0]) / 2])
    print(temp)


def func(density, overdensity, _a, _b, NSIDE, sky_masked_fraction, sky_surveyed_fraction, masked_clusters,
         unmasked_clusters):
    if NSIDE != 0:
        pixel_area = (4 / 3) * np.pi / (12 * NSIDE ** 2)
    else:
        pixel_area = (4 / 3) * np.pi
    masked_cluster_expectation = (1 + overdensity) * sky_masked_fraction * density * sky_surveyed_fraction * \
                                 pixel_area
    unmasked_cluster_expectation = (1 - sky_masked_fraction) * density * sky_surveyed_fraction * pixel_area
    ln_masked_prob = np.zeros(np.shape(masked_cluster_expectation))
    ln_unmasked_prob = np.zeros(np.shape(unmasked_cluster_expectation))
    expectation_cutoff = 0.01
    valid = masked_cluster_expectation > expectation_cutoff
    ln_masked_prob[valid] = (np.log(masked_cluster_expectation[valid]) * masked_clusters[valid] -
                             masked_cluster_expectation[valid])
    valid = unmasked_cluster_expectation > expectation_cutoff
    ln_unmasked_prob[valid] = (np.log(unmasked_cluster_expectation[valid]) * unmasked_clusters[valid] -
                               unmasked_cluster_expectation[valid])
    ln_prob = np.sum(ln_masked_prob) + np.sum(ln_unmasked_prob)
    return ln_prob


if args.catalogue == "sdss":
    random_mask = data.load_mask("sdss_mask", raise_dir=args.raise_dir)
else:
    raise ValueError

if args.data_mask == "sdss_act":
    point_mask = data.load_mask("act_point", raise_dir=args.raise_dir)
    temp1 = data.load_mask("sdss_mask", raise_dir=args.raise_dir)
    temp1.map = np.int_(temp1.map)
    #temp2 = data.load_mask("act_point", raise_dir=args.raise_dir)
    temp1.map = 1 - temp1.map
    overdensity_mask = toolkit.CombinationMask(temp1, point_mask, invert=True, use_and=False)
    sky_mask_frac = 0.009859934289099422
    data_set = np.float64(np.array((
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_1.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_2.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_3.fits").map,
        toolkit.HealpyMask("../" * args.raise_dir + f"code/binned_results/sdss_mask_act_point_256_4.fits").map
    )))
    filter_set = "n_only"
elif args.data_mask == "sdss_planck":
    point_mask = data.load_mask("planck_modified_point", raise_dir=args.raise_dir)
    temp1 = data.load_mask("sdss_mask", raise_dir=args.raise_dir)
    temp1.map = np.int_(temp1.map)
    # temp2 = data.load_mask("act_point", raise_dir=args.raise_dir)
    temp1.map = 1 - temp1.map
    overdensity_mask = toolkit.CombinationMask(temp1, point_mask, invert=True, use_and=False)
    sky_mask_frac = 0.0138443493342983
    sdss_mask = data.load_mask("sdss_mask", args.raise_dir)
    temp = point_mask.map + 1j * sdss_mask.map
    data_set = np.float64(np.array((
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
pool.close()
pool.join()
results = []
for thread in threads:
    results.append(thread.get())
