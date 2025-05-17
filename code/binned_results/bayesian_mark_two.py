import multiprocessing.pool, multiprocessing
import matplotlib.pyplot as plt
from toolkit import toolkit, data
import numpy as np
import healpy as hp
import warnings
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--overdensity", type=float, help="Overdensity", default=0.2)
parser.add_argument("-p", "--path", type=str, help="Output path", default="./temp.npy")
parser.add_argument("-t", "--threads", type=int, help="Number of threads", default=5)
parser.add_argument("-n", "--processes", type=int, help="Number of processes", default=5)
parser.add_argument("-r", "--realisations", type=int, help="Number of realisations", default=1)
parser.add_argument("-a", "--target", type=int, default=80000)
parser.add_argument("-i", "--invert_bias", default=False, type=lambda x: (str(x).lower() == 'true'))
parser.add_argument("-f", "--flat_prior", default=False, type=lambda x: (str(x).lower() == 'true'))
parser.add_argument("-d", "--debug", type=float, help="Debug", default=10.0)
parser.add_argument("-s", "--raise_dir", type=int, help="Number of folders to raise", default=2)
# TODO: Debug args below this
parser.add_argument("--data_mask", default="sdss_spt")
parser.add_argument("--catalogue", default="random")
parser.add_argument("--lon_shift", type=float, default=0.0)
to_print = 20
lock = multiprocessing.Lock()

args = parser.parse_args()

#NSIDES = [0, 1]
NSIDES = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
#NSIDES = [32]

def data_filter(redshift, richness):
    global args
    return (args.min_z < redshift < args.max_z) and (args.min_r < richness < args.max_r)

raise_dir = args.raise_dir
#cat_name = "sdss"


data_mask = args.data_mask
if data_mask == "sdss_act":
    mask_name = "act_point"
elif data_mask == "sdss_planck":
    mask_name = "planck_modified_point"
elif data_mask == "sdss_spt":
    mask_name = "spt_point"
else:
    raise ValueError

#overdensity_min = -0.95
#overdensity_max = 1 / (1 + overdensity_min)
#overdensity_steps = 1001
#overdensity = np.linspace(overdensity_min, overdensity_max, overdensity_steps)

def map_to_overdensity(prior):
    final = np.zeros(prior.shape)
    temp = prior >= 0
    final[temp] = prior[temp]
    temp = np.bitwise_not(temp)
    final[temp] = -1 * prior[temp] / (prior[temp] - 1)
    return final

prior_max = 0.8
prior_min = - prior_max
overdensity_steps = 10001
#prior = np.linspace(prior_min, prior_max, overdensity_steps)
overdensity = np.linspace(0.5, 1.5, overdensity_steps) - 1 #map_to_overdensity(prior)
density_steps = 1000
lon_shift = 0.0

density_range = 0.03

toolkit.plt_use_tex()

def data_filter(z, r):
    return (z > 0) and (r > 20)

print("Loading mask")
sdss_mask = data.load_mask("sdss_mask", raise_dir, lon_shift=args.lon_shift)
mask = data.load_mask(mask_name, raise_dir)
#Lock = multiprocessing.Lock()
to_write = []
hashmap_cache = {}
if lon_shift != 0.0:
    print("Generating rotated mask")
    data_set = toolkit.gen_mask_comparison_map(sdss_mask, mask, write=False, NSIDE=256, NSIDE_internal=4096)
else:
    print("Loading masked fractions")
    data_set = np.float64(np.array((
        toolkit.HealpyMask("../" * raise_dir + f"data/sdss_mask_{mask_name}_512_1.fits").map,
        toolkit.HealpyMask("../" * raise_dir + f"data/sdss_mask_{mask_name}_512_2.fits").map,
        toolkit.HealpyMask("../" * raise_dir + f"data/sdss_mask_{mask_name}_512_3.fits").map,
        toolkit.HealpyMask("../" * raise_dir + f"data/sdss_mask_{mask_name}_512_4.fits").map
    )))

for NSIDE in NSIDES:
    if NSIDE != 0:
        data_array = np.array((
            hp.ud_grade(data_set[0], NSIDE),
            hp.ud_grade(data_set[1], NSIDE),
            hp.ud_grade(data_set[2], NSIDE),
            hp.ud_grade(data_set[3], NSIDE)
        ))
    else:
        data_array = np.array((
            (np.mean(data_set[0]),),
            (np.mean(data_set[1]),),
            (np.mean(data_set[2]),),
            (np.mean(data_set[3]),),
        ))
    hashmap_cache[NSIDE] = data_array

def to_thread():
    np.random.seed()
    array_to_return = []
    if args.catalogue == "random":
        #print("Generating catalogue")
        cat = toolkit.ClusterCatalogue()
        random_points = toolkit.gen_random_coords(args.target, sdss_mask)[::-1].transpose()
        if args.overdensity != 0.0:
            bias_points = toolkit.gen_random_coords(args.target * args.overdensity, sdss_mask)[::-1].transpose()
            if not args.invert_bias:
                bias_points = bias_points[mask.lookup_point(*bias_points.transpose()) == 0.0]
            else:
                bias_points = bias_points[mask.lookup_point(*bias_points.transpose()) != 0.0]
            cat.lon_lat = np.append(arr=random_points, values=bias_points, axis=0)
        else:
            cat.lon_lat = random_points
    elif args.catalogue == "sdss":
        cat = data.load_catalogue("sdss", raise_dir)
    elif args.catalogue == "sdss_filtered":
        cat = data.load_catalogue("sdss", raise_dir)
        cat.load_data(selection_function=data_filter, requested_fields=["ZRED", "LAMBDA_CHISQ"], lon_lat=True)
        cat.lon_lat[:, 0] -= args.lon_shift
        cat.lon_lat[cat.lon_lat[:, 0] < 0, 0] += 360
        cat.lon_lat[cat.lon_lat[:, 0] > 360, 0] -= 360
    else:
        raise ValueError

    for NSIDE in NSIDES:
        results = np.zeros(overdensity_steps)  # replace with mutex
        #print("Resizing sky fractions")
        data_array = hashmap_cache[NSIDE]
        if NSIDE != 0:
            binmap = toolkit.HealpixBinMap(NSIDE)
        else:
            binmap = toolkit.ConstantBinMap()
        #print("Creating binmap")
        binmap.set_mask(mask)
        binmap.bin_catalogue(cat)
        #print("Dividing catalogue")
        output = binmap.divide_sample(mask, data_array, filter_fully_masked=False, filter_empty=False)

        #cat = toolkit.ClusterCatalogue()
        #cat.lon_lat = np.append(random_points, bias_points[:int(len(random_points) * overdensity * sky_mask_frac)], axis=0)
        #print(output)

        sky_masked_fraction = output[1]
        cluster_masked_fraction = output[2]
        n = output[3]
        #sky_surveyed_fraction = (data_array[1] + data_array[3])[binmap.final_filter]

        masked_clusters = np.int_(np.round(n * cluster_masked_fraction))
        unmasked_clusters = np.int_(np.round(n * (1 - cluster_masked_fraction)))
        #shared_arr = multiprocessing.Array("d", np.zeros(overdensity_steps))

        #if NSIDE != 0:
        #    pixel_area = 4 * np.pi / (12 * NSIDE ** 2)
        #else:
        #    pixel_area = 4 * np.pi

        def func(i):
            #print(sky_masked_fraction[i], masked_clusters[i], unmasked_clusters[i])
            expectation_cutoff = 0
            #expectation_cutoff = 1000 * (NSIDE / 32768) ** 2  # At this scale, quantization errors caused by the calculations of
            # f_s become significant. If f_s or 1 - f_s is less than this value, the pixel is rejected.
            # Resolution of the calculation was 8192, for reference.
            #expectation_cutoff = args.debug
            if not (expectation_cutoff < sky_masked_fraction[i] < 1 - expectation_cutoff):
                return 0
            #if (masked_clusters[i] + unmasked_clusters[i]) < 1:
            #    return 0
            #if (unmasked_clusters[i] + masked_clusters[i]) < 5:
            #    return 0
            #if (unmasked_clusters[i] + masked_clusters[i]) * min(sky_masked_fraction[i], 1 - sky_surveyed_fraction[i]) < 0.001:
            #    return 0
            #if (unmasked_clusters[i] + masked_clusters[i]) < 5 or \
            #        (unmasked_clusters[i] + masked_clusters[i]) * min(sky_surveyed_fraction[i], 1 - sky_surveyed_fraction[i]) \
            #                                                                                    < expectation_cutoff:
            #    return 0
            #if unmasked_clusters[i] < 1:
            #    return 0  # quicker than applying a NaN filter later
            """ln_prob = np.zeros(overdensity_steps)
            density_mid = (unmasked_clusters[i] + masked_clusters[i]) / (pixel_area * sky_surveyed_fraction[i])
            density_mid = unmasked_clusters[i] / (pixel_area * sky_surveyed_fraction[i] * (1 - sky_masked_fraction[i]))
            density_min = density_mid * max(0, 1 - 10/np.sqrt(unmasked_clusters[i]))
            density_max = density_mid * max(0, 1 + 10/np.sqrt(unmasked_clusters[i]))
            print(density_min, density_max)
            density_min = 1e4
            density_delta = min(density_mid - density_min, density_mid * 100 /
                               np.sqrt((unmasked_clusters[i] + masked_clusters[i])))
            density_max = 1e5
            density = np.outer(np.linspace(density_min, density_max, density_steps), np.ones(np.shape(overdensity)))
            density = np.outer(np.linspace(density_mid - density_delta, density_mid + density_delta, density_steps)
                               , np.ones(np.shape(overdensity)))
            alpha = np.outer(np.ones(np.shape(density[:, 0])), overdensity)
            masked_cluster_expectation = (1 + alpha) * sky_masked_fraction[i] * density * sky_surveyed_fraction[i] * \
                                        pixel_area
            unmasked_cluster_expectation = (1 - sky_masked_fraction[i]) * density * sky_surveyed_fraction[i] * pixel_area
            if np.min((masked_cluster_expectation, unmasked_cluster_expectation)) < expectation_cutoff:  # prevents NaN forming later - only triggers if
                # the masked fraction of the pixel is small
                return 0
            temp = ((np.log(masked_cluster_expectation) * masked_clusters[i] - masked_cluster_expectation)
                             + (np.log(unmasked_cluster_expectation) * unmasked_clusters[i] - unmasked_cluster_expectation))
            temp -= np.max(temp)
            temp = np.exp(temp)
            #plt.imshow(temp)
            #plt.title("prob 2")
            #plt.show()
            with np.errstate(all="ignore"):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    temp = ((np.log(masked_cluster_expectation) * masked_clusters[i] - masked_cluster_expectation)
                                     + (np.log(unmasked_cluster_expectation) * unmasked_clusters[i] -
                                        unmasked_cluster_expectation))
                    temp_min = np.min(temp)
                    if np.isnan(temp_min):
                        return 0
            temp[np.isnan(temp)] = temp_min
            temp = temp - np.max(temp)
            temp = np.exp(temp)
            temp = temp / np.sum(temp)
            ln_prob = np.log(np.sum(temp, axis=0))
            #ln_prob = np.sum((np.log(masked_cluster_expectation) * masked_clusters[i] - masked_cluster_expectation)
            #                 + (np.log(unmasked_cluster_expectation) * unmasked_clusters[i] - unmasked_cluster_expectation)
            #                  , axis=0)
            ln_prob -= np.max(ln_prob)"""
            """global to_print
            lock.acquire()
            if to_print > 0:
                print(to_print)
                to_print -= 1
                #plt.imshow(temp, aspect=0.01, interpolation='none', cmap="plasma", extent=(prior_min, prior_max, float(density[0][0]), float(density[1][1])))
                #plt.plot((args.overdensity + 1, args.overdensity + 1), (density[0][0], density[1][1]), color="k")
                #plt.title(str(to_print + 1))
                #plt.ylabel("Density")
                #plt.xlabel("Overdensity")
                #plt.colorbar()
                #plt.show()
                plt.plot(overdensity + 1, np.exp(ln_prob) / np.sum(np.exp(ln_prob)), label="Numerical")
                #plt.plot((args.overdensity + 1, args.overdensity + 1), (np.exp(np.min(ln_prob)), np.exp(np.max(ln_prob))), color="k")
                debug = np.exp(masked_clusters[i] * np.log(overdensity + 1) - (masked_clusters[i] + unmasked_clusters[i])
                               * np.log(1 + (sky_masked_fraction[i] / (1 - sky_masked_fraction[i])) * (1 + overdensity)))
                debug = debug / np.sum(debug)
                plt.plot(overdensity + 1, debug, linestyle="dotted", color="r", label="Analytic")
                plt.legend()
                #plt.plot(np.array((1, 1)) * (((masked_clusters[i] - (masked_clusters[i] + unmasked_clusters[i]) * sky_masked_fraction[i]) /
                #         (unmasked_clusters[i] * sky_masked_fraction[i])) + 1), (np.exp(np.min(ln_prob)), np.exp(np.max(ln_prob))), color="k", linestyle="dashed")
                plt.xlabel("Overdensity")
                plt.ylabel("LDF")
                plt.title(rf"{to_print + 1}, $n_{{m}}$: {masked_clusters[i]}, $n_{{nm}}$: {unmasked_clusters[i]}, $f_{{m}}$: {sky_masked_fraction[i]}, E: \
                    {1 + (masked_clusters[i] - (masked_clusters[i] + unmasked_clusters[i]) * sky_masked_fraction[i]) / (unmasked_clusters[i] * sky_masked_fraction[i])}")
                plt.xscale("log")
                plt.show()
            else:
                #exit()
                thread.interrupt_main(KeyboardInterrupt)
            lock.release()"""
            with (np.errstate(all="ignore")):
                with warnings.catch_warnings():
                    #debug = masked_clusters[i] * np.log(overdensity + 1) - (masked_clusters[i] + unmasked_clusters[i]) \
                    #               * np.log(1 + (sky_masked_fraction[i] / (1 - sky_masked_fraction[i])) * (1 + overdensity))
                    #debug = masked_clusters[i] * np.log(1 + overdensity) - (masked_clusters[i] + unmasked_clusters[i]) \
                    #        * np.log(sky_masked_fraction[i] * (1 + overdensity) + (1 - sky_masked_fraction[i]))
                    if args.flat_prior:
                        debug = masked_clusters[i] * np.log(1 + overdensity) - (masked_clusters[i] + unmasked_clusters[i] + 1) \
                                * np.log(sky_masked_fraction[i] * (1 + overdensity) + (1 - sky_masked_fraction[i]))
                    else:
                        debug = masked_clusters[i] * np.log(1 + overdensity) - (masked_clusters[i] + unmasked_clusters[i]) \
                            * np.log(sky_masked_fraction[i] * (1 + overdensity) + (1 - sky_masked_fraction[i]))
                    if np.isnan(debug).all():
                        print(masked_clusters[i], unmasked_clusters[i], sky_masked_fraction[i], sky_masked_fraction[i])
                    debug[np.isnan(debug)] = np.nanmin(debug)
                    debug = debug - np.nanmax(debug)
                    ln_prob = debug
            #return ln_prob
            nonlocal lock
            lock.acquire()
            nonlocal results
            results += ln_prob
            lock.release()
            return None

        pool = multiprocessing.pool.ThreadPool(args.threads)
        thread_objects = []
        lock = multiprocessing.Lock()
        #print("Initiating threads")
        for pixel_num in range(len(masked_clusters)):
            #TODO: Thread this
            thread_objects.append(pool.apply_async(func, args=(pixel_num,)))
            # TODO: Constrain range of density - calculate density of unmasked regions
            # This is a problem when no unmasked clusters exist - introduce a prior on density?

        #for i in range(len(masked_clusters)):
            #results += thread_objects[i].get()
            #print(thread_objects[i])
            #print(f"r: {np.shape(results)}")
            #print(f"t: {np.shape(thread_objects[i].get())}")
            #results += thread_objects[i].get()
            #print(thread_objects[i].get())
        pool.close()
        pool.join()

        results = results - np.max(results)
        results = np.exp(results)
        results = (results / np.sum(results))
        y_cum = np.cumsum(results) * 100  # convert to percentiles
        search_percentiles = [0.15, 2.5, 16, 25, 50, 75, 84, 97.5, 99.85]
        colours = ["m", "c", "g", "b", "r", "b", "g", "c", "m"]
        labels = [r"$3 \sigma$", r"$2 \sigma$", r"$1 \sigma$", r"$25 \%$", "Median", None, None, None, None]
        x_vals = []

        plt.plot(overdensity, results)
        plt.show()

        for percentile in search_percentiles:
            lower = 0
            upper = len(y_cum) - 1
            while True:
                #print(lower, upper)
                if y_cum[int((upper + lower) / 2)] < percentile:
                    lower = int((upper + lower) / 2)
                else:
                    upper = int((upper + lower) / 2)
                if (upper - lower) in (0, 1):
                    #x_vals.append(overdensity[lower])
                    #i = search_percentiles.index(percentile)
                    #plt.plot((overdensity[lower], overdensity[lower]), (0, results[lower]), color=colours[i],
                    #         label=labels[i])
                    x_vals.append(overdensity[lower])
                    i = search_percentiles.index(percentile)
                    plt.plot((overdensity[lower], overdensity[lower]), (0, results[lower]), color=colours[i],
                             label=labels[i])
                    break
        #to_write.append([NSIDE, *x_vals])
        print(f"NSIDE {NSIDE}: {x_vals[4]} +/- {x_vals[6] / 2 - x_vals[2] / 2}")
        #array_to_return.append([NSIDE, *x_vals])
        array_to_return.append([NSIDE, *y_cum])
        #print(f"NSIDE {NSIDE}: {x_vals[4]} +/- {x_vals[6] / 2 - x_vals[2] / 2}")
        plt.clf()
    return array_to_return

globalPool = multiprocessing.Pool(args.processes)
globalThreadObjects = []
for j in range(args.realisations):
    globalThreadObjects.append(globalPool.apply_async(to_thread))
for thread in globalThreadObjects:
    to_write.append(thread.get())

to_write = np.array(to_write)
#print(to_write)
np.save(args.path, to_write)