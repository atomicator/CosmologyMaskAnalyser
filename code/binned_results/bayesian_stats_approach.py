import multiprocessing.pool
import matplotlib.pyplot as plt
from toolkit import toolkit, data
import numpy as np
import healpy as hp

#NSIDEs = [4, 8, 16, 32]
NSIDEs = [0]

for NSIDE in NSIDEs:
    raise_dir = 2
    cat_name = "sdss"
    mask_name = "act_point"
    lon_shift = 0

    #overdensity_min = [0.0, 0.0, 0.0][param_set]
    #overdensity_max = [0.3, 0.3, 0.3][param_set]
    overdensity_min = -0.0  # 0.0
    overdensity_max = 0.4  # 0.4
    overdensity_steps = 1000
    density_min = 8.0e4  # 8e4
    density_max = 8.2e4  # 8.2e4
    #density_min = [8.0e4, 8.1e4, 8.2e4][param_set]
    #density_max = [8.3e4, 8.4e4, 8.5e4][param_set]
    density_steps = 1000

    toolkit.plt_use_tex()


    def data_filter(z, r):
        return (z > 0) and (r > 20)


    cat = data.load_catalogue(cat_name, raise_dir)
    cat.load_data(selection_function=data_filter, requested_fields=["ZRED", "LAMBDA_CHISQ"], lon_lat=True)
    cat.lon_lat[:, 0] -= lon_shift
    cat.lon_lat[cat.lon_lat[:, 0] < 0, 0] += 360
    cat.lon_lat[cat.lon_lat[:, 0] > 360, 0] -= 360
    mask = data.load_mask(mask_name, raise_dir)
    sdss_mask = data.load_mask("sdss_mask", raise_dir, lon_shift=lon_shift)
    if lon_shift != 0.0:
        data_set = toolkit.gen_mask_comparison_map(sdss_mask, mask, write=False, NSIDE=256, NSIDE_internal=4096)
    else:
        data_set = np.float64(np.array((
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_1.fits").map,
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_2.fits").map,
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_3.fits").map,
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_4.fits").map
        )))

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
    binmap.set_mask(mask)
    binmap.bin_catalogue(cat)
    output = binmap.divide_sample(mask, data_array, filter_fully_masked=False, filter_empty=False)

    #print(output)

    sky_masked_fraction = output[1]
    cluster_masked_fraction = output[2]
    n = output[3]
    sky_surveyed_fraction = (data_array[1] + data_array[3])[binmap.final_filter]

    #print(f"data_array: {data_array}")

    masked_clusters = np.int_(np.round(n * cluster_masked_fraction))
    unmasked_clusters = np.int_(np.round(n * (1 - cluster_masked_fraction)))
    results = np.zeros((overdensity_steps, density_steps))


    def func(density, overdensity, a, b):
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
        results[a][b] = ln_prob
        # print(f"""DEBUG DATA")
        # density: {density}, overdensity: {overdensity}
        # sky masked fraction: {sky_masked_fraction}
        # sky surveyed fraction: {sky_surveyed_fraction}
        # masked_cluster_expectation: {masked_cluster_expectation}
        # masked_clusters: {masked_clusters}
        # unmasked_cluster_expectation: {unmasked_cluster_expectation}
        # unmasked_clusters: {unmasked_clusters}
        # ln_prob: {ln_prob}
        # END""")
        return ln_prob


    pool = multiprocessing.pool.Pool()
    thread_objects = []
    for i in range(overdensity_steps):
        thread_objects.append([])
        for j in range(density_steps):
            thread_objects[-1].append(pool.apply_async(func, args=(density_min + ((j + 0.5) / density_steps) *
                                                                   (density_max - density_min), overdensity_min +
                                                                   ((i + 0.5) / overdensity_steps) *
                                                                   (overdensity_max - overdensity_min), i, j)))
    for i in range(overdensity_steps):
        for j in range(density_steps):
            results[i][j] = thread_objects[i][j].get()

    #results = results - np.mean(results)
    print(results)
    results = results - np.max(results)
    print(results)
    results = np.exp(results)
    results = ((results / np.sum(results)) / ((overdensity_max - overdensity_min) * (density_max - density_min) /
                                              (overdensity_steps * density_steps)))
    #print(results)
    x, y = np.meshgrid(np.linspace(density_min, density_max, density_steps + 1),
                       np.linspace(overdensity_min, overdensity_max, overdensity_steps + 1))
    #print(np.shape(results))
    #print(np.shape(x))
    plt.pcolormesh(x, y, results)
    cbar = plt.colorbar(label=r"Likelihood Density Function $(\text{Clusters}^{-1} . \Omega)$")
    plt.clim(0, np.max(results))
    plt.title(f"NSIDE: {NSIDE}")
    plt.xlabel(r"Unmasked density $(\text{Clusters} . \Omega^{-1})$")
    plt.ylabel("Overdensity")
    plt.savefig(f"{NSIDE}.png")
    plt.show()

    #print(np.sum(results) * (overdensity_max - overdensity_min) * (density_max - density_min) /
    #      (overdensity_steps * density_steps))
    x = np.linspace(overdensity_min, overdensity_max, overdensity_steps)
    y = np.sum(results, axis=1)
    y = y * (density_max - density_min) / density_steps
    plt.plot(x, y, color="k")
    plt.ylim(0, np.max(y) * 1.1)
    #plt.show()

    #print(np.sum(y) * (overdensity_max - overdensity_min) / overdensity_steps)
    # Percentiles - need 68-95-99.7
    y_cum = np.cumsum(y / np.sum(y)) * 100  # convert to percentiles
    #print(y_cum)

    search_percentiles = [0.15, 2.5, 16, 25, 50, 75, 84, 97.5, 99.85]
    colours = ["m", "c", "g", "b", "r", "b", "g", "c", "m"]
    labels = [r"$3 \sigma$", r"$2 \sigma$", r"$1 \sigma$", r"$25 \%$", "Median", None, None, None, None]
    x_vals = []

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
                x_vals.append(x[lower])
                i = search_percentiles.index(percentile)
                plt.plot((x[lower], x[lower]), (0, y[lower]), color=colours[i], label=labels[i])
                break
    print(NSIDE, search_percentiles, x_vals)
    plt.legend()
    plt.title(f"NSIDE: {NSIDE}")
    plt.ylabel("Likelihood Density")
    plt.savefig(f"{NSIDE}_likelihood.png")
    plt.xlabel("Overdensity")
    plt.show()
    np.save(f"{NSIDE}.npy", x_vals)
