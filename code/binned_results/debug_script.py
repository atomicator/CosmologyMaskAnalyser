from toolkit import toolkit_new as toolkit, weights, data
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import argparse
import multiprocessing.pool

parser = argparse.ArgumentParser()

parser.add_argument("--catalogue", default="sdss_filtered")
parser.add_argument("--save_path", default="test.png")
parser.add_argument("--raise_path", type=int, default=2)
parser.add_argument("--weight_function", default="overdensity")
parser.add_argument("--min_z", type=float, default=0.0)
parser.add_argument("--min_r", type=float, default=20.0)
parser.add_argument("--max_z", type=float, default=20.0)
parser.add_argument("--max_r", type=float, default=10000.0)
parser.add_argument("--data_mask", default="sdss_act")
parser.add_argument("--mask_set", default="both")
parser.add_argument("--lon_shift", type=float, default=0.0)
args = parser.parse_args()


def data_filter(redshift, richness):
    global args
    return (args.min_z < redshift < args.max_z) and (args.min_r < richness < args.max_r)


raise_dir = args.raise_path
if args.catalogue == "full_sky":
    cat = toolkit.ClusterCatalogue("./" + raise_dir * "../" + "code/binned_results/random_full_sky.fits", table=True)
    data_name = "random full sky"
    cat.load_data()
elif args.catalogue == "sdss_random":
    cat = toolkit.ClusterCatalogue("./" + raise_dir * "../" + "code/binned_results/test.fits", table=True)
    cat.load_data()
    data_name = "sdss random 400k"
elif args.catalogue == "sdss":
    cat = data.load_catalogue("sdss", raise_dir)
    data_name = "sdss actual"
    data_mask = "sdss"
elif args.catalogue == "sdss_filtered":
    cat = data.load_catalogue("sdss", raise_dir)
    cat.load_with_selection(data_filter, ["ZRED", "LAMBDA_CHISQ"], lon_lat=True)
    cat.lon_lat[:, 0] -= args.lon_shift
    cat.lon_lat[cat.lon_lat[:, 0] < 0, 0] += 360
    cat.lon_lat[cat.lon_lat[:, 0] > 360, 0] -= 360
    data_name = "\n" + rf"Filtered: ${args.min_z} < z < {args.max_z}$, ${args.min_r} < r < {args.max_r}$"
elif args.catalogue in ("10m", "400k", "80k"):
    cat = toolkit.ClusterCatalogue("./" + raise_dir * "../" + f"code/binned_results/random_sdss_{args.catalogue}.fits",
                                   table=True)
    cat.load_data()
    data_name = f"sdss random {args.catalogue}"
elif args.catalogue == "planck_point_biased":
    cat = toolkit.ClusterCatalogue("./" + raise_dir * "../" + "code/binned_results/bias.fits", table=True)
    cat.load_data()
    cat1 = toolkit.ClusterCatalogue("./" + raise_dir * "../" + "code/binned_results/test.fits", table=True)
    cat1.load_data()
    cat.lon_lat = np.append(cat1.lon_lat, cat.lon_lat, axis=0)
    data_name = "sdss biased"
elif args.catalogue == "full_sky_inigo":
    cat = toolkit.ClusterCatalogue("../binned_results/test.fits", hdu=1, table=True)
    cat.load_data()
    data = np.load("../../data/random_catalogue_400k_inigo.npy")
    cat.lon_lat = np.array([data[0] * 180 / np.pi, 90 - data[1] * 180 / np.pi]).transpose()
    data_name = "full sky 400k from Inigo"
elif args.catalogue == "act_bias":
    cat = toolkit.ClusterCatalogue("./" + raise_dir * "../" + "code/binned_results/random_sdss_400k.fits", table=True)
    cat.load_data()
    cat1 = toolkit.ClusterCatalogue("./" + raise_dir * "../" + "code/binned_results/act_bias2.fits", table=True)
    cat1.load_data()
    cat.lon_lat = np.append(cat1.lon_lat[0:259 - 61], cat.lon_lat, axis=0)
    data_name = "sdss biased"
else:
    raise ValueError

data_mask = args.data_mask

N = len(cat.lon_lat)
print(N)
a = 0
print(a)

if args.weight_function == "excess":
    y_axis_label = r"Excess"
    weight_function = weights.excess_measurement
    set_f_zero = True
    convert_to_mask_frac = False
    filter_set = "n_only"
elif args.weight_function == "density":
    y_axis_label = r"Error ($\%$)"
    weight_function = weights.density_weighting
    set_f_zero = False
    convert_to_mask_frac = True
    filter_set = "n_only"
elif args.weight_function == "scatter":
    y_axis_label = ""
    weight_function = weights.regression_weighting
    set_f_zero = True
    convert_to_mask_frac = False
    filter_set = "all"
elif args.weight_function == "ratio":
    y_axis_label = ""
    weight_function = weights.ratio
    set_f_zero = True
    convert_to_mask_frac = False
    filter_set = "n_only"
elif args.weight_function == "overdensity":
    y_axis_label = r"Overdensity"
    weight_function = weights.overdensity
    set_f_zero = True
    convert_to_mask_frac = False
    filter_set = "n_only"
else:
    raise ValueError

if data_mask == "sdss_planck":
    if args.mask_set == "both":
        mask_names = ["planck_modified_point", "planck_modified_galactic"]
    elif args.mask_set == "point":
        mask_names = ["planck_modified_point"]
    else:
        raise ValueError
    labels = ["Point", "Galactic", "Total", "Old Galactic"]
elif data_mask == "sdss_act":
    if args.mask_set == "both":
        mask_names = ["act_point", "act_galactic"]
    elif args.mask_set == "point":
        mask_names = ["act_point"]
    else:
        raise ValueError
    labels = ["Point", "Galactic", "Total", "Old Galactic"]
elif data_mask == "full_sky":
    mask_names = ["planck_point_test", "planck_galactic_test"]
    labels = ["Planck Total", "Planck Galactic"]
elif data_mask == "sdss_act_lon_shift":
    mask_names = ["act_point_lon_test"]
    labels = ["ACT point shifted"]
elif data_mask == "sdss_planck_lon_shift":
    mask_names = ["planck_modified_point"]
    labels = ["Planck point shifted"]
else:
    raise ValueError

NSIDES = [1, 2, 4, 8, 16, 32, 64, 128]
#NSIDES = [8]
run_const = True
x_len = len(NSIDES)
if run_const:
    x_len += 1
save_path = args.save_path
f = []
result_set = []
toolkit.plt_use_tex()
error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:mint green", "pink"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple"]

fig = plt.figure()
ax = fig.add_subplot(111)


def run_const():  # Note: These don't work if variables are passed as parameters, I don't know why
    data = np.array((
        (np.mean(data_set[0]),),
        (np.mean(data_set[1]),),
        (np.mean(data_set[2]),),
        (np.mean(data_set[3]),),
    ))
    binmap = toolkit.ConstantBinMap()
    binmap.set_mask(mask)
    binmap.bin_catalogue(cat)
    #binmap.load_catalogue(cat)
    #output = binmap.divide_sample(mask, data, True, filter_set, a)
    output = binmap.divide_sample(mask, data)
    mixed = weight_function(*output[1:], skip_n_filter=True)
    print(f"Mixed C: {mixed[0]} +/- {mixed[1]}")
    if convert_to_mask_frac:
        final = np.array(
            [output[0][0] + (1 - output[0][0] - output[0][1]) * mixed[0], (1 - output[0][0] - output[0][1]) * mixed[1]])
    else:
        final = np.array(mixed)
    print(f"Final C: {final[0]} +/- {final[1]}")
    return final


def run_nside(n):
    try:
        data = np.array((
            hp.ud_grade(data_set[0], n),
            hp.ud_grade(data_set[1], n),
            hp.ud_grade(data_set[2], n),
            hp.ud_grade(data_set[3], n)
        ))
        binmap = toolkit.HealpixBinMap(n)
        binmap.set_mask(mask)
        binmap.bin_catalogue(cat)
        #binmap.load_catalogue(cat)
        #output = binmap.divide_sample(mask, data, False, filter_set, a)
        output = binmap.divide_sample(mask, data)
        mixed = weight_function(*output[1:])
        print(f"Mixed {n}: {mixed[0]} +/- {mixed[1]}")
        print(output[0])
        if convert_to_mask_frac:
            final = np.array([output[0][0] * 100 + (1 - output[0][0] - output[0][1]) * mixed[0],
                              (1 - output[0][0] - output[0][1]) * mixed[1]])
        else:
            final = np.array(mixed)
        print(f"Final {n}: {final[0]} +/- {final[1]}")
    except ValueError:
        final = np.array([np.nan, np.nan])
    return final


for mask_name in mask_names:
    mask = data.load_mask(mask_name, raise_dir)
    mask.set_fig_ax(fig, ax)
    if args.lon_shift != 0.0:
        #if args.lon_shift  :
        print("test2")
        sdss_mask = data.load_mask("sdss_mask", raise_dir, lon_shift=args.lon_shift)
        data_set = toolkit.gen_mask_comparison_map(sdss_mask, mask, write=False, NSIDE=256, NSIDE_internal=4096)
    elif data_mask == "full_sky":
        data_set = np.array((
            np.zeros(12 * 2048 ** 2),
            1 - mask.map,
            np.zeros(12 * 2048 ** 2),
            mask.map
        ))
    elif data_mask == "sdss_planck":
        sdss_mask = data.load_mask("sdss_mask", raise_dir)
        temp = mask.map + 1j * sdss_mask.map
        data_set = np.float64(np.array((
            temp == 0,
            temp == 1j,
            temp == 1,
            temp == 1 + 1j
        )))
    elif data_mask == "sdss_act":
        data_set = np.float64(np.array((
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_1.fits").map,
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_2.fits").map,
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_3.fits").map,
            toolkit.HealpyMask("../" * raise_dir + f"code/binned_results/sdss_mask_{mask_name}_256_4.fits").map
        )))
    elif data_mask == "sdss_act_lon_shift" or "sdss_planck_lon_shift":
        sdss_mask = data.load_mask("sdss_mask", raise_dir)
        data_set = toolkit.gen_mask_comparison_map(sdss_mask, mask, write=False, NSIDE=128, NSIDE_internal=1024)
    else:
        raise ValueError

    results = np.zeros((x_len, 2))
    pool = multiprocessing.pool.Pool(processes=x_len)
    thread_objects = [[] for i in range(x_len)]
    index = 0
    if not set_f_zero:
        f.append((1 - np.sum(mask.lookup_point(*cat.lon_lat.transpose())) / len(cat.lon_lat)) * 100)
    else:
        f.append(0)
    print("running const")
    if run_const:
        #thread_objects[index] = pool.apply_async(run_const)
        thread_objects[index] = pool.apply_async(toolkit.run_const, (data_set, mask, cat, weight_function,
                                                                     convert_to_mask_frac))
        index += 1

    print("running everything else")

    for n in NSIDES:
        print(n)
        #thread_objects[index] = pool.apply_async(run_nside, (n,))
        thread_objects[index] = pool.apply_async(toolkit.run_nside, (n, data_set, mask, cat, weight_function,
                                                                     convert_to_mask_frac))
        index += 1

    print("all running")

    for i in range(x_len):
        print(f"i: {i}")
        results[i] = thread_objects[i].get()
        print(results[i])

    results = np.array(results).transpose()
    result_set.append(results)

print(result_set)
plt.clf()
fig = plt.figure()
ax = fig.add_subplot()
for i in range(len(result_set)):
    x = NSIDES.copy()
    if run_const:
        x = [0.5] + x
    print(result_set[i][1])
    print()
    ax.errorbar(x, result_set[i][0] - f[i], result_set[i][1], marker="+", ecolor=error_bar_colors[i],
                ls="none", color=line_colors[i], capsize=3, capthick=1, label=labels[i])
    ax.plot(x, result_set[i][0] - f[i], color=line_colors[i])

ax.set_xscale("log", base=2)
ax.set_xlim(1 / 2 * np.sqrt(1 / 2), NSIDES[-1] * np.sqrt(2))

ax.plot([1 / 2, NSIDES[-1]], np.zeros(2), color="k")
ax.set_xticks([0.5] + NSIDES, ["C"] + NSIDES)

ax.legend()
ax.set_xlabel("NSIDE")
ax.set_ylabel(y_axis_label)
ax.set_title("SDSS data and the ACT masks")
plt.savefig(save_path)
np.save(save_path[:-4] + ".npy", np.array(result_set))
plt.show()
