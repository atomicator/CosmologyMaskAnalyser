from toolkit import toolkit, weights
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--catalogue", default="sdss_random")
parser.add_argument("--save_path", default="test.pdf")
parser.add_argument("--raise_path", type=int, default=2)
parser.add_argument("--weight_function", default="excess_measurement")
parser.add_argument("--min_z", type=float, default=0.0)
parser.add_argument("--min_r", type=float, default=10.0)
parser.add_argument("--max_z", type=float, default=20.0)
parser.add_argument("--max_r", type=float, default=10000.0)

args = parser.parse_args()


def filter(redshift, richness):
    global args
    return (args.min_z < redshift < args.max_z) and (args.min_r < richness < args.max_r)


raise_dir = args.raise_path
if args.catalogue == "full_sky":
    cat = toolkit.StarCatalogue("./" + raise_dir * "../" + "code/binned_results/random_full_sky.fits", table=True)
    data_name = "random full sky"
    cat.load_lon_lat()
elif args.catalogue == "sdss_random":
    #cat = toolkit.StarCatalogue("./" + raise_dir * "../" + "code/binned_results/test.fits", table=True)
    cat = toolkit.StarCatalogue("./" + raise_dir * "../" + "code/binned_results/random_sdss_10m.fits", table=True)
    cat.load_lon_lat()
    data_name = "sdss random"
elif args.catalogue == "sdss":
    cat = toolkit.load_catalogue("sdss", raise_dir)
    data_name = "sdss actual"
elif args.catalogue == "sdss_filtered":
    cat = toolkit.load_catalogue("sdss", raise_dir)
    cat.load_with_selection(filter, ["ZRED", "LAMBDA"], lon_lat=True)
    data_name = "\n" + rf"Filtered: ${args.min_z} < z < {args.max_z}$, ${args.min_r} < r < {args.max_r}$"
else:
    raise ValueError

print(len(cat.lon_lat))

if args.weight_function == "excess_measurement":
    y_axis_label = r"Excess"
    weight_function = weights.excess_measurement
    set_f_zero = True
    convert_to_mask_frac = False
    filter_set = "overkill"
elif args.weight_function == "density_weighting":
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
else:
    raise ValueError

#mask_names = ["planck_modified_point", "planck_modified_galactic", "planck_modified_total"]
#mask_names = ["planck_modified_point", "planck_modified_galactic", "planck_modified_total", "planck_galactic"]
mask_names = ["planck_modified_point", "planck_modified_galactic"]
labels = ["Point", "Galactic", "Total", "Old Galactic"]
#NSIDES = [1, 2, 4, 8, 16, 32]
#NSIDES = [2, 8, 32]
#NSIDES = [32]
#NSIDES = [1, 2, 4, 8, 16, 32, 64, 128, 256]
NSIDES = [1, 4, 8, 32, 128]
run_const = True
save_path = args.save_path
f = []
result_set = []
toolkit.plt_use_tex()
error_bar_colors = ["xkcd:aqua blue", "orange", "xkcd:mint green", "pink"]
line_colors = ["xkcd:electric blue", "red", "xkcd:grass green", "purple"]
#filter_set = "all"

fig = plt.figure()
ax = fig.add_subplot(111)

"""NSIDES = [4]
run_const = False
weight_function = weights.scatter
mask_names = ["planck_modified_galactic"]
"""
for mask_name in mask_names:
    mask = toolkit.load_mask(mask_name, raise_dir)
    if args.catalogue == "full_sky":
        data_set = np.array((
            np.zeros(12 * 2048 ** 2),
            1 - mask.map,
            np.zeros(12 * 2048 ** 2),
            mask.map
        ))
    else:
        sdss_mask = toolkit.load_mask("sdss_mask", raise_dir)
        temp = mask.map + 1j * sdss_mask.map
        data_set = np.float_(np.array((
            temp == 0,
            temp == 1j,
            temp == 1,
            temp == 1+1j
        )))

    results = []
    if not set_f_zero:
        f.append((1 - np.sum(mask.lookup_point(*cat.lon_lat.transpose())) / len(cat.lon_lat)) * 100)
    else:
        f.append(0)
    print(f"{mask_name}: f = {f[-1]}")
    if run_const:
        data = np.array((
            np.mean(data_set[0]),
            np.mean(data_set[1]),
            np.mean(data_set[2]),
            np.mean(data_set[3])
        ))
        binmap = toolkit.ConstantMap()
        binmap.bin_catalogue(cat)
        binmap.load_catalogue(cat)
        output = binmap.divide_sample(mask, data, True, filter_set)
        mixed = weight_function(*output[1:])
        print(f"Mixed C: {mixed[0]} +/- {mixed[1]}")
        if convert_to_mask_frac:
            final = np.array([output[0][0] + (1 - output[0][0] - output[0][1]) * mixed[0], (1 - output[0][0] - output[0][1]) * mixed[1]])
        else:
            final = np.array(mixed)
        print(f"Final C: {final[0]} +/- {final[1]}")
        results.append(final)

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
            output = binmap.divide_sample(mask, data, False, filter_set)
            mixed = weight_function(*output[1:])
            print(f"Mixed {n}: {mixed[0]} +/- {mixed[1]}")
            print(output[0])
            if convert_to_mask_frac:
                final = np.array([output[0][0] * 100 + (1 - output[0][0] - output[0][1]) * mixed[0], (1 - output[0][0] - output[0][1]) * mixed[1]])
            else:
                final = np.array(mixed)
            print(f"Final {n}: {final[0]} +/- {final[1]}")
        except ValueError:
            final = np.array([np.NaN, np.NaN])
        results.append(final)
    results = np.array(results).transpose()
    result_set.append(results)

print(result_set)

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
ax.set_xlim(1/2 * np.sqrt(1/2), NSIDES[-1] * np.sqrt(2))

ax.plot([1/2, NSIDES[-1]], np.zeros(2), color="k")
ax.set_xticks([0.5] + NSIDES, ["C"] + NSIDES)

ax.legend()
ax.set_xlabel("NSIDE")
ax.set_ylabel(y_axis_label)
ax.set_title(f"Binning algorithm using weight: {weight_function.__name__} and data: {data_name}")
plt.savefig(save_path)
plt.show()
