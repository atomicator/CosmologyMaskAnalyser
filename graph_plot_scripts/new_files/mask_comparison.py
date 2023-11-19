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
parser.add_argument("--resolution", type=int, default=100000)
parser.add_argument("--save_path", help="The path to save the output")

args = parser.parse_args()
mask_names = [args.mask_one, args.mask_two]

masks = [toolkit.load_mask(mask_names[0], args.path_raise), toolkit.load_mask(mask_names[1], args.path_raise)]

mask_frac = toolkit.fraction_masked_pair(masks[0], masks[1], n=args.resolution)
unmasked_fraction_exact = mask_frac[2] / (mask_frac[0] + mask_frac[2])

print(mask_frac)

if args.save_path:
    np.savetxt("../" * args.path_raise + args.save_path, mask_frac, delimiter=",")
