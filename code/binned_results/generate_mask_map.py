import matplotlib.pyplot as plt

from toolkit import toolkit
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--mask_one", help="The first mask to use", default="sdss_mask")
parser.add_argument("--mask_two", help="The second mask to use", default="planck_modified_point")
parser.add_argument("--path_raise", type=int, default=2)
parser.add_argument("--nside_internal", type=int, default=4096)
parser.add_argument("--nside", type=int, default=256)
parser.add_argument("--save_path")

args = parser.parse_args()
mask_names = [args.mask_one, args.mask_two]

mask1 = toolkit.load_mask(args.mask_one)
mask2 = toolkit.load_mask(args.mask_two)

fig = plt.figure()
ax = fig.add_subplot(111)

toolkit.gen_mask_comparison_map(mask1, mask2, NSIDE=args.nside, NSIDE_internal=args.nside_internal, name=f"{args.mask_one}_{args.mask_two}")
