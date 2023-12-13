from toolkit import toolkit
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--mask_one", choices=["sdss_mask", "planck_point", "planck_galactic"],
                    help="The first mask to use", default="sdss_mask")
parser.add_argument("--mask_two", choices=["sdss_mask", "planck_point", "planck_galactic"],
                    help="The second mask to use", default="planck_galactic")
parser.add_argument("--path_raise", type=int, default=2)
parser.add_argument("--res", type=int, default=1e4)
parser.add_argument("--save_path")

args = parser.parse_args()
mask_names = [args.mask_one, args.mask_two]

mask1 = toolkit.load_mask("sdss_mask")
mask2 = toolkit.load_mask("planck_galactic")

toolkit.gen_mask_comparison_map(mask1, mask2, NSIDE=1, res=int(1e3))

