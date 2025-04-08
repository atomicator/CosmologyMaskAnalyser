import matplotlib.pyplot as plt
from toolkit import toolkit, data
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--mask_one", help="The first mask to use", default="sdss_mask")
parser.add_argument("--mask_two", help="The second mask to use", default="planck_modified_point")
parser.add_argument("--path_raise", type=int, default=2)
parser.add_argument("--nside_internal", type=int, default=65536)
parser.add_argument("--nside", type=int, default=512)
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("--save_path")

args = parser.parse_args()
mask_names = [args.mask_one, args.mask_two]
print("Loading masks")
#mask1 = data.load_mask(args.mask_one, raise_dir=args.path_raise)
print(f"Loaded {args.mask_one}")
#mask2 = data.load_mask(args.mask_two, raise_dir=args.path_raise)
print(f"Loaded {args.mask_two}")
kwargs = {"raise_dir":args.path_raise}
toolkit.gen_mask_comparison_map(data.load_mask, args=[], kwargs=kwargs, mask1_name=args.mask_one,
                                mask2_name=args.mask_one, NSIDE=args.nside, NSIDE_internal=args.nside_internal,
                                name="../" * args.path_raise + f"./data/{args.mask_one}_{args.mask_two}",
                                num_thread=args.threads)
