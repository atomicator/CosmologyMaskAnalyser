import random
from toolkit import toolkit
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--mask_one", choices=["sdss_mask", "planck_point", "planck_galactic"],
                    help="The first mask to use", default="sdss_mask")
parser.add_argument("--mask_two", choices=["sdss_mask", "planck_point", "planck_galactic"],
                    help="The second mask to use", default="planck_galactic")
parser.add_argument("--path_raise", type=int, default=2)
parser.add_argument("--bootstrap_iterations", type=int, default=5)
parser.add_argument("--save_path")

args = parser.parse_args()
mask_names = [args.mask_one, args.mask_two]

print(mask_names)

print(1 == 1.0)

mask = [toolkit.load_mask(mask_names[0], args.path_raise), toolkit.load_mask(mask_names[1], args.path_raise)]
"""
mask[0] = toolkit.HealpyMask("../../data/redmapper_dr8_public_v6.3_zmask.fits", mask_using_latlon=False,
                               hdu=1, partial=True)
mask[1] = toolkit.HealpyMask("../../data/planck_point_mask.fits")

mask[0].map[mask[0].map > 0.4] = 1.0
mask[0].map[mask[0].map < 0.3] = 0
mask[0].map = (mask[0].map - 1) * -1
"""
data = mask[1].compare(mask[0])

bootstrap_samples = int(1e3)
k = len(data)
mean_estimates = [[np.real(np.sum(data) / k), np.imag(np.sum(data) / k),
                   np.sum(data == 1 + 1j) / k, np.sum(data == 1j) / k, np.sum(data == 1) / k,
                   np.sum(data == 0) / k]]

# TODO: implement bootstrapping to estimate errors

for i in range(bootstrap_samples):
    print(f"Bootstrap iteration: {i}")
    new_data = np.array(random.choices(data, k=k))
    estimate = [np.real(np.sum(new_data) / k), np.imag(np.sum(new_data) / k),
                           np.sum(new_data == 1 + 1j) / k, np.sum(new_data == 1j) / k,
                           np.sum(new_data == 1) / k, np.sum(new_data) / k]
    print(estimate)
    mean_estimates.append(estimate)

print("mean, std")
for i in range(len(mean_estimates[0])):
    print(f"{np.mean(mean_estimates[i])}")

mean_estimates = np.array(mean_estimates)

print(f"Planck masked fraction:         {mean_estimates[0][0]} +/- {np.std(mean_estimates[:, 0])}")
print(f"SDSS masked fraction:           {mean_estimates[0][1]} +/- {np.std(mean_estimates[:, 1])}")
print(f"Fraction masked by both:        {mean_estimates[0][2]} +/- {np.std(mean_estimates[:, 2])}")
print(f"Fraction masked by Planck only: {mean_estimates[0][3]} +/- {np.std(mean_estimates[:, 3])}")
print(f"Fraction masked by SDSS only:   {mean_estimates[0][4]} +/- {np.std(mean_estimates[:, 4])}")
print(f"Fraction masked by neither:     {mean_estimates[0][5]} +/- {np.std(mean_estimates[:, 5])}")

#TODO: Plot some of the results from above as a histogram

if args.save_path:
    np.savetxt(args.save_path, mean_estimates, delimiter=",")
