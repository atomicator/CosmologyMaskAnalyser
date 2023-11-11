import random
from toolkit import toolkit
import numpy as np

sdss_mask = toolkit.HealpyMask("../../data/redmapper_dr8_public_v6.3_zmask.fits", mask_using_latlon=False,
                               hdu=1, partial=True)
planck_mask = toolkit.HealpyMask("../../data/planck_point_mask.fits")

sdss_mask.map[sdss_mask.map > 0.4] = 1.0
sdss_mask.map[sdss_mask.map < 0.3] = 0
sdss_mask.map = (sdss_mask.map - 1) * -1

data = planck_mask.compare(sdss_mask)

bootstrap_samples = int(1e3)

mean_estimates = [[np.real(np.sum(data) / len(data)), np.imag(np.sum(data) / len(data)),
                   np.sum(data == 1 + 1j) / len(data), np.sum(data == 1j) / len(data), np.sum(data == 1) / len(data),
                   np.sum(data == 0) / len(data)]]

# TODO: implement bootstrapping to estimate errors

for i in range(bootstrap_samples):
    print(f"Bootstrap iteration: {i}")
    new_data = random.choices(data, k=len(data))
    mean_estimates.append([np.real(np.sum(new_data) / len(new_data)), np.imag(np.sum(new_data) / len(new_data)),
                           np.sum(new_data == 1 + 1j) / len(new_data), np.sum(new_data == 1j) / len(new_data),
                           np.sum(new_data == 1) / len(new_data), np.sum(new_data) / len(new_data)])

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

#TODO: Plot som eof the results from above as a histogram
