from toolkit import toolkit
import numpy as np

data = toolkit.get_file_info("../data/planck_galactic_mask.fits", hdu=0)
toolkit.get_header_info(data)

mask = toolkit.HealpyMask("../data/planck_galactic_mask.fits")

print(mask.lookup_point(np.pi / 2, 0))

print(mask.calc_exact_unmasked_fraction())

mask.plot("../graphs/planck_test.png")
