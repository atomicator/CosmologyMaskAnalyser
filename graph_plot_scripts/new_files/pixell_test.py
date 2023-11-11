from toolkit import toolkit

mask = toolkit.PixellMask("../../data/ACT_mask.fits", hdu=1, step=6)
mask.plot(".../graphs/act_mask/toolkit_test.png", title="Toolkit plot test")
print(mask.calc_exact_unmasked_fraction())
print(mask.simulate_unmasked_fraction(1e6))
