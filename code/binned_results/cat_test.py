from toolkit import toolkit, weights
import numpy as np
import astropy

print(np.array((0,)) / 2)
exit()
"""def filter(redshift):
    global results
    results.append(redshift)
    return redshift > 0.0


results = []
#cat = toolkit.load_catalogue("sdss")
#cat = toolkit.ClusterCatalogue("../../data/dr8_run_redmapper_v5.10_lgt5_catalog.fit", hdu=1)


max_z = 0.5
min_z = 0.1

cat = astropy.io.fits.open("../../data/dr8_run_redmapper_v5.10_lgt5_catalog.fit", hdu=1)

print(cat[1].columns)"""
"""
cat = toolkit.ClusterCatalogue("../../data/dr8_run_redmapper_v5.10_lgt5_catalog.fit", hdu=1)
#cat.load_with_selection(filter, ["Z_LAMBDA"], lon_lat=True)
#cat.load_with_selection(filter, ["ZRED"], lon_lat=True)
cat.load_with_selection(filter, ["LAMBDA_CHISQ"], lon_lat=True)
results = np.array(results)
print(np.min(results), np.max(results))
median = np.median(results)
upper = np.median(results[results > median])
lower = np.median(results[results < median])

print(upper, median, lower)
"""

"""data = (
    np.mean(toolkit.HealpyMask("./sdss_mask_planck_modified_point_256_1.fits").map),
    np.mean(toolkit.HealpyMask("./sdss_mask_planck_modified_point_256_2.fits").map),
    np.mean(toolkit.HealpyMask("./sdss_mask_planck_modified_point_256_3.fits").map),
    np.mean(toolkit.HealpyMask("./sdss_mask_planck_modified_point_256_4.fits").map)
)"""
"""point_mask = toolkit.load_mask("planck_modified_point", 2)
sdss_mask = toolkit.load_mask("sdss_mask", 2)
temp = point_mask.map + 1j * sdss_mask.map
data = np.float_(np.array((
    np.mean(np.float_(temp == 0)),
    np.mean(np.float_(temp == 1j)),
    np.mean(np.float_(temp == 1)),
    np.mean(np.float_(temp == 1 + 1j))
)))"""
"""
galactic_masked_fraction = data[0] + data[1]"""

mask_name = "sdss_mask_planck_modified_point"
data = (
    np.mean(toolkit.HealpyMask(f"./{mask_name}_256_1.fits").map),
    np.mean(toolkit.HealpyMask(f"./{mask_name}_256_2.fits").map),
    np.mean(toolkit.HealpyMask(f"./{mask_name}_256_3.fits").map),
    np.mean(toolkit.HealpyMask(f"./{mask_name}_256_4.fits").map)
)
#point_masked_fraction = data[0] + data[1]
#final = point_masked_fraction / (1 - galactic_masked_fraction)
final = data[1] / (data[1] + data[3])
#print((0.015567657065637092 + 0.05) * 401074 * data[1] / (data[1] + data[3]))
#print(0.05 * 401074)
# ACT: 0.015567657065637092
#data = (
#    np.mean(toolkit.load_mask("planck_modified_point").map),
#    np.mean(toolkit.load_mask("planck_modified_galactic").map)
#)
#final = (1 - data[0]) / data[1]

# ACT: 0.009859934289099422
#401074
# Planck - 0.0138443493342983
print(final)


#  Results:
# ACT point: 0.01771962453092031
# Planck point: 0.0138443493342983