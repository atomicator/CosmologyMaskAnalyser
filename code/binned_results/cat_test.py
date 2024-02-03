from toolkit import toolkit, weights
import numpy as np


def filter(redshift):
    global results
    results.append(redshift)
    return redshift > 0.45


results = []
#cat = toolkit.load_catalogue("sdss")
cat = toolkit.StarCatalogue("../../data/dr8_run_redmapper_v5.10_lgt5_catalog.fit", hdu=1)

max_z = 0.5
min_z = 0.1


#cat.load_with_selection(filter, ["Z_LAMBDA"], lon_lat=True)
cat.load_with_selection(filter, ["ZRED"], lon_lat=True)
#cat.load_with_selection(filter, ["R_LAMBDA"], lon_lat=True)
results = np.array(results)
print(np.min(results), np.max(results))
median = np.median(results)
upper = np.median(results[results > median])
lower = np.median(results[results < median])

print(upper, median, lower)
