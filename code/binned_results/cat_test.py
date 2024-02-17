from toolkit import toolkit, weights
import numpy as np
import astropy


def filter(redshift):
    global results
    results.append(redshift)
    return redshift > 0.45


results = []
#cat = toolkit.load_catalogue("sdss")
#cat = toolkit.StarCatalogue("../../data/dr8_run_redmapper_v5.10_lgt5_catalog.fit", hdu=1)


max_z = 0.5
min_z = 0.1

cat = astropy.io.fits.open("../../data/dr8_run_redmapper_v5.10_lgt5_catalog.fit", hdu=1)

print(cat[1].columns)

#cat.load_with_selection(filter, ["Z_LAMBDA"], lon_lat=True)
#cat.load_with_selection(filter, ["ZRED"], lon_lat=True)
#cat.load_with_selection(filter, ["R_LAMBDA"], lon_lat=True)
#results = np.array(results)
#print(np.min(results), np.max(results))
#median = np.median(results)
#upper = np.median(results[results > median])
#lower = np.median(results[results < median])

#print(upper, median, lower)
