# A bunch of hardcoded data files - adjust as necessary
from toolkit.toolkit import HealpyMask, PixellMask, CombinationMask, ClusterCatalogue
import numpy as np
import healpy as hp
import toolkit.filters as filters
import os

def load_mask(mask, raise_dir=2, nside=8, invert=False, lon_shift=None):
    value = None
    if mask == "planck_galactic":
        value = HealpyMask("../" * raise_dir + "data/planck_galactic_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "planck_point" or mask == "planck_modified_total":
        value = HealpyMask("../" * raise_dir + "data/HFI_PCCS_SZ-selfunc-inter-cosmo_2.02.fits", partial=False,
                           mask_using_latlon=True)
    elif mask == "planck_survey":
        value = HealpyMask("../" * raise_dir + "data/planck_survey_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "sdss_mask":
        value = HealpyMask("../" * raise_dir + "data/redmapper_dr8_public_v6.3_zmask.fits", mask_using_latlon=False,
                           hdu=1, partial=True, lon_shift=lon_shift)
        value.map[value.map > 0.4] = 1.0
        value.map[value.map < 0.3] = 0
        rotator = hp.Rotator(coord=["C", "G"])
        value.map = rotator.rotate_map_pixel(value.map)
        value.mask_using_latlon = True
        value.map[value.map > 0.5] = 1.0
        value.map[value.map < 0.5] = 0
        # value.map = (value.map - 1) * -1
    elif mask == "planck_modified_galactic":
        point_mask = load_mask("planck_point", raise_dir=raise_dir)
        galactic_mask = load_mask("planck_galactic", raise_dir=raise_dir)
        compound_map = load_mask("planck_point", raise_dir=raise_dir)
        value = load_mask("planck_point", raise_dir=raise_dir)
        compound_map.map = 2 * point_mask.map + galactic_mask.map
        value.map = np.float64(np.bitwise_not(compound_map.map == 0))
    elif mask == "planck_modified_point":
        point_mask = load_mask("planck_point", raise_dir=raise_dir)
        galactic_mask = load_mask("planck_galactic", raise_dir=raise_dir)
        compound_map = load_mask("planck_point", raise_dir=raise_dir)
        value = load_mask("planck_point", raise_dir=raise_dir)
        compound_map.map = 2 * point_mask.map + galactic_mask.map
        value.map = np.float64(np.bitwise_not(compound_map.map == 1))
    elif mask == "comparison_sdss_planck_galactic":
        planck_only = HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_2.fits")
        neither = HealpyMask("../../data/cached_results/sdss_mask_planck_point_32_4.fits")
        data = planck_only.map / (planck_only.map + neither.map + 1e-30)
        data = hp.pixelfunc.ud_grade(data, nside)
        planck_only.map = data
        planck_only.NSIDE = nside
        planck_only.NPIX = hp.nside2npix(nside)
        return planck_only
    elif mask.lower() == "act":
        value = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, mask_using_latlon=False, invert=invert)
    elif mask == "sdss_planck_point_only":
        mask1 = load_mask("sdss_mask", raise_dir)
        mask2 = load_mask("planck_modified_point", raise_dir)
        mask2.map = 1 - mask2.map
        value = CombinationMask(mask1, mask2)
    elif mask == "planck_galactic_test":
        value = HealpyMask("../" * raise_dir + "data/planck_galactic_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "planck_point_test":
        value = HealpyMask("../" * raise_dir + "data/planck_point_mask.fits", partial=True, mask_using_latlon=True)
    elif mask == "act_galactic":
        act_mask = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
        act_graph_filter = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
        act_graph_filter.map = np.load("../" * raise_dir + "data/act_galactic_mask_array.npy")
        test = filters.SquareMaskFilter((-10, 20), (40, 100), lonlat=False)
        final_filter = CombinationMask(act_mask, test, use_and=False, invert=False)
        value = CombinationMask(act_graph_filter, final_filter, use_and=True, invert=False)
    elif mask == "act_point":
        act_mask = PixellMask("../" * raise_dir + "data/ACT_mask.fits", hdu=1, invert=False, mask_using_latlon=False)
        galactic_mask = load_mask("act_galactic", raise_dir, nside, invert)
        galactic_mask.invert = True
        value = CombinationMask(act_mask, galactic_mask, use_and=False)
    elif mask == "spt_total":
        mask_files = os.listdir("../" * raise_dir + "./data/test")
        temp_path = mask_files.pop(0)
        mask = PixellMask("../" * raise_dir + "./data/test/" + temp_path,
                          mask_using_latlon=False, init_val=0, suppress_warnings=True)
        while mask_files:
            temp_path = mask_files.pop(0)
            temp_mask = PixellMask("../" * raise_dir + "./data/test/" + temp_path,
                                   mask_using_latlon=False, init_val=0, suppress_warnings=True)
            mask = CombinationMask(mask, temp_mask, lon_shift=0, use_and=False)
        value = mask
    elif mask == "spt_point":
        mask_files = os.listdir("../" * raise_dir + "./data/test")
        overwrite_path = "../" * raise_dir + "./data/test2/"
        temp_path = mask_files.pop(0)
        mask = PixellMask("../" * raise_dir + "./data/test/" + temp_path,
                                  mask_using_latlon=False, init_val=1, suppress_warnings=True)
        mask.map = np.load(overwrite_path + "point_" + temp_path + ".npy")
        while mask_files:
            temp_path = mask_files.pop(0)
            temp_mask = PixellMask("../" * raise_dir + "./data/test/" + temp_path,
                                           mask_using_latlon=False, init_val=1, suppress_warnings=True)
            temp_mask.map = np.load(overwrite_path + "point_" + temp_path + ".npy")
            mask = CombinationMask(mask, temp_mask, lon_shift=0, use_and=True)
        value = mask
    elif mask == "spt_galactic":
        mask_files = os.listdir("../" * raise_dir + "./data/test")
        overwrite_path = "../" * raise_dir + "./data/test2/"
        temp_path = mask_files.pop(0)
        mask = PixellMask("../" * raise_dir + "./data/test/" + temp_path,
                                  mask_using_latlon=False, init_val=0, suppress_warnings=True)
        mask.map = (mask.map + np.load(overwrite_path + "point_" + temp_path + ".npy") + 1) % 2
        while mask_files:
            temp_path = mask_files.pop(0)
            temp_mask = PixellMask("../" * raise_dir + "./data/test/" + temp_path,
                                           mask_using_latlon=False, init_val=0, suppress_warnings=True)
            temp_mask.map = (temp_mask.map + np.load(overwrite_path + "point_" + temp_path + ".npy") + 1) % 2
            mask = CombinationMask(mask, temp_mask, lon_shift=0, use_and=False)
        value = mask
    else:
        raise ValueError(f"{mask} is not a recognised mask")
    return value


def load_catalogue(cat, raise_dir=2):
    if cat == "sdss":
        data = ClusterCatalogue(raise_dir * "../" + "data/sdss_catalogue.fits", hdu=1)
        data.load_data()
    else:
        raise ValueError(f"{cat} is not a recognised catalogue")
    return data
