#!/bin/bash
# Generate random distributions
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=2 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_2.fits" &
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=8 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_8.fits" &
#python3 "./../code/binned_results/weighted_random_distribution.py" --raise_path=1 --nside=32 --target=3000000 --save_path="./../data/binned_random_catalogues/sdss_nside_32.fits" &
wait
# Write new graphs

#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="full_sky" --weight_function="excess_measurement" --save_path="./full_sky_excess.pdf" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="full_sky" --weight_function="density_weighting" --save_path="./full_sky_density.pdf" &

#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_random" --weight_function="excess_measurement" --save_path="./sdss_random_excess.pdf" &
#python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_random" --weight_function="density_weighting" --save_path="./sdss_random_density.pdf" &

#nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss" --weight_function="excess_measurement" --save_path="./sdss_data_excess.pdf"
#nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss" --weight_function="density_weighting" --save_path="./sdss_data_density.pdf"

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./min_z_3.pdf" --min_z=0.3 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./min_z_4.pdf" --min_z=0.4 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./min_z_5.pdf" --min_z=0.5
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./z_bin_1.pdf" --max_z=0.3 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./z_bin_2.pdf" --min_z=0.3 --max_z=0.4 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./z_bin_3.pdf" --min_z=0.4 --max_z=0.5
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./min_r_6.pdf" --min_r=0.6 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./min_r_65.pdf" --min_r=0.65 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./min_r_7.pdf" --min_r=0.7
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./r_bin_1.pdf" --max_r=0.6 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./r_bin_2.pdf" --min_r=0.6 --max_r=0.65 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="density_weighting" --save_path="./r_bin_3.pdf" --min_r=0.65 --max_r=0.7 &
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_min_z_3.pdf" --min_z=0.3 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_min_z_4.pdf" --min_z=0.4 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_min_z_5.pdf" --min_z=0.5
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_z_bin_1.pdf" --max_z=0.3 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_z_bin_2.pdf" --min_z=0.3 --max_z=0.4 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_z_bin_3.pdf" --min_z=0.4 --max_z=0.5
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_min_r_6.pdf" --min_r=0.6 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_min_r_65.pdf" --min_r=0.65 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_min_r_7.pdf" --min_r=0.7
wait

nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_r_bin_1.pdf" --max_r=0.6 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_r_bin_2.pdf" --min_r=0.6 --max_r=0.65 &
nice -20 python3 "./../code/binned_results/debug_script.py" --raise_path=1 --catalogue="sdss_filtered"  \
--weight_function="excess_measurement" --save_path="./excess_r_bin_3.pdf" --min_r=0.65 --max_r=0.7 &
wait
